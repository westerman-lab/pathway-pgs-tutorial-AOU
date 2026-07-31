#!/usr/bin/env python3
"""Core utilities for the interactive pathway-PGS tutorial.

The notebook validates inputs, builds a reproducible scoring command, executes
the configured calculation engine, and summarizes aggregate outputs.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


REQUIRED_GWAS_KEYS = ("chr", "bp", "snp", "a1", "a2", "stat", "p")
AUTOSOMES = tuple(str(i) for i in range(1, 23))


@dataclass
class WorkflowConfig:
    """Configuration shared by the notebook UI and scoring command builder."""

    project_name: str = "t2d_reactome_pathway_pgs"
    genome_build: str = "GRCh38"
    base_gwas: str = ""
    base_separator: str = "auto"
    gwas_columns: dict[str, str] = field(
        default_factory=lambda: {
            "chr": "CHR",
            "bp": "BP",
            "snp": "SNP",
            "a1": "A1",
            "a2": "A2",
            "stat": "BETA",
            "p": "P",
        }
    )
    statistic_type: str = "BETA"
    target_type: str = "bed"
    target_prefix: str = ""
    target_list: str = ""
    target_keep_file: str = ""
    pathway_input_mode: str = "gtf_gmt"
    gtf_file: str = ""
    gmt_file: str = ""
    gmt_second_column_is_description: bool = True
    snp_set_file: str = ""
    window_5_bp: int = 35_000
    window_3_bp: int = 10_000
    pvalue_thresholds: list[float] = field(default_factory=lambda: [1.0])
    clump_kb: int = 1_000
    clump_r2: float = 0.1
    base_maf_column: str = ""
    base_maf_min: float | None = None
    base_info_column: str = ""
    base_info_min: float | None = None
    phenotype_file: str = ""
    phenotype_column: str = ""
    binary_target: bool = False
    covariate_file: str = ""
    covariate_columns: str = ""
    no_regression: bool = True
    set_permutations: int = 0
    score_method: str = "sum"
    threads: int = 4
    prsice_r: str = ""
    prsice_binary: str = ""
    rscript: str = "Rscript"
    output_dir: str = "analysis/results/pathway_prs_tutorial"
    output_prefix: str = "pathway_pgs"
    controlled_workspace_acknowledged: bool = False

    @classmethod
    def from_json(cls, path: str | Path) -> "WorkflowConfig":
        with Path(path).open() as handle:
            return cls(**json.load(handle))

    def to_json(self, path: str | Path) -> Path:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w") as handle:
            json.dump(asdict(self), handle, indent=2, sort_keys=True)
        return destination


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r")


def _detect_separator(path: Path, configured: str = "auto") -> str | None:
    if configured != "auto":
        return {"tab": "\t", "comma": ",", "whitespace": None}.get(
            configured, configured
        )
    with _open_text(path) as handle:
        first = handle.readline()
    if "\t" in first:
        return "\t"
    if "," in first:
        return ","
    return None


def _metadata_skiprows(path: Path) -> int:
    """Skip VCF-style metadata while preserving a possible '#CHROM' header."""
    count = 0
    with _open_text(path) as handle:
        for raw in handle:
            if raw.startswith("##"):
                count += 1
                continue
            break
    return count


def _read_preview(path: Path, separator: str | None, nrows: int = 5) -> pd.DataFrame:
    options = {"nrows": nrows, "skiprows": _metadata_skiprows(path)}
    if separator is None:
        return pd.read_csv(path, sep=r"\s+", **options)
    return pd.read_csv(path, sep=separator, **options)


def infer_gwas_schema(path: str | Path) -> dict[str, Any]:
    """Infer common GWAS column names and the effect-statistic representation."""
    source = Path(path).expanduser()
    separator = _detect_separator(source, "auto")
    preview = _read_preview(source, separator, nrows=3)
    columns = [str(column) for column in preview.columns]

    def normalized(value: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", value.lower().lstrip("#"))

    lookup: dict[str, str] = {}
    for column in columns:
        lookup.setdefault(normalized(column), column)

    def choose(*candidates: str) -> str:
        for candidate in candidates:
            match = lookup.get(normalized(candidate))
            if match is not None:
                return match
        return ""

    inferred = {
        "chr": choose("CHR", "CHROM", "CHROMOSOME"),
        "bp": choose("BP", "POS", "POSITION"),
        "snp": choose("SNP", "RSID", "RSIDS", "MARKERNAME", "ID", "VARIANT_ID"),
        "a1": choose(
            "A1", "EFFECT_ALLELE", "EA", "ALT", "ALLELE1", "RISK_ALLELE"
        ),
        "a2": choose(
            "A2", "OTHER_ALLELE", "NEA", "REF", "ALLELE2", "NON_EFFECT_ALLELE"
        ),
        "p": choose(
            "P", "PVAL", "PVALUE", "P_VALUE", "FIXED_EFFECTS_P_VALUE"
        ),
    }
    beta_column = choose(
        "BETA", "EFFECT", "EFFECT_SIZE", "LOG_OR", "FIXED_EFFECTS_BETA"
    )
    odds_ratio_column = choose("OR", "ODDS_RATIO")
    if beta_column:
        inferred["stat"] = beta_column
        statistic_type = "BETA"
    elif odds_ratio_column:
        inferred["stat"] = odds_ratio_column
        statistic_type = "OR"
    else:
        inferred["stat"] = ""
        statistic_type = "BETA"

    missing = [key for key in REQUIRED_GWAS_KEYS if not inferred.get(key)]
    return {
        "separator": (
            "tab" if separator == "\t" else "comma" if separator == "," else "whitespace"
        ),
        "statistic_type": statistic_type,
        "gwas_columns": inferred,
        "maf_column": choose("MAF", "MINOR_ALLELE_FREQUENCY"),
        "info_column": choose("INFO", "INFO_SCORE", "IMPUTATION_INFO"),
        "available_columns": columns,
        "missing_required_fields": missing,
    }


def _path_status(path_text: str, label: str, required: bool = True) -> dict[str, Any]:
    path = Path(path_text).expanduser() if path_text else None
    exists = bool(path and path.exists())
    return {
        "check": label,
        "status": "PASS" if exists or not required else "FAIL",
        "detail": str(path) if path else ("optional" if not required else "not provided"),
    }


def _resolve_executable(value: str) -> str | None:
    expanded = str(Path(value).expanduser()) if value else ""
    if expanded and Path(expanded).exists():
        return expanded
    return shutil.which(value) if value else None


def read_gmt_summary(
    path: str | Path,
    second_column_is_description: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return pathway-level sizes and a small pathway/gene preview."""
    source = Path(path).expanduser()
    rows: list[dict[str, Any]] = []
    preview: list[dict[str, str]] = []
    with _open_text(source) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"Invalid GMT row {line_number}: fewer than 2 fields")
            pathway = fields[0]
            gene_start = 2 if second_column_is_description else 1
            genes = [g for g in fields[gene_start:] if g]
            rows.append(
                {
                    "pathway": pathway,
                    "n_genes": len(set(genes)),
                    "line_number": line_number,
                }
            )
            if len(preview) < 20:
                preview.extend(
                    {"pathway": pathway, "gene": gene}
                    for gene in genes[: min(5, len(genes))]
                )
    summary = pd.DataFrame(rows)
    return summary, pd.DataFrame(preview)


def read_gtf_summary(path: str | Path, maximum_rows: int = 250_000) -> dict[str, Any]:
    """Inspect a GTF without loading the full annotation into memory."""
    source = Path(path).expanduser()
    features: dict[str, int] = {}
    chromosomes: set[str] = set()
    rows = 0
    with _open_text(source) as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            rows += 1
            chromosomes.add(fields[0])
            features[fields[2]] = features.get(fields[2], 0) + 1
            if rows >= maximum_rows:
                break
    return {
        "rows_inspected": rows,
        "inspection_truncated": rows >= maximum_rows,
        "n_chromosome_labels": len(chromosomes),
        "top_features": sorted(features.items(), key=lambda item: -item[1])[:10],
    }


def _parse_gtf_attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in value.split(";"):
        item = item.strip()
        if not item:
            continue
        key, separator, raw_value = item.partition(" ")
        if separator:
            attributes[key] = raw_value.strip().strip('"')
    return attributes


def _normalize_chromosome(value: object) -> str:
    chromosome = str(value).strip()
    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]
    return chromosome


def _pathway_member_tokens(
    gmt_file: str | Path,
    second_column_is_description: bool,
) -> set[str]:
    tokens: set[str] = set()
    with _open_text(Path(gmt_file).expanduser()) as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            gene_start = 2 if second_column_is_description else 1
            if len(fields) > gene_start:
                tokens.update(token for token in fields[gene_start:] if token)
    return tokens


def _pathway_gene_intervals(config: WorkflowConfig) -> pd.DataFrame:
    members = _pathway_member_tokens(
        config.gmt_file,
        config.gmt_second_column_is_description,
    )
    rows: list[dict[str, Any]] = []
    with _open_text(Path(config.gtf_file).expanduser()) as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = _parse_gtf_attributes(fields[8])
            gene_id = attributes.get("gene_id", "").split(".")[0]
            gene_name = attributes.get("gene_name", "")
            if gene_id not in members and gene_name not in members:
                continue
            chromosome = _normalize_chromosome(fields[0])
            if chromosome not in AUTOSOMES:
                continue
            start = int(fields[3])
            end = int(fields[4])
            if fields[6] == "-":
                start = max(1, start - config.window_3_bp)
                end += config.window_5_bp
            else:
                start = max(1, start - config.window_5_bp)
                end += config.window_3_bp
            rows.append(
                {
                    "chromosome": chromosome,
                    "start": start,
                    "end": end,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                }
            )
    intervals = pd.DataFrame(rows)
    if intervals.empty:
        raise ValueError(
            "No GMT members matched GTF gene_id or gene_name values. "
            "Check the genome build and gene identifier type."
        )
    return intervals


def _merge_intervals(intervals: pd.DataFrame) -> pd.DataFrame:
    merged: list[dict[str, int | str]] = []
    for chromosome, group in intervals.groupby("chromosome", sort=False):
        ordered = group.sort_values(["start", "end"])
        current_start: int | None = None
        current_end: int | None = None
        for start, end in ordered[["start", "end"]].itertuples(index=False):
            start = int(start)
            end = int(end)
            if current_start is None:
                current_start, current_end = start, end
            elif start <= int(current_end) + 1:
                current_end = max(int(current_end), end)
            else:
                merged.append(
                    {
                        "chromosome": chromosome,
                        "start": int(current_start),
                        "end": int(current_end),
                    }
                )
                current_start, current_end = start, end
        if current_start is not None:
            merged.append(
                {
                    "chromosome": chromosome,
                    "start": int(current_start),
                    "end": int(current_end),
                }
            )
    return pd.DataFrame(merged)


def build_pathway_variant_union_bed(
    config: WorkflowConfig,
    output_file: str | Path,
    *,
    chunksize: int = 250_000,
) -> dict[str, Any]:
    """Build an aggregate-only BED1 union for a targeted PGEN conversion.

    Base-GWAS positions are retained when they fall in any pathway-gene window.
    The resulting three-column file can be passed to PLINK 2 with
    ``--extract bed1``. Definitive allele harmonization, pathway-specific
    clumping, and scoring occur after conversion.
    """
    if config.pathway_input_mode != "gtf_gmt":
        raise ValueError("The BED1 union builder requires GTF + GMT pathway inputs.")
    intervals = _pathway_gene_intervals(config)
    merged = _merge_intervals(intervals)
    interval_lookup: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chromosome, group in merged.groupby("chromosome"):
        ordered = group.sort_values("start")
        interval_lookup[str(chromosome)] = (
            ordered["start"].to_numpy(dtype=np.int64),
            ordered["end"].to_numpy(dtype=np.int64),
        )

    source = Path(config.base_gwas).expanduser()
    separator = _detect_separator(source, config.base_separator)
    read_kwargs: dict[str, Any] = {
        "usecols": [
            config.gwas_columns["chr"],
            config.gwas_columns["bp"],
        ],
        "chunksize": chunksize,
        "skiprows": _metadata_skiprows(source),
    }
    if separator is None:
        read_kwargs.update({"sep": r"\s+", "engine": "python"})
    else:
        read_kwargs["sep"] = separator

    selected: list[pd.DataFrame] = []
    n_gwas_rows = 0
    for chunk in pd.read_csv(source, **read_kwargs):
        n_gwas_rows += len(chunk)
        chromosome = chunk[config.gwas_columns["chr"]].map(_normalize_chromosome)
        position = pd.to_numeric(
            chunk[config.gwas_columns["bp"]], errors="coerce"
        )
        for chrom, indices in chromosome.groupby(chromosome).groups.items():
            if chrom not in interval_lookup:
                continue
            positions = position.loc[indices].to_numpy(dtype=float)
            valid = np.isfinite(positions)
            integer_positions = np.zeros(len(positions), dtype=np.int64)
            integer_positions[valid] = positions[valid].astype(np.int64)
            starts, ends = interval_lookup[chrom]
            interval_index = np.searchsorted(
                starts, integer_positions, side="right"
            ) - 1
            inside = valid & (interval_index >= 0)
            valid_index = interval_index.clip(min=0)
            inside &= integer_positions <= ends[valid_index]
            if inside.any():
                selected.append(
                    pd.DataFrame(
                        {
                            "chromosome": chrom,
                            "start": integer_positions[inside],
                            "end": integer_positions[inside],
                        }
                    )
                )

    if not selected:
        raise ValueError(
            "No GWAS positions overlapped pathway-gene windows. "
            "Check chromosome labels and genome builds."
        )
    union = (
        pd.concat(selected, ignore_index=True)
        .drop_duplicates(["chromosome", "start", "end"])
        .sort_values(["chromosome", "start"], key=lambda column: (
            pd.to_numeric(column, errors="coerce")
            if column.name == "chromosome"
            else column
        ))
    )
    destination = Path(output_file).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    union.to_csv(destination, sep="\t", header=False, index=False)
    audit = {
        "status": "PASS",
        "pathway_member_tokens": len(
            _pathway_member_tokens(
                config.gmt_file,
                config.gmt_second_column_is_description,
            )
        ),
        "matched_gtf_genes": len(intervals),
        "merged_gene_windows": len(merged),
        "gwas_rows_inspected": n_gwas_rows,
        "unique_gwas_positions_in_pathway_windows": len(union),
        "output_bed1": str(destination),
        "person_level_data_read": False,
    }
    with destination.with_suffix(destination.suffix + ".audit.json").open("w") as handle:
        json.dump(audit, handle, indent=2)
    return audit


def build_variant_mapping_union_bed(
    config: WorkflowConfig,
    mapping_file: str | Path,
    output_file: str | Path,
    *,
    chunksize: int = 250_000,
) -> dict[str, Any]:
    """Build a BED1 union from a user-provided SNP-to-gene mapping.

    The mapping may contain chromosome and position columns directly, or a
    variant-ID column whose values match the configured GWAS SNP column.
    Only aggregate variant coordinates are written.
    """
    source = Path(mapping_file).expanduser()
    separator = _detect_separator(source, "auto")
    mapping = pd.read_csv(source, sep=r"\s+" if separator is None else separator)
    normalized = {
        re.sub(r"[^a-z0-9]", "", str(column).lower()): str(column)
        for column in mapping.columns
    }

    def choose(*names: str) -> str:
        for name in names:
            key = re.sub(r"[^a-z0-9]", "", name.lower())
            if key in normalized:
                return normalized[key]
        return ""

    variant_column = choose("variant", "variant_id", "snp", "rsid", "rsids")
    chromosome_column = choose("chr", "chrom", "chromosome", "contig")
    position_column = choose("bp", "pos", "position", "base_pair_location")
    selected: list[pd.DataFrame] = []
    n_gwas_rows = 0

    if chromosome_column and position_column:
        chromosome = mapping[chromosome_column].map(_normalize_chromosome)
        position = pd.to_numeric(mapping[position_column], errors="coerce")
        valid = chromosome.ne("") & position.notna()
        selected.append(
            pd.DataFrame(
                {
                    "chromosome": chromosome.loc[valid],
                    "start": position.loc[valid].astype(np.int64),
                    "end": position.loc[valid].astype(np.int64),
                }
            )
        )
        mapping_mode = "chromosome_position"
    elif variant_column:
        variant_ids = set(mapping[variant_column].dropna().astype(str))
        gwas_path = Path(config.base_gwas).expanduser()
        gwas_separator = _detect_separator(gwas_path, config.base_separator)
        read_kwargs: dict[str, Any] = {
            "usecols": [
                config.gwas_columns["chr"],
                config.gwas_columns["bp"],
                config.gwas_columns["snp"],
            ],
            "chunksize": chunksize,
            "skiprows": _metadata_skiprows(gwas_path),
        }
        if gwas_separator is None:
            read_kwargs.update({"sep": r"\s+", "engine": "python"})
        else:
            read_kwargs["sep"] = gwas_separator
        for chunk in pd.read_csv(gwas_path, **read_kwargs):
            n_gwas_rows += len(chunk)
            keep = chunk[config.gwas_columns["snp"]].astype(str).isin(variant_ids)
            if not keep.any():
                continue
            chromosome = chunk.loc[keep, config.gwas_columns["chr"]].map(
                _normalize_chromosome
            )
            position = pd.to_numeric(
                chunk.loc[keep, config.gwas_columns["bp"]], errors="coerce"
            )
            valid = chromosome.ne("") & position.notna()
            selected.append(
                pd.DataFrame(
                    {
                        "chromosome": chromosome.loc[valid],
                        "start": position.loc[valid].astype(np.int64),
                        "end": position.loc[valid].astype(np.int64),
                    }
                )
            )
        mapping_mode = "variant_id_matched_to_gwas"
    else:
        raise ValueError(
            "Custom mapping must contain either recognizable chromosome and "
            "position columns or a recognizable variant-ID column."
        )

    if not selected:
        raise ValueError(
            "No custom-mapping variants could be located. Check the mapping "
            "columns, variant IDs, and genome build."
        )
    union = (
        pd.concat(selected, ignore_index=True)
        .drop_duplicates(["chromosome", "start", "end"])
        .sort_values(
            ["chromosome", "start"],
            key=lambda column: (
                pd.to_numeric(column, errors="coerce")
                if column.name == "chromosome"
                else column
            ),
        )
    )
    destination = Path(output_file).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    union.to_csv(destination, sep="\t", header=False, index=False)
    audit = {
        "status": "PASS",
        "mapping_mode": mapping_mode,
        "mapping_rows": len(mapping),
        "gwas_rows_inspected": n_gwas_rows,
        "unique_mapping_positions": len(union),
        "output_bed1": str(destination),
        "person_level_data_read": False,
    }
    with destination.with_suffix(destination.suffix + ".audit.json").open("w") as handle:
        json.dump(audit, handle, indent=2)
    return audit


def _read_candidate_bed1(path: str | Path) -> dict[str, set[int]]:
    table = pd.read_csv(
        Path(path).expanduser(),
        sep=r"\s+",
        header=None,
        usecols=[0, 1],
        names=["chromosome", "position"],
    )
    table["chromosome"] = table["chromosome"].map(_normalize_chromosome)
    table["position"] = pd.to_numeric(table["position"], errors="coerce")
    table = table.dropna().astype({"position": np.int64})
    return {
        chromosome: set(group["position"].tolist())
        for chromosome, group in table.groupby("chromosome")
    }


def _pvar_path(prefix: str) -> Path:
    candidates = [
        Path(prefix + ".pvar"),
        Path(prefix + ".pvar.gz"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"No .pvar or .pvar.gz file found for prefix: {prefix}")


def _pvar_candidate_alleles(
    *,
    pgen_pattern: str,
    chromosomes: Iterable[int],
    candidate_positions: dict[str, set[int]],
) -> tuple[dict[str, str], dict[str, Any]]:
    allele_map: dict[str, str] = {}
    n_position_rows = 0
    n_non_biallelic_or_non_snp = 0
    for chromosome_number in chromosomes:
        chromosome = str(chromosome_number)
        positions = candidate_positions.get(chromosome, set())
        if not positions:
            continue
        prefix = (
            pgen_pattern.replace("{chr}", chromosome).replace("#", chromosome)
        )
        source = _pvar_path(prefix)
        with _open_text(source) as handle:
            header: list[str] | None = None
            for raw in handle:
                if raw.startswith("##"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) == 1:
                    fields = raw.rstrip("\n").split()
                if raw.startswith("#"):
                    header = [value.lstrip("#") for value in fields]
                    continue
                if header is None:
                    raise ValueError(f"PVAR header not found: {source}")
                row = dict(zip(header, fields))
                position = int(row["POS"])
                if position not in positions:
                    continue
                n_position_rows += 1
                ref = row["REF"].upper()
                alt = row["ALT"].upper()
                if (
                    len(ref) != 1
                    or len(alt) != 1
                    or ref not in "ACGT"
                    or alt not in "ACGT"
                ):
                    n_non_biallelic_or_non_snp += 1
                    continue
                raw_chromosome = row["CHROM"]
                first, second = sorted((ref, alt))
                key = f"{chromosome}:{position}:{first}:{second}"
                target_id = f"{raw_chromosome}:{position}:{ref}:{alt}"
                if key in allele_map and allele_map[key] != target_id:
                    raise ValueError(
                        f"Ambiguous biallelic target records for {key}: "
                        f"{allele_map[key]} and {target_id}"
                    )
                allele_map[key] = target_id
    audit = {
        "candidate_pvar_rows": n_position_rows,
        "biallelic_snp_target_records": len(allele_map),
        "non_biallelic_or_non_snp_pvar_rows": n_non_biallelic_or_non_snp,
    }
    return allele_map, audit


def _canonicalize_harmonized_gwas(
    frame: pd.DataFrame,
    config: WorkflowConfig,
) -> pd.DataFrame:
    """Return a compact GWAS table with PRSice-safe canonical headers."""
    columns = config.gwas_columns
    statistic_name = "OR" if config.statistic_type.upper() == "OR" else "BETA"
    canonical = pd.DataFrame(
        {
            "CHR": frame[columns["chr"]].map(_normalize_chromosome),
            "BP": pd.to_numeric(frame[columns["bp"]], errors="coerce"),
            "SNP": frame[columns["snp"]].astype(str),
            "A1": frame[columns["a1"]].astype(str).str.upper(),
            "A2": frame[columns["a2"]].astype(str).str.upper(),
            statistic_name: pd.to_numeric(frame[columns["stat"]], errors="coerce"),
            "P": pd.to_numeric(frame[columns["p"]], errors="coerce"),
        }
    )
    if "ORIGINAL_SNP_ID" in frame.columns:
        canonical["ORIGINAL_SNP_ID"] = frame["ORIGINAL_SNP_ID"].astype(str)
    if config.base_maf_column and config.base_maf_column in frame.columns:
        canonical["MAF"] = pd.to_numeric(
            frame[config.base_maf_column], errors="coerce"
        )
    if config.base_info_column and config.base_info_column in frame.columns:
        canonical["INFO"] = pd.to_numeric(
            frame[config.base_info_column], errors="coerce"
        )
    return canonical


def harmonize_gwas_ids_to_pgen(
    config: WorkflowConfig,
    *,
    pgen_pattern: str,
    chromosomes: Iterable[int],
    candidate_bed1_file: str | Path,
    output_gwas: str | Path,
    chunksize: int = 250_000,
) -> dict[str, Any]:
    """Rewrite base-GWAS SNP IDs to match coordinate/allele PGEN IDs.

    Matching uses chromosome, position, and the unordered A1/A2 allele pair.
    Effect alleles and effect sizes are not changed. The adapter later assigns
    the corresponding ``CHROM:POS:REF:ALT`` IDs to converted target variants.
    """
    candidate_positions = _read_candidate_bed1(candidate_bed1_file)
    chromosomes = list(chromosomes)
    allele_map, pvar_audit = _pvar_candidate_alleles(
        pgen_pattern=pgen_pattern,
        chromosomes=chromosomes,
        candidate_positions=candidate_positions,
    )
    source = Path(config.base_gwas).expanduser()
    separator = _detect_separator(source, config.base_separator)
    read_kwargs: dict[str, Any] = {
        "chunksize": chunksize,
        "skiprows": _metadata_skiprows(source),
    }
    if separator is None:
        read_kwargs.update({"sep": r"\s+", "engine": "python"})
    else:
        read_kwargs["sep"] = separator

    chromosome_column = config.gwas_columns["chr"]
    position_column = config.gwas_columns["bp"]
    snp_column = config.gwas_columns["snp"]
    a1_column = config.gwas_columns["a1"]
    a2_column = config.gwas_columns["a2"]
    p_column = config.gwas_columns["p"]
    matched_chunks: list[pd.DataFrame] = []
    n_gwas_rows = 0
    n_valid_snp_rows = 0
    for chunk in pd.read_csv(source, **read_kwargs):
        n_gwas_rows += len(chunk)
        chromosome = chunk[chromosome_column].map(_normalize_chromosome)
        position = pd.to_numeric(chunk[position_column], errors="coerce")
        a1 = chunk[a1_column].astype(str).str.upper()
        a2 = chunk[a2_column].astype(str).str.upper()
        valid = (
            position.notna()
            & a1.str.match(r"^[ACGT]$")
            & a2.str.match(r"^[ACGT]$")
            & (a1 != a2)
        )
        n_valid_snp_rows += int(valid.sum())
        first = a1.where(a1 <= a2, a2)
        second = a2.where(a1 <= a2, a1)
        keys = (
            chromosome
            + ":"
            + position.fillna(-1).astype(np.int64).astype(str)
            + ":"
            + first
            + ":"
            + second
        )
        target_ids = keys.map(allele_map)
        keep = valid & target_ids.notna()
        if keep.any():
            selected = chunk.loc[keep].copy()
            if "ORIGINAL_SNP_ID" not in selected.columns:
                selected.insert(
                    selected.columns.get_loc(snp_column) + 1,
                    "ORIGINAL_SNP_ID",
                    selected[snp_column].astype(str),
                )
            selected[snp_column] = target_ids.loc[keep].to_numpy()
            matched_chunks.append(selected)

    if not matched_chunks:
        raise ValueError(
            "No GWAS variants matched candidate PVAR records by position and alleles."
        )
    harmonized = pd.concat(matched_chunks, ignore_index=True)
    harmonized["_numeric_p"] = pd.to_numeric(
        harmonized[p_column], errors="coerce"
    ).fillna(np.inf)
    n_before_deduplication = len(harmonized)
    harmonized = (
        harmonized.sort_values("_numeric_p")
        .drop_duplicates(snp_column, keep="first")
        .drop(columns="_numeric_p")
    )
    harmonized = _canonicalize_harmonized_gwas(harmonized, config)
    destination = Path(output_gwas).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    harmonized.to_csv(
        destination,
        sep="\t",
        index=False,
        compression="gzip" if destination.suffix == ".gz" else None,
    )
    audit = {
        "status": "PASS",
        "gwas_rows_inspected": n_gwas_rows,
        "valid_biallelic_gwas_snp_rows": n_valid_snp_rows,
        "matched_gwas_rows_before_deduplication": n_before_deduplication,
        "unique_harmonized_target_ids": len(harmonized),
        "duplicate_target_ids_removed": n_before_deduplication - len(harmonized),
        "harmonized_gwas": str(destination),
        "effect_alleles_or_effect_sizes_changed": False,
        "person_level_data_read": False,
        **pvar_audit,
    }
    with destination.with_suffix(destination.suffix + ".audit.json").open("w") as handle:
        json.dump(audit, handle, indent=2)
    return audit


def _bim_candidate_alleles(
    *,
    bed_pattern: str,
    chromosomes: Iterable[int],
    candidate_positions: dict[str, set[int]],
) -> tuple[dict[str, str], dict[str, Any]]:
    """Read only BIM metadata and index target IDs by position and allele pair."""
    allele_map: dict[str, str] = {}
    n_position_rows = 0
    n_non_biallelic_or_non_snp = 0
    for chromosome_number in chromosomes:
        chromosome = str(chromosome_number)
        positions = candidate_positions.get(chromosome, set())
        if not positions:
            continue
        prefix = bed_pattern.replace("{chr}", chromosome).replace("#", chromosome)
        source = Path(prefix + ".bim")
        if not source.exists():
            raise FileNotFoundError(f"Target BIM not found: {source}")
        table = pd.read_csv(
            source,
            sep=r"\s+",
            header=None,
            names=["chromosome", "variant_id", "cm", "position", "a1", "a2"],
            dtype={"chromosome": str, "variant_id": str, "a1": str, "a2": str},
        )
        table["chromosome"] = table["chromosome"].map(_normalize_chromosome)
        table = table[
            (table["chromosome"] == chromosome)
            & table["position"].isin(positions)
        ]
        n_position_rows += len(table)
        for row in table.itertuples(index=False):
            first, second = sorted((str(row.a1).upper(), str(row.a2).upper()))
            if (
                len(first) != 1
                or len(second) != 1
                or first not in "ACGT"
                or second not in "ACGT"
            ):
                n_non_biallelic_or_non_snp += 1
                continue
            key = f"{chromosome}:{int(row.position)}:{first}:{second}"
            if key in allele_map and allele_map[key] != row.variant_id:
                continue
            allele_map[key] = str(row.variant_id)
    return allele_map, {
        "candidate_bim_rows": n_position_rows,
        "biallelic_snp_target_records": len(allele_map),
        "non_biallelic_or_non_snp_bim_rows": n_non_biallelic_or_non_snp,
    }


def harmonize_gwas_ids_to_bed(
    config: WorkflowConfig,
    *,
    bed_pattern: str | None = None,
    chromosomes: Iterable[int] | None = None,
    candidate_bed1_file: str | Path,
    output_gwas: str | Path,
    chunksize: int = 250_000,
) -> dict[str, Any]:
    """Rewrite GWAS SNP IDs to the IDs already present in AoU BIM files."""
    candidate_positions = _read_candidate_bed1(candidate_bed1_file)
    bed_pattern = bed_pattern or config.target_prefix
    if not bed_pattern:
        raise ValueError("A chromosome-wise target BED prefix is required.")
    chromosomes = list(chromosomes or sorted(map(int, candidate_positions)))
    allele_map, bim_audit = _bim_candidate_alleles(
        bed_pattern=bed_pattern,
        chromosomes=chromosomes,
        candidate_positions=candidate_positions,
    )
    source = Path(config.base_gwas).expanduser()
    separator = _detect_separator(source, config.base_separator)
    read_kwargs: dict[str, Any] = {
        "chunksize": chunksize,
        "skiprows": _metadata_skiprows(source),
    }
    if separator is None:
        read_kwargs.update({"sep": r"\s+", "engine": "python"})
    else:
        read_kwargs["sep"] = separator

    destination = Path(output_gwas).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    matched_rows = 0
    inspected_rows = 0
    wrote_header = False
    for chunk in pd.read_csv(source, **read_kwargs):
        inspected_rows += len(chunk)
        chrom = chunk[config.gwas_columns["chr"]].map(_normalize_chromosome)
        pos = pd.to_numeric(chunk[config.gwas_columns["bp"]], errors="coerce")
        a1 = chunk[config.gwas_columns["a1"]].astype(str).str.upper()
        a2 = chunk[config.gwas_columns["a2"]].astype(str).str.upper()
        first = np.minimum(a1, a2)
        second = np.maximum(a1, a2)
        keys = (
            chrom.astype(str)
            + ":"
            + pos.fillna(-1).astype(np.int64).astype(str)
            + ":"
            + first
            + ":"
            + second
        )
        target_ids = keys.map(allele_map)
        keep = target_ids.notna()
        if not keep.any():
            continue
        retained = chunk.loc[keep].copy()
        retained["ORIGINAL_SNP_ID"] = retained[config.gwas_columns["snp"]].values
        retained[config.gwas_columns["snp"]] = target_ids.loc[keep].values
        retained = _canonicalize_harmonized_gwas(retained, config)
        matched_rows += len(retained)
        retained.to_csv(
            destination,
            sep="\t",
            index=False,
            mode="wt" if not wrote_header else "at",
            header=not wrote_header,
            compression="gzip" if destination.suffix == ".gz" else None,
        )
        wrote_header = True
    if not wrote_header:
        raise ValueError("No GWAS variants matched the selected AoU BED files.")
    audit = {
        "status": "PASS",
        "gwas_rows_inspected": inspected_rows,
        "harmonized_gwas_rows": matched_rows,
        "unique_harmonized_target_ids": len(allele_map),
        "harmonized_gwas": str(destination),
        "effect_alleles_or_effect_sizes_changed": False,
        "person_level_data_read": False,
        **bim_audit,
    }
    with destination.with_suffix(destination.suffix + ".audit.json").open("w") as handle:
        json.dump(audit, handle, indent=2)
    return audit


def build_snp_set_from_variant_gene_mapping(
    *,
    mapping_file: str | Path,
    gmt_file: str | Path,
    output_file: str | Path,
    harmonized_gwas_file: str | Path,
    harmonized_snp_column: str = "SNP",
) -> dict[str, Any]:
    """Create pathway SNP sets from a user-provided variant-to-gene table."""
    mapping_path = Path(mapping_file).expanduser()
    separator = _detect_separator(mapping_path, "auto")
    mapping = pd.read_csv(mapping_path, sep=r"\s+" if separator is None else separator)
    normalized = {
        re.sub(r"[^a-z0-9]", "", str(column).lower()): str(column)
        for column in mapping.columns
    }

    def choose(*names: str) -> str:
        for name in names:
            key = re.sub(r"[^a-z0-9]", "", name.lower())
            if key in normalized:
                return normalized[key]
        return ""

    variant_column = choose("variant", "variant_id", "snp", "rsid", "rsids")
    gene_column = choose("gene", "gene_id", "gene_name", "symbol", "ensembl_gene_id")
    chromosome_column = choose("chr", "chrom", "chromosome", "contig")
    position_column = choose("bp", "pos", "position", "base_pair_location")
    if not gene_column or not (
        variant_column or (chromosome_column and position_column)
    ):
        raise ValueError(
            "Custom mapping must contain a gene column plus either a variant-ID "
            "column or chromosome and position columns."
        )

    harmonized_path = Path(harmonized_gwas_file).expanduser()
    separator = _detect_separator(harmonized_path, "auto")
    harmonized = pd.read_csv(
        harmonized_path,
        sep=r"\s+" if separator is None else separator,
    )
    if "ORIGINAL_SNP_ID" not in harmonized.columns:
        harmonized["ORIGINAL_SNP_ID"] = harmonized[harmonized_snp_column]
    id_map = dict(
        zip(
            harmonized["ORIGINAL_SNP_ID"].astype(str),
            harmonized[harmonized_snp_column].astype(str),
        )
    )
    selected_columns = [gene_column]
    selected_columns.extend(
        [
            column
            for column in (variant_column, chromosome_column, position_column)
            if column
        ]
    )
    mapping = mapping[selected_columns].dropna(subset=[gene_column]).copy()
    mapping["gene_key"] = mapping[gene_column].astype(str).str.replace(r"\.\d+$", "", regex=True)
    mapping["target_snp"] = np.nan
    if variant_column:
        mapping["target_snp"] = mapping[variant_column].astype(str).map(id_map)
    if chromosome_column and position_column:
        harmonized_normalized = {
            re.sub(r"[^a-z0-9]", "", str(column).lower()): str(column)
            for column in harmonized.columns
        }

        def choose_harmonized(*names: str) -> str:
            for name in names:
                key = re.sub(r"[^a-z0-9]", "", name.lower())
                if key in harmonized_normalized:
                    return harmonized_normalized[key]
            return ""

        harmonized_chr = choose_harmonized("chr", "chrom", "chromosome", "contig")
        harmonized_pos = choose_harmonized("bp", "pos", "position", "base_pair_location")
        if harmonized_chr and harmonized_pos:
            coordinate_map = dict(
                zip(
                    harmonized[harmonized_chr].map(_normalize_chromosome)
                    + ":"
                    + pd.to_numeric(
                        harmonized[harmonized_pos], errors="coerce"
                    ).fillna(-1).astype(np.int64).astype(str),
                    harmonized[harmonized_snp_column].astype(str),
                )
            )
            mapping_coordinates = (
                mapping[chromosome_column].map(_normalize_chromosome)
                + ":"
                + pd.to_numeric(mapping[position_column], errors="coerce")
                .fillna(-1)
                .astype(np.int64)
                .astype(str)
            )
            mapping["target_snp"] = mapping["target_snp"].fillna(
                mapping_coordinates.map(coordinate_map)
            )
    mapping = mapping.dropna(subset=["target_snp"])
    gene_to_snps = {
        gene: sorted(set(group["target_snp"].astype(str)))
        for gene, group in mapping.groupby("gene_key")
    }

    destination = Path(output_file).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    n_pathways = 0
    n_memberships = 0
    with _open_text(Path(gmt_file).expanduser()) as source, destination.open("w") as handle:
        for raw in source:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            pathway, description, *genes = fields
            variants = sorted(
                {
                    snp
                    for gene in genes
                    for snp in gene_to_snps.get(re.sub(r"\.\d+$", "", gene), [])
                }
            )
            if not variants:
                continue
            handle.write("\t".join([pathway, description, *variants]) + "\n")
            n_pathways += 1
            n_memberships += len(variants)
    audit = {
        "status": "PASS",
        "pathways_with_variants": n_pathways,
        "pathway_variant_memberships": n_memberships,
        "output_snp_set": str(destination),
    }
    return audit


def target_prefixes(config: WorkflowConfig) -> list[str]:
    """Resolve one target prefix or a chromosome-specific target-list file."""
    if config.target_list:
        source = Path(config.target_list).expanduser()
        prefixes = []
        with source.open() as handle:
            for raw in handle:
                value = raw.strip().split()[0] if raw.strip() else ""
                if value:
                    prefixes.append(value)
        return prefixes
    if config.target_prefix:
        if "#" in config.target_prefix:
            return [config.target_prefix.replace("#", chrom) for chrom in AUTOSOMES]
        return [config.target_prefix]
    return []


def validate_inputs(config: WorkflowConfig) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Validate files, tools, columns, and basic parameter consistency."""
    checks: list[dict[str, Any]] = []
    details: dict[str, Any] = {}

    checks.append(_path_status(config.base_gwas, "GWAS summary statistics"))
    checks.append(_path_status(config.prsice_r, "Scoring R wrapper"))
    checks.append(_path_status(config.prsice_binary, "Scoring executable"))

    rscript_path = _resolve_executable(config.rscript)
    checks.append(
        {
            "check": "Rscript executable",
            "status": "PASS" if rscript_path else "FAIL",
            "detail": rscript_path or config.rscript,
        }
    )

    if config.pathway_input_mode == "gtf_gmt":
        checks.append(_path_status(config.gtf_file, "Gene annotation GTF"))
        checks.append(_path_status(config.gmt_file, "Pathway GMT"))
    elif config.pathway_input_mode == "snp_set":
        checks.append(_path_status(config.snp_set_file, "Pathway SNP-set file"))
    else:
        checks.append(
            {
                "check": "Pathway input mode",
                "status": "FAIL",
                "detail": f"Unsupported mode: {config.pathway_input_mode}",
            }
        )

    prefixes = target_prefixes(config)
    if not prefixes:
        checks.append(
            {"check": "Target genotype", "status": "FAIL", "detail": "not provided"}
        )
    elif config.target_type == "bed":
        missing: list[str] = []
        for prefix in prefixes:
            for suffix in (".bed", ".bim", ".fam"):
                if not Path(prefix + suffix).exists():
                    missing.append(prefix + suffix)
        checks.append(
            {
                "check": "Target PLINK BED files",
                "status": "PASS" if not missing else "FAIL",
                "detail": (
                    f"{len(prefixes)} prefix(es)"
                    if not missing
                    else f"{len(missing)} missing; first: {missing[0]}"
                ),
            }
        )
    elif config.target_type == "bgen":
        missing = [prefix + ".bgen" for prefix in prefixes if not Path(prefix + ".bgen").exists()]
        checks.append(
            {
                "check": "Target BGEN files",
                "status": "PASS" if not missing else "FAIL",
                "detail": (
                    f"{len(prefixes)} prefix(es)"
                    if not missing
                    else f"{len(missing)} missing; first: {missing[0]}"
                ),
            }
        )
    else:
        checks.append(
            {
                "check": "Target genotype type",
                "status": "FAIL",
                "detail": "The scoring engine accepts BED or BGEN; PGEN needs the adapter",
            }
        )

    base_path = Path(config.base_gwas).expanduser() if config.base_gwas else None
    if base_path and base_path.exists():
        try:
            separator = _detect_separator(base_path, config.base_separator)
            preview = _read_preview(base_path, separator)
            details["gwas_preview"] = preview
            missing_columns = [
                config.gwas_columns[key]
                for key in REQUIRED_GWAS_KEYS
                if config.gwas_columns.get(key) not in preview.columns
            ]
            checks.append(
                {
                    "check": "GWAS column mapping",
                    "status": "PASS" if not missing_columns else "FAIL",
                    "detail": (
                        "all required columns found"
                        if not missing_columns
                        else "missing: " + ", ".join(missing_columns)
                    ),
                }
            )
        except Exception as error:
            checks.append(
                {
                    "check": "GWAS preview",
                    "status": "FAIL",
                    "detail": f"{type(error).__name__}: {error}",
                }
            )

    thresholds = np.asarray(config.pvalue_thresholds, dtype=float)
    valid_thresholds = (
        len(thresholds) > 0
        and np.isfinite(thresholds).all()
        and (thresholds > 0).all()
        and (thresholds <= 1).all()
    )
    checks.append(
        {
            "check": "GWAS p-value thresholds",
            "status": "PASS" if valid_thresholds else "FAIL",
            "detail": ", ".join(f"{x:g}" for x in thresholds),
        }
    )

    checks.append(
        {
            "check": "Clumping parameters",
            "status": "PASS"
            if config.clump_kb > 0 and 0 < config.clump_r2 <= 1
            else "FAIL",
            "detail": f"{config.clump_kb} kb; r2={config.clump_r2:g}",
        }
    )

    if config.no_regression:
        checks.append(
            {
                "check": "Analysis mode",
                "status": "PASS",
                "detail": "score construction only; no phenotype regression",
            }
        )
    else:
        checks.append(_path_status(config.phenotype_file, "Phenotype file"))
        checks.append(
            {
                "check": "Phenotype column",
                "status": "PASS" if config.phenotype_column else "FAIL",
                "detail": config.phenotype_column or "not provided",
            }
        )
        checks.append(_path_status(config.covariate_file, "Covariate file", required=False))

    if not config.controlled_workspace_acknowledged:
        checks.append(
            {
                "check": "Controlled-workspace acknowledgement",
                "status": "WARN",
                "detail": "required before running person-level scoring",
            }
        )

    if config.pathway_input_mode == "gtf_gmt":
        gmt = Path(config.gmt_file).expanduser() if config.gmt_file else None
        if gmt and gmt.exists():
            try:
                pathway_summary, gene_preview = read_gmt_summary(
                    gmt,
                    config.gmt_second_column_is_description,
                )
                details["pathway_summary"] = pathway_summary
                details["pathway_gene_preview"] = gene_preview
                checks.append(
                    {
                        "check": "Pathway definitions",
                        "status": "PASS" if len(pathway_summary) else "FAIL",
                        "detail": (
                            f"{len(pathway_summary):,} pathways; "
                            f"median {pathway_summary['n_genes'].median():.0f} genes"
                            if len(pathway_summary)
                            else "no pathways found"
                        ),
                    }
                )
            except Exception as error:
                checks.append(
                    {
                        "check": "Pathway definitions",
                        "status": "FAIL",
                        "detail": f"{type(error).__name__}: {error}",
                    }
                )
        gtf = Path(config.gtf_file).expanduser() if config.gtf_file else None
        if gtf and gtf.exists():
            try:
                details["gtf_summary"] = read_gtf_summary(gtf)
            except Exception as error:
                checks.append(
                    {
                        "check": "GTF preview",
                        "status": "FAIL",
                        "detail": f"{type(error).__name__}: {error}",
                    }
                )

    return pd.DataFrame(checks), details


def build_prset_command(config: WorkflowConfig) -> list[str]:
    """Build the pathway-scoring command without executing it."""
    columns = config.gwas_columns
    command = [
        config.rscript,
        config.prsice_r,
        "--prsice",
        config.prsice_binary,
        "--base",
        config.base_gwas,
    ]

    if config.target_list:
        command.extend(["--target-list", config.target_list])
    else:
        command.extend(["--target", config.target_prefix])

    if config.target_type == "bgen":
        command.extend(["--type", "bgen"])

    if config.target_keep_file:
        command.extend(["--keep", config.target_keep_file])

    column_flags = {
        "chr": "--chr",
        "bp": "--bp",
        "snp": "--snp",
        "a1": "--A1",
        "a2": "--A2",
        "stat": "--stat",
        "p": "--pvalue",
    }
    for key, flag in column_flags.items():
        value = str(columns.get(key, "")).strip()
        if not value:
            raise ValueError(f"GWAS column '{key}' is missing.")
        if value.startswith("#"):
            raise ValueError(
                f"GWAS column '{value}' is not safe for the scoring wrapper. "
                "Run Step 2 to create the canonical harmonized GWAS file."
            )
        command.extend([flag, value])

    if config.statistic_type.upper() == "OR":
        command.append("--or")
    else:
        command.append("--beta")

    if config.pathway_input_mode == "gtf_gmt":
        command.extend(
            [
                "--gtf",
                config.gtf_file,
                "--msigdb",
                config.gmt_file,
                "--wind-5",
                str(config.window_5_bp),
                "--wind-3",
                str(config.window_3_bp),
            ]
        )
    else:
        command.extend(["--snp-set", config.snp_set_file])

    command.extend(
        [
            "--clump-kb",
            str(config.clump_kb),
            "--clump-r2",
            str(config.clump_r2),
            "--bar-levels",
            ",".join(f"{value:g}" for value in config.pvalue_thresholds),
            "--fastscore",
            "--score",
            config.score_method,
            "--thread",
            str(config.threads),
            "--print-snp",
        ]
    )
    if 1.0 not in set(float(value) for value in config.pvalue_thresholds):
        command.append("--no-full")

    if config.base_maf_column and config.base_maf_min is not None:
        command.extend(
            ["--base-maf", f"{config.base_maf_column}:{config.base_maf_min:g}"]
        )
    if config.base_info_column and config.base_info_min is not None:
        command.extend(
            ["--base-info", f"{config.base_info_column}:{config.base_info_min:g}"]
        )

    if config.no_regression:
        command.extend(["--no-regress", "--all-score"])
    else:
        command.extend(
            [
                "--pheno",
                config.phenotype_file,
                "--pheno-col",
                config.phenotype_column,
                "--binary-target",
                "T" if config.binary_target else "F",
            ]
        )
        if config.covariate_file:
            command.extend(["--cov", config.covariate_file])
        if config.covariate_columns:
            command.extend(["--cov-col", config.covariate_columns])
        if config.set_permutations > 0:
            command.extend(["--set-perm", str(config.set_permutations)])

    out_prefix = Path(config.output_dir) / config.output_prefix
    command.extend(["--out", str(out_prefix)])
    return command


def command_as_shell(command: Iterable[str]) -> str:
    """Render a command safely for display and provenance."""
    return " \\\n  ".join(shlex.quote(str(value)) for value in command)


def run_prset(
    config: WorkflowConfig,
    *,
    execute: bool = False,
    timeout_seconds: int | None = None,
) -> dict[str, Any]:
    """Validate and optionally run pathway scoring.

    ``execute=False`` is the default so opening the tutorial cannot start an
    expensive analysis accidentally.
    """
    checks, _ = validate_inputs(config)
    failures = checks.loc[checks["status"] == "FAIL"]
    if len(failures):
        raise ValueError(
            "Input validation failed:\n"
            + failures[["check", "detail"]].to_string(index=False)
        )
    command = build_prset_command(config)
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = config.to_json(output_dir / "run_config.json")
    manifest = {
        "project_name": config.project_name,
        "genome_build": config.genome_build,
        "execute": execute,
        "command": command,
        "command_shell": command_as_shell(command),
        "config_path": str(config_path),
        "started_at_unix": time.time(),
    }
    if not execute:
        manifest["status"] = "DRY_RUN"
        with (output_dir / "run_manifest.json").open("w") as handle:
            json.dump(manifest, handle, indent=2)
        return manifest

    if not config.controlled_workspace_acknowledged:
        raise PermissionError(
            "Confirm that person-level outputs will remain in the controlled workspace."
        )

    completed = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout_seconds,
    )
    manifest.update(
        {
            "status": "PASS" if completed.returncode == 0 else "FAIL",
            "return_code": completed.returncode,
            "finished_at_unix": time.time(),
            "stdout_tail": completed.stdout[-10_000:],
            "stderr_tail": completed.stderr[-10_000:],
        }
    )
    with (output_dir / "run_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2)
    if completed.returncode != 0:
        raise RuntimeError(
            "Pathway scoring failed. Review run_manifest.json and the run log.\n"
            + completed.stderr[-3000:]
        )
    return manifest


def _read_result_table(path: Path) -> pd.DataFrame:
    for separator in ("\t", None, ","):
        try:
            frame = pd.read_csv(
                path,
                sep=separator if separator is not None else r"\s+",
                engine="python",
            )
            if len(frame.columns) > 1:
                return frame
        except Exception:
            continue
    raise ValueError(f"Could not parse result table: {path}")


def discover_outputs(config: WorkflowConfig) -> pd.DataFrame:
    output_dir = Path(config.output_dir)
    prefix = config.output_prefix
    rows = []
    if output_dir.exists():
        for path in sorted(output_dir.glob(prefix + "*")):
            rows.append(
                {
                    "file": path.name,
                    "size_mb": path.stat().st_size / (1024**2),
                    "person_level": (
                        path.name.endswith(".best")
                        or path.name.endswith(".all.score")
                        or path.name.endswith(".all_score")
                    ),
                }
            )
    return pd.DataFrame(rows)


def load_aggregate_results(config: WorkflowConfig) -> dict[str, pd.DataFrame]:
    """Load only aggregate pathway-scoring result tables for visualization."""
    output_dir = Path(config.output_dir)
    prefix = config.output_prefix
    result: dict[str, pd.DataFrame] = {}
    candidates = {
        "summary": output_dir / f"{prefix}.summary",
        "model_fit": output_dir / f"{prefix}.prsice",
    }
    for label, path in candidates.items():
        if path.exists():
            result[label] = _read_result_table(path)
    return result


def summarize_aggregate_results(config: WorkflowConfig) -> pd.DataFrame:
    tables = load_aggregate_results(config)
    if "summary" not in tables:
        output_dir = Path(config.output_dir)
        candidates = [
            output_dir / f"{config.output_prefix}.all.score",
            output_dir / f"{config.output_prefix}.all_score",
            output_dir / f"{config.output_prefix}.best",
        ]
        score_file = next((path for path in candidates if path.exists()), None)
        if score_file is None:
            return pd.DataFrame()

        separator = _detect_separator(score_file, "auto")
        read_separator = r"\s+" if separator is None else separator
        preview = pd.read_csv(score_file, sep=read_separator, nrows=5)
        excluded = {
            "FID", "IID", "ID", "PHENO", "PHENO1", "CNT", "CNT2",
            "#FID", "#IID",
        }
        score_columns = [column for column in preview.columns if str(column) not in excluded]
        if not score_columns:
            return pd.DataFrame()

        count = pd.Series(0, index=score_columns, dtype="int64")
        total = pd.Series(0.0, index=score_columns)
        total_squared = pd.Series(0.0, index=score_columns)
        minimum = pd.Series(np.inf, index=score_columns)
        maximum = pd.Series(-np.inf, index=score_columns)
        for chunk in pd.read_csv(
            score_file,
            sep=read_separator,
            usecols=score_columns,
            chunksize=2_000,
        ):
            numeric = chunk.apply(pd.to_numeric, errors="coerce")
            count += numeric.count()
            total += numeric.sum(skipna=True)
            total_squared += numeric.pow(2).sum(skipna=True)
            minimum = pd.concat([minimum, numeric.min(skipna=True)], axis=1).min(axis=1)
            maximum = pd.concat([maximum, numeric.max(skipna=True)], axis=1).max(axis=1)
        denominator = count.clip(lower=1)
        mean = total / denominator
        variance = (total_squared / denominator - mean.pow(2)).clip(lower=0)
        summary = pd.DataFrame(
            {
                "score_name": score_columns,
                "n": count.values,
                "score_mean": mean.values,
                "score_sd": np.sqrt(variance.values),
                "score_min": minimum.values,
                "score_max": maximum.values,
            }
        )
        summary["all_finite"] = np.isfinite(
            summary[["score_mean", "score_sd", "score_min", "score_max"]]
        ).all(axis=1)
        aggregate_path = output_dir / f"{config.output_prefix}.aggregate_score_summary.tsv"
        summary.to_csv(aggregate_path, sep="\t", index=False)
        return summary
    summary = tables["summary"].copy()
    rename = {
        "Set": "pathway",
        "Threshold": "threshold",
        "PRS.R2": "prs_r2",
        "P": "p_value",
        "Num_SNP": "n_snps",
        "Coefficient": "coefficient",
    }
    summary = summary.rename(columns={key: value for key, value in rename.items() if key in summary})
    if "p_value" in summary:
        summary["minus_log10_p"] = -np.log10(
            pd.to_numeric(summary["p_value"], errors="coerce").clip(lower=1e-300)
        )
    return summary


def plot_pathway_definition_qc(pathway_summary: pd.DataFrame):
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].hist(pathway_summary["n_genes"], bins=35, color="#287271", edgecolor="white")
    axes[0].set_xlabel("Genes per pathway")
    axes[0].set_ylabel("Number of pathways")
    axes[0].set_title("Pathway-size distribution")

    top = pathway_summary.nlargest(15, "n_genes").sort_values("n_genes")
    axes[1].barh(top["pathway"], top["n_genes"], color="#E07A5F")
    axes[1].set_xlabel("Genes")
    axes[1].set_title("Largest pathways")
    axes[1].tick_params(axis="y", labelsize=8)
    figure.tight_layout()
    return figure


def plot_prset_results(summary: pd.DataFrame):
    import matplotlib.pyplot as plt

    if summary.empty:
        raise ValueError("No aggregate pathway-scoring summary is available.")
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    if "score_sd" in summary:
        values = pd.to_numeric(summary["score_sd"], errors="coerce").dropna()
        axes[0].hist(values, bins=min(30, max(5, len(values))), color="#277DA1", edgecolor="white")
        axes[0].set_xlabel("Score standard deviation")
        axes[0].set_ylabel("Number of pathway scores")
        axes[0].set_title("Pathway-score variability")
    elif "n_snps" in summary:
        values = pd.to_numeric(summary["n_snps"], errors="coerce").dropna()
        axes[0].hist(values, bins=35, color="#3D7EA6", edgecolor="white")
        axes[0].set_xlabel("SNPs retained per pathway")
        axes[0].set_ylabel("Number of pathways")
        axes[0].set_title("Post-clumping pathway sizes")
    else:
        axes[0].text(0.5, 0.5, "Num_SNP not found", ha="center", va="center")
        axes[0].set_axis_off()

    if {"score_name", "score_sd"}.issubset(summary.columns):
        top = summary.nlargest(min(15, len(summary)), "score_sd").sort_values("score_sd")
        axes[1].barh(top["score_name"], top["score_sd"], color="#43AA8B")
        axes[1].set_xlabel("Score standard deviation")
        axes[1].set_title("Most variable pathway scores")
        axes[1].tick_params(axis="y", labelsize=8)
    elif {"pathway", "minus_log10_p"}.issubset(summary.columns):
        top = summary.nlargest(20, "minus_log10_p").sort_values("minus_log10_p")
        axes[1].barh(top["pathway"], top["minus_log10_p"], color="#8F5D90")
        axes[1].set_xlabel(r"$-\log_{10}(p)$")
        axes[1].set_title("Top pathway associations")
        axes[1].tick_params(axis="y", labelsize=8)
    elif {"pathway", "prs_r2"}.issubset(summary.columns):
        top = summary.nlargest(20, "prs_r2").sort_values("prs_r2")
        axes[1].barh(top["pathway"], top["prs_r2"], color="#8F5D90")
        axes[1].set_xlabel(r"Pathway PRS $R^2$")
        axes[1].set_title("Top pathway models")
        axes[1].tick_params(axis="y", labelsize=8)
    else:
        axes[1].text(0.5, 0.5, "Association fields not found", ha="center", va="center")
        axes[1].set_axis_off()

    figure.tight_layout()
    return figure


def generate_pgen_adapter_script(
    *,
    pgen_pattern: str,
    chromosomes: Iterable[int],
    extract_file: str,
    keep_file: str,
    output_dir: str,
    extract_format: str = "id",
    set_coordinate_allele_ids: bool = False,
    plink2: str = "plink2",
    threads: int = 4,
    memory_mb: int = 16_000,
) -> str:
    """Generate, but do not execute, an AoU PGEN-to-BED adapter script.

    The extract file must contain the de-duplicated union of target variants.
    Converting the entire AoU WGS dataset is intentionally not supported.
    """
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# The scoring engine does not read PGEN directly. Convert only the",
        "# union of pathway variants and only the selected participant cohort.",
        f"mkdir -p {shlex.quote(output_dir)}",
        f": > {shlex.quote(str(Path(output_dir) / 'target_prefixes.txt'))}",
        "",
    ]
    for chrom in chromosomes:
        prefix = pgen_pattern.replace("{chr}", str(chrom)).replace("#", str(chrom))
        out_prefix = str(Path(output_dir) / f"target_chr{chrom}")
        extract_arguments = (
            ["--extract", "bed1", extract_file]
            if extract_format == "bed1"
            else ["--extract", extract_file]
        )
        id_arguments = (
            [
                "--set-all-var-ids",
                "@:#:$r:$a",
                "--new-id-max-allele-len",
                "100",
                "missing",
                "--snps-only",
                "just-acgt",
                "--max-alleles",
                "2",
                "--rm-dup",
                "force-first",
            ]
            if set_coordinate_allele_ids
            else []
        )
        command = [
            plink2,
            "--pfile",
            prefix,
            *id_arguments,
            *extract_arguments,
            "--keep",
            keep_file,
            "--make-bed",
            "--threads",
            str(threads),
            "--memory",
            str(memory_mb),
            "--out",
            out_prefix,
        ]
        lines.append(" ".join(shlex.quote(value) for value in command))
        lines.append(
            f"printf '%s\\n' {shlex.quote(out_prefix)} >> "
            f"{shlex.quote(str(Path(output_dir) / 'target_prefixes.txt'))}"
        )
        lines.append("")
    return "\n".join(lines)


def write_markdown_report(
    config: WorkflowConfig,
    summary: pd.DataFrame | None = None,
) -> Path:
    """Create a concise aggregate-only report after a dry run or completed run."""
    checks, details = validate_inputs(config)
    outputs = discover_outputs(config)
    if summary is None:
        summary = summarize_aggregate_results(config)
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / "PATHWAY_PRS_REPORT.md"

    lines = [
        f"# Pathway PRS report: {config.project_name}",
        "",
        "## Definition",
        "",
        f"- Genome build: `{config.genome_build}`",
        f"- Pathway input: `{config.pathway_input_mode}`",
        f"- P-value threshold(s): `{config.pvalue_thresholds}`",
        f"- LD clumping: `{config.clump_kb} kb`, `r2={config.clump_r2}`",
        f"- Score method: `{config.score_method}`",
        "",
        "## Input validation",
        "",
        checks.to_markdown(index=False),
        "",
    ]
    if "pathway_summary" in details:
        pathway_summary = details["pathway_summary"]
        lines.extend(
            [
                "## Pathway definitions",
                "",
                f"- Pathways: {len(pathway_summary):,}",
                f"- Unique pathway-gene memberships: {int(pathway_summary['n_genes'].sum()):,}",
                f"- Median genes per pathway: {pathway_summary['n_genes'].median():.0f}",
                "",
            ]
        )
    if len(summary):
        lines.extend(
            [
                "## Aggregate pathway PGS results",
                "",
                summary.head(30).to_markdown(index=False),
                "",
            ]
        )
    if len(outputs):
        public_outputs = outputs.loc[~outputs["person_level"]]
        lines.extend(
            [
                "## Aggregate output files",
                "",
                public_outputs.to_markdown(index=False),
                "",
            ]
        )
    lines.extend(
        [
            "## Privacy",
            "",
            "- Individual-level PRS files must remain inside the controlled workspace.",
            "- Only aggregate QC tables and approved figures should be exported.",
            "",
        ]
    )
    destination.write_text("\n".join(lines))
    return destination


# Public notebook-facing names. The legacy function names remain available so
# existing configurations and tests continue to work.
build_pathway_pgs_command = build_prset_command
run_pathway_pgs = run_prset
plot_pathway_pgs_results = plot_prset_results
