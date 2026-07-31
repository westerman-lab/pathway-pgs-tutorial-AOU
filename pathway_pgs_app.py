#!/usr/bin/env python3
"""Compact ipywidgets interface for the Pathway PGS tutorial."""

from __future__ import annotations

import os
import shutil
import urllib.request
import zipfile
from pathlib import Path

import ipywidgets as widgets
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import HTML, Markdown, clear_output, display

from pathway_prs_core import (
    WorkflowConfig,
    build_pathway_pgs_command,
    build_snp_set_from_variant_gene_mapping,
    build_pathway_variant_union_bed,
    build_variant_mapping_union_bed,
    command_as_shell,
    harmonize_gwas_ids_to_bed,
    infer_gwas_schema,
    plot_pathway_definition_qc,
    plot_pathway_pgs_results,
    run_pathway_pgs,
    summarize_aggregate_results,
    validate_inputs,
    write_markdown_report,
)


HOME = Path.home()
TUTORIAL_DIR = Path(__file__).resolve().parent
DEMO_DIR = TUTORIAL_DIR / "demo_data"
ANALYSIS = HOME / "analysis"
DATA = ANALYSIS / "data"
RESULTS = ANALYSIS / "results" / "pathway_prs_tutorial"
ADAPTER = DATA / "pathway_prs_adapter"
WORK = ADAPTER
LABEL = {"description_width": "155px"}
STEP_BUTTON = widgets.Layout(
    width="auto",
    min_width="145px",
    height="40px",
    flex="1 1 155px",
    margin="0 8px 8px 0",
)
CARD = widgets.Layout(
    width="100%",
    border_top="1px solid #d8e0e5",
    border_right="1px solid #d8e0e5",
    border_bottom="1px solid #d8e0e5",
    border_left="1px solid #d8e0e5",
    padding="14px 16px",
    margin="0 0 12px 0",
)


def _full_layout() -> widgets.Layout:
    """Return an independent full-width layout for one widget."""
    return widgets.Layout(width="100%")


def _text(description: str, value: str = "", placeholder: str = ""):
    return widgets.Text(
        description=description,
        value=value,
        placeholder=placeholder,
        style=LABEL,
        layout=_full_layout(),
    )


def _workflow_status(step: int, label: str) -> str:
    return f"""
    <div style="
        background:#f3f7f8;
        border-left:4px solid #177e89;
        padding:8px 12px;
        margin:4px 0 8px 0;
        color:#24323d;
    ">
      <b>Step {step} of 5</b>
      <span style="margin-left:8px">{label}</span>
    </div>
    """


def _walk_candidates(root: Path, predicate, maximum_depth: int = 6) -> list[Path]:
    if not root.exists():
        return []
    found: list[Path] = []
    root_depth = len(root.parts)
    for directory, subdirectories, filenames in os.walk(root):
        current = Path(directory)
        if len(current.parts) - root_depth >= maximum_depth:
            subdirectories[:] = []
        subdirectories[:] = [
            name
            for name in subdirectories
            if name not in {"results", "__pycache__", ".git", ".ipynb_checkpoints"}
        ]
        for filename in filenames:
            path = current / filename
            if predicate(path):
                found.append(path)
    return found


def _ranked(candidates: list[Path], preferred: tuple[str, ...]) -> str:
    if not candidates:
        return ""

    def rank(path: Path):
        name = str(path).lower()
        score = sum(10 ** (len(preferred) - index) for index, token in enumerate(preferred) if token in name)
        return (-score, len(str(path)), str(path))

    return str(sorted(candidates, key=rank)[0])


def discover_local_inputs() -> dict[str, str]:
    """Find likely local tools and reference files without reading person-level data."""
    known_gwas = DATA / "t2d_gwas_sumstats" / "t2d_gwas_sumstats.tsv.gz"
    gwas_candidates = [known_gwas] if known_gwas.exists() else _walk_candidates(
        DATA,
        lambda path: any(token in path.name.lower() for token in ("gwas", "sumstat"))
        and path.suffix.lower() in {".gz", ".tsv", ".txt", ".csv"},
    )
    gtf_candidates = _walk_candidates(
        DATA,
        lambda path: path.name.lower().endswith((".gtf", ".gtf.gz")),
    )
    gmt_candidates = _walk_candidates(
        DATA,
        lambda path: path.name.lower().endswith(".gmt"),
    )
    wrapper_candidates = [
        ANALYSIS / "tools" / "PRSice" / "PRSice.R",
        HOME / "PRSice" / "PRSice.R",
    ]
    executable_candidates = [
        ANALYSIS / "tools" / "PRSice" / "PRSice_linux",
        HOME / "PRSice" / "PRSice_linux",
    ]
    return {
        "gwas": _ranked(gwas_candidates, ("t2d", "gwas", "sumstat")),
        "gtf": _ranked(gtf_candidates, ("grch38", "gencode", "homo_sapiens")),
        "gmt": _ranked(gmt_candidates, ("reactome", "pathway")),
        "wrapper": next((str(path) for path in wrapper_candidates if path.exists()), ""),
        "executable": next(
            (str(path) for path in executable_candidates if path.exists()), ""
        ),
        "rscript": shutil.which("Rscript") or "Rscript",
        "plink2": shutil.which("plink2")
        or "/opt/workbench-tools/binaries/bin/plink2",
    }


def launch_pathway_pgs_app():
    """Return the complete compact Pathway PGS application."""
    discovered = discover_local_inputs()

    title = widgets.HTML(
        """
        <div style="border-left:5px solid #177e89;padding:8px 14px;margin:2px 0 12px 0">
          <div style="font-size:22px;font-weight:700">Pathway PGS workspace</div>
          <div style="color:#177e89;font-weight:600;margin-top:2px">
            Made by Sijia Zhu, Westerman's Lab
          </div>
          <div style="color:#4b5563;margin-top:3px">
            Run a short demonstration first, then switch to a complete or
            custom pathway analysis.
          </div>
        </div>
        """
    )

    run_mode = widgets.ToggleButtons(
        description="Run mode",
        options=[("Quick demo", "demo"), ("Full analysis", "full")],
        value="demo",
        style=LABEL,
    )
    project_name = _text("Run name", "pathway_pgs_quick_demo")
    genome_build = widgets.Dropdown(
        description="Genome build",
        options=["GRCh38", "GRCh37"],
        value="GRCh38",
        style=LABEL,
        layout=_full_layout(),
    )
    mapping_method = widgets.Dropdown(
        description="SNP-to-gene method",
        options=[
            ("Physical distance", "physical"),
            ("Custom mapping table", "custom"),
        ],
        value="physical",
        style=LABEL,
        layout=_full_layout(),
    )
    annotation_source = widgets.Dropdown(
        description="Gene annotation",
        options=[
            ("GRCh38 / Ensembl", "grch38"),
            ("GRCh37 / Ensembl", "grch37"),
            ("Custom GTF", "custom"),
        ],
        value="grch38",
        style=LABEL,
        layout=_full_layout(),
    )
    pathway_source = widgets.Dropdown(
        description="Pathway database",
        options=[
            ("Reactome", "reactome"),
            ("GO Biological Process", "go_bp"),
            ("WikiPathways", "wikipathways"),
            ("Custom GMT", "custom"),
        ],
        value="reactome",
        style=LABEL,
        layout=_full_layout(),
    )
    base_gwas = _text("GWAS file", str(DEMO_DIR / "demo_gwas_chr21.tsv"))
    gtf_file = _text(
        "GTF file",
        str(DEMO_DIR / "demo_chr21.gtf"),
        "GRCh38 .gtf or .gtf.gz",
    )
    gmt_file = _text(
        "GMT file",
        str(DEMO_DIR / "demo_pathways.gmt"),
        "Reactome .gmt",
    )
    variant_gene_file = _text(
        "SNP-to-gene file",
        "",
        "Table with variant + gene, or chr + position + gene",
    )
    keep_file = _text(
        "Participant keep file",
        str(ADAPTER / "quick_demo_100.keep"),
        "Controlled-workspace participant IDs",
    )
    demo_summary = widgets.HTML(
        """
        <div style="background:#eef7f6;border:1px solid #b9d9d5;padding:10px 12px">
          <b>Quick demo</b><br>
          100 AoU participants &nbsp;|&nbsp; chromosome 21 &nbsp;|&nbsp;
          2 example variants &nbsp;|&nbsp; 2 example pathways<br>
          <span style="color:#53636c">Technical validation only; not a scientific result.</span>
        </div>
        """
    )

    engine_status = widgets.HTML()
    install_engine = widgets.Button(
        description="Install scoring engine",
        icon="download",
        layout=widgets.Layout(width="190px", height="36px"),
    )
    setup_message = widgets.Output()

    auto_find = widgets.Button(
        description="Find local inputs", button_style="info", icon="search"
    )
    detect_gwas = widgets.Button(description="Read GWAS header", icon="table")
    input_message = widgets.Output(
        layout=widgets.Layout(
            border_top="1px solid #d7dde5",
            border_right="1px solid #d7dde5",
            border_bottom="1px solid #d7dde5",
            border_left="1px solid #d7dde5",
            padding="8px",
            margin="6px 0",
            max_height="130px",
            overflow="auto",
        )
    )

    upload_help = widgets.HTML(
        "<div style='color:#53636c;padding:4px 0'>"
        "To use custom files, upload them with the JupyterLab file browser, "
        "then enter or paste their paths above. The interface can detect common "
        "local files automatically."
        "</div>"
    )

    base_separator = widgets.Dropdown(
        description="GWAS separator",
        options=["auto", "tab", "whitespace", "comma"],
        value="auto",
        style=LABEL,
        layout=_full_layout(),
    )
    statistic_type = widgets.ToggleButtons(
        description="Effect statistic",
        options=["BETA", "OR"],
        value="BETA",
        style=LABEL,
    )
    column_widgets = {
        key: _text(f"GWAS {key}", default)
        for key, default in {
            "chr": "CHR",
            "bp": "BP",
            "snp": "SNP",
            "a1": "A1",
            "a2": "A2",
            "stat": "BETA",
            "p": "P",
        }.items()
    }
    base_maf_column = _text("MAF column", "")
    base_info_column = _text("INFO column", "")

    prsice_r = _text("Scoring R wrapper", discovered["wrapper"])
    prsice_binary = _text("Scoring executable", discovered["executable"])
    rscript = _text("Rscript", discovered["rscript"])
    output_dir = _text("Output folder", str(RESULTS))
    output_prefix = _text("Output prefix", project_name.value)
    bed_pattern = _text(
        "AoU BED pattern",
        (
            "/home/jupyter/workspace/vwb-aou-datasets-controlled-v9/v9/"
            "wgs/short_read/snpindel/acaf_threshold/plink_bed/"
            "acaf_threshold.chr#"
        ),
    )
    plink2_path = _text("PLINK2", discovered["plink2"])

    gmt_description = widgets.Checkbox(
        description="GMT field 2 is a description or URL",
        value=True,
        indent=False,
        layout=_full_layout(),
    )
    window_5 = widgets.IntText(
        description="5-prime window (bp)",
        value=35_000,
        style=LABEL,
        layout=_full_layout(),
    )
    window_3 = widgets.IntText(
        description="3-prime window (bp)",
        value=10_000,
        style=LABEL,
        layout=_full_layout(),
    )
    thresholds = _text("P-value threshold(s)", "1")
    clump_kb = widgets.IntText(
        description="Clump window (kb)",
        value=1_000,
        style=LABEL,
        layout=_full_layout(),
    )
    clump_r2 = widgets.FloatText(
        description="Clump r-squared",
        value=0.1,
        style=LABEL,
        layout=_full_layout(),
    )
    score_method = widgets.Dropdown(
        description="Score method",
        options=["sum", "avg", "std", "con-std"],
        value="sum",
        style=LABEL,
        layout=_full_layout(),
    )
    threads = widgets.IntSlider(
        description="Threads",
        min=1,
        max=32,
        value=4,
        style=LABEL,
        continuous_update=False,
        layout=_full_layout(),
    )
    chromosomes = _text(
        "Chromosomes", "21"
    )

    phenotype_file = _text("Phenotype file", "")
    phenotype_column = _text("Phenotype column", "")
    binary_target = widgets.Checkbox(
        description="Binary phenotype", value=False, indent=False, layout=_full_layout()
    )
    covariate_file = _text("Covariate file", "")
    covariate_columns = _text("Covariate columns", "")
    run_association = widgets.Checkbox(
        description="Run phenotype association after scoring",
        value=False,
        indent=False,
        layout=_full_layout(),
    )

    advanced = widgets.Accordion(
        children=[
            widgets.VBox(
                [
                    genome_build,
                    base_separator,
                    statistic_type,
                    *column_widgets.values(),
                    base_maf_column,
                    base_info_column,
                ]
            ),
            widgets.VBox(
                [
                    window_5,
                    window_3,
                    gmt_description,
                    thresholds,
                    clump_kb,
                    clump_r2,
                    score_method,
                    threads,
                    chromosomes,
                ]
            ),
            widgets.VBox(
                [
                    prsice_r,
                    prsice_binary,
                    rscript,
                    plink2_path,
                    bed_pattern,
                    output_dir,
                    output_prefix,
                ]
            ),
            widgets.VBox(
                [
                    run_association,
                    phenotype_file,
                    phenotype_column,
                    binary_target,
                    covariate_file,
                    covariate_columns,
                ]
            ),
        ],
        selected_index=None,
        layout=_full_layout(),
    )
    for index, label in enumerate(
        ["GWAS mapping", "Scoring parameters", "Tools and output", "Optional association"]
    ):
        advanced.set_title(index, label)

    workspace_ack = widgets.Checkbox(
        description="Keep person-level files inside the controlled workspace",
        value=False,
        indent=False,
        layout=_full_layout(),
    )
    execute_ack = widgets.Checkbox(
        description="I reviewed the inputs and approve execution",
        value=False,
        indent=False,
        layout=_full_layout(),
    )

    check_button = widgets.Button(
        description="1. Check", button_style="info", icon="check", layout=STEP_BUTTON
    )
    prepare_button = widgets.Button(
        description="2. Prepare data",
        icon="cogs",
        disabled=True,
        layout=STEP_BUTTON,
    )
    preview_button = widgets.Button(
        description="3. Review",
        icon="search",
        disabled=True,
        layout=STEP_BUTTON,
    )
    run_button = widgets.Button(
        description="4. Calculate",
        button_style="danger",
        icon="play",
        disabled=True,
        layout=STEP_BUTTON,
    )
    results_button = widgets.Button(
        description="5. Results",
        button_style="success",
        icon="bar-chart",
        disabled=True,
        layout=STEP_BUTTON,
    )
    workflow_status = widgets.HTML(
        value=_workflow_status(1, "Check files and detected columns")
    )
    action_output = widgets.Output(
        layout=widgets.Layout(
            border_top="1px solid #d7dde5",
            border_right="1px solid #d7dde5",
            border_bottom="1px solid #d7dde5",
            border_left="1px solid #d7dde5",
            padding="10px",
            margin="10px 0 0 0",
            max_height="430px",
            overflow="auto",
        )
    )

    harmonized_gwas = ADAPTER / "gwas_for_aou_bed.tsv.gz"
    union_bed = ADAPTER / "pathway_variant_union.bed1.tsv"
    custom_snp_set = ADAPTER / "custom_pathway_snp_sets.gmt"
    target_list_file = ADAPTER / "aou_bed_target_prefixes.txt"
    source_gwas = {"path": ""}

    def refresh_engine_status() -> None:
        ready = Path(prsice_r.value).exists() and Path(prsice_binary.value).exists()
        color = "#eaf6ef" if ready else "#fff7e6"
        border = "#9bc9ac" if ready else "#e5c172"
        label = "Scoring engine ready" if ready else "Scoring engine not found"
        detail = (
            "No setup is needed."
            if ready
            else "Install it once, or select an existing copy under Advanced settings."
        )
        engine_status.value = (
            f"<div style='background:{color};border:1px solid {border};"
            "padding:8px 12px'>"
            f"<b>{label}</b><br><span style='color:#53636c'>{detail}</span></div>"
        )
        install_engine.disabled = ready
        if ready:
            install_engine.description = "Scoring engine ready"
        else:
            install_engine.description = "Install scoring engine"

    def handle_install_engine(_):
        install_engine.disabled = True
        with setup_message:
            clear_output()
            try:
                destination = ANALYSIS / "tools" / "PRSice"
                destination.mkdir(parents=True, exist_ok=True)
                archive = destination / "PRSice_linux_2.3.5.zip"
                url = (
                    "https://github.com/choishingwan/PRSice/releases/download/"
                    "2.3.5/PRSice_linux.zip"
                )
                print("Downloading the Linux scoring engine from its official release...")
                urllib.request.urlretrieve(url, archive)
                with zipfile.ZipFile(archive) as handle:
                    handle.extractall(destination)
                wrapper = next(destination.rglob("PRSice.R"), None)
                executable = next(destination.rglob("PRSice_linux"), None)
                if wrapper is None or executable is None:
                    raise FileNotFoundError(
                        "The downloaded release did not contain the expected files."
                    )
                executable.chmod(executable.stat().st_mode | 0o111)
                prsice_r.value = str(wrapper)
                prsice_binary.value = str(executable)
                refresh_engine_status()
                print("Setup complete. Continue with Step 1.")
            except Exception as error:
                install_engine.disabled = False
                print(f"{type(error).__name__}: {error}")

    def selected_chromosomes() -> list[int]:
        return [
            int(value.strip())
            for value in chromosomes.value.split(",")
            if value.strip()
        ]

    def apply_schema():
        schema = infer_gwas_schema(base_gwas.value.strip())
        base_separator.value = schema["separator"]
        statistic_type.value = schema["statistic_type"]
        for key, value in schema["gwas_columns"].items():
            column_widgets[key].value = value
        base_maf_column.value = schema["maf_column"]
        base_info_column.value = schema["info_column"]
        return schema

    def current_config() -> WorkflowConfig:
        p_values = [
            float(value.strip())
            for value in thresholds.value.split(",")
            if value.strip()
        ]
        custom_mapping = mapping_method.value == "custom"
        return WorkflowConfig(
            project_name=project_name.value.strip(),
            genome_build=genome_build.value,
            base_gwas=base_gwas.value.strip(),
            base_separator=base_separator.value,
            gwas_columns={
                key: widget.value.strip() for key, widget in column_widgets.items()
            },
            statistic_type=statistic_type.value,
            target_type="bed",
            target_prefix="",
            target_list=str(target_list_file),
            target_keep_file=keep_file.value.strip(),
            pathway_input_mode="snp_set" if custom_mapping else "gtf_gmt",
            gtf_file=gtf_file.value.strip(),
            gmt_file=gmt_file.value.strip(),
            gmt_second_column_is_description=gmt_description.value,
            window_5_bp=int(window_5.value),
            window_3_bp=int(window_3.value),
            snp_set_file=str(custom_snp_set) if custom_mapping else "",
            pvalue_thresholds=p_values,
            clump_kb=int(clump_kb.value),
            clump_r2=float(clump_r2.value),
            base_maf_column=base_maf_column.value.strip(),
            base_maf_min=0.01 if base_maf_column.value.strip() else None,
            base_info_column=base_info_column.value.strip(),
            base_info_min=0.8 if base_info_column.value.strip() else None,
            phenotype_file=phenotype_file.value.strip(),
            phenotype_column=phenotype_column.value.strip(),
            binary_target=binary_target.value,
            covariate_file=covariate_file.value.strip(),
            covariate_columns=covariate_columns.value.strip(),
            no_regression=not run_association.value,
            score_method=score_method.value,
            threads=int(threads.value),
            prsice_r=prsice_r.value.strip(),
            prsice_binary=prsice_binary.value.strip(),
            rscript=rscript.value.strip(),
            output_dir=output_dir.value.strip(),
            output_prefix=output_prefix.value.strip(),
            controlled_workspace_acknowledged=workspace_ack.value,
        )

    def refresh_local_inputs(_=None):
        found = discover_local_inputs()
        if found["gwas"]:
            base_gwas.value = found["gwas"]
        if found["gtf"]:
            gtf_file.value = found["gtf"]
        if found["gmt"]:
            gmt_file.value = found["gmt"]
        if found["wrapper"]:
            prsice_r.value = found["wrapper"]
        if found["executable"]:
            prsice_binary.value = found["executable"]
        rscript.value = found["rscript"]
        plink2_path.value = found["plink2"]
        refresh_engine_status()
        with input_message:
            clear_output()
            print("Local input search complete.")
            for label, widget in (
                ("GWAS", base_gwas),
                ("GTF", gtf_file),
                ("GMT", gmt_file),
            ):
                name = Path(widget.value).name if widget.value else "not found"
                print(f"{label}: {name}")

    def select_reference_files(_=None):
        if run_mode.value == "demo":
            return
        found = discover_local_inputs()
        if annotation_source.value == "grch38":
            genome_build.value = "GRCh38"
        elif annotation_source.value == "grch37":
            genome_build.value = "GRCh37"
        if annotation_source.value != "custom":
            candidates = _walk_candidates(
                DATA,
                lambda path: path.name.lower().endswith((".gtf", ".gtf.gz"))
                and genome_build.value.lower() in str(path).lower(),
            )
            gtf_file.value = _ranked(candidates, (genome_build.value.lower(), "homo_sapiens"))
        if pathway_source.value != "custom":
            tokens = {
                "reactome": ("reactome",),
                "go_bp": ("go", "biological", "process"),
                "wikipathways": ("wikipathway",),
            }[pathway_source.value]
            candidates = _walk_candidates(
                DATA,
                lambda path: path.name.lower().endswith(".gmt")
                and any(token in path.name.lower() for token in tokens),
            )
            gmt_file.value = _ranked(candidates, tokens)
        if found["gwas"] and not base_gwas.value:
            base_gwas.value = found["gwas"]
        gtf_file.disabled = annotation_source.value != "custom"
        gmt_file.disabled = pathway_source.value != "custom"

    def refresh_mapping_controls(_=None):
        if run_mode.value == "demo":
            return
        physical = mapping_method.value == "physical"
        annotation_source.layout.display = "" if physical else "none"
        gtf_file.layout.display = "" if physical else "none"
        variant_gene_file.layout.display = "none" if physical else ""

    def apply_run_mode(change=None):
        demo = run_mode.value == "demo"
        workspace_ack.value = False
        execute_ack.value = False
        prepare_button.disabled = True
        preview_button.disabled = True
        run_button.disabled = True
        results_button.disabled = True
        workflow_status.value = _workflow_status(1, "Check files and detected columns")
        action_output.clear_output()
        annotation_source.disabled = demo
        pathway_source.disabled = demo
        hidden_in_demo = [
            mapping_method,
            annotation_source,
            pathway_source,
            base_gwas,
            gtf_file,
            gmt_file,
            variant_gene_file,
            keep_file,
            auto_find,
            detect_gwas,
            upload_help,
        ]
        for widget in hidden_in_demo:
            widget.layout.display = "none" if demo else ""
        demo_summary.layout.display = "" if demo else "none"
        if demo:
            project_name.value = "pathway_pgs_quick_demo"
            base_gwas.value = str(DEMO_DIR / "demo_gwas_chr21.tsv")
            gtf_file.value = str(DEMO_DIR / "demo_chr21.gtf")
            gmt_file.value = str(DEMO_DIR / "demo_pathways.gmt")
            keep_file.value = str(ADAPTER / "quick_demo_100.keep")
            chromosomes.value = "21"
            thresholds.value = "1"
            mapping_method.value = "physical"
            input_message.clear_output()
            with input_message:
                print("Quick demo: 100 AoU participants, chromosome 21, two example pathways.")
                print("Expected runtime after setup: usually a few minutes.")
                print("Demo outputs verify the workflow and are not scientific results.")
        else:
            project_name.value = "t2d_reactome_pathway_pgs"
            chromosomes.value = ",".join(str(number) for number in range(1, 23))
            keep_file.value = str(ADAPTER / "participants.keep")
            base_gwas.value = discovered["gwas"]
            select_reference_files()
            refresh_mapping_controls()
            input_message.clear_output()
            with input_message:
                print("Full analysis: review the detected files or enter custom paths.")
                print("Complete Step 1 before preparing the AoU genotype target.")

    def ensure_demo_keep_file() -> None:
        if run_mode.value != "demo" or Path(keep_file.value).exists():
            return
        prefix = bed_pattern.value.replace("#", "21").replace("{chr}", "21")
        fam = Path(prefix + ".fam")
        if not fam.exists():
            raise FileNotFoundError(f"AoU chromosome 21 FAM not found: {fam}")
        ids = pd.read_csv(fam, sep=r"\s+", header=None, usecols=[0, 1], nrows=100)
        destination = Path(keep_file.value)
        destination.parent.mkdir(parents=True, exist_ok=True)
        ids.to_csv(destination, sep="\t", header=False, index=False)

    def handle_detect(_):
        with input_message:
            clear_output()
            try:
                ensure_demo_keep_file()
                schema = apply_schema()
                if schema["missing_required_fields"]:
                    print(
                        "Could not identify: "
                        + ", ".join(schema["missing_required_fields"])
                    )
                else:
                    print("GWAS header detected. No manual column entry is needed.")
                print(", ".join(schema["available_columns"]))
            except Exception as error:
                print(f"{type(error).__name__}: {error}")

    def handle_check(_):
        workflow_status.value = _workflow_status(1, "Checking inputs")
        with action_output:
            clear_output()
            try:
                ensure_demo_keep_file()
                schema = apply_schema()
                rows = []
                required = [
                    ("GWAS", base_gwas.value),
                    ("GMT", gmt_file.value),
                    ("Participant keep file", keep_file.value),
                    ("Scoring wrapper", prsice_r.value),
                    ("Scoring executable", prsice_binary.value),
                    ("Rscript", rscript.value),
                    ("PLINK2", plink2_path.value),
                ]
                if mapping_method.value == "physical":
                    required.insert(1, ("GTF", gtf_file.value))
                else:
                    required.insert(1, ("SNP-to-gene table", variant_gene_file.value))
                for label, value in required:
                    resolved = shutil.which(value) if label in {"Rscript", "PLINK2"} else None
                    exists = bool(value and (Path(value).expanduser().exists() or resolved))
                    rows.append(
                        {
                            "input": label,
                            "status": "PASS" if exists else "MISSING",
                            "location": (
                                Path(value).name
                                if value and "/" in value
                                else value or "not found"
                            ),
                        }
                    )
                rows.append(
                    {
                        "input": "GWAS column detection",
                        "status": (
                            "PASS"
                            if not schema["missing_required_fields"]
                            else "REVIEW"
                        ),
                        "location": (
                            "all required columns identified"
                            if not schema["missing_required_fields"]
                            else ", ".join(schema["missing_required_fields"])
                        ),
                    }
                )
                prefixes = [
                    bed_pattern.value.replace("#", str(chrom)).replace(
                        "{chr}", str(chrom)
                    )
                    for chrom in selected_chromosomes()
                ]
                missing_target_files = [
                    prefix + suffix
                    for prefix in prefixes
                    for suffix in (".bed", ".bim", ".fam")
                    if not Path(prefix + suffix).exists()
                ]
                rows.append(
                    {
                        "input": "AoU target genotypes",
                        "status": "PASS" if not missing_target_files else "MISSING",
                        "location": (
                            f"{len(prefixes)} chromosome BED prefix(es)"
                            if not missing_target_files
                            else Path(missing_target_files[0]).name
                        ),
                    }
                )
                if (
                    genome_build.value != "GRCh38"
                    and "vwb-aou-datasets-controlled" in bed_pattern.value
                ):
                    rows.append(
                        {
                            "input": "Genome-build compatibility",
                            "status": "REVIEW",
                            "location": "AoU v9 WGS target is GRCh38",
                        }
                    )
                if base_gwas.value != str(harmonized_gwas):
                    source_gwas["path"] = base_gwas.value
                display(pd.DataFrame(rows))
                if all(row["status"] == "PASS" for row in rows):
                    prepare_button.disabled = False
                    workflow_status.value = _workflow_status(
                        2, "Inputs ready; prepare the AoU target"
                    )
                    print("Inputs ready. Continue to Step 2.")
                else:
                    prepare_button.disabled = True
                    preview_button.disabled = True
                    run_button.disabled = True
                    results_button.disabled = True
                    workflow_status.value = _workflow_status(
                        1, "Resolve missing inputs"
                    )
                    print("Resolve MISSING or REVIEW rows before continuing.")
            except Exception as error:
                print(f"{type(error).__name__}: {error}")

    def handle_prepare(_):
        workflow_status.value = _workflow_status(
            2, "Preparing selected AoU variants and participants"
        )
        prepare_button.disabled = True
        with action_output:
            clear_output()
            if not workspace_ack.value:
                print("Confirm controlled-workspace handling before preparing genotypes.")
                prepare_button.disabled = False
                return
            try:
                schema = apply_schema()
                if schema["missing_required_fields"]:
                    raise ValueError(
                        "GWAS fields unresolved: "
                        + ", ".join(schema["missing_required_fields"])
                    )
                preparation_config = current_config()
                if source_gwas["path"]:
                    preparation_config.base_gwas = source_gwas["path"]
                custom_mapping = mapping_method.value == "custom"
                step_total = 4 if custom_mapping else 3
                if custom_mapping:
                    print(f"1/{step_total} Locating variants from the custom SNP-to-gene table...")
                    union_audit = build_variant_mapping_union_bed(
                        preparation_config,
                        variant_gene_file.value.strip(),
                        union_bed,
                    )
                else:
                    print(f"1/{step_total} Building the deduplicated pathway-variant union...")
                    union_audit = build_pathway_variant_union_bed(
                        preparation_config, union_bed
                    )
                display(pd.Series(union_audit, name="value").to_frame())
                print(f"2/{step_total} Harmonizing GWAS alleles with AoU BIM metadata...")
                display(
                    pd.Series(
                        harmonize_gwas_ids_to_bed(
                            preparation_config,
                            bed_pattern=bed_pattern.value.strip(),
                            chromosomes=selected_chromosomes(),
                            candidate_bed1_file=union_bed,
                            output_gwas=harmonized_gwas,
                        ),
                        name="value",
                    ).to_frame()
                )
                base_gwas.value = str(harmonized_gwas)
                apply_schema()
                if custom_mapping:
                    print(f"3/{step_total} Building pathway SNP sets from the custom mapping...")
                    display(
                        pd.Series(
                            build_snp_set_from_variant_gene_mapping(
                                mapping_file=variant_gene_file.value.strip(),
                                gmt_file=gmt_file.value.strip(),
                                output_file=custom_snp_set,
                                harmonized_gwas_file=harmonized_gwas,
                                harmonized_snp_column=column_widgets["snp"].value,
                            ),
                            name="value",
                        ).to_frame()
                    )
                print(f"{step_total}/{step_total} Registering existing AoU chromosome-wise BED files...")
                target_list_file.parent.mkdir(parents=True, exist_ok=True)
                prefixes = [
                    bed_pattern.value.replace("#", str(chrom)).replace("{chr}", str(chrom))
                    for chrom in selected_chromosomes()
                ]
                missing = [prefix + ".bed" for prefix in prefixes if not Path(prefix + ".bed").exists()]
                if missing:
                    raise FileNotFoundError(f"AoU BED file not found: {missing[0]}")
                target_list_file.write_text("\n".join(prefixes) + "\n")
                preview_button.disabled = False
                workflow_status.value = _workflow_status(
                    3, "AoU target ready; review the analysis"
                )
                print("AoU target metadata are ready. No genotype conversion was performed.")
            except Exception as error:
                prepare_button.disabled = False
                print(f"{type(error).__name__}: {error}")

    def handle_preview(_):
        workflow_status.value = _workflow_status(3, "Reviewing QC and command")
        with action_output:
            clear_output()
            try:
                config = current_config()
                checks, details = validate_inputs(config)
                display(checks)
                if (checks["status"] == "FAIL").any():
                    print("Preview stopped because required inputs failed validation.")
                    return
                display(
                    HTML(
                        "<pre style='white-space:pre-wrap'>"
                        + command_as_shell(build_pathway_pgs_command(config))
                        + "</pre>"
                    )
                )
                manifest = run_pathway_pgs(config, execute=False)
                display(manifest)
                if "pathway_summary" in details:
                    display(plot_pathway_definition_qc(details["pathway_summary"]))
                    plt.show()
                run_button.disabled = False
                workflow_status.value = _workflow_status(
                    4, "Review complete; calculate pathway scores"
                )
                print("Preview complete. Nothing was executed.")
            except Exception as error:
                print(f"{type(error).__name__}: {error}")

    def handle_run(_):
        workflow_status.value = _workflow_status(4, "Calculating pathway scores")
        run_button.disabled = True
        with action_output:
            clear_output()
            if not workspace_ack.value or not execute_ack.value:
                print("Confirm both safety boxes before running.")
                run_button.disabled = False
                return
            try:
                manifest = run_pathway_pgs(current_config(), execute=True)
                display(manifest)
                results_button.disabled = False
                workflow_status.value = _workflow_status(
                    5, "Scoring complete; view aggregate results"
                )
                print("Pathway scoring complete. Continue to Step 5.")
            except Exception as error:
                run_button.disabled = False
                print(f"{type(error).__name__}: {error}")

    def handle_results(_):
        workflow_status.value = _workflow_status(5, "Aggregate results")
        with action_output:
            clear_output()
            try:
                config = current_config()
                summary = summarize_aggregate_results(config)
                if summary.empty:
                    print("No aggregate result file was found for this run.")
                    return
                if "score_sd" in summary:
                    finite = int(summary["all_finite"].sum()) if "all_finite" in summary else len(summary)
                    display(
                        HTML(
                            "<div style='display:flex;gap:10px;flex-wrap:wrap;margin-bottom:12px'>"
                            f"<div style='border:1px solid #d8e0e5;padding:10px 16px'><b>{len(summary)}</b><br>scores</div>"
                            f"<div style='border:1px solid #d8e0e5;padding:10px 16px'><b>{finite}</b><br>finite scores</div>"
                            f"<div style='border:1px solid #d8e0e5;padding:10px 16px'><b>{summary['n'].max():,}</b><br>participants</div>"
                            "</div>"
                        )
                    )
                display(summary.head(30).style.hide(axis="index"))
                figure = plot_pathway_pgs_results(summary)
                figure_dir = Path(config.output_dir) / "figures"
                figure_dir.mkdir(parents=True, exist_ok=True)
                figure.savefig(figure_dir / "pathway_pgs_results.png", dpi=300, bbox_inches="tight")
                figure.savefig(figure_dir / "pathway_pgs_results.pdf", bbox_inches="tight")
                display(figure)
                plt.close(figure)
                report = write_markdown_report(config, summary)
                print(f"Report: {report}")
                print(f"Figures: {figure_dir}")
            except Exception as error:
                print(f"{type(error).__name__}: {error}")

    project_name.observe(
        lambda change: setattr(output_prefix, "value", change["new"].strip())
        if change["new"].strip()
        else None,
        names="value",
    )
    run_mode.observe(apply_run_mode, names="value")
    mapping_method.observe(refresh_mapping_controls, names="value")
    annotation_source.observe(select_reference_files, names="value")
    pathway_source.observe(select_reference_files, names="value")
    auto_find.on_click(refresh_local_inputs)
    install_engine.on_click(handle_install_engine)
    detect_gwas.on_click(handle_detect)
    check_button.on_click(handle_check)
    prepare_button.on_click(handle_prepare)
    preview_button.on_click(handle_preview)
    run_button.on_click(handle_run)
    results_button.on_click(handle_results)

    required_inputs = widgets.VBox(
        [
            widgets.HTML("<b>Required inputs</b>"),
            run_mode,
            project_name,
            widgets.HBox(
                [engine_status, install_engine],
                layout=widgets.Layout(
                    width="100%",
                    flex_flow="row wrap",
                    align_items="center",
                ),
            ),
            setup_message,
            demo_summary,
            mapping_method,
            annotation_source,
            pathway_source,
            base_gwas,
            gtf_file,
            gmt_file,
            variant_gene_file,
            keep_file,
            widgets.HBox(
                [auto_find, detect_gwas],
                layout=widgets.Layout(flex_flow="row wrap"),
            ),
            input_message,
            upload_help,
        ]
    )
    required_inputs.layout = CARD
    workflow_buttons = widgets.HBox(
        [check_button, prepare_button, preview_button, run_button, results_button],
        layout=widgets.Layout(
            width="100%",
            flex_flow="row wrap",
        ),
    )
    apply_run_mode()
    refresh_engine_status()
    workflow_section = widgets.VBox(
        [workflow_status, workflow_buttons, action_output],
        layout=CARD,
    )
    return widgets.VBox(
        [
            title,
            required_inputs,
            widgets.VBox(
                [widgets.HTML("<b>Advanced settings</b>"), advanced], layout=CARD
            ),
            widgets.VBox(
                [
                    widgets.HTML("<b>Safety confirmation</b>"),
                    workspace_ack,
                    execute_ack,
                ],
                layout=CARD,
            ),
            widgets.HTML("<b>Workflow</b>"),
            workflow_section,
        ],
        layout=_full_layout(),
    )
