#!/usr/bin/env python3
"""Build the code-first All of Us Pathway PGS tutorial notebook."""

import json
from pathlib import Path
from textwrap import dedent


ROOT = Path(__file__).resolve().parent
DESTINATION = ROOT / "Pathway_PGS_AoU_Tutorial.ipynb"


def markdown(source: str):
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": dedent(source).strip().splitlines(keepends=True),
    }


def code(source: str, *, tags: list[str] | None = None):
    """Create a visible code cell.

    The tutorial intentionally exposes executable code so readers can inspect
    each method and parameter choice instead of treating the workflow as a
    black box.
    """
    cell = {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": dedent(source).strip().splitlines(keepends=True),
    }
    if tags:
        cell["metadata"]["tags"] = tags
    return cell


notebook = {
    "cells": [],
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

notebook["cells"] = [
    markdown(
        r"""
        # Interactive Pathway PGS Tutorial for All of Us

        **Made by Sijia Zhu, Westerman Lab**

        This code-first notebook explains how pathway polygenic scores are
        constructed in the All of Us controlled workspace. Each section shows
        the methodological choice, the Python code, the generated command, and
        the relevant quality-control output.

        The optional interface at the end is a convenience layer. It does not
        replace the step-by-step analysis shown in the notebook.
        """
    ),
    markdown(
        r"""
        ## Learning goals

        After completing the tutorial, a reader should be able to:

        1. identify the required GWAS, gene, pathway, and genotype inputs;
        2. explain and modify SNP-to-gene and gene-to-pathway mapping;
        3. inspect allele harmonization with All of Us genotypes;
        4. choose LD-clumping and GWAS p-value thresholds;
        5. reproduce the pathway-score command and interpret its QC outputs.

        Run the cells **in order**. Expensive operations are disabled by
        default and require an explicit change from `False` to `True`.
        """
    ),
    markdown(
        r"""
        ## Statistical definition

        For participant $i$, pathway $k$, and GWAS p-value threshold $\tau$,

        $$
        P_{ik}^{(\tau)}
        =
        \sum_j
        G_{ij}\,\widehat{\beta}_j\,
        A_{jk}\,L_{jk}\,
        \mathbf{1}(p_j\leq\tau).
        $$

        - $G_{ij}$: effect-allele dosage for variant $j$.
        - $\widehat{\beta}_j$: GWAS effect estimate.
        - $A_{jk}$: variant-to-pathway membership after SNP-to-gene mapping.
        - $L_{jk}$: retention after pathway-specific LD clumping.
        - $\mathbf{1}(p_j\leq\tau)$: GWAS p-value threshold indicator.

        A variant may contribute to several pathways, but it is counted only
        once within any one pathway after clumping.
        """
    ),
    markdown(
        r"""
        ## 0. Load the tutorial code

        This cell imports the exact functions used below. Their complete
        implementations are available in `pathway_prs_core.py` and can be
        opened beside the notebook at any time.
        """
    ),
    code(
        r"""
        from pathlib import Path
        import json
        import shutil
        import sys

        import matplotlib.pyplot as plt
        import pandas as pd
        from IPython.display import display

        candidates = [
            Path.cwd(),
            Path.cwd() / "pathway-pgs-tutorial-AOU",
            Path.cwd() / "pathway-pgs-tutorial",
            Path.cwd() / "pathway_prs_tutorial",
            Path.home() / "pathway-pgs-tutorial-AOU",
            Path.home() / "pathway-pgs-tutorial",
            Path.home() / "pathway_prs_tutorial",
        ]
        tutorial_dir = next(
            (path for path in candidates if (path / "pathway_prs_core.py").exists()),
            None,
        )
        if tutorial_dir is None:
            raise FileNotFoundError("Open this notebook from the tutorial folder.")

        sys.path.insert(0, str(tutorial_dir.resolve()))

        from pathway_pgs_app import discover_local_inputs, launch_pathway_pgs_app
        from pathway_prs_core import (
            WorkflowConfig,
            build_pathway_pgs_command,
            build_pathway_variant_union_bed,
            build_snp_set_from_variant_gene_mapping,
            build_variant_mapping_union_bed,
            command_as_shell,
            harmonize_gwas_ids_to_bed,
            infer_gwas_schema,
            plot_pathway_definition_qc,
            plot_pathway_pgs_results,
            read_gmt_summary,
            read_gtf_summary,
            run_pathway_pgs,
            summarize_aggregate_results,
            validate_inputs,
            write_markdown_report,
        )

        print(f"Tutorial code: {tutorial_dir.resolve()}")
        print("Core implementation: pathway_prs_core.py")
        """
    ),
    markdown(
        r"""
        ## 1. Make the analysis choices

        Edit this cell before a scientific run. The demonstration uses two
        synthetic chromosome 21 variants and 100 All of Us participants. It
        verifies the workflow only and is not a scientific result.

        Important choices are kept together so the analysis can be reviewed
        and reproduced. `PREPARE_AOU_DATA` and `RUN_SCORING` remain `False`
        until the corresponding sections are reviewed.
        """
    ),
    code(
        r"""
        # --------------------------- User choices ---------------------------
        RUN_MODE = "demo"                 # "demo" or "full"
        GENOME_BUILD = "GRCh38"           # AoU v9 WGS is GRCh38
        MAPPING_METHOD = "physical"        # "physical" or "custom"
        PATHWAY_DATABASE = "Reactome"      # descriptive provenance label

        WINDOW_5_BP = 35_000               # upstream, relative to gene strand
        WINDOW_3_BP = 10_000               # downstream, relative to gene strand
        P_VALUE_THRESHOLDS = [1.0]
        CLUMP_KB = 1_000
        CLUMP_R2 = 0.10
        SCORE_METHOD = "sum"
        THREADS = 4

        PREPARE_AOU_DATA = False           # change only in Section 4
        RUN_SCORING = False                # change only in Section 6
        KEEP_PERSON_LEVEL_DATA_IN_AOU = True
        # -------------------------------------------------------------------

        home = Path.home()
        demo_dir = tutorial_dir / "demo_data"
        work_dir = home / "analysis" / "data" / "pathway_prs_adapter"
        output_dir = home / "analysis" / "results" / "pathway_prs_tutorial"
        discovered = discover_local_inputs()

        if RUN_MODE == "demo":
            run_name = "pathway_pgs_quick_demo"
            gwas_file = demo_dir / "demo_gwas_chr21.tsv"
            gtf_file = demo_dir / "demo_chr21.gtf"
            gmt_file = demo_dir / "demo_pathways.gmt"
            custom_mapping_file = Path("")
            keep_file = work_dir / "quick_demo_100.keep"
            chromosomes = [21]
        else:
            run_name = "pathway_pgs_full_analysis"
            gwas_file = Path(discovered["gwas"] or "/path/to/gwas.tsv.gz")
            gtf_file = Path(discovered["gtf"] or "/path/to/genes.gtf.gz")
            gmt_file = Path(discovered["gmt"] or "/path/to/pathways.gmt")
            custom_mapping_file = Path("/path/to/custom_snp_to_gene.tsv")
            keep_file = work_dir / "participants.keep"
            chromosomes = list(range(1, 23))

        bed_pattern = (
            "/home/jupyter/workspace/vwb-aou-datasets-controlled-v9/v9/"
            "wgs/short_read/snpindel/acaf_threshold/plink_bed/"
            "acaf_threshold.chr#"
        )
        union_bed = work_dir / "pathway_variant_union.bed1.tsv"
        harmonized_gwas = work_dir / "gwas_for_aou_bed.tsv.gz"
        target_list = work_dir / "aou_bed_target_prefixes.txt"
        custom_snp_set = work_dir / "custom_pathway_snp_sets.gmt"

        choices = pd.Series({
            "run mode": RUN_MODE,
            "genome build": GENOME_BUILD,
            "SNP-to-gene mapping": MAPPING_METHOD,
            "pathway database": PATHWAY_DATABASE,
            "5-prime window (bp)": WINDOW_5_BP,
            "3-prime window (bp)": WINDOW_3_BP,
            "p-value thresholds": P_VALUE_THRESHOLDS,
            "clump window (kb)": CLUMP_KB,
            "clump r-squared": CLUMP_R2,
            "score method": SCORE_METHOD,
            "chromosomes": chromosomes,
        }, name="selected value")
        display(choices.to_frame())
        """,
        tags=["parameters"],
    ),
    markdown(
        r"""
        ## 2. Inspect the scientific inputs

        The GWAS header is detected rather than assumed. The code also checks
        the pathway file and gene annotation before any participant genotype
        is read. Review the displayed column mapping carefully, especially the
        effect allele (`a1`) and effect statistic (`stat`).
        """
    ),
    code(
        r"""
        schema = infer_gwas_schema(gwas_file)

        detected_columns = pd.DataFrame({
            "required role": list(schema["gwas_columns"].keys()),
            "detected GWAS column": list(schema["gwas_columns"].values()),
        })
        display(detected_columns)
        print("Effect statistic:", schema["statistic_type"])
        print("Separator:", schema["separator"])
        print("Unresolved fields:", schema["missing_required_fields"] or "none")

        pathway_summary, pathway_gene_preview = read_gmt_summary(gmt_file)
        display(pathway_summary.head())
        display(pathway_gene_preview.head(10))

        if MAPPING_METHOD == "physical":
            display(pd.Series(read_gtf_summary(gtf_file), name="GTF QC").to_frame())

        figure = plot_pathway_definition_qc(pathway_summary)
        display(figure)
        plt.close(figure)
        """
    ),
    markdown(
        r"""
        ## 3. Define variant-to-pathway membership

        ### Physical-distance mapping

        A variant maps to a gene when its position falls inside the gene body
        or the selected strand-aware flanking windows. Genes are then assigned
        to pathways using the GMT membership file.

        ### Custom mapping

        A user may instead supply an eQTL, chromatin, regulatory, or other
        variant-to-gene table. It must contain a gene column and either a
        variant-ID column or chromosome and position columns.

        The code below creates one deduplicated union of required variants.
        This union limits later genotype access; it does not yet calculate any
        participant score.
        """
    ),
    code(
        r"""
        def make_config(base_path):
            current_schema = infer_gwas_schema(base_path)
            return WorkflowConfig(
                project_name=run_name,
                genome_build=GENOME_BUILD,
                base_gwas=str(base_path),
                base_separator=current_schema["separator"],
                gwas_columns=current_schema["gwas_columns"],
                statistic_type=current_schema["statistic_type"],
                target_type="bed",
                target_list=str(target_list),
                target_keep_file=str(keep_file),
                pathway_input_mode=(
                    "gtf_gmt" if MAPPING_METHOD == "physical" else "snp_set"
                ),
                gtf_file=str(gtf_file),
                gmt_file=str(gmt_file),
                snp_set_file=str(custom_snp_set),
                window_5_bp=WINDOW_5_BP,
                window_3_bp=WINDOW_3_BP,
                pvalue_thresholds=P_VALUE_THRESHOLDS,
                clump_kb=CLUMP_KB,
                clump_r2=CLUMP_R2,
                score_method=SCORE_METHOD,
                threads=THREADS,
                prsice_r=discovered["wrapper"],
                prsice_binary=discovered["executable"],
                rscript=discovered["rscript"],
                output_dir=str(output_dir),
                output_prefix=run_name,
                controlled_workspace_acknowledged=KEEP_PERSON_LEVEL_DATA_IN_AOU,
            )

        source_config = make_config(gwas_file)

        print("The preparation cell will execute one of these two branches:")
        print("  physical -> build_pathway_variant_union_bed(...)" )
        print("  custom   -> build_variant_mapping_union_bed(...)" )
        print("Selected branch:", MAPPING_METHOD)
        """
    ),
    markdown(
        r"""
        ## 4. Prepare the All of Us target

        This step performs three transparent operations:

        1. create the deduplicated pathway-variant union;
        2. match chromosome, position, REF, and ALT against AoU BED metadata;
        3. register the existing chromosome-wise AoU files and selected people.

        It does **not** copy the full WGS dataset. Set `PREPARE_AOU_DATA = True`
        in Section 1 only after reviewing the inputs and this code.
        """
    ),
    code(
        r"""
        def ensure_demo_keep_file():
            # Select the first 100 demo participants from the AoU FAM file.
            if RUN_MODE != "demo" or keep_file.exists():
                return
            prefix = bed_pattern.replace("#", "21")
            fam = Path(prefix + ".fam")
            if not fam.exists():
                raise FileNotFoundError(f"AoU chromosome 21 FAM not found: {fam}")
            ids = pd.read_csv(fam, sep=r"\s+", header=None, usecols=[0, 1], nrows=100)
            keep_file.parent.mkdir(parents=True, exist_ok=True)
            ids.to_csv(keep_file, sep="\t", header=False, index=False)

        if not PREPARE_AOU_DATA:
            print("Preparation is OFF. Review this cell, then set PREPARE_AOU_DATA = True.")
        else:
            ensure_demo_keep_file()
            work_dir.mkdir(parents=True, exist_ok=True)

            if MAPPING_METHOD == "physical":
                mapping_qc = build_pathway_variant_union_bed(source_config, union_bed)
            else:
                mapping_qc = build_variant_mapping_union_bed(
                    source_config, custom_mapping_file, union_bed
                )
            display(pd.Series(mapping_qc, name="mapping QC").to_frame())

            harmonization_qc = harmonize_gwas_ids_to_bed(
                source_config,
                bed_pattern=bed_pattern,
                chromosomes=chromosomes,
                candidate_bed1_file=union_bed,
                output_gwas=harmonized_gwas,
            )
            display(pd.Series(harmonization_qc, name="harmonization QC").to_frame())

            if MAPPING_METHOD == "custom":
                snp_set_qc = build_snp_set_from_variant_gene_mapping(
                    mapping_file=custom_mapping_file,
                    gmt_file=gmt_file,
                    output_file=custom_snp_set,
                    harmonized_gwas_file=harmonized_gwas,
                    harmonized_snp_column="SNP",
                )
                display(pd.Series(snp_set_qc, name="custom SNP-set QC").to_frame())

            prefixes = [bed_pattern.replace("#", str(chrom)) for chrom in chromosomes]
            target_list.write_text("\n".join(prefixes) + "\n")
            print("Prepared GWAS:", harmonized_gwas)
            print("AoU target list:", target_list)
        """
    ),
    markdown(
        r"""
        ## 5. Review the exact scoring model

        After preparation, this section validates all inputs and prints the
        exact command. Nothing is scored here. Check the effect statistic,
        allele columns, pathway mapping, p-value thresholds, LD parameters,
        participant file, and output path before continuing.
        """
    ),
    code(
        r"""
        if not harmonized_gwas.exists():
            print("Prepared GWAS not found. Complete Section 4 first.")
            scoring_config = None
        else:
            scoring_config = make_config(harmonized_gwas)
            checks, details = validate_inputs(scoring_config)
            display(checks)

            scoring_command = build_pathway_pgs_command(scoring_config)
            print("\nExact command (review only):\n")
            print(command_as_shell(scoring_command))

            dry_run_manifest = run_pathway_pgs(scoring_config, execute=False)
            print("\nDry-run status:", dry_run_manifest["status"])
            print("Saved configuration:", dry_run_manifest["config_path"])
        """
    ),
    markdown(
        r"""
        ## 6. Calculate pathway scores

        This is the only cell that launches the scoring engine. It applies
        pathway-specific LD clumping and computes weighted dosage sums for the
        selected participants. Set `RUN_SCORING = True` in Section 1 only after
        the validation table and exact command in Section 5 are correct.
        """
    ),
    code(
        r"""
        if not RUN_SCORING:
            print("Scoring is OFF. Review Section 5, then set RUN_SCORING = True.")
        elif scoring_config is None:
            raise RuntimeError("Complete data preparation before scoring.")
        elif not KEEP_PERSON_LEVEL_DATA_IN_AOU:
            raise PermissionError("Confirm that person-level outputs remain inside AoU.")
        else:
            run_manifest = run_pathway_pgs(scoring_config, execute=True)
            display(pd.Series(run_manifest, name="run manifest").to_frame())
        """
    ),
    markdown(
        r"""
        ## 7. Inspect QC and results

        The final section reads aggregate summaries, checks whether scores are
        finite and variable, creates PNG/PDF figures, and writes a concise
        report. Participant-level score files remain in the controlled
        workspace and are not displayed or exported by this section.
        """
    ),
    code(
        r"""
        if scoring_config is None:
            print("No prepared scoring configuration is available.")
        else:
            summary = summarize_aggregate_results(scoring_config)
            if summary.empty:
                print("No completed score output found. Run Section 6 first.")
            else:
                display(summary.head(30))
                figure = plot_pathway_pgs_results(summary)
                figure_dir = Path(scoring_config.output_dir) / "figures"
                figure_dir.mkdir(parents=True, exist_ok=True)
                figure.savefig(
                    figure_dir / "pathway_pgs_results.png",
                    dpi=300,
                    bbox_inches="tight",
                )
                figure.savefig(
                    figure_dir / "pathway_pgs_results.pdf",
                    bbox_inches="tight",
                )
                display(figure)
                plt.close(figure)
                report = write_markdown_report(scoring_config, summary)
                print("Report:", report)
        """
    ),
    markdown(
        r"""
        ## Optional guided interface

        The interface below runs the same functions used in the visible cells
        above. It is useful for demonstrations and parameter exploration, but
        readers should first understand Sections 1-7 and review every generated
        command before execution.
        """
    ),
    code(
        r"""
        app = launch_pathway_pgs_app()
        display(app)
        """
    ),
]

DESTINATION.write_text(json.dumps(notebook, indent=1) + "\n")
print(DESTINATION)
