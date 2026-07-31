#!/usr/bin/env python3
"""Build the compact lab-facing Pathway PGS notebook."""

from pathlib import Path
from textwrap import dedent

import nbformat as nbf


ROOT = Path(__file__).resolve().parent
DESTINATION = ROOT / "Pathway_PGS_AoU_Tutorial.ipynb"


def markdown(source: str):
    return nbf.v4.new_markdown_cell(dedent(source).strip())


def code(source: str):
    cell = nbf.v4.new_code_cell(dedent(source).strip())
    cell.metadata["tags"] = ["hide-input"]
    cell.metadata["jupyter"] = {"source_hidden": True}
    return cell


notebook = nbf.v4.new_notebook()
notebook.metadata.update(
    {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3"},
    }
)
notebook.cells = [
    markdown(
        r"""
        # Pathway Polygenic Score Tutorial

        **Made by Sijia Zhu, Westerman's Lab**

        This notebook provides a guided workflow for constructing pathway-level
        polygenic scores inside the All of Us controlled workspace. Start with
        **Quick demo** to produce example scores, QC, and figures using 100
        participants. Then switch to **Full analysis** for scientific data.

        The interface guides the user through five steps: check, prepare,
        review, calculate, and view results. It can locate or install the
        supported scoring engine and detects standard GWAS columns
        automatically.

        The recommended AoU workflow uses the existing chromosome-wise PLINK
        BED files directly. It does not create a second genotype copy.
        """
    ),
    markdown(
        r"""
        ## Statistical workflow

        | Stage | Input | Output |
        |---|---|---|
        | SNP weights | GWAS summary statistics | Effect allele, beta, and p-value |
        | SNP to gene | Physical windows or a custom mapping | Gene-level SNP sets |
        | Gene to pathway | Pathway GMT | Pathway-level SNP sets |
        | C+T | AoU genotypes | LD-clumped pathway variants |
        | Scoring | Dosage and GWAS beta | Participant-by-pathway PGS matrix |

        For participant $i$ and pathway $k$,

        $$
        P_{ik}
        =
        \sum_j
        G_{ij}\,\widehat{\beta}_j\,
        A_{jk}\,L_{jk}\,
        \mathbf{1}(p_j\leq\tau).
        $$

        - $G_{ij}$ is effect-allele dosage.
        - $\widehat{\beta}_j$ is the GWAS effect estimate.
        - $A_{jk}$ indicates SNP-to-pathway membership.
        - $L_{jk}$ indicates retention after pathway-specific LD clumping.
        - $\tau$ is the selected GWAS p-value threshold.
        """
    ),
    markdown(
        r"""
        ## Start

        Run the single cell below. The interface finds common local inputs and
        detects GWAS columns automatically. Optional settings remain collapsed.

        Quick demo requires no scientific input files. It uses two synthetic
        chromosome 21 variants and two synthetic pathways solely to verify that
        the software and AoU genotype connection work. Results include an
        aggregate QC table, PNG and PDF figures, and a concise report.

        Full analysis requires:

        1. GWAS summary statistics;
        2. a genome-matched GTF or custom SNP-to-gene mapping;
        3. Reactome, GO Biological Process, WikiPathways, or a custom GMT;
        4. participant keep file.
        """
    ),
    code(
        r"""
        from pathlib import Path
        import sys
        from IPython.display import display

        tutorial_dir = Path.cwd()
        if not (tutorial_dir / "pathway_pgs_app.py").exists():
            tutorial_dir = tutorial_dir / "pathway_prs_tutorial"
        if not (tutorial_dir / "pathway_pgs_app.py").exists():
            raise FileNotFoundError("Open this notebook from the tutorial folder.")

        sys.path.insert(0, str(tutorial_dir.resolve()))
        from pathway_pgs_app import launch_pathway_pgs_app

        print("Pathway PGS interface initialized. Controls should appear below.")
        display(launch_pathway_pgs_app())
        """
    ),
]

nbf.write(notebook, DESTINATION)
print(DESTINATION)
