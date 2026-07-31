# Interactive Pathway PGS Tutorial for All of Us

A code-first Jupyter tutorial for constructing pathway-level polygenic scores
inside the All of Us Researcher Workbench.

## Purpose

This repository is designed for teaching and reproducible analysis. The
notebook walks through the code and scientific choices rather than treating
pathway scoring as a one-click black box.

Readers can inspect and modify:

1. GWAS column detection and effect-allele definitions;
2. physical or custom SNP-to-gene mapping;
3. gene-to-pathway membership;
4. allele harmonization with All of Us genotypes;
5. LD-clumping and p-value thresholds;
6. the exact scoring command, QC, figures, and output report.

An optional interactive helper is included at the end of the notebook. It
calls the same functions shown in the visible code cells.

## Score definition

For participant `i`, pathway `k`, and GWAS p-value threshold `tau`:

```text
P[i,k] = sum_j G[i,j] * beta[j] * A[j,k] * L[j,k] * I(p[j] <= tau)
```

`G` is effect-allele dosage, `beta` is the GWAS effect estimate, `A` is
variant-to-pathway membership, and `L` indicates retention after
pathway-specific LD clumping.

## Quick start in All of Us

Clone the repository in a Workbench terminal:

```bash
cd ~
git clone https://github.com/westerman-lab/pathway-pgs-tutorial-AOU.git
cd pathway-pgs-tutorial-AOU
pip install --user -r requirements.txt
```

Open `Pathway_PGS_AoU_Tutorial.ipynb` in JupyterLab and run the cells in order.
Do not use **Run All** for a full analysis.

The notebook defaults to a short technical demonstration with 100 AoU
participants, chromosome 21, two synthetic variants, and two synthetic
pathways. It is not a scientific result.

## Notebook sections

| Section | What the reader sees |
|---|---|
| 0. Load code | Exact imported functions and source-file location |
| 1. Choices | Editable genome build, mapping, clumping, and scoring parameters |
| 2. Inputs | Detected GWAS columns, pathway contents, and GTF QC |
| 3. Mapping | Physical-distance and custom SNP-to-gene alternatives |
| 4. AoU preparation | Deduplicated variant union and allele harmonization code |
| 5. Review | Validation table and complete generated command |
| 6. Scoring | Explicitly enabled score calculation |
| 7. Results | Aggregate QC, figures, and report |

Expensive operations are off by default:

```python
PREPARE_AOU_DATA = False
RUN_SCORING = False
```

Set each value to `True` only after reviewing its section.

## Scientific inputs

| Input | Required content |
|---|---|
| GWAS summary statistics | Chromosome, position, variant ID, effect allele, other allele, beta or OR, and p-value |
| Gene annotation | Genome-build-matched GTF for physical mapping, or a custom SNP-to-gene table |
| Pathway definition | Reactome, GO Biological Process, WikiPathways, or a custom GMT |
| Participant keep file | Intended PLINK FID and IID values inside the controlled workspace |

All of Us v9 WGS uses GRCh38. GRCh37 inputs must be lifted or otherwise
harmonized before use with the default target.

## Mapping choices

**Physical distance** assigns variants within gene bodies and configurable,
strand-aware upstream and downstream windows.

**Custom mapping** accepts eQTL, regulatory, chromatin, or other maps. The
table must include a gene column and either a variant-ID column or chromosome
and position columns.

## Outputs

The default output directory is:

```text
/home/jupyter/analysis/results/pathway_prs_tutorial/
```

It contains the saved configuration, exact command manifest, aggregate QC,
PNG/PDF figures, and a Markdown report. Person-level scores must remain inside
the controlled All of Us workspace and must not be committed to this
repository.

## Repository contents

- `Pathway_PGS_AoU_Tutorial.ipynb`: code-visible teaching notebook;
- `pathway_prs_core.py`: mapping, harmonization, scoring, QC, and plotting;
- `pathway_pgs_app.py`: optional guided interface;
- `demo_data/`: synthetic technical demonstration inputs;
- `build_notebook.py`: reproducible notebook generator;
- `example_config.json`: configuration example;
- `test_pathway_prs_core.py`: workflow tests;
- `test_notebook_structure.py`: educational-structure tests.
- `CHANGELOG.md`: release notes.

## Testing

```bash
python3 build_notebook.py
python3 -m pytest -q
```

Current release: `v2.3.0`.
