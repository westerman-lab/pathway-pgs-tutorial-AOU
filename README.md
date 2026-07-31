# Pathway PGS Tutorial for All of Us

A Jupyter tutorial for constructing pathway-level polygenic scores
(pPGS) in the All of Us Researcher Workbench controlled workspace.

The notebook teaches the analysis rather than hiding it behind a one-click
tool. Each stage exposes the scientific choice, executable code, exact scoring
command, and quality-control output. An optional guided interface calls the
same functions after the reader has reviewed the method.

## What You Will Learn

The tutorial walks through the complete pathway-scoring workflow:

1. read GWAS effect alleles, effect sizes, and p-values;
2. map variants to genes using physical distance or a custom mapping table;
3. map genes to biological pathways;
4. harmonize GWAS variants and alleles with All of Us genotypes;
5. apply pathway-specific LD clumping and p-value thresholds;
6. calculate participant-level pathway scores;
7. inspect aggregate QC, figures, and a reproducible report.

For participant `i`, pathway `k`, and p-value threshold `tau`,

```text
P[i,k,tau] = sum_j G[i,j] * beta[j] * A[j,k] * L[j,k] * I(p[j] <= tau)
```

- `G[i,j]`: effect-allele dosage for variant `j`.
- `beta[j]`: GWAS effect estimate.
- `A[j,k]`: 1 when variant `j` maps to a gene in pathway `k`.
- `L[j,k]`: 1 when the variant remains after pathway-specific LD clumping.
- `I(p[j] <= tau)`: GWAS p-value threshold indicator.

## Start With the Quick Demo

The quickest way to understand the workflow is the built-in technical demo.
It uses 100 All of Us participants, chromosome 21, two synthetic variants,
and two synthetic pathways. It verifies the software and genotype connection;
it is not a scientific result.

In an All of Us Workbench terminal:

```bash
cd ~
git clone https://github.com/westerman-lab/pathway-pgs-tutorial-AOU.git
cd pathway-pgs-tutorial-AOU
pip install --user -r requirements.txt
```

Open `Pathway_PGS_AoU_Tutorial.ipynb` in JupyterLab. Run the cells in order.
Do not use **Run All** for a full analysis because preparation and scoring are
intentional checkpoints.

For the first demo:

1. Leave `RUN_MODE = "demo"`.
2. Review Sections 0-3.
3. Set `PREPARE_AOU_DATA = True` and run Section 4.
4. Review the validation table and exact command in Section 5.
5. Set `RUN_SCORING = True` and run Section 6.
6. Run Section 7 to create aggregate QC, figures, and the report.

Both expensive switches are `False` by default:

```python
PREPARE_AOU_DATA = False
RUN_SCORING = False
```

## Move From Demo to Scientific Analysis

Change `RUN_MODE` to `"full"`, then provide the four scientific inputs below.

| Input | Required content | Main decision |
|---|---|---|
| GWAS summary statistics | Chromosome, position, variant ID, effect and other alleles, beta or OR, and p-value | Trait, ancestry, genome build, and effect allele |
| SNP-to-gene definition | Genome-build-matched GTF or a custom regulatory mapping table | Physical distance, eQTL, chromatin, or another justified map |
| Pathway definition | Reactome, GO Biological Process, WikiPathways, or custom GMT | Database version and pathway scope |
| Participant keep file | PLINK FID and IID values inside the controlled workspace | Analysis cohort and exclusions |

All of Us v9 WGS uses GRCh38. GRCh37 files must be lifted and harmonized before
they are combined with the default target.

### Scientific Choices to Record

| Choice | Why it matters |
|---|---|
| Effect statistic | `BETA` is used directly; an odds ratio requires the corresponding OR setting |
| Effect allele | Dosage must refer to the same allele as the GWAS effect estimate |
| SNP-to-gene mapping | Changes which variants enter each pathway |
| Gene windows | Changes physical-distance pathway membership |
| LD clumping | Reduces correlated variant redundancy within each pathway |
| P-value thresholds | Controls how much GWAS signal enters each pathway score |
| Score method | Defines how retained weighted dosages are aggregated |

## Notebook Map

| Section | Purpose | Expected output |
|---|---|---|
| 0. Load code | Import the visible workflow functions | Source locations |
| 1. Choices | Set mapping, clumping, threshold, and execution controls | Reproducible configuration |
| 2. Inputs | Detect GWAS columns and inspect pathway and GTF files | Input QC |
| 3. Mapping | Compare physical and custom SNP-to-gene routes | Selected mapping branch |
| 4. AoU preparation | Build a deduplicated variant union and harmonize alleles | Canonical GWAS and target list |
| 5. Review | Validate inputs and print the exact command | Dry-run manifest |
| 6. Scoring | Run pathway-specific clumping and scoring | Participant-by-pathway scores |
| 7. Results | Summarize scores without exporting participant IDs | Aggregate table, figures, report |

The optional interface at the end mirrors these same stages:
**Check inputs**, **Map & harmonize**, **Review command**, **Calculate scores**,
and **Inspect results**.

## Understanding the Outputs

The default output directory is:

```text
/home/jupyter/analysis/results/pathway_prs_tutorial/
```

Important outputs include:

- `run_config.json`: all analysis parameters;
- `run_manifest.json`: exact command, status, and log tail;
- participant-by-pathway score output: controlled-workspace data;
- `*.aggregate_score_summary.tsv`: score count, mean, SD, minimum, and maximum;
- `figures/pathway_pgs_results.png` and `.pdf`: aggregate QC figures;
- Markdown report: concise methods and output inventory.

A successful technical run should have finite scores and nonzero variability
for pathways containing informative variants. Biological interpretation still
depends on appropriate GWAS, mapping, pathway, ancestry, and phenotype choices.

## Data Safety

Participant-level genotype and score files must remain inside the All of Us
controlled workspace. Do not commit, download, or display participant IDs or
person-level scores. The tutorial's result view reports aggregate statistics
only.

## Troubleshooting

**The scoring engine is missing**

Use the setup control once, or provide existing paths to the R wrapper and
Linux executable. Then rerun input validation.

**A GWAS column is not detected**

Inspect the header in Section 2. Confirm chromosome, position, SNP, alleles,
effect statistic, and p-value. The harmonized GWAS uses canonical column names.

**No variants match the AoU target**

Check genome build, chromosome naming, position, REF/ALT alleles, and whether
the selected variants exist in the target genotype release.

**A score has zero variance**

Check the number of retained variants, LD settings, p-value threshold, allele
harmonization, and minor-allele frequency in the selected participants.

**The notebook appears to stop before scoring**

This is expected while the safety switches are `False`. Review the preceding
section before enabling each expensive step.

## Repository Contents

- `Pathway_PGS_AoU_Tutorial.ipynb`: code-visible teaching notebook;
- `pathway_prs_core.py`: mapping, harmonization, command, QC, and plotting code;
- `pathway_pgs_app.py`: optional guided interface;
- `demo_data/`: synthetic demonstration inputs;
- `build_notebook.py`: reproducible notebook generator;
- `example_config.json`: configuration example;
- `test_pathway_prs_core.py`: unit and synthetic workflow tests;
- `test_notebook_execution.py`: safe top-to-bottom notebook test;
- `test_notebook_structure.py`: educational-structure tests;
- `CHANGELOG.md`: release history.

## Development and Verification

```bash
python3 build_notebook.py
python3 -m pytest -q
```

The automated suite checks mapping, harmonization, command construction,
synthetic score execution, aggregate result loading, interface launch, and safe
top-to-bottom notebook execution. Real All of Us WGS execution must be tested
inside an authorized controlled workspace.

Current release: `v2.4.0`.
