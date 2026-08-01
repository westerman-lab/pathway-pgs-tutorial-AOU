# From Variant Effects to Pathway PGS in All of Us

An end-to-end Jupyter workflow to construct the interpretable pathway-level polygenic
scores (pPGS) in All of Us workbench.

The tutorial connects every stage of the analysis in one reproducible workflow:

**Variants -> genes -> biological pathways -> individual-level pPGS**

Users can:

- connect GWAS summary statistics with All of Us genotype data;
- map variants to genes and biological pathways;
- calculate pathway-level polygenic scores for selected participants;
- customize SNP-to-gene mapping, pathway definitions, LD clumping, and p-value thresholds;
- automatically generate QC summaries, figures, and reproducible reports.

All participant-level data remain within the All of Us controlled workspace.

**Made by Westerman Lab**

## Pathway Score Definition

For participant $i$, pathway $k$, and GWAS p-value threshold $\tau$,

$$
P_{ik}^{(\tau)}=
\sum_j G_{ij}\widehat{\beta}_jA_{jk}L_{jk}
\mathbf{1}(p_j\leq\tau).
$$

- $G_{ij}$ is the effect-allele dosage for variant $j$.
- $\widehat{\beta}_j$ is the GWAS effect estimate.
- $A_{jk}=1$ when variant $j$ maps to a gene in pathway $k$.
- $L_{jk}=1$ when variant $j$ remains after pathway-specific LD clumping.
- $\mathbf{1}(p_j\leq\tau)$ applies the selected GWAS p-value threshold.

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

Open `Pathway_PGS_AoU_Tutorial.ipynb` in JupyterLab and run the cells in order.

For the first demo:

1. Leave `RUN_MODE = "demo"`.
2. Run the notebook from top to bottom.
3. Review the explanation, code, and QC output at each step.
4. Inspect the final aggregate table, figures, and report.

The demo prepares and scores its two synthetic variants automatically:

```python
PREPARE_AOU_DATA = RUN_MODE == "demo"
RUN_SCORING = RUN_MODE == "demo"
```

Changing `RUN_MODE` to `"full"` turns both actions off. Enable them only after
reviewing the real files, cohort, mapping method, and generated command.

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

## Tutorial Roadmap

| Step | Scientific question | Expected output |
|---|---|---|
| 1. Choose evidence and participants | Which variant effects, thresholds, pathways, and people will be analyzed? | Reproducible settings |
| 2. Inspect weights and definitions | Do the GWAS columns, genes, and pathways mean what we expect? | Input QC and pathway preview |
| 3. Link variants to pathways | Which variants belong to each biological pathway? | Variant-to-pathway membership |
| 4. Match variants to AoU | Are the same variants and alleles available in the AoU genotypes? | Harmonized scoring input |
| 5. Control LD and calculate pPGS | Which variants remain, and how are participant scores calculated? | pPGS, aggregate QC, figures, and report |

The guided interface at the end repeats the same visible workflow through:
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

**A full analysis stops before preparation or scoring**

This is expected after changing `RUN_MODE` to `"full"`. Review the scientific
inputs and exact command before enabling each action.

## Repository Contents

- `Pathway_PGS_AoU_Tutorial.ipynb`: code-visible teaching notebook;
- `pathway_prs_core.py`: mapping, harmonization, command, QC, and plotting code;
- `pathway_pgs_app.py`: optional guided interface;
- `demo_data/`: synthetic demonstration inputs;
- `build_notebook.py`: notebook syntax, structure, and release-output validator;
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

Current release: `v2.5.0`.
