# Changelog

## 2.5.0

- Reorganized the tutorial into five biological and statistical pPGS steps.
- Added beginner-focused explanations of variant weights, pathway membership,
  allele harmonization, LD clumping, p-value thresholds, and score interpretation.
- Added concise QC tables and a demo-to-full-analysis replacement guide.
- Added a completion message and a clearly separated guided interface.
- Removed saved AoU runtime output and widget state from the public notebook.

## 2.4.0

- Added beginner-focused explanations of pPGS inputs, decisions, and outputs.
- Added a live summary of the selected mapping, pathway, LD, and threshold settings.
- Renamed workflow actions to describe the underlying analytical step.
- Added safe top-to-bottom notebook execution and synthetic end-to-end scoring tests.
- Reorganized the README around a first demo, scientific analysis, QC, and troubleshooting.

## 2.3.0

- Reframed the project as an All of Us-specific, code-first tutorial.
- Split the notebook into seven visible analysis sections.
- Exposed scientific choices, preparation code, the generated command, and QC.
- Disabled expensive preparation and scoring by default.
- Retained the interactive interface as an optional helper at the end.
- Added notebook-structure tests and fixed empty scoring-engine path detection.
