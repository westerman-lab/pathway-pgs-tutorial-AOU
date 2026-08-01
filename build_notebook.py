#!/usr/bin/env python3
"""Validate and normalize the canonical tutorial notebook for release."""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK = ROOT / "Pathway_PGS_AoU_Tutorial.ipynb"
REQUIRED_HEADINGS = [
    "# From Variant Effects to Pathway PGS in All of Us",
    "## Step 1: Choose the variant evidence",
    "## Step 2: Inspect the variant weights",
    "## Step 3: Link variants to genes and pathways",
    "## Step 4: Match pathway variants to All of Us genotypes",
    "## Step 5: Control LD and calculate pathway scores",
    "## The guided interface",
]


def main() -> None:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    if notebook.get("nbformat") != 4:
        raise ValueError("The tutorial must use notebook format 4.")

    clean_cells = []
    for cell in notebook["cells"]:
        source = cell.get("source", "")
        source_text = "".join(source) if isinstance(source, list) else source

        if cell["cell_type"] == "code" and not source_text.strip():
            continue

        if cell["cell_type"] == "code":
            compile(source_text, f"notebook-cell-{len(clean_cells) + 1}", "exec")
            cell["execution_count"] = None
            cell["outputs"] = []

        clean_cells.append(cell)

    combined_source = "\n".join(
        "".join(cell.get("source", []))
        if isinstance(cell.get("source", []), list)
        else cell.get("source", "")
        for cell in clean_cells
    )
    missing = [heading for heading in REQUIRED_HEADINGS if heading not in combined_source]
    if missing:
        raise ValueError(f"Missing required tutorial headings: {missing}")

    notebook["cells"] = clean_cells
    notebook.setdefault("metadata", {}).pop("widgets", None)
    notebook["metadata"]["tutorial_version"] = "2.5.0"
    for index, cell in enumerate(clean_cells, start=1):
        cell["id"] = f"pathway-pgs-{index:02d}"

    NOTEBOOK.write_text(
        json.dumps(notebook, ensure_ascii=False, indent=1) + "\n",
        encoding="utf-8",
    )
    print(f"Validated and normalized {NOTEBOOK.name} ({len(clean_cells)} cells)")


if __name__ == "__main__":
    main()
