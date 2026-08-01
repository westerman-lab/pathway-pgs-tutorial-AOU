import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK = ROOT / "Pathway_PGS_AoU_Tutorial.ipynb"


def load_notebook():
    return json.loads(NOTEBOOK.read_text())


def test_notebook_is_code_first():
    notebook = load_notebook()
    code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]

    assert len(code_cells) >= 8
    assert all("hide-input" not in cell.get("metadata", {}).get("tags", []) for cell in code_cells)
    assert all(
        not cell.get("metadata", {}).get("jupyter", {}).get("source_hidden", False)
        for cell in code_cells
    )


def test_notebook_exposes_methods_and_safety_switches():
    notebook = load_notebook()
    source = "\n".join(
        "".join(cell.get("source", [])) for cell in notebook["cells"]
    )

    required_text = [
        "From Variant Effects to Pathway PGS in All of Us",
        "Step 1: Choose the variant evidence",
        "Step 2: Inspect the variant weights",
        "Step 3: Link variants to genes and pathways",
        "Step 4: Match pathway variants to All of Us genotypes",
        "Step 5: Control LD and calculate pathway scores",
        'PREPARE_AOU_DATA = RUN_MODE == "demo"',
        'RUN_SCORING = RUN_MODE == "demo"',
        "build_pathway_variant_union_bed",
        "harmonize_gwas_ids_to_bed",
        "build_pathway_pgs_command",
        "Exact scoring command (review before running)",
        "The guided interface",
        "Move from the demo to a complete analysis",
    ]
    for text in required_text:
        assert text in source


def test_notebook_does_not_hide_the_primary_workflow():
    notebook = load_notebook()
    code_source = "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )
    assert code_source.index("infer_gwas_schema") < code_source.index(
        "launch_pathway_pgs_app()"
    )


def test_public_notebook_contains_no_saved_runtime_output():
    notebook = load_notebook()
    code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]

    assert all(cell.get("execution_count") is None for cell in code_cells)
    assert all(not cell.get("outputs") for cell in code_cells)
    assert all("".join(cell.get("source", [])).strip() for cell in code_cells)
