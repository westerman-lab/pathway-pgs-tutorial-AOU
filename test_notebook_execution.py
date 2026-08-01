import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK = ROOT / "Pathway_PGS_AoU_Tutorial.ipynb"


def test_notebook_executes_safely_from_top_to_bottom(tmp_path: Path):
    """Run every code cell while expensive AoU operations remain disabled."""
    runner = r'''
import json
from pathlib import Path

notebook = json.loads(Path("Pathway_PGS_AoU_Tutorial.ipynb").read_text())
namespace = {"__name__": "__main__"}
for number, cell in enumerate(notebook["cells"], start=1):
    if cell["cell_type"] != "code":
        continue
    source = "".join(cell.get("source", []))
    source = source.replace(
        'PREPARE_AOU_DATA = RUN_MODE == "demo"',
        'PREPARE_AOU_DATA = False',
    )
    source = source.replace(
        'RUN_SCORING = RUN_MODE == "demo"',
        'RUN_SCORING = False',
    )
    exec(compile(source, f"notebook-cell-{number}", "exec"), namespace)

assert namespace["PREPARE_AOU_DATA"] is False
assert namespace["RUN_SCORING"] is False
assert namespace["app"].__class__.__name__ == "VBox"
print("SAFE_NOTEBOOK_EXECUTION_PASS")
'''
    environment = os.environ.copy()
    environment["HOME"] = str(tmp_path)
    environment["MPLBACKEND"] = "Agg"
    completed = subprocess.run(
        [sys.executable, "-c", runner],
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        timeout=90,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "SAFE_NOTEBOOK_EXECUTION_PASS" in completed.stdout
    assert not (tmp_path / "analysis").exists()
