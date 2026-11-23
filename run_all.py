"""Run the whole analysis from the raw files:   python run_all.py
Each step can also be run on its own from the scripts/ folder (python 01_expression.py, ...), in order."""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = ["01_expression.py", "02_mtdna_quantities.py", "03_mtdna_statistics.py",
         "04_kr_replicate_20180314.py", "05_summary.py", "06_supplementary_tables.py", "07_figures.py",
         "08_copy_number_vs_damage.py", "09_coe3_2014_vs_2017.py"]
for s in STEPS:
    print(f"\n========== {s} ==========", flush=True)
    r = subprocess.run([sys.executable, s], cwd=HERE / "scripts")
    if r.returncode:
        sys.exit(f"{s} failed (exit {r.returncode})")
print("\nAll steps finished. Outputs are in results/ (start with results/SUMMARY.md).")
