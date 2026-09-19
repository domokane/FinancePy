from pathlib import Path
import sys

# Add the FinancePy repository root to Python's import path.
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

print("Adding", financepy_path, "to path")
