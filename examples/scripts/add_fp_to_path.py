import sys
import os

# ============================================================================
# FINANCEPY EXAMPLES - Add Fp To Path
# ============================================================================

financepy_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

if financepy_path not in sys.path:
    sys.path.append(financepy_path)

print("Adding", financepy_path, "to path")
