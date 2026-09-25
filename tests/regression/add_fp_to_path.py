import sys
import os

if 1 == 0:
    financepy_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if financepy_path not in sys.path:
        sys.path.append(financepy_path)

    print("Adding", financepy_path, "to path")
