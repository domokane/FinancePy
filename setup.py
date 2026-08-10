from datetime import datetime
from pathlib import Path
import sys
import tomllib

from setuptools import setup


ROOT = Path(__file__).resolve().parent


def generate_init():

    with open(ROOT / "pyproject.toml", "rb") as f:
        pyproject = tomllib.load(f)

    version = pyproject["project"]["version"]

    build_dt = datetime.now().strftime("%d %b %Y at %H:%M")

    contents = f'''\
__version__ = "{version}"

CR = "\\n"
S =  "#############################################################" + CR
S += "#  FINANCEPY Version {version} - Built on {build_dt}  #" + CR
S += "#  This software is distributed FREE AND WITHOUT WARRANTY   #" + CR
S += "#  Report issues at https://github.com/domokane/FinancePy   #" + CR
S += "#############################################################"
S += CR

print(S)
'''

    init_file = ROOT / "financepy" / "__init__.py"
    init_file.write_text(contents, encoding="utf-8")

    print(f"Generated {init_file}")
    print(f"FinancePy version: {version}")
    print(f"Build date/time:   {build_dt}")


if __name__ == "__main__":

    generate_init()

    # Only invoke setuptools when a command has been supplied.
    if len(sys.argv) > 1:
        setup()
    else:
        print("No setuptools command supplied - generated __init__.py only.")