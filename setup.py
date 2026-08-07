from datetime import datetime
from pathlib import Path

from setuptools import setup


def generate_init():
    version = "1.1.0"
    build_dt = datetime.now().strftime("%d %b %Y at %H:%M")

    contents = f'''\
__version__ = "{version}"

CR = "\\n"
S = "##################################################################" + CR
S += "#   FINANCEPY Version {version} - This build: {build_dt}   #" + CR
S += "#   This software is distributed FREE AND WITHOUT ANY WARRANTY   #" + CR
S += "# Report bugs as issues at https://github.com/domokane/FinancePy #" + CR
S += "##################################################################"
S += CR

print(S)
'''

    Path("financepy/__init__.py").write_text(contents, encoding="utf-8")


if __name__ == "__main__":
    generate_init()
    setup()