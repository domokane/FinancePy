from pathlib import Path
import re
from setuptools import setup, find_packages

# Paths
this_dir = Path(__file__).parent

# Long description from README
long_description = (this_dir / "README.md").read_text(encoding="utf-8")

# Read version from financepy/__init__.py (single source of truth)
init_text = (this_dir / "financepy" / "__init__.py").read_text(encoding="utf-8")
m = re.search(r'^__version__\s*=\s*[\'"]([^\'"]+)[\'"]', init_text, re.MULTILINE)
if not m:
    raise RuntimeError("Unable to find __version__ in financepy/__init__.py")
version = m.group(1)

setup(
    name="financepy",
    version=version,
    author="Dominic O'Kane",
    description="A Finance Securities Valuation Library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/domokane/FinancePy",
    license="GPLv3",
    keywords=["FINANCE", "OPTIONS", "BONDS", "VALUATION", "DERIVATIVES"],
    packages=find_packages(exclude=("tests", "docs")),
    include_package_data=True,  # <-- fixed typo
    package_data={"": ["*.npz"]},  # or be explicit: {"financepy": ["data/*.npz"]}
    python_requires=">=3.8",  # consider modern floor; 3.6 is EOL
    install_requires=[
        "numpy",
        "numba",
        "scipy",
        "llvmlite",  # keep in sync with numba versions if you pin later
        "ipython",  # consider moving to extras if not required at runtime
        "matplotlib",
        "pandas",
    ],
    extras_require={
        "dev": ["pytest", "build", "twine"],
    },
    entry_points={
        "console_scripts": [
            # example: "financepy=financepy.__main__:main"
        ]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    project_urls={
        "Source": "https://github.com/domokane/FinancePy",
        "Tracker": "https://github.com/domokane/FinancePy/issues",
    },
)
