from pathlib import Path
import re
from setuptools import setup, find_packages

ROOT = Path(__file__).parent
long_description = (ROOT / "README.md").read_text(encoding="utf-8")

# --- Single source of truth: read __version__ from financepy/__init__.py ---
init_text = (ROOT / "financepy" / "__init__.py").read_text(encoding="utf-8")
m = re.search(r'^__version__\s*=\s*[\'"]([^\'"]+)[\'"]', init_text, flags=re.M)
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
    license="GPL-3.0-or-later",  # SPDX identifier
    keywords=["FINANCE", "OPTIONS", "BONDS", "VALUATION", "DERIVATIVES"],
    packages=find_packages(exclude=("tests*", "docs*", "notebooks*", "examples*")),
    include_package_data=True,
    # Add curated data files (e.g. under financepy/data/)
    package_data={"financepy": ["data/*.npz"]},
    python_requires=">=3.8",
    install_requires=[
        "numpy",
        "scipy",
        "numba",  # you can add bounds like "numba>=0.59,<0.61"
        "llvmlite",  # keep aligned with numba version
    ],
    extras_require={
        "viz": ["matplotlib", "pandas"],
        "dev": ["pytest", "build", "twine", "ipython"],
    },
    entry_points={
        "console_scripts": [
            # e.g. "financepy=financepy.__main__:main",
        ]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        # If you keep `license="GPL-3.0-or-later"` above,
        # you don’t need a Trove license classifier here.
    ],
    project_urls={
        "Source": "https://github.com/domokane/FinancePy",
        "Tracker": "https://github.com/domokane/FinancePy/issues",
    },
)
