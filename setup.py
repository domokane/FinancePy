import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="financepy",
    version="0.151",
    author="Dominic O'Kane",
    author_email="okane.dominic@gmail.com",
    description="A Finance Library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/domokane/FinancePy",
    keywords=['FINANCE','OPTIONS','BONDS','VALUATION','DERIVATIVES'],
    install_requires=['numpy','numba','scipy'],
    packages=setuptools.find_packages(),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)