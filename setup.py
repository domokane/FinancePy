from datetime import datetime
from version import __version__
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("version.py", "r") as fh:
    version_number = fh.read()
    start = version_number.find("\"")
    end = version_number[start+1:].find("\"")
    version_number_str = str(version_number[start+1:start+end+1])
    version_number_str = version_number_str.replace('\n', '')

print(">>>" + version_number_str + "<<<")

###############################################################################
cr = "\n"

with open('financepy//__init__.template', 'r') as file:
    filedata = file.read()

# Replace the target string
filedata = filedata.replace('__version__', "'" + str(__version__) + "'")

now = datetime.now()
dt_string = now.strftime("%d %b %Y at %H:%M")

# Replace the target string
filedata = filedata.replace('__dateandtime__', dt_string)

# Write the file out again
with open('./financepy//__init__.py', 'w') as file:
    file.write(filedata)

###############################################################################

setuptools.setup(
    name="financepy",
    version=version_number_str,
    author="Dominic O'Kane",
    author_email="quant@financepy.com",
    description="A Finance Library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/domokane/FinancePy",
    keywords=['FINANCE', 'OPTIONS', 'BONDS', 'VALUATION', 'DERIVATIVES'],
    install_requires=['numpy', 'numba', 'scipy'],
    package_data={'': ['*.npz'], },
    include_package_date=True,
    packages=setuptools.find_packages(),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
