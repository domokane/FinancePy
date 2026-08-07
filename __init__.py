# financepy/__init__.py

__version__ = "1.1.0"

from build_info import __build_datetime__

CR = "\n"
S = "##################################################################" + CR
S += "# FINANCEPY Version " + str(__version__) + " - This build: " + str(__build_datetime__) + " # " + CR
S += "# This software is distributed FREE AND WITHOUT ANY WARRANTY     #" + CR
S += "# Report bugs as issues at https://github.com/domokane/FinancePy #" + CR
S += "##################################################################"
S += CR

print(S)