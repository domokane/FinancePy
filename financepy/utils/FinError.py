##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

###############################################################################
# Suppress error traceback messages in Jupyter Notebook
###############################################################################

import traceback
import sys

from IPython import get_ipython

ipython = get_ipython()

###############################################################################

def _hide_traceback(exc_tuple=None, filename=None, tb_offset=None,
                    exception_only=False, running_compiled_code=False):
    etype, value, _ = sys.exc_info()

    msg = ipython._showtraceback(etype, value,
                                 ipython.InteractiveTB.get_exception_only(
                                     etype, value))

    return msg


def func_name():
    return traceback.extract_stack(None, 2)[0][2]


def suppressTraceback():
#    print(sys.tracebacklimit)
#    print(ipython.showtrackeback)

    sys.tracebacklimit = 0
    ipython.showtraceback = _hide_traceback

###############################################################################


class FinError(Exception):
    ''' Simple error class specific to FinPy. Need to decide how to handle
    FinancePy errors. Work in progress. '''

    def __init__(self,
                 message: str):
        ''' Create FinError object by passing a message string. '''
        self._message = message

    def _print(self):
        print("FinError:", self._message)

###############################################################################
