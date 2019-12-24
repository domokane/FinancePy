# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 20:49:10 2016

@author: Dominic O'Kane
"""
import traceback


def func_name():
    import traceback
    return traceback.extract_stack(None, 2)[0][2]


def isNotEqual(x, y, tol=1e-6):
    if abs(x - y) > tol:
        return True
    return False


class FinError(Exception):
    ''' Simple error class specific to FinPy. Need to decide how to handle
    FinancePy errors. Work in progress. '''

    def __init__(self, message):
        ''' Create FinError object by passing a message string. '''
        self._message = message

    def print(self):
        print("FinError:", self._message)
