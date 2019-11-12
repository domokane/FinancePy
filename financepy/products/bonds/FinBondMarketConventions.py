# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:16:14 2019

@author: Dominic
"""
#TODO Add Japan 

import enum as Enum

from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.products.bonds.FinBond import FinBondAccruedTypes


class FinCountry(Enum):
    AUSTRIA = 1,
    BELGIUM = 2,
    FINLAND = 3,
    FRANCE = 4,
    GERMANY = 5,
    IRELAND = 6,
    ITALY = 7,
    LUXEMBOURG = 8,
    NETHERLANDS = 9,
    PORTUGAL = 10,
    SPAIN = 11,
    UNITED_KINGDOM = 12,
    UNITED_STATES = 13


def getBondMarketConventions(country):
    ''' Returns the day count convention for accrued interest, the frequency
    and the number of days from trade date to settlement date. '''

    annual = FinFrequencyTypes.ANNUAL
    semi_annual = FinFrequencyTypes.SEMI_ANNUAL
    act_act = FinBondAccruedTypes.ACT_ACT_ICMA

    if country == FinCountry.AUSTRIA:
        return (act_act, annual, 3)
    elif country == FinCountry.BELGIUM:
        return (act_act, annual, 3)
    elif country == FinCountry.FINLAND:
        return (act_act, annual, 3)
    elif country == FinCountry.FRANCE:
        return (act_act, annual, 3)
    elif country == FinCountry.GERMANY:
        return (act_act, annual, 3)
    elif country == FinCountry.IRELAND:
        return (act_act, annual, 3)
    elif country == FinCountry.ITALY:
        return (act_act, semi_annual, 3)
    elif country == FinCountry.LUXEMBOURG:
        return (act_act, annual, 3)
    elif country == FinCountry.NETHERLANDS:
        return (act_act, annual, 3)
    elif country == FinCountry.PORTUGAL:
        return (act_act, annual, 3)
    elif country == FinCountry.SPAIN:
        return (act_act, annual, 3)
    elif country == FinCountry.UNITED_STATES:
        return (act_act, annual, 1)
    elif country == FinCountry.UNITED_KINGDOM:
        return (act_act, annual, 1)
