##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

#TODO  Add Japan 


from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinError import FinError
from ...finutils.FinHelperFunctions import labelToString

from enum import Enum


class FinBondMarkets(Enum):
    AUSTRIA = 1,
    BELGIUM = 2,
    CYPRUS = 3,
    ESTONIA = 4,
    FINLAND = 5,
    FRANCE = 6,
    GERMANY = 7,
    GREECE = 8,
    IRELAND = 9,
    ITALY = 10,
    LATVIA = 11,
    LITHUANIA = 12,
    LUXEMBOURG = 13,
    MALTA = 14,
    NETHERLANDS = 15,
    PORTUGAL = 16,
    SLOVAKIA = 17,
    SLOVENIA = 18,
    SPAIN = 19,
    ESM = 20,
    EFSF = 21,
    BULGARIA = 22,
    CROATIA = 23,
    CZECH_REPUBLIC = 24,
    DENMARK = 25,
    HUNGARY = 26,
    POLAND = 27,
    ROMANIA = 28,
    SWEDEN = 29,
    JAPAN = 30,
    SWITZERLAND = 31,
    UNITED_KINGDOM = 32,
    UNITED_STATES = 33


def getTreasuryBondMarketConventions(country):
    ''' Returns the day count convention for accrued interest, the frequency
    and the number of days from trade date to settlement date.
    This is for Treasury markets. And for secondary bond markets. '''

    annual = FinFrequencyTypes.ANNUAL
    semi_annual = FinFrequencyTypes.SEMI_ANNUAL
    act_act = FinDayCountTypes.ACT_ACT_ICMA
    thirtye360 = FinDayCountTypes.THIRTY_E_360

    if country == FinBondMarkets.AUSTRIA:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.BELGIUM:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.CYPRUS:
        return (act_act, semi_annual, 2)
    elif country == FinBondMarkets.ESTONIA:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.FINLAND:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.FRANCE:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.GERMANY:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.GREECE:
        return (act_act, annual, 3)
    elif country == FinBondMarkets.IRELAND:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.ITALY:
        return (act_act, semi_annual, 2)
    elif country == FinBondMarkets.LATVIA:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.LITHUANIA:
        return (act_act, annual, 1)
    elif country == FinBondMarkets.LUXEMBOURG:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.MALTA:
        return (act_act, semi_annual, 2)
    elif country == FinBondMarkets.NETHERLANDS:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.PORTUGAL:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.SLOVAKIA:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.SLOVENIA:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.SPAIN:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.ESM:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.EFSF:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.BULGARIA:
        return (act_act, semi_annual, 0)
    elif country == FinBondMarkets.CROATIA:
        return (act_act, semi_annual, 3)
    elif country == FinBondMarkets.CZECH_REPUBLIC:
        return (act_act, semi_annual, 2)
    elif country == FinBondMarkets.DENMARK:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.HUNGARY:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.POLAND:
        return (act_act, semi_annual, 2)
    elif country == FinBondMarkets.ROMANIA:
        return (act_act, semi_annual, 2)
    elif country == FinBondMarkets.SWEDEN:
        return (thirtye360, annual, 2)
    elif country == FinBondMarkets.JAPAN:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.SWITZERLAND:
        return (act_act, annual, 2)
    elif country == FinBondMarkets.UNITED_STATES:
        return (act_act, semi_annual, 2)
    elif country == FinBondMarkets.UNITED_KINGDOM:
        return (act_act, semi_annual, 1)
    else:
        raise FinError("Unknown country.")
