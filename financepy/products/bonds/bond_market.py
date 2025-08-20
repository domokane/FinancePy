##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO  Add Japan


from enum import Enum

from ...utils.frequency import FrequencyTypes
from ...utils.day_count import DayCountTypes
from ...utils.calendar import CalendarTypes


class BondMarkets(Enum):
    AUSTRIA = 1
    BELGIUM = 2
    CYPRUS = 3
    ESTONIA = 4
    FINLAND = 5
    FRANCE = 6
    GERMANY = 7
    GREECE = 8
    IRELAND = 9
    ITALY = 10
    LATVIA = 11
    LITHUANIA = 12
    LUXEMBOURG = 13
    MALTA = 14
    NETHERLANDS = 15
    PORTUGAL = 16
    SLOVAKIA = 17
    SLOVENIA = 18
    SPAIN = 19
    ESM = 20
    EFSF = 21
    BULGARIA = 22
    CROATIA = 23
    CZECH_REPUBLIC = 24
    DENMARK = 25
    HUNGARY = 26
    POLAND = 27
    ROMANIA = 28
    SWEDEN = 29
    JAPAN = 30
    SWITZERLAND = 31
    UNITED_KINGDOM = 32
    UNITED_STATES = 33
    AUSTRALIA = 34
    NEW_ZEALAND = 35
    NORWAY = 36
    SOUTH_AFRICA = 37


########################################################################################


def get_bond_market_conventions(country):
    """Returns the day count convention for accrued interest, the frequency
    and the number of days from trade date to settlement date.
    This is for Treasury markets. And for secondary bond markets."""

    annual = FrequencyTypes.ANNUAL
    semi_annual = FrequencyTypes.SEMI_ANNUAL
    act_act = DayCountTypes.ACT_ACT_ICMA
    thirtye360 = DayCountTypes.THIRTY_E_360
    thirty360 = DayCountTypes.THIRTY_360_BOND

    # TODO: CHECK CONVENTIONS
    # RETURNS
    # ACCRUAL CONVENTION
    # COUPON FREQUENCY
    # SETTLEMENT DAYS
    # NUM EX DIVIDEND DAYS AND CALENDAR TO USE

    if country == BondMarkets.AUSTRIA:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.AUSTRALIA:
        return (act_act, annual, 2, 7, CalendarTypes.NONE)
    if country == BondMarkets.BELGIUM:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.CYPRUS:
        return (act_act, semi_annual, 2, 0, None)
    if country == BondMarkets.ESTONIA:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.FINLAND:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.FRANCE:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.GERMANY:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.GREECE:
        return (act_act, annual, 3, 0, None)
    if country == BondMarkets.IRELAND:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.ITALY:
        return (act_act, semi_annual, 2, 0, None)
    if country == BondMarkets.LATVIA:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.LITHUANIA:
        return (act_act, annual, 1, 0, None)
    if country == BondMarkets.LUXEMBOURG:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.MALTA:
        return (act_act, semi_annual, 2, 0, None)
    if country == BondMarkets.NETHERLANDS:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.NEW_ZEALAND:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.NORWAY:
        return (thirty360, annual, 2, 2, None)
    if country == BondMarkets.PORTUGAL:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.SLOVAKIA:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.SLOVENIA:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.SPAIN:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.ESM:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.EFSF:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.BULGARIA:
        return (act_act, semi_annual, 0, 0, None)
    if country == BondMarkets.CROATIA:
        return (act_act, semi_annual, 3, 0, None)
    if country == BondMarkets.CZECH_REPUBLIC:
        return (act_act, semi_annual, 2, 0, None)
    if country == BondMarkets.DENMARK:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.HUNGARY:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.POLAND:
        return (act_act, semi_annual, 2, 0, None)
    if country == BondMarkets.ROMANIA:
        return (act_act, semi_annual, 2, 0, None)
    if country == BondMarkets.SOUTH_AFRICA:
        return (act_act, annual, 2, 10, CalendarTypes.NONE)  # CHECK
    if country == BondMarkets.SWEDEN:
        return (thirtye360, annual, 2, 0, None)
    if country == BondMarkets.JAPAN:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.SWITZERLAND:
        return (act_act, annual, 2, 0, None)
    if country == BondMarkets.UNITED_STATES:
        return (act_act, semi_annual, 2, 0, None)
    if country == BondMarkets.UNITED_KINGDOM:
        # OR 7 DAYS ?
        return (act_act, semi_annual, 1, 6, CalendarTypes.UNITED_KINGDOM)

    print("Unknown Country:", country)
    return (None, None, None, None, None)
