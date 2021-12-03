##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.products.inflation.FinInflationBond import FinInflationBond
from financepy.products.bonds import YTMCalcType
from financepy.products.inflation.FinInflationIndexCurve import FinInflationIndexCurve
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinInflationBondBBG():

    ##########################################################################
    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
    # Look for CPI Bond example
    ##########################################################################

    testCases.banner("BLOOMBERG US TIPS EXAMPLE")
    settlement_date = Date(21, 7, 2017)
    issue_date = Date(15, 7, 2010)
    maturity_date = Date(15, 7, 2020)
    coupon = 0.0125
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0
    baseCPIValue = 218.08532

    bond = FinInflationBond(issue_date,
                            maturity_date,
                            coupon,
                            freq_type,
                            accrual_type,
                            face,
                            baseCPIValue)

    testCases.header("FIELD", "VALUE")
    clean_price = 104.03502

    yld = bond.current_yield(clean_price)
    testCases.print("Current Yield = ", yld)

    ###########################################################################
    # Inherited functions that just calculate real yield without CPI adjustments
    ###########################################################################

    ytm = bond.yield_to_maturity(settlement_date,
                                 clean_price,
                                 YTMCalcType.UK_DMO)

    testCases.print("UK DMO REAL Yield To Maturity = ", ytm)

    ytm = bond.yield_to_maturity(settlement_date,
                                 clean_price,
                                 YTMCalcType.US_STREET)

    testCases.print("US STREET REAL Yield To Maturity = ", ytm)

    ytm = bond.yield_to_maturity(settlement_date,
                                 clean_price,
                                 YTMCalcType.US_TREASURY)

    testCases.print("US TREASURY REAL Yield To Maturity = ", ytm)

    full_price = bond.full_price_from_ytm(settlement_date, ytm)
    testCases.print("Full Price from REAL YTM = ", full_price)

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    testCases.print("Clean Price from Real YTM = ", clean_price)

    accddays = bond._accrued_days
    testCases.print("Accrued Days = ", accddays)

    accrued_interest = bond._accrued_interest
    testCases.print("REAL Accrued Interest = ", accrued_interest)

    ###########################################################################
    # Inflation functions that calculate nominal yield with CPI adjustment
    ###########################################################################

    refCPIValue = 244.65884

    ###########################################################################

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    testCases.print("Clean Price from Real YTM = ", clean_price)

    inflationAccd = bond.calc_inflation_accrued_interest(settlement_date,
                                                         refCPIValue)

    testCases.print("Inflation Accrued = ", inflationAccd)

    lastCpnCPIValue = 244.61839

    clean_price = bond.flat_price_from_yield_to_maturity(settlement_date, ytm,
                                                         lastCpnCPIValue,
                                                         YTMCalcType.US_TREASURY)

    testCases.print("Flat Price from Real YTM = ", clean_price)

    principal = bond.inflation_principal(settlement_date,
                                         ytm,
                                         refCPIValue,
                                         YTMCalcType.US_TREASURY)

    testCases.print("Inflation Principal = ", principal)

    ###########################################################################

    duration = bond.dollar_duration(settlement_date, ytm)
    testCases.print("Dollar Duration = ", duration)

    modified_duration = bond.modified_duration(settlement_date, ytm)
    testCases.print("Modified Duration = ", modified_duration)

    macauley_duration = bond.macauley_duration(settlement_date, ytm)
    testCases.print("Macauley Duration = ", macauley_duration)

    conv = bond.convexity_from_ytm(settlement_date, ytm)
    testCases.print("Convexity = ", conv)


###############################################################################
###############################################################################

def test_FinInflationBondStack():

    ##########################################################################
    # https://stackoverflow.com/questions/57676724/failing-to-obtain-correct-accrued-interest-with-quantlib-inflation-bond-pricer-i
    ##########################################################################

    testCases.banner("=============================")
    testCases.banner("QUANT FINANCE US TIPS EXAMPLE")
    testCases.banner("=============================")
    settlement_date = Date(23, 8, 2019)
    issue_date = Date(25, 9, 2013)
    maturity_date = Date(22, 3, 2068)
    coupon = 0.00125
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0
    baseCPIValue = 249.70

    ###########################################################################
    # Discount curve
    discount_curve = DiscountCurveFlat(settlement_date,
                                       0.01033692,
                                       FrequencyTypes.ANNUAL,
                                       DayCountTypes.ACT_ACT_ISDA)

    lag = 3
    fixingCPI = 244.65884
    fixingDate = settlement_date.add_months(-lag)

    ###########################################################################
    # Create Index Curve
    months = range(0, 12, 1)
    fixingDates = Date(31, 8, 2018).add_months(months)
    fixingRates = [284.2, 284.1, 284.5, 284.6, 285.6, 283.0, 285.0,
                   285.1, 288.2, 289.2, 289.6, 289.5]
    inflationIndex = FinInflationIndexCurve(fixingDates, fixingRates, lag)
#    print(inflationIndex)
    ###########################################################################

    zciisData = [(Date(31, 7, 2020), 3.1500000000137085),
                 (Date(31, 7, 2021), 3.547500000013759),
                 (Date(31, 7, 2022), 3.675000000013573),
                 (Date(31, 7, 2023), 3.7250000000134342),
                 (Date(31, 7, 2024), 3.750000000013265),
                 (Date(31, 7, 2025), 3.7430000000129526),
                 (Date(31, 7, 2026), 3.741200000012679),
                 (Date(31, 7, 2027), 3.7337000000123632),
                 (Date(31, 7, 2028), 3.725000000011902),
                 (Date(31, 7, 2029), 3.720000000011603),
                 (Date(31, 7, 2030), 3.712517289063011),
                 (Date(31, 7, 2031), 3.7013000000108764),
                 (Date(31, 7, 2032), 3.686986039205209),
                 (Date(31, 7, 2033), 3.671102614032895),
                 (Date(31, 7, 2034), 3.655000000009778),
                 (Date(31, 7, 2035), 3.6394715951305834),
                 (Date(31, 7, 2036), 3.624362044800966),
                 (Date(31, 7, 2037), 3.6093619727979087),
                 (Date(31, 7, 2038), 3.59421438364369),
                 (Date(31, 7, 2039), 3.5787000000081948),
                 (Date(31, 7, 2040), 3.5626192748395624),
                 (Date(31, 7, 2041), 3.545765016376823),
                 (Date(31, 7, 2042), 3.527943521613608),
                 (Date(31, 7, 2043), 3.508977137925462),
                 (Date(31, 7, 2044), 3.48870000000685),
                 (Date(31, 7, 2045), 3.467083068721011),
                 (Date(31, 7, 2046), 3.4445738220594935),
                 (Date(31, 7, 2047), 3.4216470902302065),
                 (Date(31, 7, 2048), 3.3986861494999188),
                 (Date(31, 7, 2049), 3.376000000005752),
                 (Date(31, 7, 2050), 3.3538412080641233),
                 (Date(31, 7, 2051), 3.3324275806807746),
                 (Date(31, 7, 2052), 3.311938788306623),
                 (Date(31, 7, 2053), 3.2925208131865835),
                 (Date(31, 7, 2054), 3.274293040759302),
                 (Date(31, 7, 2055), 3.2573541974782794),
                 (Date(31, 7, 2056), 3.241787355503245),
                 (Date(31, 7, 2057), 3.227664186159851),
                 (Date(31, 7, 2058), 3.2150486140060774),
                 (Date(31, 7, 2059), 3.204000000004159),
                 (Date(31, 7, 2060), 3.1945334946674064),
                 (Date(31, 7, 2061), 3.1865047145143377),
                 (Date(31, 7, 2062), 3.179753073456304),
                 (Date(31, 7, 2063), 3.1741427790361154),
                 (Date(31, 7, 2064), 3.1695593261025223),
                 (Date(31, 7, 2065), 3.1659065919088736),
                 (Date(31, 7, 2066), 3.163104428386987),
                 (Date(31, 7, 2067), 3.1610866681252903),
                 (Date(31, 7, 2068), 3.1597994770515836),
                 (Date(31, 7, 2069), 3.159200000003204),
                 (Date(31, 7, 2070), 3.159242349440139),
                 (Date(31, 7, 2071), 3.1598400898057433),
                 (Date(31, 7, 2072), 3.16090721831932),
                 (Date(31, 7, 2073), 3.162369676612098),
                 (Date(31, 7, 2074), 3.1641636543027207)]

    zcDates = []
    zcRates = []
    for i in range(0, len(zciisData)):
        zcDates.append(zciisData[i][0])
        zcRates.append(zciisData[i][1]/100.0)

    inflationZeroCurve = DiscountCurveZeros(settlement_date,
                                            zcDates,
                                            zcRates,
                                            FrequencyTypes.ANNUAL,
                                            DayCountTypes.ACT_ACT_ISDA)

#    print(inflationZeroCurve)

    ###########################################################################

    bond = FinInflationBond(issue_date,
                            maturity_date,
                            coupon,
                            freq_type,
                            accrual_type,
                            face,
                            baseCPIValue)

    testCases.header("FIELD", "VALUE")
    clean_price = 104.03502

    yld = bond.current_yield(clean_price)
    testCases.print("Current Yield = ", yld)

    return

    ###########################################################################
    # Inherited functions that just calculate real yield without CPI adjustments
    ###########################################################################

    ytm = bond.yield_to_maturity(settlement_date,
                                 clean_price,
                                 YTMCalcType.UK_DMO)

    testCases.print("UK DMO REAL Yield To Maturity = ", ytm)

    ytm = bond.yield_to_maturity(settlement_date,
                                 clean_price,
                                 YTMCalcType.US_STREET)

    testCases.print("US STREET REAL Yield To Maturity = ", ytm)

    ytm = bond.yield_to_maturity(settlement_date,
                                 clean_price,
                                 YTMCalcType.US_TREASURY)

    testCases.print("US TREASURY REAL Yield To Maturity = ", ytm)

    full_price = bond.full_price_from_ytm(settlement_date, ytm)
    testCases.print("Full Price from REAL YTM = ", full_price)

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    testCases.print("Clean Price from Real YTM = ", clean_price)

    accddays = bond._accrued_days
    testCases.print("Accrued Days = ", accddays)

    accrued_interest = bond._accrued_interest
    testCases.print("REAL Accrued Interest = ", accrued_interest)

    ###########################################################################
    # Inflation functions that calculate nominal yield with CPI adjustment
    ###########################################################################

    ###########################################################################

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    testCases.print("Clean Price from Real YTM = ", clean_price)

    inflationAccd = bond.calc_inflation_accrued_interest(settlement_date,
                                                         refCPIValue)

    testCases.print("Inflation Accrued = ", inflationAccd)

    lastCpnCPIValue = 244.61839

    clean_price = bond.flat_price_from_yield_to_maturity(settlement_date, ytm,
                                                         lastCpnCPIValue,
                                                         YTMCalcType.US_TREASURY)

    testCases.print("Flat Price from Real YTM = ", clean_price)

    principal = bond.inflation_principal(settlement_date,
                                         ytm,
                                         refCPIValue,
                                         YTMCalcType.US_TREASURY)

    testCases.print("Inflation Principal = ", principal)

    ###########################################################################

    duration = bond.dollar_duration(settlement_date, ytm)
    testCases.print("Dollar Duration = ", duration)

    modified_duration = bond.modified_duration(settlement_date, ytm)
    testCases.print("Modified Duration = ", modified_duration)

    macauley_duration = bond.macauley_duration(settlement_date, ytm)
    testCases.print("Macauley Duration = ", macauley_duration)

    conv = bond.convexity_from_ytm(settlement_date, ytm)
    testCases.print("Convexity = ", conv)

###############################################################################


test_FinInflationBondBBG()
test_FinInflationBondStack()
testCases.compareTestCases()
