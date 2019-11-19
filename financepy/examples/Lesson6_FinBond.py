###########################################################################

from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
from financepy.products.bonds.FinBond import FinBond, FinYieldConventions
from financepy.market.curves.FinFlatCurve import FinFlatCurve
from financepy.market.curves.FinFlatCurve import FinCompoundingMethods

###############################################################################

print("BLOOMBERG US TREASURY EXAMPLE")
settlementDate = FinDate(2017, 7, 21)
maturityDate = FinDate(2027, 5, 15)
coupon = 0.02375
freqType = FinFrequencyTypes.SEMI_ANNUAL
accrualType = FinDayCountTypes.ACT_ACT_ICMA

# By setting the face to 100 we expect a price of par to be 100.0
face = 100.0

bond = FinBond(maturityDate,
               coupon,
               freqType,
               accrualType,
               face)

cleanPrice = 99.780842  # if face is 1 then this must be 0.99780842

yld = bond.currentYield(cleanPrice)
print("Current Yield = ", yld)

ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                           FinYieldConventions.UK_DMO)
print("UK DMO Yield To Maturity = ", ytm)

ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                           FinYieldConventions.US_STREET)
print("US STREET Yield To Maturity = ", ytm)

ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                           FinYieldConventions.US_TREASURY)
print("US TREASURY Yield To Maturity = ", ytm)

fullPrice = bond.fullPriceFromYield(settlementDate, ytm)
print("Full Price = ", fullPrice)

cleanPrice = bond.cleanPriceFromYield(settlementDate, ytm)
print("Clean Price = ", cleanPrice)

accd = bond._accrued
print("Accrued = ", accd)

accddays = bond._accruedDays
print("Accrued Days = ", accddays)

duration = bond.dollarDuration(settlementDate, ytm)
print("Dollar Duration = ", duration)

modifiedDuration = bond.modifiedDuration(settlementDate, ytm)
print("Modified Duration = ", modifiedDuration)

macauleyDuration = bond.macauleyDuration(settlementDate, ytm)
print("Macauley Duration = ", macauleyDuration)

conv = bond.convexityFromYield(settlementDate, ytm)
print("Convexity = ", conv)

discountRate = 0.0275
flatCurve = FinFlatCurve(settlementDate,
                         discountRate,
                         FinCompoundingMethods.SEMI_ANNUAL)

asw = bond.assetSwapSpread(settlementDate, cleanPrice, flatCurve)
print("Discounted on LIBOR Curve ASW (bp):", asw * 10000)

oas = bond.optionAdjustedSpread(settlementDate, cleanPrice, flatCurve)
print("Discounted on LIBOR Curve OAS (bp):", oas * 10000)
