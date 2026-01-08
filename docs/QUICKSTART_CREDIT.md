## Quickstart Guide to Credit Derivatives
We wish to set up and value a CDS contract.

To do this we need the Date and other utility classes such as Frequency and DayCounts. The easiest way to do this is to use a wildcard import.
```
from financepy.utils import *
```
As we will need interest rates to build the CDS curve we load in the rates modules
```
from financepy.products.rates import *
```
And of course we need the credit modules
```
from financepy.products.credit import *
```
In each case we load everything. As we become more expert in the library we can learn to just import exactly what we need.

We now define our CDS contract. It was traded first on 1st Dec 2024. We value it today 10th Jan 2026, and it matures on the 20th March 2029. The contractual coupon is set at 100 basis points. And we are long protection with a notional of $1,000,000. We therefore define
```
effective_dt = Date(1, 12, 2024)
maturity_dt = Date(20, 3, 2028)
cds_coupon = 0.010
notional = ONE_MILLION
long_protection = True
```
We can then define the CDS contract as follows
```
cds_contract = CDS(value_dt, maturity_dt, cds_coupon, notional, long_protection)
```
We can print the contract to ensure that it is as we wish
```
***Output***
OBJECT TYPE: CDS
STEP-IN DATE: 01-DEC-2024
MATURITY: 20-MAR-2028
NOTIONAL: 1000000
LONG PROTECTION: True
RUN COUPON: 100.0bp
DAYCOUNT: DayCountTypes.ACT_360
FREQUENCY: FrequencyTypes.QUARTERLY
CALENDAR: CalendarTypes.WEEKEND
BUSDAYRULE: BusDayAdjustTypes.FOLLOWING
DATEGENRULE: DateGenRuleTypes.BACKWARD
PAYMENT_dt, YEAR_FRAC, ACCRUAL_START, ACCRUAL_END, FLOW
20-DEC-2024,     0.252778, 20-SEP-2024, 19-DEC-2024,  2527.777778
20-MAR-2025,     0.250000, 20-DEC-2024, 19-MAR-2025,  2500.000000
20-JUN-2025,     0.255556, 20-MAR-2025, 19-JUN-2025,  2555.555556
22-SEP-2025,     0.261111, 20-JUN-2025, 21-SEP-2025,  2611.111111
22-DEC-2025,     0.252778, 22-SEP-2025, 21-DEC-2025,  2527.777778
20-MAR-2026,     0.244444, 22-DEC-2025, 19-MAR-2026,  2444.444444
22-JUN-2026,     0.261111, 20-MAR-2026, 21-JUN-2026,  2611.111111
21-SEP-2026,     0.252778, 22-JUN-2026, 20-SEP-2026,  2527.777778
21-DEC-2026,     0.252778, 21-SEP-2026, 20-DEC-2026,  2527.777778
22-MAR-2027,     0.252778, 21-DEC-2026, 21-MAR-2027,  2527.777778
21-JUN-2027,     0.252778, 22-MAR-2027, 20-JUN-2027,  2527.777778
20-SEP-2027,     0.252778, 21-JUN-2027, 19-SEP-2027,  2527.777778
20-DEC-2027,     0.252778, 20-SEP-2027, 19-DEC-2027,  2527.777778
20-MAR-2028,     0.255556, 20-DEC-2027, 20-MAR-2028,  2555.555556
```
We see the details of the contract that would be found on its term sheet. The premium leg starts paying in December 2024.

To value this CDS, we first need to load the rates curve. This is usually based on IBOR instruments, specifically Deposits, Futures and Swaps. However to keep things simple for the quick start, we will assume a flat discount curve. We therefore define
```
interest_rate = 0.03
```
We also pull in the required curve type
```
from financepy.market.curves import DiscountCurveFlat
```
We then construct this flat curve
```
discount_curve = DiscountCurveFlat(value_dt, interest_rate, FrequencyTypes.CONTINUOUS)
```
Next we need to construct the CDS curve. For this we need to specify a term structure of CDS contracts with their corresponding par CDS spread. We do not need to use CDS maturity dates - we can simply specify their tenor i.e. how many years they each last. We first create the settlement date which is value date plus one
```
settle_dt = value_dt.add_days(1)
```
and then we specify the 1Y, 2Y, 3Y and 5Y CDS contracts with CDS par spreads of 80, 85, 90, 95 bps
```
cds1 = CDS(settle_dt, "1Y", 0.0065)
cds2 = CDS(settle_dt, "2Y", 0.0070)
cds3 = CDS(settle_dt, "3Y", 0.0075)
cds5 = CDS(settle_dt, "5Y", 0.0080)
```
We store these in a list
```
cdss = [cds1, cds2, cds3, cds5]
```
We also need to specify a recovery rate (market standard is to assume 40%)
```
recovery_rate = 0.40
```
and we create our CDS curve
```
cds_curve = CDSCurve(value_dt, cds_list, discount_curve, recovery_rate)
```
We can print this to see what is stored
```
print(cds_curve)
```
***Output:***
```
OBJECT TYPE: CDSCurve
TIME,SURVIVAL_PROBABILITY
 0.0000000,  1.0000000
 1.2136986,  0.9865973
 2.2164384,  0.9739514
 3.2164384,  0.9597968
 5.2164384,  0.9314422
 ```
 We see the term structure of years and risk-neutral survival probabilities.

 We can now value the CDS contract. First we calculate the par CDS spread for this maturity. Recall that our CDS contract protection costs us 100bp a year.

```
spd = cds_contract.par_spread(value_dt, cds_curve, recovery_rate) * 10000.0
print("FAIR CDS SPREAD %10.5f bp"% spd)
```
Which gives
```
FAIR CDS SPREAD   69.999 bp
```
This makes sense. It agrees with the current market par spread of 70bps for a 2Y contract. Note that a 2Y contract matures in 2 years on the next IMM date which is the 20 March 2028.

The value of the contract is given by the difference between the par spread of 70bps and the 100bps that we are paying for the 2.25 years of protection remaining. It should be negative and equal to 30bps times 2.25 years. Let us see.
```
v = cds_contract.value(settle_dt, cds_curve, recovery_rate)
```
This has two components:
```
dirty_pv = v['dirty_pv']
clean_pv = v['clean_pv']
```
We print these
```
print("DIRTY VALUE %12.2f"% dirty_pv)
print("CLEAN VALUE %12.2f"% clean_pv)
```
which gives
```
DIRTY VALUE      -6752.97
CLEAN VALUE      -6475.19
```
We see that the dirty price (that the contract is worth) is approximately equal to -30bps x 2.25 x $1m = -$6,750.

We can get this as a clean price
```
cleanp = cds_contract.clean_price(settle_dt, cds_curve, recovery_rate)
print("CLEAN PRICE %12.6f"% cleanp)
```
which gives
```
CLEAN PRICE   100.646919
```
