## Quickstart Guide to Rate Derivatives
To analyse a bond first load up the Date and other utility classes such as Frequency and DayCounts. The easiest way to do this is to use a wildcard import. The downside is that it imports a lot of unnecessary classes and will take a second or so. But it is probably the simplest way to start until you become familiar with the library structure.

```
from financepy.utils import *
```

Then load up all of the bond classes from the bond products module. Once again you will load classes you will not use but this is the simplest way to begin.
```
from financepy.products.bonds import *
```

Now use the Date class to define the bond issue and maturity dates
```
issue_dt = Date(13, 5, 2010)
maturity_dt = Date(13, 5, 2022)
```
Now set the 2.7% coupon - this is the annualised coupon which is paid in discrete chunks depending on the payment frequency
```
coupon = 0.027
```
Here the payment frequency for the coupon is SEMI ANNUAL so there are two payments of 0.0135 paid every six months
```
freq_type = FrequencyTypes.SEMI_ANNUAL
```
The Day Count convention used for accrued interest is defined as
```
dc_type = DayCountTypes.THIRTY_E_360
```
Finally create the bond object
```
bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)
```
We can print the bond object to see what we stored
```
print(bond)
```
**Output:**
```
OBJECT TYPE: Bond
ISSUE DATE: 13-MAY-2010
MATURITY DATE: 13-MAY-2022
COUPON (%): 2.7
FREQUENCY: FrequencyTypes.SEMI_ANNUAL
DAY COUNT TYPE: DayCountTypes.THIRTY_E_360
EX_DIV DAYS: 0
```
---

We then get the bond payments for a specific settlement date
```
settle_dt = Date(21, 7, 2017)
```
The function to print the bond payments is
```
bond.print_payments(settle_dt)
```
**Output:**
```
 13-NOV-2017      1.35000
 13-MAY-2018      1.35000
 13-NOV-2018      1.35000
 13-MAY-2019      1.35000
 13-NOV-2019      1.35000
 13-MAY-2020      1.35000
 13-NOV-2020      1.35000
 13-MAY-2021      1.35000
 13-NOV-2021      1.35000
 13-MAY-2022    101.35000
---
```
To get the accrued interest use
```
print("Accrued = %12.2f"% bond.accrued_int)
```
**Output:**
```
Accrued =         0.51
```
We can also calculate the number of accrued days
```
print("Accrued Days = ", bond.accrued_days)
```
**Output:**
```
Accrued Days =  68
```
Now we set the bond (clean) price
```
clean_price = 101.581564
```
From this we can extract the current yield
```
bond.current_yield(clean_price)*100.0
```
**Output:**
```
2.657962620067555
```
We can get the yield to maturity using US STREET convention
```
bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.US_STREET)*100
```
**Output:**
```
2.3499999794774293
```
We can also get the yield to maturity using the US TREASURY convention
```
bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.US_TREASURY)*100
```
**Output:**
```
2.3496418129692254
```

