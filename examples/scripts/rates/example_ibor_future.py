# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

from financepy.products.rates.ibor_future import IborFuture
from financepy.utils.date_format import set_date_format, DateFormatTypes
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - IborFuture
# ============================================================================


set_date_format(DateFormatTypes.UK_LONG)

########################################################################################




########################################################################################

# ============================================================================
# 1. FIN IBOR FUTURE
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. FIN IBOR FUTURE")
print("=" * 78)

today_date = Date(5, 5, 2020)

print("VALUES")

for i in range(1, 12):
    fut = IborFuture(today_date, i, "3M")
    print(fut)

    fra = fut.to_fra(0.020, 0.0)
    print(fra)

