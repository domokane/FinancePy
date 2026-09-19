# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCount, DayCountTypes
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - DayCount
# ============================================================================



########################################################################################




########################################################################################

# ============================================================================
# 1. FIN DAY COUNT
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. FIN DAY COUNT")
print("=" * 78)

print("DAY_COUNT_METHOD", "START", "END", "ALPHA")

fin_freq = FrequencyTypes.ANNUAL

for day_count_method in DayCountTypes:

    start_dt = Date(1, 1, 2019)
    next_dt = start_dt
    num_days = 20
    day_count = DayCount(day_count_method)

    for _ in range(0, num_days):
        next_dt = next_dt.add_days(7)
        dcf = day_count.year_frac(start_dt, next_dt, next_dt, fin_freq)

        print(str(day_count_method), str(start_dt), str(next_dt), dcf[0])

