# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.



from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.schedule import Schedule
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - Date
# ============================================================================

########################################################################################




########################################################################################

# ============================================================================
# 1. DT ADJUST
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. DT ADJUST")
print("=" * 78)

start_dt = Date(28, 2, 2008)
end_dt = Date(28, 2, 2011)

freq_type = FrequencyTypes.SEMI_ANNUAL
cal_type = CalendarTypes.NONE
bd_type = BusDayAdjustTypes.FOLLOWING
dg_type = DateGenRuleTypes.BACKWARD

print("NO ADJUSTMENTS", "DATE")
schedule = Schedule(start_dt, end_dt, freq_type, cal_type, bd_type, dg_type)

for dt in schedule.adjusted_dts:
    print("Date:", dt)

print("")
print("NO WEEKENDS AND FOLLOWING", "DATE")
freq_type = FrequencyTypes.SEMI_ANNUAL
cal_type = CalendarTypes.WEEKEND
bd_type = BusDayAdjustTypes.FOLLOWING
dg_type = DateGenRuleTypes.BACKWARD

schedule = Schedule(start_dt, end_dt, freq_type, cal_type, bd_type, dg_type)

for dt in schedule.adjusted_dts:
    print("Date:", dt)

print("")
print("NO WEEKENDS AND MODIFIED FOLLOWING", "DATE")
freq_type = FrequencyTypes.SEMI_ANNUAL
cal_type = CalendarTypes.WEEKEND
bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
dg_type = DateGenRuleTypes.BACKWARD

schedule = Schedule(start_dt, end_dt, freq_type, cal_type, bd_type, dg_type)

for dt in schedule.adjusted_dts:
    print("Date:", dt)

print("")
print("NO WEEKENDS AND US HOLIDAYS AND MODIFIED FOLLOWING", "DATE")

freq_type = FrequencyTypes.SEMI_ANNUAL
cal_type = CalendarTypes.UNITED_STATES
bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
dg_type = DateGenRuleTypes.BACKWARD

start_dt = Date(4, 7, 2008)
end_dt = Date(4, 7, 2011)

schedule = Schedule(start_dt, end_dt, freq_type, cal_type, bd_type, dg_type)

for dt in schedule.adjusted_dts:
    print("Date:", dt)

