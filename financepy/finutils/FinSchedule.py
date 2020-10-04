##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from .FinError import FinError
from .FinDate import FinDate
from .FinCalendar import (FinCalendar, FinCalendarTypes)
from .FinCalendar import (FinBusDayAdjustTypes, FinDateGenRuleTypes)
from .FinFrequency import (FinFrequency, FinFrequencyTypes)
from .FinHelperFunctions import labelToString
from .FinHelperFunctions import checkArgumentTypes

###############################################################################


class FinSchedule(object):
    ''' A Schedule is a vector of dates generated according to ISDA standard
    rules which starts on the next date after the start date and runs up to
    an end date. Dates are adjusted to a provided calendar. The zeroth
    element is the previous coupon date (PCD) and the first element is the
    Next Coupon Date (NCD). '''

    def __init__(self,
                 startDate: FinDate,   # Also known as the effective date
                 endDate: FinDate,  # Also known as the termination date
                 frequencyType: FinFrequencyTypes = FinFrequencyTypes.ANNUAL,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD,
                 adjustTerminationDate:bool = False):
        ''' Create FinSchedule object which calculates a sequence of dates in
        line with market convention for fixed income products. '''

        checkArgumentTypes(self.__init__, locals())

        # validation complete
        self._startDate = startDate
        self._endDate = endDate
        self._frequencyType = frequencyType
        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType
        
        # TODO USE THIS TO ADJUST OR NOT
        self._adjustTerminationDate = adjustTerminationDate
        
        self._adjustedDates = None
        
        self._generate()

###############################################################################

    def scheduleDates(self):
        ''' Returns a list of the schedule of FinDates. '''

        if self._adjustedDates is None:
            raise FinError("Schedule dates have not been generated.")

        return self._adjustedDates

###############################################################################

    def _generate(self):
        ''' Generate schedule of dates according to specified date generation
        rules and also adjust these dates for holidays according to the
        specified business day convention and the specified calendar. '''

        # NEED TO INCORPORATE ADJUST TERMINATION DATE FLAG

        calendar = FinCalendar(self._calendarType)
        frequency = FinFrequency(self._frequencyType)
        numMonths = int(12 / frequency)

        unadjustedScheduleDates = []
        self._adjustedDates = []

        if self._dateGenRuleType == FinDateGenRuleTypes.BACKWARD:

            nextDate = self._endDate
            flowNum = 0

            while nextDate > self._startDate:
                unadjustedScheduleDates.append(nextDate)
                nextDate = nextDate.addMonths(-numMonths)
                flowNum += 1

            # Add on the Previous Coupon Date
            unadjustedScheduleDates.append(nextDate)
            flowNum += 1

            # reverse order and holiday adjust dates
            # the first date is not adjusted as this was provided
            dt = unadjustedScheduleDates[flowNum - 1]
            self._adjustedDates.append(dt)

            for i in range(1, flowNum):

                dt = calendar.adjust(unadjustedScheduleDates[flowNum - i - 1],
                                     self._busDayAdjustType)

                self._adjustedDates.append(dt)

        elif self._dateGenRuleType == FinDateGenRuleTypes.FORWARD:

            # This needs checking
            nextDate = self._startDate
            flowNum = 0

            unadjustedScheduleDates.append(nextDate)
            flowNum = 1

            while nextDate < self._endDate:
                unadjustedScheduleDates.append(nextDate)
                nextDate = nextDate.addMonths(numMonths)
                flowNum = flowNum + 1

            # The first date is not adjusted as it is given
            for i in range(1, flowNum):

                dt = calendar.adjust(unadjustedScheduleDates[i],
                                     self._busDayAdjustType)

                self._adjustedDates.append(dt)

            self._adjustedDates.append(self._endDate)

        return self._adjustedDates

###############################################################################

    def _generate_alternative(self):
        ''' This adjusts each date BEFORE generating the next date.
        Generate schedule of dates according to specified date generation
        rules and also adjust these dates for holidays according to
        the business day convention and the specified calendar. '''

        # print("======= SCHEDULE HAS CHANGED - MUST TEST =============")

        self._adjustedDates = []
        calendar = FinCalendar(self._calendarType)
        frequency = FinFrequency(self._frequencyType)
        numMonths = int(12 / frequency)

        unadjustedScheduleDates = []

        if self._dateGenRuleType == FinDateGenRuleTypes.BACKWARD:

            nextDate = self._endDate
            print("END:", nextDate)
            flowNum = 0

            while nextDate > self._startDate:
                unadjustedScheduleDates.append(nextDate)
                nextDate = nextDate.addMonths(-numMonths)
                nextDate = calendar.adjust(nextDate, self._busDayAdjustType)
                flowNum += 1

            # Add on the Previous Coupon Date
            unadjustedScheduleDates.append(nextDate)
            flowNum += 1

            # reverse order and holiday adjust dates
            for i in range(0, flowNum):
                dt = unadjustedScheduleDates[flowNum - i - 1]
                self._adjustedDates.append(dt)

        elif self._dateGenRuleType == FinDateGenRuleTypes.FORWARD:

            # This needs checking
            nextDate = self._startDate
            flowNum = 0

            unadjustedScheduleDates.append(nextDate)
            flowNum = 1

            while nextDate < self._endDate:
                unadjustedScheduleDates.append(nextDate)
                nextDate = nextDate.addMonths(numMonths)
                nextDate = calendar.adjust(nextDate, self._busDayAdjustType)
                flowNum = flowNum + 1

            for i in range(1, flowNum):

                dt = unadjustedScheduleDates[i]
                self._adjustedDates.append(dt)

            self._adjustedDates.append(self._endDate)

        return self._adjustedDates

##############################################################################

    def __repr__(self):
        ''' Print out the details of the schedule and the actual dates. This
        can be used for providing transparency on schedule calculations. '''

        s = labelToString("START DATE", self._startDate)
        s += labelToString("END DATE", self._endDate)
        s += labelToString("FREQUENCY", self._frequencyType)
        s += labelToString("CALENDAR", self._calendarType)
        s += labelToString("BUSDAYRULE", self._busDayAdjustType)
        s += labelToString("DATEGENRULE", self._dateGenRuleType, "")

        if len(self._adjustedDates) > 0:
            s += "\n\n"
            s += labelToString("PCD", self._adjustedDates[0], "")

        if len(self._adjustedDates) > 1:
            s += "\n"
            s += labelToString("NCD", self._adjustedDates[1:], "",
                               listFormat=True)

        return s

###############################################################################

    def _print(self):
        ''' Print out the details of the schedule and the actual dates. This
        can be used for providing transparency on schedule calculations. '''
        print(self)

###############################################################################
