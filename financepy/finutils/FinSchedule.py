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
    rules which starts on the next date after the effective date and runs up to
    a termination date. Dates are adjusted to a provided calendar. The zeroth
    element is the previous coupon date (PCD) and the first element is the
    Next Coupon Date (NCD). '''

    def __init__(self,
                 effectiveDate: FinDate,   # Also known as the start date
                 terminationDate: FinDate,  # Also known as the termination date
                 freqType: FinFrequencyTypes = FinFrequencyTypes.ANNUAL,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD,
                 adjustTerminationDate:bool = False, # MAY BE REMOVED
                 endOfMonthFlag:bool = True) : # ONLY USED IF TERMINATION DATE IS EOM
        ''' Create FinSchedule object which calculates a sequence of dates in
        line with market convention for fixed income products. 
        
        First we get the unadjusted dates by either stepping forward from the 
        effective date or stepping back from the termination date in steps 
        determined by the period - i.e. the number of months between payments.

        We then adjust all of these dates if they fall on a weekend or holiday 
        according to the calendar specified. These dates are then adjusted in 
        accordance with the business date adjustment.
        
        The End of Month (EOM) logic needs to be understood. Some cases:
        - All unadjusted dates will be EOM if the termination date (TD) is EOM 
        and the EOM flag is true. 
        - If the termination date is EOM and the EOM flag is false, all of
        the unadjusted dates will fall on the 30 if a 30 is a previous date. 
        - If the termination date is the 28 Feb then all dates will be 
        28th if the EOM flag is false. 
        - It is rare for a TD to be EOM and not to have all coupons EOM so
        this is the default.'''

        checkArgumentTypes(self.__init__, locals())

        # validation complete
        self._effectiveDate = effectiveDate
        self._terminationDate = terminationDate
        self._freqType = freqType
        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType
        
        self._adjustTerminationDate = adjustTerminationDate

        if self._terminationDate.isEOM() and endOfMonthFlag is True:
            self._endOfMonthFlag = True
        else:
            self._endOfMonthFlag = False

        self._adjustedDates = None

        self._generate()

###############################################################################

    def scheduleDates(self):
        ''' Returns a list of the schedule of FinDates. '''

        if self._adjustedDates is None:
            self._generate()

        return self._adjustedDates

###############################################################################

    def _generate(self):
        ''' Generate schedule of dates according to specified date generation
        rules and also adjust these dates for holidays according to the
        specified business day convention and the specified calendar. '''

        # NEED TO INCORPORATE ADJUST TERMINATION DATE FLAG

        calendar = FinCalendar(self._calendarType)
        frequency = FinFrequency(self._freqType)
        numMonths = int(12 / frequency)

        unadjustedScheduleDates = []
        self._adjustedDates = []

        if self._dateGenRuleType == FinDateGenRuleTypes.BACKWARD:

            nextDate = self._terminationDate
            flowNum = 0

            while nextDate > self._effectiveDate:

                unadjustedScheduleDates.append(nextDate)
                nextDate = nextDate.addMonths(-numMonths)
                
                if self._endOfMonthFlag is True:
                    nextDate = nextDate.EOM()

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
            nextDate = self._effectiveDate
            flowNum = 0

            unadjustedScheduleDates.append(nextDate)
            flowNum = 1

            while nextDate < self._terminationDate:
                unadjustedScheduleDates.append(nextDate)
                nextDate = nextDate.addMonths(numMonths)
                flowNum = flowNum + 1

            # The first date is not adjusted as it is given
            for i in range(1, flowNum):

                dt = calendar.adjust(unadjustedScheduleDates[i],
                                     self._busDayAdjustType)

                self._adjustedDates.append(dt)

            self._adjustedDates.append(self._terminationDate)

        if self._adjustedDates[0] < self._effectiveDate:
            self._adjustedDates[0] = self._effectiveDate

        return self._adjustedDates

##############################################################################

    def __repr__(self):
        ''' Print out the details of the schedule and the actual dates. This
        can be used for providing transparency on schedule calculations. '''

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EFFECTIVE DATE", self._effectiveDate)
        s += labelToString("END DATE", self._terminationDate)
        s += labelToString("FREQUENCY", self._freqType)
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


###############################################################################

    # def _generate_alternative(self):
    #     ''' This holiday adjusts each date BEFORE generating the next date.
    #     Generate schedule of dates according to specified date generation
    #     rules and also adjust these dates for holidays according to
    #     the business day convention and the specified calendar. 
        
    #     I believe this is incorrect. DO NOT USE. '''

    #     print("======= DO NOT USE THIS =============")

    #     self._adjustedDates = []
    #     calendar = FinCalendar(self._calendarType)
    #     frequency = FinFrequency(self._freqType)
    #     numMonths = int(12 / frequency)

    #     unadjustedScheduleDates = []

    #     if self._dateGenRuleType == FinDateGenRuleTypes.BACKWARD:

    #         nextDate = self._terminationDate
    #         print("END:", nextDate)
    #         flowNum = 0

    #         while nextDate > self._effectiveDate:
    #             unadjustedScheduleDates.append(nextDate)
    #             nextDate = nextDate.addMonths(-numMonths)
    #             nextDate = calendar.adjust(nextDate, self._busDayAdjustType)
    #             flowNum += 1

    #         # Add on the Previous Coupon Date
    #         unadjustedScheduleDates.append(nextDate)
    #         flowNum += 1

    #         # reverse order and holiday adjust dates
    #         for i in range(0, flowNum):
    #             dt = unadjustedScheduleDates[flowNum - i - 1]
    #             self._adjustedDates.append(dt)

    #     elif self._dateGenRuleType == FinDateGenRuleTypes.FORWARD:

    #         # This needs checking
    #         nextDate = self._startDate
    #         flowNum = 0

    #         unadjustedScheduleDates.append(nextDate)
    #         flowNum = 1

    #         while nextDate < self._terminationDate:
    #             unadjustedScheduleDates.append(nextDate)
    #             nextDate = nextDate.addMonths(numMonths)
    #             nextDate = calendar.adjust(nextDate, self._busDayAdjustType)
    #             flowNum = flowNum + 1

    #         for i in range(1, flowNum):

    #             dt = unadjustedScheduleDates[i]
    #             self._adjustedDates.append(dt)

    #         self._adjustedDates.append(self._terminationDate)

    #     if self._adjustedDates[0] < self._effectiveDate:
    #         self._adjustedDates[0] = self._effectiveDate

    #     return self._adjustedDates