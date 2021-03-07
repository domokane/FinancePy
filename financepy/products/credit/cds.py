##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from numba import njit, float64, int64
from math import exp, log
from copy import deepcopy

from ...utils.date import Date
from ...utils.FinError import FinError
from ...utils.calendar import Calendar, CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.global_vars import gDaysInYear
from ...utils.math import ONE_MILLION
from ...utils.helpers import labelToString, tableToString
from ...market.discount.interpolator import FinInterpTypes, _uinterpolate

from ...utils.helpers import check_argument_types

useFlatHazardRateIntegral = True
standard_recovery_rate = 0.40


###############################################################################
# TODO: Perform protection leg pv analytically using fact that hazard rate and
#       interest rates are flat between their combined node points. Right now I
#       do not find the protection leg PV calculations to be a bottleneck,
#       especially given the speedup benefits of using NUMBA.
###############################################################################


@njit(float64[:](float64, float64, float64[:], float64[:], float64[:],
                 float64[:], float64[:], float64[:], int64),
      fastmath=True, cache=True)
def _riskyPV01_NUMBA(teff,
                     accrual_factorPCDToNow,
                     paymentTimes,
                     year_fracs,
                     npIborTimes,
                     npIborValues,
                     npSurvTimes,
                     npSurvValues,
                     pv01_method):
    """ Fast calculation of the risky PV01 of a CDS using NUMBA.
    The output is a numpy array of the full and clean risky PV01."""

    method = FinInterpTypes.FLAT_FWD_RATES.value

    couponAccruedIndicator = 1

    # Method 0 : This is the market standard which assumes that the coupon
    # accrued is treated as though on average default occurs roughly midway
    # through a coupon period.

    tncd = paymentTimes[1]

    # The first coupon is a special case which needs to be handled carefully
    # taking into account what coupon has already accrued and what has not
    qeff = _uinterpolate(teff, npSurvTimes, npSurvValues, method)
    q1 = _uinterpolate(tncd, npSurvTimes, npSurvValues, method)
    z1 = _uinterpolate(tncd, npIborTimes, npIborValues, method)

    # this is the part of the coupon accrued from previous coupon date to now
    # accrual_factorPCDToNow = dayCount.year_frac(pcd,teff)

    # reference credit survives to the premium payment date
    fullRPV01 = q1 * z1 * year_fracs[1]

    # coupon accrued from previous coupon to today paid in full at default
    # before coupon payment
    fullRPV01 = fullRPV01 + z1 * \
                (qeff - q1) * accrual_factorPCDToNow * couponAccruedIndicator

    # future accrued from now to coupon payment date assuming default roughly
    # midway
    fullRPV01 += 0.5 * z1 * \
                 (qeff - q1) * (year_fracs[1] - accrual_factorPCDToNow) \
                 * couponAccruedIndicator

    for it in range(2, len(paymentTimes)):

        t2 = paymentTimes[it]

        q2 = _uinterpolate(t2, npSurvTimes, npSurvValues, method)
        z2 = _uinterpolate(t2, npIborTimes, npIborValues, method)

        accrual_factor = year_fracs[it]

        # full coupon is paid at the end of the current period if survives to
        # payment date
        fullRPV01 += q2 * z2 * accrual_factor

        #######################################################################

        if couponAccruedIndicator == 1:

            if useFlatHazardRateIntegral:
                # This needs to be updated to handle small h+r
                tau = accrual_factor
                h12 = -log(q2 / q1) / tau
                r12 = -log(z2 / z1) / tau
                alpha = h12 + r12
                expTerm = 1.0 - exp(-alpha * tau) - alpha * \
                          tau * exp(-alpha * tau)
                dfullRPV01 = q1 * z1 * h12 * \
                             expTerm / abs(alpha * alpha + 1e-20)
            else:
                dfullRPV01 = 0.50 * (q1 - q2) * z2 * accrual_factor

            fullRPV01 = fullRPV01 + dfullRPV01

        #######################################################################

        q1 = q2

    cleanRPV01 = fullRPV01 - accrual_factorPCDToNow

    return np.array([fullRPV01, cleanRPV01])


###############################################################################

@njit(float64(float64, float64, float64[:], float64[:], float64[:], float64[:],
              float64, int64, int64), fastmath=True, cache=True)
def _protectionLegPV_NUMBA(teff,
                           tmat,
                           npIborTimes,
                           npIborValues,
                           npSurvTimes,
                           npSurvValues,
                           contract_recovery_rate,
                           num_steps_per_year,
                           protMethod):
    """ Fast calculation of the CDS protection leg PV using NUMBA to speed up
    the numerical integration over time. """

    method = FinInterpTypes.FLAT_FWD_RATES.value
    dt = (tmat - teff) / num_steps_per_year
    t = teff
    z1 = _uinterpolate(t, npIborTimes, npIborValues, method)
    q1 = _uinterpolate(t, npSurvTimes, npSurvValues, method)

    prot_pv = 0.0
    small = 1e-8

    if useFlatHazardRateIntegral is True:

        for _ in range(0, num_steps_per_year):
            t = t + dt
            z2 = _uinterpolate(t, npIborTimes, npIborValues, method)
            q2 = _uinterpolate(t, npSurvTimes, npSurvValues, method)
            # This needs to be updated to handle small h+r
            h12 = -log(q2 / q1) / dt
            r12 = -log(z2 / z1) / dt
            expTerm = exp(-(r12 + h12) * dt)
            dprot_pv = h12 * (1.0 - expTerm) * q1 * z1 / \
                       (abs(h12 + r12) + small)
            prot_pv += dprot_pv
            q1 = q2
            z1 = z2

    else:

        for _ in range(0, num_steps_per_year):
            t += dt
            z2 = _uinterpolate(t, npIborTimes, npIborValues, method)
            q2 = _uinterpolate(t, npSurvTimes, npSurvValues, method)
            dq = q1 - q2
            dprot_pv = 0.5 * (z1 + z2) * dq
            prot_pv += dprot_pv
            q1 = q2
            z1 = z2

    prot_pv = prot_pv * (1.0 - contract_recovery_rate)
    return prot_pv


###############################################################################
###############################################################################
###############################################################################


class FinCDS(object):
    """ A class which manages a Credit Default Swap. It performs schedule
    generation and the valuation and risk management of CDS. """

    def __init__(self,
                 step_in_date: Date,  # Date protection starts
                 maturity_date_or_tenor: (Date, str),  # Date or tenor
                 running_coupon: float,  # Annualised coupon on premium fee leg
                 notional: float = ONE_MILLION,
                 long_protection: bool = True,
                 freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 day_count_type: DayCountTypes = DayCountTypes.ACT_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create a CDS from the step-in date, maturity date and coupon """

        check_argument_types(self.__init__, locals())

        if type(maturity_date_or_tenor) == Date:
            maturity_date = maturity_date_or_tenor
        else:
            # To get the next CDS date we move on by the tenor and then roll to
            # the next CDS date after that. We do not holiday adjust it. That
            # is handled in the schedule generation.
            maturity_date = step_in_date.addTenor(maturity_date_or_tenor)
            maturity_date = maturity_date.nextCDSDate()

        if step_in_date > maturity_date:
            raise FinError("Step in date after maturity date")

        self._step_in_date = step_in_date
        self._maturity_date = maturity_date
        self._running_coupon = running_coupon
        self._notional = notional
        self._long_protection = long_protection
        self._day_count_type = day_count_type
        self._date_gen_rule_type = date_gen_rule_type
        self._calendar_type = calendar_type
        self._freq_type = freq_type
        self._bus_day_adjust_type = bus_day_adjust_type

        self._generateAdjustedCDSPaymentDates()
        self._calcFlows()

    ###############################################################################

    def _generateAdjustedCDSPaymentDates(self):
        """ Generate CDS payment dates which have been holiday adjusted."""
        frequency = annual_frequency(self._freq_type)
        calendar = Calendar(self._calendar_type)
        start_date = self._step_in_date
        end_date = self._maturity_date

        self._adjusted_dates = []
        numMonths = int(12.0 / frequency)

        unadjustedScheduleDates = []

        if self._date_gen_rule_type == DateGenRuleTypes.BACKWARD:

            next_date = end_date
            flowNum = 0

            while next_date > start_date:
                unadjustedScheduleDates.append(next_date)
                next_date = next_date.addMonths(-numMonths)
                flowNum += 1

            # Add on the Previous Coupon Date
            unadjustedScheduleDates.append(next_date)
            flowNum += 1

            # reverse order
            for i in range(0, flowNum):
                dt = unadjustedScheduleDates[flowNum - i - 1]
                self._adjusted_dates.append(dt)

            # holiday adjust dates except last one
            for i in range(0, flowNum - 1):
                dt = calendar.adjust(self._adjusted_dates[i],
                                     self._bus_day_adjust_type)

                self._adjusted_dates[i] = dt

            finalDate = self._adjusted_dates[flowNum - 1]

            # Final date is moved forward by one day
            self._adjusted_dates[flowNum - 1] = finalDate.addDays(1)

        elif self._date_gen_rule_type == DateGenRuleTypes.FORWARD:

            next_date = start_date
            flowNum = 0

            unadjustedScheduleDates.append(next_date)
            flowNum = 1

            while next_date < end_date:
                unadjustedScheduleDates.append(next_date)
                next_date = next_date.addMonths(numMonths)
                flowNum = flowNum + 1

            for i in range(1, flowNum):
                dt = calendar.adjust(unadjustedScheduleDates[i],
                                     self._bus_day_adjust_type)

                self._adjusted_dates.append(dt)

            finalDate = end_date.addDays(1)
            self._adjusted_dates.append(finalDate)

        else:
            raise FinError("Unknown DateGenRuleType:" +
                           str(self._date_gen_rule_type))

    ###############################################################################

    def _calcFlows(self):
        """ Calculate cash flow amounts on premium leg. """
        payment_dates = self._adjusted_dates
        dayCount = DayCount(self._day_count_type)

        self._accrual_factors = []
        self._flows = []

        self._accrual_factors.append(0.0)
        self._flows.append(0.0)

        num_flows = len(payment_dates)

        for it in range(1, num_flows):
            t0 = payment_dates[it - 1]
            t1 = payment_dates[it]
            accrual_factor = dayCount.year_frac(t0, t1)[0]
            flow = accrual_factor * self._running_coupon * self._notional

            self._accrual_factors.append(accrual_factor)
            self._flows.append(flow)

    ###############################################################################

    def value(self,
              valuation_date,
              issuer_curve,
              contract_recovery_rate=standard_recovery_rate,
              pv01_method=0,
              prot_method=0,
              num_steps_per_year=25):
        """ Valuation of a CDS contract on a specific valuation date given
        an issuer curve and a contract recovery rate."""

        rpv01 = self.riskyPV01(valuation_date,
                               issuer_curve,
                               pv01_method)

        fullRPV01 = rpv01['full_rpv01']
        cleanRPV01 = rpv01['clean_rpv01']

        prot_pv = self.protectionLegPV(valuation_date,
                                       issuer_curve,
                                       contract_recovery_rate,
                                       num_steps_per_year,
                                       prot_method)

        fwdDf = 1.0

        if self._long_protection:
            longProt = +1
        else:
            longProt = -1

        fullPV = fwdDf * longProt * \
                 (prot_pv - self._running_coupon * fullRPV01 * self._notional)
        cleanPV = fwdDf * longProt * \
                  (prot_pv - self._running_coupon * cleanRPV01 * self._notional)

        #        print("protLeg", prot_pv, "cleanRPV01", cleanRPV01, "value", cleanPV)

        return {'full_pv': fullPV, 'clean_pv': cleanPV}

    ###############################################################################

    def creditDV01(self,
                   valuation_date,
                   issuer_curve,
                   contract_recovery_rate=standard_recovery_rate,
                   pv01_method=0,
                   prot_method=0,
                   num_steps_per_year=25):
        """ Calculation of the change in the value of the CDS contract for a
        one basis point change in the level of the CDS curve."""

        v0 = self.value(valuation_date,
                        issuer_curve,
                        contract_recovery_rate,
                        pv01_method,
                        prot_method,
                        num_steps_per_year)

        bump = 0.0001  # 1 basis point

        # we create a deep copy to avoid state issues
        bumpedIssuerCurve = deepcopy(issuer_curve)
        for cds in bumpedIssuerCurve._cds_contracts:
            cds._running_coupon += bump

        bumpedIssuerCurve._buildCurve()

        v1 = self.value(valuation_date,
                        bumpedIssuerCurve,
                        contract_recovery_rate,
                        pv01_method,
                        prot_method,
                        num_steps_per_year)

        creditDV01 = (v1['full_pv'] - v0['full_pv'])
        return creditDV01

    ###############################################################################

    def interestDV01(self,
                     valuation_date: Date,
                     issuer_curve,
                     contract_recovery_rate=standard_recovery_rate,
                     pv01_method: int = 0,
                     prot_method: int = 0,
                     num_steps_per_year: int = 25):
        """ Calculation of the interest DV01 based on a simple bump of
        the discount factors and reconstruction of the CDS curve. """

        v0 = self.value(valuation_date,
                        issuer_curve,
                        contract_recovery_rate,
                        pv01_method,
                        prot_method,
                        num_steps_per_year)

        # we create a deep copy to avoid state issues
        new_issuer_curve = deepcopy(issuer_curve)

        bump = 0.0001  # 1 basis point

        for depo in new_issuer_curve._libor_curve._usedDeposits:
            depo._deposit_rate += bump
        for fra in new_issuer_curve._libor_curve._usedFRAs:
            fra._fraRate += bump
        for swap in new_issuer_curve._libor_curve._usedSwaps:
            swap._fixed_leg._coupon += bump

        new_issuer_curve._libor_curve._buildCurve()

        new_issuer_curve._buildCurve()

        v1 = self.value(valuation_date,
                        new_issuer_curve,
                        contract_recovery_rate,
                        pv01_method,
                        prot_method,
                        num_steps_per_year)

        interestDV01 = (v1['full_pv'] - v0['full_pv'])
        return interestDV01

    ###############################################################################

    def cashSettlementAmount(self,
                             valuation_date,
                             settlement_date,
                             issuer_curve,
                             contract_recovery_rate=standard_recovery_rate,
                             pv01_method=0,
                             prot_method=0,
                             num_steps_per_year=25):
        """ Value of the contract on the settlement date including accrued
        interest. """

        v = self.value(valuation_date,
                       issuer_curve,
                       contract_recovery_rate,
                       pv01_method,
                       prot_method,
                       num_steps_per_year)

        libor_curve = issuer_curve._libor_curve
        df = libor_curve.df(settlement_date)
        v = v / df
        return v

    ###############################################################################

    def clean_price(self,
                    valuation_date,
                    issuer_curve,
                    contract_recovery_rate=standard_recovery_rate,
                    pv01_method=0,
                    prot_method=0,
                    num_steps_per_year=52):
        """ Value of the CDS contract excluding accrued interest. """

        riskyPV01 = self.riskyPV01(valuation_date, issuer_curve, pv01_method)

        cleanRPV01 = riskyPV01['clean_rpv01']

        prot_pv = self.protectionLegPV(valuation_date,
                                       issuer_curve,
                                       contract_recovery_rate,
                                       num_steps_per_year,
                                       prot_method)

        fwdDf = 1.0

        cleanPV = fwdDf * (prot_pv - self._running_coupon * cleanRPV01
                           * self._notional)
        clean_price = (self._notional - cleanPV) / self._notional * 100.0
        return clean_price

    ###############################################################################

    def riskyPV01OLD(self,
                     valuation_date,
                     issuer_curve,
                     pv01_method=0):
        """ RiskyPV01 of the contract using the OLD method. """

        payment_dates = self._adjusted_dates
        dayCount = DayCount(self._day_count_type)

        couponAccruedIndicator = 1

        # Method 0 : This is the market standard which assumes that the coupon
        # accrued is treated as though on average default occurs roughly midway
        # through a coupon period.

        teff = self._step_in_date
        pcd = payment_dates[0]  # PCD
        ncd = payment_dates[1]  # NCD

        # The first coupon is a special case which must be handled carefully
        # taking into account what coupon has already accrued and what has not
        qeff = issuer_curve.survivalProbability(teff)
        q1 = issuer_curve.survivalProbability(ncd)
        z1 = issuer_curve.df(ncd)

        # this is the part of the coupon accrued from the previous coupon date
        # to now
        accrual_factorPCDToNow = dayCount.year_frac(pcd, teff)[0]

        # full first coupon is paid at the end of the current period if the
        year_frac = dayCount.year_frac(pcd, ncd)[0]

        # reference credit survives to the premium payment date
        fullRPV01 = q1 * z1 * year_frac

        # coupon accrued from previous coupon to today paid in full at default
        # before coupon payment
        fullRPV01 = fullRPV01 + z1 * \
                    (qeff - q1) * accrual_factorPCDToNow * couponAccruedIndicator

        # future accrued from now to coupon payment date assuming default
        # roughly midway
        fullRPV01 = fullRPV01 + 0.5 * z1 * \
                    (qeff - q1) * (year_frac - accrual_factorPCDToNow) \
                    * couponAccruedIndicator

        for it in range(2, len(payment_dates)):

            t1 = payment_dates[it - 1]
            t2 = payment_dates[it]
            q2 = issuer_curve.survivalProbability(t2)
            z2 = issuer_curve.df(t2)

            accrual_factor = dayCount.year_frac(t1, t2)[0]

            # full coupon is paid at the end of the current period if survives
            # to payment date
            fullRPV01 += q2 * z2 * accrual_factor

            ###################################################################

            if couponAccruedIndicator == 1:

                if useFlatHazardRateIntegral:
                    # This needs to be updated to handle small h+r
                    tau = accrual_factor
                    h12 = -log(q2 / q1) / tau
                    r12 = -log(z2 / z1) / tau
                    alpha = h12 + r12
                    expTerm = 1.0 - exp(-alpha * tau) - \
                              alpha * tau * exp(-alpha * tau)
                    dfullRPV01 = q1 * z1 * h12 * \
                                 expTerm / abs(alpha * alpha + 1e-20)
                else:
                    dfullRPV01 = 0.50 * (q1 - q2) * z2 * accrual_factor

                fullRPV01 = fullRPV01 + dfullRPV01

            ###################################################################

            q1 = q2

        cleanRPV01 = fullRPV01 - accrual_factorPCDToNow

        #        print("OLD PV01",fullRPV01, cleanRPV01)

        return {'full_rpv01': fullRPV01, 'clean_rpv01': cleanRPV01}

    ###############################################################################

    def accrued_days(self):
        """ Number of days between the previous coupon and the currrent step
        in date. """

        # I assume accrued runs to the effective date
        payment_dates = self._adjusted_dates
        pcd = payment_dates[0]
        accrued_days = (self._step_in_date - pcd)
        return accrued_days

    ###############################################################################

    def accruedInterest(self):
        """ Calculate the amount of accrued interest that has accrued from the
        previous coupon date (PCD) to the step_in_date of the CDS contract. """

        dayCount = DayCount(self._day_count_type)
        payment_dates = self._adjusted_dates
        pcd = payment_dates[0]
        accrual_factor = dayCount.year_frac(pcd, self._step_in_date)[0]
        accruedInterest = accrual_factor * self._notional * self._running_coupon

        if self._long_protection:
            accruedInterest *= -1.0

        return accruedInterest

    ###############################################################################

    def protectionLegPV(self,
                        valuation_date,
                        issuer_curve,
                        contract_recovery_rate=standard_recovery_rate,
                        num_steps_per_year=25,
                        protMethod=0):
        """ Calculates the protection leg PV of the CDS by calling into the
        fast NUMBA code that has been defined above. """

        teff = (self._step_in_date - valuation_date) / gDaysInYear
        tmat = (self._maturity_date - valuation_date) / gDaysInYear

        libor_curve = issuer_curve._libor_curve

        v = _protectionLegPV_NUMBA(teff,
                                   tmat,
                                   libor_curve._times,
                                   libor_curve._dfs,
                                   issuer_curve._times,
                                   issuer_curve._values,
                                   contract_recovery_rate,
                                   num_steps_per_year,
                                   protMethod)

        return v * self._notional

    ###############################################################################

    def riskyPV01(self,
                  valuation_date,
                  issuer_curve,
                  pv01_method=0):
        """ The riskyPV01 is the present value of a risky one dollar paid on
        the premium leg of a CDS contract. """

        libor_curve = issuer_curve._libor_curve

        paymentTimes = []
        for it in range(0, len(self._adjusted_dates)):
            t = (self._adjusted_dates[it] - valuation_date) / gDaysInYear
            paymentTimes.append(t)

        # this is the part of the coupon accrued from the previous coupon date
        # to now
        pcd = self._adjusted_dates[0]
        eff = self._step_in_date
        dayCount = DayCount(self._day_count_type)

        accrual_factorPCDToNow = dayCount.year_frac(pcd, eff)[0]

        year_fracs = self._accrual_factors
        teff = (eff - valuation_date) / gDaysInYear

        valueRPV01 = _riskyPV01_NUMBA(teff,
                                      accrual_factorPCDToNow,
                                      np.array(paymentTimes),
                                      np.array(year_fracs),
                                      libor_curve._times,
                                      libor_curve._dfs,
                                      issuer_curve._times,
                                      issuer_curve._values,
                                      pv01_method)

        fullRPV01 = valueRPV01[0]
        cleanRPV01 = valueRPV01[1]

        #        print("NEW PV01",fullRPV01, cleanRPV01)
        return {'full_rpv01': fullRPV01, 'clean_rpv01': cleanRPV01}

    ###############################################################################

    def premiumLegPV(self,
                     valuation_date,
                     issuer_curve,
                     pv01_method=0):
        """ Value of the premium leg of a CDS. """

        fullRPV01 = self.riskyPV01(valuation_date,
                                   issuer_curve,
                                   pv01_method)['full_rpv01']

        v = fullRPV01 * self._notional * self._running_coupon
        return v

    ###############################################################################

    def parSpread(self,
                  valuation_date,
                  issuer_curve,
                  contract_recovery_rate=standard_recovery_rate,
                  num_steps_per_year=25,
                  pv01_method=0,
                  protMethod=0):
        """ Breakeven CDS coupon that would make the value of the CDS contract
        equal to zero. """

        cleanRPV01 = self.riskyPV01(valuation_date,
                                    issuer_curve,
                                    pv01_method)['clean_rpv01']

        prot = self.protectionLegPV(valuation_date,
                                    issuer_curve,
                                    contract_recovery_rate,
                                    num_steps_per_year,
                                    protMethod)

        # By convention this is calculated using the clean RPV01
        spd = prot / cleanRPV01 / self._notional
        return spd

    ###############################################################################

    def valueFastApprox(self,
                        valuation_date,
                        flatContinuousInterestRate,
                        flatCDSCurveSpread,
                        curveRecovery=standard_recovery_rate,
                        contract_recovery_rate=standard_recovery_rate):
        """ Implementation of fast valuation of the CDS contract using an
        accurate approximation that avoids curve building. """

        if type(valuation_date) is not Date:
            raise FinError("Valuation date must be a Date and not " +
                           str(valuation_date))

        t_mat = (self._maturity_date - valuation_date) / gDaysInYear
        t_eff = (self._step_in_date - valuation_date) / gDaysInYear

        h = flatCDSCurveSpread / (1.0 - curveRecovery)
        r = flatContinuousInterestRate
        fwdDf = 1.0
        bump_size = 0.0001

        if self._long_protection:
            long_protection = +1
        else:
            long_protection = -1

        # The sign of he accrued has already been sign adjusted for direction
        accrued = self.accruedInterest()

        # This is the clean RPV01 as it treats the PV01 stream as though it
        # pays just the accrued for the time between 0 and the maturity
        # It therefore omits the part that has accrued

        w = r + h
        z = np.exp(-w * t_eff) - np.exp(-w * t_mat)
        cleanRPV01 = (z / w) * 365.0 / 360.0
        prot_pv = h * (1.0 - contract_recovery_rate) * (z / w) * self._notional
        cleanPV = fwdDf * long_protection * \
                  (prot_pv - self._running_coupon * cleanRPV01 * self._notional)
        fullPV = cleanPV + fwdDf * accrued

        #######################################################################
        # bump CDS spread and calculate
        #######################################################################

        h = (flatCDSCurveSpread + bump_size) / (1.0 - contract_recovery_rate)
        r = flatContinuousInterestRate
        w = r + h
        z = np.exp(-w * t_eff) - np.exp(-w * t_mat)
        cleanRPV01 = (z / w) * 365.0 / 360.0
        prot_pv = h * (1.0 - contract_recovery_rate) * (z / w) * self._notional
        cleanPV_credit_bumped = fwdDf * long_protection * \
                                (prot_pv - self._running_coupon * cleanRPV01 * self._notional)
        fullPV_credit_bumped = cleanPV_credit_bumped \
                               + fwdDf * long_protection * accrued
        credit01 = fullPV_credit_bumped - fullPV

        #######################################################################
        # bump Rate and calculate
        #######################################################################

        h = flatCDSCurveSpread / (1.0 - contract_recovery_rate)
        r = flatContinuousInterestRate + bump_size

        w = r + h
        z = np.exp(-w * t_eff) - np.exp(-w * t_mat)
        cleanRPV01 = (z / w) * 365.0 / 360.0
        prot_pv = h * (1.0 - contract_recovery_rate) * (z / w) * self._notional
        cleanPV_ir_bumped = fwdDf * long_protection * \
                            (prot_pv - self._running_coupon * cleanRPV01 * self._notional)
        fullPV_ir_bumped = cleanPV_ir_bumped + fwdDf * long_protection * accrued
        ir01 = fullPV_ir_bumped - fullPV

        return (fullPV, cleanPV, credit01, ir01)

    ###############################################################################

    def printFlows(self, issuer_curve):

        num_flows = len(self._adjusted_dates)

        print("PAYMENT_DATE      YEAR_FRAC      FLOW           DF       SURV_PROB      NPV")

        for it in range(1, num_flows):
            dt = self._adjusted_dates[it]
            acc_factor = self._accrual_factors[it]
            flow = self._flows[it]
            z = issuer_curve.df(dt)
            q = issuer_curve.survProb(dt)
            print("%15s %10.6f %12.2f %12.6f %12.6f %12.2f" %
                  (dt, acc_factor, flow, z, q, flow * z * q))

    ###############################################################################

    def __repr__(self):
        """ print out details of the CDS contract and all of the calculated
        cash flows """
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("STEP-IN DATE", self._step_in_date)
        s += labelToString("MATURITY", self._maturity_date)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("RUNNING COUPON", self._running_coupon * 10000, "bp\n")
        s += labelToString("DAYCOUNT", self._day_count_type)
        s += labelToString("FREQUENCY", self._freq_type)
        s += labelToString("CALENDAR", self._calendar_type)
        s += labelToString("BUSDAYRULE", self._bus_day_adjust_type)
        s += labelToString("DATEGENRULE", self._date_gen_rule_type)
        s += labelToString("ACCRUED DAYS", self.accrued_days())

        header = "PAYMENT_DATE, YEAR_FRAC, FLOW"
        valueTable = [self._adjusted_dates, self._accrual_factors, self._flows]
        precision = "12.6f"

        s += tableToString(header, valueTable, precision)

        return s

    ###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
