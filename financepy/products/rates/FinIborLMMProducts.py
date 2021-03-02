##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Extend to allow term structure of volatility
# TODO: Extend to allow two fixed legs in underlying swap
# TODO: Cash settled swaptions

""" This module implements the LMM in the spot measure. It combines both model
and product specific code - I am not sure if it is better to separate these. At
the moment this seems to work ok.

THIS IS STILL IN PROTOPTYPE MODE. DO NOT USE. """

import numpy as np

from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.day_count import DayCount
from ...utils.schedule import Schedule
from ...utils.FinError import FinError
from ...utils.helper_functions import check_argument_types
from ...utils.date import Date

from ...models.rates_libor_market_model import LMMSimulateFwds1F
from ...models.rates_libor_market_model import LMMSimulateFwdsMF
from ...models.rates_libor_market_model import LMMSimulateFwdsNF
from ...models.rates_libor_market_model import FinRateModelLMMModelTypes
from ...models.rates_libor_market_model import LMMCapFlrPricer

from ...utils.global_variables import gDaysInYear
from ...utils.fin_math import ONE_MILLION

from ...utils.FinGlobalTypes import FinSwapTypes
from ...utils.FinGlobalTypes import FinCapFloorTypes

from financepy.market.volatility.FinIborCapVolCurve import FinIborCapVolCurve

###############################################################################


class FinIborLMMProducts():
    """ This is the class for pricing Ibor products using the LMM. """

    def __init__(self,
                 settlement_date: Date,
                 maturity_date: Date,
                 floatFrequencyType: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 floatDayCountType: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create a European-style swaption by defining the exercise date of
        the swaption, and all of the details of the underlying interest rate
        swap including the fixed coupon and the details of the fixed and the
        floating leg payment schedules. """

        check_argument_types(self.__init__, locals())

        if settlement_date > maturity_date:
            raise FinError("Settlement date must be before maturity date")

        """ Set up the grid for the Ibor rates that are to be simulated. These
        must be consistent with the floating rate leg of the product that is to
        be priced. """

        self._start_date = settlement_date
        self._gridDates = Schedule(settlement_date,
                                   maturity_date,
                                   floatFrequencyType,
                                   calendar_type,
                                   bus_day_adjust_type,
                                   date_gen_rule_type)._generate()

        self._accrualFactors = []
        self._floatDayCountType = floatDayCountType

        basis = DayCount(self._floatDayCountType)
        prevDt = self._gridDates[0]

        self._gridTimes = [0.0]

        for nextDt in self._gridDates[1:]:
            tau = basis.year_frac(prevDt, nextDt)[0]
            t = (nextDt - self._gridDates[0]) / gDaysInYear
            self._accrualFactors.append(tau)
            self._gridTimes.append(t)
            prevDt = nextDt

#        print(self._gridTimes)
        self._accrualFactors = np.array(self._accrualFactors)
        self._numForwards = len(self._accrualFactors)
        self._fwds = None

#        print("Num FORWARDS", self._numForwards)

###############################################################################

    def simulate1F(self,
                   discount_curve,
                   volCurve: FinIborCapVolCurve,
                   num_paths: int = 1000,
                   numeraireIndex: int = 0,
                   useSobol: bool = True,
                   seed: int = 42):
        """ Run the one-factor simulation of the evolution of the forward
        Ibors to generate and store all of the Ibor forward rate paths. """

        if num_paths < 2 or num_paths > 1000000:
            raise FinError("NumPaths must be between 2 and 1 million")

        if discount_curve._valuation_date != self._start_date:
            raise FinError("Curve anchor date not the same as LMM start date.")

        self._num_paths = num_paths
        self._numeraireIndex = numeraireIndex
        self._useSobol = useSobol

        numGridPoints = len(self._gridDates)

        self._numForwards = numGridPoints
        self._forwardCurve = []

        for i in range(1, numGridPoints):
            start_date = self._gridDates[i-1]
            end_date = self._gridDates[i]
            fwd_rate = discount_curve.fwd_rate(start_date,
                                            end_date,
                                            self._floatDayCountType)
            self._forwardCurve.append(fwd_rate)

        self._forwardCurve = np.array(self._forwardCurve)

        gammas = np.zeros(numGridPoints)
        for ix in range(1, numGridPoints):
            dt = self._gridDates[ix]
            gammas[ix] = volCurve.capletVol(dt)

        self._fwds = LMMSimulateFwds1F(self._numForwards,
                                       num_paths,
                                       numeraireIndex,
                                       self._forwardCurve,
                                       gammas,
                                       self._accrualFactors,
                                       useSobol,
                                       seed)

###############################################################################

    def simulateMF(self,
                   discount_curve,
                   numFactors: int,
                   lambdas: np.ndarray,
                   num_paths: int = 10000,
                   numeraireIndex: int = 0,
                   useSobol: bool = True,
                   seed: int = 42):
        """ Run the simulation to generate and store all of the Ibor forward
        rate paths. This is a multi-factorial version so the user must input
        a numpy array consisting of a column for each factor and the number of
        rows must equal the number of grid times on the underlying simulation
        grid. CHECK THIS. """

#        check_argument_types(self.__init__, locals())

        if num_paths < 2 or num_paths > 1000000:
            raise FinError("NumPaths must be between 2 and 1 million")

        if discount_curve._curveDate != self._start_date:
            raise FinError("Curve anchor date not the same as LMM start date.")

        print("LEN LAMBDAS", len(lambdas))
        print("LEN", len(lambdas[0]))
        # We pass a vector of vol curves, one for each factor
        if numFactors != len(lambdas):
            raise FinError("Lambda doesn't have specified number of factors.")

        numRows = len(lambdas[0])
        if numRows != self._numForwards+1:
            raise FinError("Vol Components needs same number of rows as grid")

        self._num_paths = num_paths
        self._numeraireIndex = numeraireIndex
        self._useSobol = useSobol

        self._numForwards = len(self._gridDates) - 1
        self._forwardCurve = []

        for i in range(1, self._numForwards):
            start_date = self._gridDates[i-1]
            end_date = self._gridDates[i]
            fwd_rate = discount_curve.fwd_rate(start_date, end_date,
                                            self._floatDayCountType)
            self._forwardCurve.append(fwd_rate)

        self._forwardCurve = np.array(self._forwardCurve)

        self._fwds = LMMSimulateFwdsMF(self._numForwards,
                                       numFactors,
                                       num_paths,
                                       numeraireIndex,
                                       self._forwardCurve,
                                       lambdas,
                                       self._accrualFactors,
                                       useSobol,
                                       seed)

###############################################################################

    def simulateNF(self,
                   discount_curve,
                   volCurve: FinIborCapVolCurve,
                   correlationMatrix: np.ndarray,
                   modelType: FinRateModelLMMModelTypes,
                   num_paths: int = 1000,
                   numeraireIndex: int = 0,
                   useSobol: bool = True,
                   seed: int = 42):
        """ Run the simulation to generate and store all of the Ibor forward
        rate paths using a full factor reduction of the fwd-fwd correlation
        matrix using Cholesky decomposition."""

        check_argument_types(self.__init__, locals())

        if num_paths < 2 or num_paths > 1000000:
            raise FinError("NumPaths must be between 2 and 1 million")

        if isinstance(modelType, FinRateModelLMMModelTypes) is False:
            raise FinError("Model type must be type FinRateModelLMMModelTypes")

        if discount_curve.curveDate != self._start_date:
            raise FinError("Curve anchor date not the same as LMM start date.")

        self._num_paths = num_paths
        self._volCurves = volCurve
        self._correlationMatrix = correlationMatrix
        self._modelType = modelType
        self._numeraireIndex = numeraireIndex
        self._useSobol = useSobol

        numGridPoints = len(self._gridTimes)

        self._numForwards = numGridPoints - 1
        self._forwardCurve = []

        for i in range(1, numGridPoints):
            start_date = self._gridDates[i-1]
            end_date = self._gridDates[i]
            fwd_rate = discount_curve.forwardRate(start_date,
                                                end_date,
                                                self._floatDayCountType)
            self._forwardCurve.append(fwd_rate)

        self._forwardCurve = np.array(self._forwardCurve)

        zetas = np.zeros(numGridPoints)
        for ix in range(1, numGridPoints):
            dt = self._gridDates[ix]
            zetas[ix] = volCurve.capletVol(dt)

        # This function does not use Sobol - TODO
        self._fwds = LMMSimulateFwdsNF(self._numForwards,
                                       num_paths,
                                       self._forwardCurve,
                                       zetas,
                                       correlationMatrix,
                                       self._accrualFactors,
                                       seed)

###############################################################################

    def valueSwaption(self,
                      settlement_date: Date,
                      exerciseDate: Date,
                      maturity_date: Date,
                      swaptionType: FinSwapTypes,
                      fixedCoupon: float,
                      fixedFrequencyType: FrequencyTypes,
                      fixedDayCountType: DayCountTypes,
                      notional: float = ONE_MILLION,
                      floatFrequencyType: FrequencyTypes = FrequencyTypes.QUARTERLY,
                      floatDayCountType: DayCountTypes = DayCountTypes.THIRTY_E_360,
                      calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                      bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                      date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Value a swaption in the LMM model using simulated paths of the
        forward curve. This relies on pricing the fixed leg of the swap and
        assuming that the floating leg will be worth par. As a result we only
        need simulate Ibors with the frequency of the fixed leg. """

        # Note that the simulation time steps run all the way out to the last
        # forward rate. However we only really need the forward rates at the
        # expiry date of the option. It may be worth amending the simulate
        # code to impose a limit on the time steps in order to speed up the
        # overall pricing if it requires a new run every time. However once
        # generated, the speed of pricing is not affected so this is not
        # strictly an urgent issue.

        swaptionFloatDates = Schedule(settlement_date,
                                      maturity_date,
                                      floatFrequencyType,
                                      calendar_type,
                                      bus_day_adjust_type,
                                      date_gen_rule_type)._generate()

        for swaptionDt in swaptionFloatDates:
            foundDt = False
            for gridDt in self._gridDates:
                if swaptionDt == gridDt:
                    foundDt = True
                    break
            if foundDt is False:
                raise FinError("Swaption float leg not on grid.")

        swaptionFixedDates = Schedule(settlement_date,
                                      maturity_date,
                                      fixedFrequencyType,
                                      calendar_type,
                                      bus_day_adjust_type,
                                      date_gen_rule_type)._generate()

        for swaptionDt in swaptionFixedDates:
            foundDt = False
            for gridDt in self._gridDates:
                if swaptionDt == gridDt:
                    foundDt = True
                    break
            if foundDt is False:
                raise FinError("Swaption fixed leg not on grid.")

        a = 0
        b = 0

        for gridDt in self._gridDates:
            if gridDt == exerciseDate:
                break
            else:
                a += 1

        for gridDt in self._gridDates:
            if gridDt == maturity_date:
                break
            else:
                b += 1

        if b == 0:
            raise FinError("Swaption swap maturity date is today.")

#        num_paths = 1000
#        v = LMMSwaptionPricer(fixedCoupon, a, b, num_paths,
#                              fwd0, fwds, taus, isPayer)
        v = 0.0
        return v

###############################################################################

    def valueCapFloor(self,
                      settlement_date: Date,
                      maturity_date: Date,
                      capFloorType: FinCapFloorTypes,
                      capFloorRate: float,
                      frequencyType: FrequencyTypes = FrequencyTypes.QUARTERLY,
                      day_count_type: DayCountTypes = DayCountTypes.ACT_360,
                      notional: float = ONE_MILLION,
                      calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                      bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                      date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Value a cap or floor in the LMM. """

        capFloorDates = Schedule(settlement_date,
                                 maturity_date,
                                 frequencyType,
                                 calendar_type,
                                 bus_day_adjust_type,
                                 date_gen_rule_type)._generate()

        for capFloorletDt in capFloorDates:
            foundDt = False
            for gridDt in self._gridDates:
                if capFloorletDt == gridDt:
                    foundDt = True
                    break
            if foundDt is False:
                raise FinError("CapFloor date not on grid.")

        numFowards = len(capFloorDates)
        num_paths = self._num_paths
        K = capFloorRate
        isCap = 0
        if capFloorType == FinCapFloorTypes.CAP:
            isCap = 1

        fwd0 = self._forwardCurve
        fwds = self._fwds
        taus = self._accrualFactors

        v = LMMCapFlrPricer(numFowards, num_paths, K, fwd0, fwds, taus, isCap)

        # Sum the cap/floorlets to get cap/floor value
        v_capFloor = 0.0
        for v_capFloorLet in v:
            v_capFloor += v_capFloorLet * notional

        return v_capFloor

###############################################################################

    def __repr__(self):
        """ Function to allow us to print the LMM Products details. """

        s = "Function not written"
        return s

###############################################################################

    def _print(self):
        """ Alternative print method. """

        print(self)

###############################################################################
