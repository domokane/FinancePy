###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import SwapTypes
from financepy.utils.global_types import FinExerciseTypes
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.bermudan_swaption import IborBermudanSwaption
from financepy.models.black import Black
from financepy.models.bk_tree import BKTree
from financepy.models.hw_tree import HWTree
from financepy.models.bdt_tree import BDTTree
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from FinTestCases import FinTestCases, globalTestCaseMode


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_IborBermudanSwaptionBKModel():
    """ Replicate examples in paper by Leif Andersen which can be found at
    file:///C:/Users/Dominic/Downloads/SSRN-id155208.pdf """

    value_date = Date(1, 1, 2011)
    settle_date = value_date
    exercise_date = settle_date.add_years(1)
    swap_maturity_date = settle_date.add_years(4)

    swap_fixed_coupon = 0.060
    swap_fixed_frequency_type = FrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = DayCountTypes.ACT_365F

    libor_curve = DiscountCurveFlat(value_date,
                                    0.0625,
                                    FrequencyTypes.SEMI_ANNUAL,
                                    DayCountTypes.ACT_365F)

    fwdPAYSwap = IborSwap(exercise_date,
                          swap_maturity_date,
                          SwapTypes.PAY,
                          swap_fixed_coupon,
                          swap_fixed_frequency_type,
                          swapFixedDayCountType)

    fwdSwapValue = fwdPAYSwap.value(settle_date, libor_curve, libor_curve)

    testCases.header("LABEL", "VALUE")
    testCases.print("FWD SWAP VALUE", fwdSwapValue)

    # fwdPAYSwap.print_fixed_leg_pv()

    # Now we create the European swaptions
    fixed_leg_type = SwapTypes.PAY
    europeanSwaptionPay = IborSwaption(settle_date,
                                       exercise_date,
                                       swap_maturity_date,
                                       fixed_leg_type,
                                       swap_fixed_coupon,
                                       swap_fixed_frequency_type,
                                       swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    europeanSwaptionRec = IborSwaption(settle_date,
                                       exercise_date,
                                       swap_maturity_date,
                                       fixed_leg_type,
                                       swap_fixed_coupon,
                                       swap_fixed_frequency_type,
                                       swapFixedDayCountType)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BLACK'S MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    testCases.banner("======= ZERO VOLATILITY ========")
    model = Black(0.0000001)
    testCases.print("Black Model", model._volatility)

    valuePay = europeanSwaptionPay.value(settle_date, libor_curve, model)
    testCases.print("EUROPEAN BLACK PAY VALUE ZERO VOL:", valuePay)

    valueRec = europeanSwaptionRec.value(settle_date, libor_curve, model)
    testCases.print("EUROPEAN BLACK REC VALUE ZERO VOL:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    testCases.banner("======= 20%% BLACK VOLATILITY ========")

    model = Black(0.20)
    testCases.print("Black Model", model._volatility)

    valuePay = europeanSwaptionPay.value(settle_date, libor_curve, model)
    testCases.print("EUROPEAN BLACK PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(settle_date, libor_curve, model)
    testCases.print("EUROPEAN BLACK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BK MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    testCases.banner("=======================================================")
    testCases.banner("=======================================================")
    testCases.banner("==================== BK MODEL =========================")
    testCases.banner("=======================================================")
    testCases.banner("=======================================================")

    testCases.banner("======= 0% VOLATILITY EUROPEAN SWAPTION BK MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.000000001
    a = 0.01
    num_time_steps = 100
    model = BKTree(sigma, a, num_time_steps)

    valuePay = europeanSwaptionPay.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN BK PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    testCases.banner(
        "======= 20% VOLATILITY EUROPEAN SWAPTION BK MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    testCases.banner("BK MODEL SWAPTION CLASS EUROPEAN EXERCISE")

    valuePay = europeanSwaptionPay.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN BK PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    ###########################################################################

    # Now we create the Bermudan swaptions but only allow European exercise
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_rec = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    testCases.banner(
        "======= 0% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BK MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BK PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    testCases.banner(
        "======= 20% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BK MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.2
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BK PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    ###########################################################################
    # Now we create the Bermudan swaptions but allow Bermudan exercise
    ###########################################################################

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_rec = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    testCases.banner(
        "======= ZERO VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BK MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BK PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    testCases.banner(
        "======= 20% VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BK MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BK PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BDT MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    testCases.banner("=======================================================")
    testCases.banner("=======================================================")
    testCases.banner("======================= BDT MODEL =====================")
    testCases.banner("=======================================================")
    testCases.banner("=======================================================")

    testCases.banner("====== 0% VOLATILITY EUROPEAN SWAPTION BDT MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.00001
    num_time_steps = 200
    model = BDTTree(sigma, num_time_steps)

    valuePay = europeanSwaptionPay.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN BDT PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    testCases.banner("===== 20% VOLATILITY EUROPEAN SWAPTION BDT MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = BDTTree(sigma, num_time_steps)

    testCases.banner("BDT MODEL SWAPTION CLASS EUROPEAN EXERCISE")

    valuePay = europeanSwaptionPay.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN BDT PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    ###########################################################################

    # Now we create the Bermudan swaptions but only allow European exercise
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    bermudan_swaption_rec = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    testCases.banner(
        "======= 0% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BDT MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    model = BDTTree(sigma, num_time_steps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    testCases.banner(
        "======= 20% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BDT MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.2
    model = BDTTree(sigma, num_time_steps)

    testCases.banner("BDT MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    ###########################################################################
    # Now we create the Bermudan swaptions but allow Bermudan exercise
    ###########################################################################

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    bermudan_swaption_rec = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    testCases.banner(
        "======= ZERO VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BDT MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = BDTTree(sigma, num_time_steps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    testCases.banner(
        "======= 20% VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BDT MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = BDTTree(sigma, num_time_steps)

#    print("BDT MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BDT MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    testCases.banner("=======================================================")
    testCases.banner("=======================================================")
    testCases.banner("======================= HW MODEL ======================")
    testCases.banner("=======================================================")
    testCases.banner("=======================================================")

    testCases.banner("====== 0% VOLATILITY EUROPEAN SWAPTION HW MODEL ======")

    sigma = 0.0000001
    a = 0.1
    num_time_steps = 200
    model = HWTree(sigma, a, num_time_steps)

    valuePay = europeanSwaptionPay.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN HW PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN HW REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    testCases.banner("===== 20% VOLATILITY EUROPEAN SWAPTION BDT MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.01
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    testCases.banner("HW MODEL SWAPTION CLASS EUROPEAN EXERCISE")

    valuePay = europeanSwaptionPay.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN HW PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(value_date, libor_curve, model)
    testCases.print("EUROPEAN HW REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    ###########################################################################

    # Now we create the Bermudan swaptions but only allow European exercise
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    bermudan_swaption_rec = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    testCases.banner(
        "======= 0% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE HW MODEL ========")

    sigma = 0.000001
    model = HWTree(sigma, a, num_time_steps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    testCases.banner(
        "======= 100bp VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE HW MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.01
    model = HWTree(sigma, a, num_time_steps)

    testCases.banner("BDT MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    ###########################################################################
    # Now we create the Bermudan swaptions but allow Bermudan exercise
    ###########################################################################

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    bermudan_swaption_rec = IborBermudanSwaption(settle_date,
                                                 exercise_date,
                                                 swap_maturity_date,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swap_fixed_coupon,
                                                 swap_fixed_frequency_type,
                                                 swapFixedDayCountType)

    testCases.banner(
        "======= ZERO VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE HW MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    testCases.banner("HW MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN HW PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN HW REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)

    testCases.banner(
        "======= 100bps VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE HW MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.01
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    testCases.banner("HW MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudan_swaption_pay.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN HW PAY VALUE:", valuePay)

    valueRec = bermudan_swaption_rec.value(value_date, libor_curve, model)
    testCases.print("BERMUDAN HW REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAY MINUS RECEIVER :", payRec)
##########################################################################


test_IborBermudanSwaptionBKModel()

testCases.compareTestCases()
