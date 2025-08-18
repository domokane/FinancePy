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
from financepy.products.rates.ibor_bermudan_swaption import IborBermudanSwaption
from financepy.models.black import Black
from financepy.models.bk_tree import BKTree
from financepy.models.hw_tree import HWTree
from financepy.models.bdt_tree import BDTTree
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from FinTestCases import FinTestCases, globalTestCaseMode


test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_IborBermudanSwaptionBKModel():
    """Replicate examples in paper by Leif Andersen which can be found at
    file:///C:/Users/Dominic/Downloads/SSRN-id155208.pdf"""

    value_dt = Date(1, 1, 2011)
    settle_dt = value_dt
    exercise_dt = settle_dt.add_years(1)
    swap_maturity_dt = settle_dt.add_years(4)

    swap_fixed_cpn = 0.060
    swap_fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    swap_fixed_day_count_type = DayCountTypes.ACT_365F

    libor_curve = DiscountCurveFlat(
        value_dt, 0.0625, FrequencyTypes.SEMI_ANNUAL, DayCountTypes.ACT_365F
    )

    fwd_pay_swap = IborSwap(
        exercise_dt,
        swap_maturity_dt,
        SwapTypes.PAY,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fwd_swap_value = fwd_pay_swap.value(settle_dt, libor_curve, libor_curve)

    test_cases.header("LABEL", "VALUE")
    test_cases.print("FWD SWAP VALUE", fwd_swap_value)

    # fwd_pay_swap.print_fixed_leg_pv()

    # Now we create the European swaptions
    fixed_leg_type = SwapTypes.PAY
    european_swaption_pay = IborSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    european_swaption_rec = IborSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BLACK'S MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    test_cases.banner("======= ZERO VOLATILITY ========")
    model = Black(0.0000001)
    test_cases.print("Black Model", model.volatility)

    value_pay = european_swaption_pay.value(settle_dt, libor_curve, model)
    test_cases.print("EUROPEAN BLACK PAY VALUE ZERO VOL:", value_pay)

    value_rec = european_swaption_rec.value(settle_dt, libor_curve, model)
    test_cases.print("EUROPEAN BLACK REC VALUE ZERO VOL:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    test_cases.banner("======= 20%% BLACK VOLATILITY ========")

    model = Black(0.20)
    test_cases.print("Black Model", model.volatility)

    value_pay = european_swaption_pay.value(settle_dt, libor_curve, model)
    test_cases.print("EUROPEAN BLACK PAY VALUE:", value_pay)

    value_rec = european_swaption_rec.value(settle_dt, libor_curve, model)
    test_cases.print("EUROPEAN BLACK REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BK MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    test_cases.banner("=======================================================")
    test_cases.banner("=======================================================")
    test_cases.banner("==================== BK MODEL =========================")
    test_cases.banner("=======================================================")
    test_cases.banner("=======================================================")

    test_cases.banner("======= 0% VOLATILITY EUROPEAN SWAPTION BK MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.000000001
    a = 0.01
    num_time_steps = 100
    model = BKTree(sigma, a, num_time_steps)

    value_pay = european_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN BK PAY VALUE:", value_pay)

    value_rec = european_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN BK REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    test_cases.banner("======= 20% VOLATILITY EUROPEAN SWAPTION BK MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    test_cases.banner("BK MODEL SWAPTION CLASS EUROPEAN EXERCISE")

    value_pay = european_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN BK PAY VALUE:", value_pay)

    value_rec = european_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN BK REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    ###########################################################################

    # Now we create the Bermudan swaptions but only allow European exercise
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    test_cases.banner(
        "======= 0% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BK MODEL ========"
    )

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    test_cases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BK PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BK REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    test_cases.banner(
        "======= 20% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BK MODEL ========"
    )

    # Used BK with constant short-rate volatility
    sigma = 0.2
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    test_cases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BK PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BK REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    ###########################################################################
    # Now we create the Bermudan swaptions but allow Bermudan exercise
    ###########################################################################

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    test_cases.banner(
        "======= ZERO VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BK MODEL ========"
    )

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    test_cases.banner("BK MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BK PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BK REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    test_cases.banner(
        "======= 20% VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BK MODEL ========"
    )

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    test_cases.banner("BK MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BK PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BK REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BDT MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    test_cases.banner("=======================================================")
    test_cases.banner("=======================================================")
    test_cases.banner("======================= BDT MODEL =====================")
    test_cases.banner("=======================================================")
    test_cases.banner("=======================================================")

    test_cases.banner("====== 0% VOLATILITY EUROPEAN SWAPTION BDT MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.00001
    num_time_steps = 200
    model = BDTTree(sigma, num_time_steps)

    value_pay = european_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN BDT PAY VALUE:", value_pay)

    value_rec = european_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN BDT REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    test_cases.banner("===== 20% VOLATILITY EUROPEAN SWAPTION BDT MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = BDTTree(sigma, num_time_steps)

    test_cases.banner("BDT MODEL SWAPTION CLASS EUROPEAN EXERCISE")

    value_pay = european_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN BDT PAY VALUE:", value_pay)

    value_rec = european_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN BDT REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    ###########################################################################

    # Now we create the Bermudan swaptions but only allow European exercise
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    test_cases.banner(
        "======= 0% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BDT MODEL ========"
    )

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    model = BDTTree(sigma, num_time_steps)

    test_cases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    test_cases.banner(
        "======= 20% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BDT MODEL ========"
    )

    # Used BK with constant short-rate volatility
    sigma = 0.2
    model = BDTTree(sigma, num_time_steps)

    test_cases.banner("BDT MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    ###########################################################################
    # Now we create the Bermudan swaptions but allow Bermudan exercise
    ###########################################################################

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    test_cases.banner(
        "======= ZERO VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BDT MODEL ========"
    )

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = BDTTree(sigma, num_time_steps)

    test_cases.banner("BK MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    test_cases.banner(
        "======= 20% VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BDT MODEL ========"
    )

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = BDTTree(sigma, num_time_steps)

    #    print("BDT MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BDT MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    test_cases.banner("=======================================================")
    test_cases.banner("=======================================================")
    test_cases.banner("======================= HW MODEL ======================")
    test_cases.banner("=======================================================")
    test_cases.banner("=======================================================")

    test_cases.banner("====== 0% VOLATILITY EUROPEAN SWAPTION HW MODEL ======")

    sigma = 0.0000001
    a = 0.1
    num_time_steps = 200
    model = HWTree(sigma, a, num_time_steps)

    value_pay = european_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN HW PAY VALUE:", value_pay)

    value_rec = european_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN HW REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    test_cases.banner("===== 20% VOLATILITY EUROPEAN SWAPTION BDT MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.01
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    test_cases.banner("HW MODEL SWAPTION CLASS EUROPEAN EXERCISE")

    value_pay = european_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN HW PAY VALUE:", value_pay)

    value_rec = european_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("EUROPEAN HW REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    ###########################################################################

    # Now we create the Bermudan swaptions but only allow European exercise
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    test_cases.banner(
        "======= 0% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE HW MODEL ========"
    )

    sigma = 0.000001
    model = HWTree(sigma, a, num_time_steps)

    test_cases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    test_cases.banner(
        "======= 100bp VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE HW MODEL ========"
    )

    # Used BK with constant short-rate volatility
    sigma = 0.01
    model = HWTree(sigma, a, num_time_steps)

    test_cases.banner("BDT MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN BDT REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    ###########################################################################
    # Now we create the Bermudan swaptions but allow Bermudan exercise
    ###########################################################################

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    test_cases.banner(
        "======= ZERO VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE HW MODEL ========"
    )

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    test_cases.banner("HW MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN HW PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN HW REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)

    test_cases.banner(
        "======= 100bps VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE HW MODEL ========"
    )

    # Used BK with constant short-rate volatility
    sigma = 0.01
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    test_cases.banner("HW MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN HW PAY VALUE:", value_pay)

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    test_cases.print("BERMUDAN HW REC VALUE:", value_rec)

    pay_rec = value_pay - value_rec
    test_cases.print("PAY MINUS RECEIVER :", pay_rec)


##########################################################################


test_IborBermudanSwaptionBKModel()

test_cases.compareTestCases()
