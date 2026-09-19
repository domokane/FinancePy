from financepy.market.curves import InterpTypes
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.utils.math import ONE_MILLION
from financepy.products.rates.ibor_single_curve import IborSingleCurve

from .helpers import build_ibor_single_curve

########################################################################################


def test_ibor_swap_end_of_month_aligns_coupon_grid():
    valuation_date = Date(31, 5, 2023)

    swaps = []
    for tenor in ["6M", "1Y"]:
        swaps.append(
            IborSwap(
                effective_dt=valuation_date,
                term_dt_or_tenor=tenor,
                fixed_leg_type=SwapTypes.PAY,
                fixed_cpn=0.01,
                fixed_freq_type=FrequencyTypes.QUARTERLY,
                fixed_dc_type=DayCountTypes.ACT_360,
                float_dc_type=DayCountTypes.ACT_360,
                bd_type=BusDayAdjustTypes.NONE,
                end_of_month=True,
            )
        )

    assert swaps[0].fixed_leg.payment_dts == [
        Date(31, 8, 2023),
        Date(30, 11, 2023),
    ]
    assert swaps[1].fixed_leg.payment_dts[:2] == swaps[0].fixed_leg.payment_dts

    IborSingleCurve(
        value_dt=valuation_date,
        ibor_deposits=[],
        ibor_fras=[],
        ibor_swaps=swaps,
        interp_type=InterpTypes.FLAT_FWD_RATES,
    )


def test_ibor_swap_end_of_month_handles_non_eom_effective_date():
    swap = IborSwap(
        start_dt,
        end_dt,
        fixed_leg_type,
        fixed_coupon,
        fixed_freq_type,
        fixed_dc_type,
        notional,
        float_spread,
        float_freq_type,
        float_dc_type,
        swap_calendar_type,
        bus_day_adjust_type,
        date_gen_rule_type,
    )

    """ Now perform a valuation after the swap has seasoned but with the
    same curve being used for discounting and working out the implied
    future Libor rates. """

    value_date = Date(30, 11, 2018)
    settle_dt = value_date.add_days(2)
    libor_curve = build_ibor_single_curve(value_date)
    return first_fixing, swap, settle_dt, libor_curve


########################################################################################


def test_libor_swap():

    # I have tried to reproduce the example from the blog by Ioannis Rigopoulos
    # https://blog.deriscope.com/index.php/en/excel-interest-rate-swap-price-dual-bootstrapping-curve

    start_dt = Date(27, 12, 2017)
    end_dt = Date(27, 12, 2067)

    # The swap is a long dated fixed vs Euribor interest rate swap that runs
    # until Dec 27, 2067. The exact details are shown below, along with its
    # Bloomberg valuation of 388,147.49 EUR as of Nov 30, 2018:

    first_fixing, swap, settle_dt, libor_curve = _load_test_swap_and_curve(start_dt, end_dt)

    v = swap.value(settle_dt, libor_curve, libor_curve, first_fixing)

    assert round(v, 0) == 392684.0


########################################################################################


def test_libor_swap_cashflow_report():

    # as test_LiborSwap but with extra output
    start_dt = Date(27, 12, 2017)
    end_dt = Date(27, 12, 2067)

    first_fixing, swap, settle_dt, libor_curve = _load_test_swap_and_curve(start_dt, end_dt)

    print(swap)

    v = swap.value(settle_dt, libor_curve, libor_curve, first_fixing, pv_only=False)

    sum_of_cashflows = v[1]["payment_pv"].sum()
    assert round(sum_of_cashflows - v[0], 4) == 0


########################################################################################


def test_dp_example():

    #  http://www.derivativepricing.com/blogpage.asp?id=8

    start_dt = Date(14, 11, 2011)
    end_dt = Date(14, 11, 2016)
    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    swap_cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    fixed_dc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_leg_type = SwapTypes.PAY
    fixed_cpn = 0.0124
    notional = ONE_MILLION

    swap = IborSwap(
        start_dt,
        end_dt,
        fixed_leg_type,
        fixed_cpn=fixed_cpn,
        fixed_freq_type=fixed_freq_type,
        fixed_dc_type=fixed_dc_type,
        float_freq_type=FrequencyTypes.SEMI_ANNUAL,
        float_dc_type=DayCountTypes.ACT_360,
        bd_type=BusDayAdjustTypes.NONE,
        end_of_month=True,
    )

    assert swap.fixed_leg.payment_dts == [
        Date(31, 5, 2023),
        Date(31, 8, 2023),
        Date(30, 11, 2023),
        Date(29, 2, 2024),
        Date(15, 5, 2024),
    ]


def test_ibor_swap_end_of_month_applies_to_float_leg_schedule():
    swap = IborSwap(
        effective_dt=Date(31, 5, 2023),
        term_dt_or_tenor="1Y",
        fixed_leg_type=SwapTypes.PAY,
        fixed_cpn=0.01,
        fixed_freq_type=FrequencyTypes.ANNUAL,
        fixed_dc_type=DayCountTypes.ACT_360,
        float_freq_type=FrequencyTypes.QUARTERLY,
        float_dc_type=DayCountTypes.ACT_360,
        bd_type=BusDayAdjustTypes.NONE,
        end_of_month=True,
    )

    assert swap.float_leg.payment_dts == [
        Date(31, 8, 2023),
        Date(30, 11, 2023),
        Date(29, 2, 2024),
        Date(31, 5, 2024),
    ]
    assert swap.float_leg.start_accrued_dts == [
        Date(31, 5, 2023),
        Date(31, 8, 2023),
        Date(30, 11, 2023),
        Date(29, 2, 2024),
    ]

    value_dt = start_dt

    curve = DiscountCurve(value_dt, dts, np.array(dfs), InterpTypes.FLAT_FWD_RATES)

    v = swap.value(value_dt, curve, curve)

    # This is essentially zero
    assert round(v, 1) == 15.7


###################################################################################


def test_ibor_swap_end_of_month_aligns_coupon_grid():
    valuation_date = Date(31, 5, 2023)

    swaps = []
    for tenor in ["6M", "1Y"]:
        swaps.append(
            IborSwap(
                effective_dt=valuation_date,
                term_dt_or_tenor=tenor,
                fixed_leg_type=SwapTypes.PAY,
                fixed_cpn=0.01,
                fixed_freq_type=FrequencyTypes.QUARTERLY,
                fixed_dc_type=DayCountTypes.ACT_360,
                float_dc_type=DayCountTypes.ACT_360,
                bd_type=BusDayAdjustTypes.NONE,
                end_of_month=True,
            )
        )

    assert swaps[0].fixed_leg.payment_dts == [
        Date(31, 8, 2023),
        Date(30, 11, 2023),
    ]
    assert swaps[1].fixed_leg.payment_dts[:2] == swaps[0].fixed_leg.payment_dts

    IborSingleCurve(
        value_dt=valuation_date,
        ibor_deposits=[],
        ibor_fras=[],
        ibor_swaps=swaps,
        interp_type=InterpTypes.FLAT_FWD_RATES,
    )


def test_ibor_swap_end_of_month_handles_non_eom_effective_date():
    swap = IborSwap(
        effective_dt=Date(15, 5, 2023),
        term_dt_or_tenor="1Y",
        fixed_leg_type=SwapTypes.PAY,
        fixed_cpn=0.01,
        fixed_freq_type=FrequencyTypes.QUARTERLY,
        fixed_dc_type=DayCountTypes.ACT_360,
        float_dc_type=DayCountTypes.ACT_360,
        bd_type=BusDayAdjustTypes.NONE,
        end_of_month=True,
    )

    assert swap.fixed_leg.payment_dts == [
        Date(31, 5, 2023),
        Date(31, 8, 2023),
        Date(30, 11, 2023),
        Date(29, 2, 2024),
        Date(15, 5, 2024),
    ]


def test_ibor_swap_end_of_month_applies_to_float_leg_schedule():
    swap = IborSwap(
        effective_dt=Date(31, 5, 2023),
        term_dt_or_tenor="1Y",
        fixed_leg_type=SwapTypes.PAY,
        fixed_cpn=0.01,
        fixed_freq_type=FrequencyTypes.ANNUAL,
        fixed_dc_type=DayCountTypes.ACT_360,
        float_freq_type=FrequencyTypes.QUARTERLY,
        float_dc_type=DayCountTypes.ACT_360,
        bd_type=BusDayAdjustTypes.NONE,
        end_of_month=True,
    )

    assert swap.float_leg.payment_dts == [
        Date(31, 8, 2023),
        Date(30, 11, 2023),
        Date(29, 2, 2024),
        Date(31, 5, 2024),
    ]
    assert swap.float_leg.start_accrued_dts == [
        Date(31, 5, 2023),
        Date(31, 8, 2023),
        Date(30, 11, 2023),
        Date(29, 2, 2024),
    ]


########################################################################################

# if __name__ == '__main__':
#     test_libor_swap_cashflow_report()

test_libor_swap()
