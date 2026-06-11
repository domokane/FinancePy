from financepy.market.curves import InterpTypes
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.utils import (
    BusDayAdjustTypes,
    Date,
    DayCountTypes,
    FrequencyTypes,
    SwapTypes,
)


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
