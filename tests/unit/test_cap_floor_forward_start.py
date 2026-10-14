import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.bachelier import Bachelier
from financepy.models.black import Black
from financepy.products.rates.ibor_cap_floor import IborCapFloor
from financepy.utils.calendar import CalendarTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import CapFloorTypes

# A cap or floor starting a year after the valuation date: every optionlet, the first included,
# expires when its rate fixes, measured from the valuation date. Reference values are from
# QuantLib's BlackCapFloorEngine and BachelierCapFloorEngine on the same flat curve, with an
# index fixing on the start of each period.

VALUE_DT = Date(1, 1, 2025)
START_DT = Date(1, 1, 2026)
MATURITY_DT = Date(1, 1, 2028)


@pytest.mark.parametrize(
    "opt_type, strike, model, expected",
    [
        (CapFloorTypes.CAP, 0.04, Black(0.20), 13322.88),
        (CapFloorTypes.CAP, 0.04, Bachelier(0.01), 14772.78),
        (CapFloorTypes.FLOOR, 0.05, Black(0.20), 15122.23),
        (CapFloorTypes.FLOOR, 0.05, Bachelier(0.01), 15647.21),
    ],
)
def test_forward_starting_cap_floor_expiries_from_value_date(opt_type, strike, model, expected):
    curve = FlatDiscountCurve(VALUE_DT, 0.045)
    cap_floor = IborCapFloor(
        START_DT,
        MATURITY_DT,
        opt_type,
        strike,
        freq_type=FrequencyTypes.QUARTERLY,
        accrual_dc_type=DayCountTypes.ACT_360,
        cal_type=CalendarTypes.NONE,
    )

    value = cap_floor.value(VALUE_DT, curve, model)

    assert value == pytest.approx(expected, abs=0.05)
