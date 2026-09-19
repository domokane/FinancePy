# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import pytest

from financepy.utils.date_format import set_date_format, DateFormatTypes
from financepy.utils.date import Date
from financepy.utils.error import FinError
from financepy.products.rates.ibor_future import IborFuture

set_date_format(DateFormatTypes.UK_LONG)

########################################################################################


def test_fin_ibor_future():

    today_date = Date(5, 5, 2020)

    expected_dates = [
        (1, Date(17, 6, 2020), Date(16, 9, 2020)),
        (4, Date(17, 3, 2021), Date(16, 6, 2021)),
        (7, Date(15, 12, 2021), Date(16, 3, 2022)),
        (10, Date(21, 9, 2022), Date(21, 12, 2022)),
    ]

    for future_number, delivery_dt, end_dt in expected_dates:
        fut = IborFuture(today_date, future_number, "3M")
        fra = fut.to_fra(97.50, 0.0)

        assert fut.future_tenor == "3M"
        assert fut.delivery_dt == delivery_dt
        assert fut.end_of_interest_period == end_dt
        assert fra.start_dt == delivery_dt
        assert fra.maturity_dt == end_dt
        assert fra.fra_rate == pytest.approx(0.0250)


def test_fin_ibor_future_one_month_contracts():

    today_date = Date(5, 5, 2020)

    expected_dates = [
        (1, Date(20, 5, 2020), Date(17, 6, 2020)),
        (2, Date(17, 6, 2020), Date(15, 7, 2020)),
        (3, Date(15, 7, 2020), Date(19, 8, 2020)),
        (4, Date(19, 8, 2020), Date(16, 9, 2020)),
    ]

    for future_number, delivery_dt, end_dt in expected_dates:
        fut = IborFuture(today_date, future_number, "1m")
        fra = fut.to_fra(97.50, 0.0)

        assert fut.future_tenor == "1M"
        assert fut.delivery_dt == delivery_dt
        assert fut.end_of_interest_period == end_dt
        assert fra.start_dt == delivery_dt
        assert fra.maturity_dt == end_dt


def test_fin_ibor_future_settlement_amount():

    today_date = Date(5, 5, 2020)
    fut = IborFuture(today_date, 1, "3M")

    assert fut.last_trading_dt == Date(15, 6, 2020)
    assert fut.accrual_factor() == pytest.approx(91.0 / 360.0)
    assert fut.basis_point_value() == pytest.approx(25.2777777778)
    assert fut.settlement_amount(97.50, 97.51) == pytest.approx(
        fut.basis_point_value()
    )
    assert fut.settlement_amount(97.51, 97.50) == pytest.approx(
        -fut.basis_point_value()
    )
    assert fut.settlement_amount(97.50, 97.51, 3.0) == pytest.approx(
        3.0 * fut.basis_point_value()
    )


def test_fin_ibor_future_invalid_inputs():

    today_date = Date(5, 5, 2020)

    with pytest.raises(FinError):
        IborFuture(today_date, 0, "3M")

    with pytest.raises(FinError):
        IborFuture(today_date, 1, "2M")
