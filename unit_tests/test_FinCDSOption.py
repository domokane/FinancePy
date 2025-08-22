########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import sys, os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .helpers import build_full_issuer_curve
from financepy.utils.date import Date
from financepy.products.credit.cds_option import CDSOption


# This reproduces example on page 38 of Open Gamma note on CDS Option
trade_dt = Date(5, 2, 2014)
_, issuer_curve = build_full_issuer_curve(trade_dt)
step_in_dt = trade_dt.add_days(1)
value_dt = step_in_dt
expiry_dt = Date(20, 3, 2014)
maturity_dt = Date(20, 6, 2019)
notional = 100.0


def test_cds_option():
    volatility = 0.3

    strike_result = [(100, 4.0007), (150, 1.5874), (200, 0.0955), (300, 0.0)]

    for strike, result in strike_result:
        cdsOption = CDSOption(expiry_dt, maturity_dt, strike / 10000.0, notional)

        v = cdsOption.value(value_dt, issuer_curve, volatility)

        vol = cdsOption.implied_volatility(value_dt, issuer_curve, v)

        assert round(v, 4) == result
        assert vol == 0.3
