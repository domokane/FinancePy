###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from helpers import buildFullIssuerCurve
from financepy.utils.date import Date
from financepy.products.credit.cds_option import CDSOption


# This reproduces example on page 38 of Open Gamma note on CDS Option
tradeDate = Date(5, 2, 2014)
_, issuer_curve = buildFullIssuerCurve(tradeDate)
step_in_date = tradeDate.add_days(1)
valuation_date = step_in_date
expiry_date = Date(20, 3, 2014)
maturity_date = Date(20, 6, 2019)
notional = 100.0


def test_cds_option():
    volatility = 0.3

    strike_result = [
        (100, 4.0007),
        (150, 1.5874),
        (200, 0.0956),
        (300, 0.0)
    ]

    for strike, result in strike_result:
        cdsOption = CDSOption(expiry_date,
                              maturity_date,
                              strike / 10000.0,
                              notional)

        v = cdsOption.value(valuation_date,
                            issuer_curve,
                            volatility)

        vol = cdsOption.implied_volatility(valuation_date,
                                           issuer_curve,
                                           v)

        assert round(v, 4) == result
        assert vol == 0.3
