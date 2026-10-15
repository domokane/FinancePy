import numpy as np
import pytest
from scipy.integrate import quad
from financepy.models.bs_power_option import PowerOption, power_option_value
from financepy.utils.global_types import OptionTypes
from financepy.utils.error import FinError


@pytest.mark.parametrize("power", [0.0, 0.5, 1.0, 2.0, -1.0])
@pytest.mark.parametrize("kind", [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT])
def test_power_option_matches_independent_gauss_hermite_expectation(power,kind):
    spot, strike, t, rate, carry, volatility = 1.14, 1.2, 1.5, 0.045, 0.022, 0.23
    def discounted_payoff(z):
        terminal = spot * np.exp((rate - carry - .5*volatility**2)*t + volatility*np.sqrt(t)*z)
        transformed = terminal**power
        payoff = (max(transformed-strike,0) if kind == OptionTypes.EUROPEAN_CALL
                  else max(strike-transformed,0))
        return np.exp(-rate*t-z*z/2)*payoff/np.sqrt(2*np.pi)
    independent = quad(discounted_payoff, -12,12,epsabs=1e-11,epsrel=1e-11,
                       limit=300)[0]
    analytic = power_option_value(spot,strike,power,t,rate,carry,volatility,kind)
    assert analytic == pytest.approx(independent, abs=1e-8)


@pytest.mark.parametrize("power", [-1.0,0.5,1.0,2.0])
@pytest.mark.parametrize("kind", [OptionTypes.EUROPEAN_CALL,OptionTypes.EUROPEAN_PUT])
def test_power_option_delta_matches_finite_difference(power,kind):
    option=PowerOption(1.2,power,kind)
    spot,h=1.14,1e-5
    numerical=(option.value(spot+h,1.5,.045,.022,.23)
               -option.value(spot-h,1.5,.045,.022,.23))/(2*h)
    assert option.delta(spot,1.5,.045,.022,.23)==pytest.approx(numerical,rel=1e-6,abs=1e-7)


def test_power_one_matches_standard_stock_forward_parity():
    call=PowerOption(100.,1.,OptionTypes.EUROPEAN_CALL).value(105,1,.04,.02,.2)
    put=PowerOption(100.,1.,OptionTypes.EUROPEAN_PUT).value(105,1,.04,.02,.2)
    assert call-put==pytest.approx(105*np.exp(-.02)-100*np.exp(-.04))


def test_zero_volatility_and_expiry():
    call=PowerOption(110.,2,OptionTypes.EUROPEAN_CALL)
    assert call.value(10,0,0,0,.2)==0
    assert call.value(12,0,0,0,.2)==34
    assert call.value(12,1,0,0,0)==34


@pytest.mark.parametrize("kwargs",[
    dict(spot=0.,strike=100.,power=2.,maturity=1.,rate=.03,carry=0.,volatility=.2),
    dict(spot=10.,strike=-1.,power=2.,maturity=1.,rate=.03,carry=0.,volatility=.2),
    dict(spot=10.,strike=1.,power=2.,maturity=-1.,rate=.03,carry=0.,volatility=.2),
])
def test_reject_invalid_parameters(kwargs):
    with pytest.raises(FinError):
        power_option_value(**kwargs,opt_type=OptionTypes.EUROPEAN_CALL)

