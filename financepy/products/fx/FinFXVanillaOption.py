##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize

from ...finutils.FinDate import FinDate
from ...finutils.FinMath import nprime
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
from ...finutils.FinGlobalTypes import FinOptionTypes
from ...products.fx.FinFXModelTypes import FinFXModel
from ...products.fx.FinFXModelTypes import FinFXModelBlackScholes
from ...products.fx.FinFXModelTypes import FinFXModelSABR
from ...models.FinModelCRRTree import crrTreeValAvg
from ...models.FinModelSABR import blackVolFromSABR
from ...finutils.FinHelperFunctions import checkArgumentTypes

from ...finutils.FinMath import N


###############################################################################


def f(volatility, *args):
    ''' This is the objective function used in the determination of the FX
    Option implied volatility which is computed in the class below. '''

    self = args[0]
    valueDate = args[1]
    spotFXRate = args[2]
    domDiscountCurve = args[3]
    forDiscountCurve = args[4]
    price = args[5]

    model = FinFXModelBlackScholes(volatility)

    vdf = self.value(valueDate,
                     spotFXRate,
                     domDiscountCurve,
                     forDiscountCurve,
                     model)['v']

    objFn = vdf - price

    return objFn

###############################################################################


def fvega(volatility, *args):
    ''' This is the derivative of the objective function with respect to the
    option volatility. It is used to speed up the determination of the FX
    Option implied volatility which is computed in the class below. '''

    self = args[0]
    valueDate = args[1]
    spotFXRate = args[2]
    domDiscountCurve = args[3]
    forDiscountCurve = args[4]

    model = FinFXModelBlackScholes(volatility)

    fprime = self.vega(valueDate,
                       spotFXRate,
                       domDiscountCurve,
                       forDiscountCurve,
                       model)

    return fprime

###############################################################################


def g(K, *args):
    ''' This is the objective function used in the determination of the FX
    Option implied strike which is computed in the class below. '''

    self = args[0]
    valueDate = args[1]
    stockPrice = args[2]
    domDiscountCurve = args[3]
    forDiscountCurve = args[4]
    delta = args[5]
    deltaType = args[6]
    volatility = args[7]

    model = FinFXModelBlackScholes(volatility)

    self._strikeFXRate = K

    deltaDict = self.delta(valueDate,
                           stockPrice,
                           domDiscountCurve,
                           forDiscountCurve,
                           model)

    deltaOut = deltaDict[deltaType]
    objFn = delta - deltaOut
    return objFn

###############################################################################


def solveForStrike(valueDate,
                   vanillaOption,
                   spotFXRate,
                   domDiscountCurve,
                   forDiscountCurve,
                   delta,
                   deltaType,
                   volatility):
    ''' This function determines the implied strike of an FX option
    given a delta and the other option details. It uses a one-dimensional
    Newton root search algorith to determine the strike that matches an
    input volatility. '''

    argtuple = (vanillaOption, valueDate, spotFXRate, domDiscountCurve,
                forDiscountCurve, delta, deltaType, volatility)

    sigma = optimize.newton(g, x0=spotFXRate, args=argtuple,
                            tol=1e-5, maxiter=100, fprime2=None)
    return sigma

###############################################################################
# ALL CCY RATES MUST BE IN NUM UNITS OF DOMESTIC PER UNIT OF FOREIGN CURRENCY
# SO EURUSD = 1.30 MEANS 1.30 DOLLARS PER EURO SO DOLLAR IS THE DOMESTIC AND
# EUR IS THE FOREIGN CURRENCY
###############################################################################


class FinFXVanillaOption():
    ''' This is a class for an FX Option trade. It permits the user to
    calculate the price of an FX Option trade which can be expressed in a
    number of ways depending on the investor or hedger's currency. It aslo
    allows the calculation of the option's delta in a number of forms as
    well as the various Greek risk sensitivies. '''

    def __init__(self,
                 expiryDate: FinDate,
                 strikeFXRate: float,  # 1 unit of foreign in domestic
                 currencyPair: str,  # FORDOM
                 optionType: FinOptionTypes,
                 notional: float,
                 premCurrency: str,
                 spotDays: int = 0):
        ''' Create the FX Vanilla Option object. Inputs include expiry date,
        strike, currency pair, option type (call or put), notional and the
        currency of the notional. And adjustment for spot days is enabled. All
        currency rates must be entered in the price in domestic currency of
        one unit of foreign. And the currency pair should be in the form FORDOM
        where FOR is the foreign currency pair currency code and DOM is the
        same for the domestic currency. '''

        checkArgumentTypes(self.__init__, locals())

        deliveryDate = expiryDate.addWeekDays(spotDays)

        ''' The FX rate the price in domestic currency ccy2 of a single unit
        of the foreign currency which is ccy1. For example EURUSD of 1.3 is the
        price in USD (CCY2) of 1 unit of EUR (CCY1)'''

        if deliveryDate < expiryDate:
            raise FinError("Delivery date must be on or after expiry date.")

        if len(currencyPair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self._expiryDate = expiryDate
        self._deliveryDate = deliveryDate

        if strikeFXRate < 0.0:
            raise FinError("Negative strike.")

        self._strikeFXRate = strikeFXRate

        self._currencyPair = currencyPair
        self._forName = self._currencyPair[0:3]
        self._domName = self._currencyPair[3:6]

        if premCurrency != self._domName and premCurrency != self._forName:
            raise FinError("Premium currency not in currency pair.")

        self._premCurrency = premCurrency

        self._notional = notional

        if optionType != FinOptionTypes.EUROPEAN_CALL and \
           optionType != FinOptionTypes.EUROPEAN_PUT and\
           optionType != FinOptionTypes.AMERICAN_CALL and \
           optionType != FinOptionTypes.AMERICAN_PUT:
            raise FinError("Unknown Option Type:" + optionType)

        self._optionType = optionType
        self._spotDays = spotDays

###############################################################################

    def value(self,
              valueDate,
              spotFXRate,  # 1 unit of foreign in domestic
              domDiscountCurve,
              forDiscountCurve,
              model):
        ''' This function calculates the value of the option using a specified
        model with the resulting value being in domestic i.e. ccy2 terms.
        Recall that Domestic = CCY2 and Foreign = CCY1 and FX rate is in
        price in domestic of one unit of foreign currency. '''

        if type(valueDate) == FinDate:
            spotDate = valueDate.addWeekDays(self._spotDays)
            tdel = (self._deliveryDate - spotDate) / gDaysInYear
            texp = (self._expiryDate - valueDate) / gDaysInYear
        else:
            tdel = valueDate
            texp = tdel

        if spotFXRate <= 0.0:
            raise FinError("spotFXRate must be greater than zero.")

#        if np.any(spotFXRate <= 0.0):
#            raise FinError("spotFXRate must be greater than zero.")

        if model._parentType != FinFXModel:
            raise FinError("Model is not inherited off type FinFXModel.")

#        if np.any(tdel < 0.0):
#            raise FinError("Time to expiry must be positive.")
    
        if tdel < 0.0:
            raise FinError("Time to expiry must be positive.")

        tdel = np.maximum(tdel, 1e-10)

        rd = domDiscountCurve.zeroRate(self._deliveryDate)
        rf = forDiscountCurve.zeroRate(self._deliveryDate)

        S0 = spotFXRate
        K = self._strikeFXRate
        F0T = S0 * np.exp((rd-rf)*tdel)

        if type(model) == FinFXModelBlackScholes or \
           type(model) == FinFXModelSABR:

            if type(model) == FinFXModelBlackScholes:
                volatility = model._volatility
            elif type(model) == FinFXModelSABR:
                volatility = blackVolFromSABR(model.alpha,
                                              model.beta,
                                              model.rho,
                                              model.nu,
                                              F0T, K, tdel)

#            if np.any(volatility) < 0.0:
#                raise FinError("Volatility should not be negative.")

            if volatility < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(S0/K)
            sqrtT = np.sqrt(texp)
            den = volatility * sqrtT
            mu = rd - rf
            v2 = volatility * volatility

            d1 = (lnS0k + mu * tdel + v2 * texp / 2.0) / den
            d2 = (lnS0k + mu * tdel - v2 * texp / 2.0) / den

            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                vdf = np.exp(-rd*tdel) * (F0T*N(d1) - K*N(d2))
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                vdf = -np.exp(-rd*tdel) * (F0T*N(-d1) - K*N(-d2))
            elif self._optionType == FinOptionTypes.AMERICAN_CALL:
                numStepsPerYear = 100
                vdf = crrTreeValAvg(S0, rd, rf, volatility, numStepsPerYear,
                          texp, FinOptionTypes.AMERICAN_CALL.value, K)['value']
            elif self._optionType == FinOptionTypes.AMERICAN_PUT:
                numStepsPerYear = 100
                vdf = crrTreeValAvg(S0, rd, rf, volatility, numStepsPerYear,
                           texp, FinOptionTypes.AMERICAN_PUT.value, K)['value']
            else:
                raise FinError("Unknown option type")

        # The option value v is in domestic currency terms but the value of
        # the option may be quoted in either currency terms and so we calculate
        # these

        if self._premCurrency == self._domName:
            notional_dom = self._notional
            notional_for = self._notional / self._strikeFXRate
        elif self._premCurrency == self._forName:
            notional_dom = self._notional * self._strikeFXRate
            notional_for = self._notional
        else:
            raise FinError("Invalid notional currency.")

        vdf = vdf
        pips_dom = vdf
        pips_for = vdf / (spotFXRate * self._strikeFXRate)

        cash_dom = vdf * notional_dom / self._strikeFXRate
        cash_for = vdf * notional_for / spotFXRate

        pct_dom = vdf / self._strikeFXRate
        pct_for = vdf / spotFXRate

        return {'v': vdf,
                "cash_dom": cash_dom,
                "cash_for": cash_for,
                "pips_dom": pips_dom,
                "pips_for": pips_for,
                "pct_dom": pct_dom,
                "pct_for": pct_for,
                "not_dom": notional_dom,
                "not_for": notional_for,
                "ccy_dom": self._domName,
                "ccy_for": self._forName}

###############################################################################

    def delta_bump(self,
                   valueDate,
                   spotFXRate,
                   ccy1DiscountCurve,
                   ccy2DiscountCurve,
                   model):
        ''' Calculation of the FX option delta by bumping the spot FX rate by
        1 cent of its value. This gives the FX spot delta. For speed we prefer
        to use the analytical calculation of the derivative given below. '''

        bump = 0.0001 * spotFXRate

        v = self.value(
            valueDate,
            spotFXRate,
            ccy1DiscountCurve,
            ccy2DiscountCurve,
            model)

        vBumped = self.value(
            valueDate,
            spotFXRate + bump,
            ccy1DiscountCurve,
            ccy2DiscountCurve,
            model)

        if type(vBumped) is dict:
            delta = (vBumped['value'] - v['value']) / bump
        else:
            delta = (vBumped - v) / bump

        return delta

###############################################################################

    def delta(self,
              valueDate,
              spotFXRate,
              domDiscountCurve,
              forDiscountCurve,
              model):
        ''' Calculation of the FX Option delta. There are several definitions
        of delta and so we are required to return a dictionary of values. The
        definitions can be found on Page 44 of Foreign Exchange Option Pricing
        by Iain Clark, published by Wiley Finance. '''

        if type(valueDate) == FinDate:
            spotDate = valueDate.addWeekDays(self._spotDays)
            tdel = (self._deliveryDate - spotDate) / gDaysInYear
            texp = (self._expiryDate - valueDate) / gDaysInYear
        else:
            tdel = valueDate
            texp = tdel

        if np.any(spotFXRate <= 0.0):
            raise FinError("Spot FX Rate must be greater than zero.")

        if model._parentType != FinFXModel:
            raise FinError("Model is not inherited off type FinFXModel.")

        if np.any(tdel < 0.0):
            raise FinError("Time to expiry must be positive.")

        tdel = np.maximum(tdel, 1e-10)

        domDf = domDiscountCurve._df(tdel)
        rd = -np.log(domDf)/tdel

        forDf = forDiscountCurve._df(tdel)
        rf = -np.log(forDf)/tdel

        S0 = spotFXRate
        K = self._strikeFXRate
        F = S0 * np.exp((rd-rf)*tdel)

        if type(model) == FinFXModelBlackScholes:

            volatility = model._volatility

            if np.any(volatility < 0.0):
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(S0/K)
            sqrtT = np.sqrt(texp)
            den = volatility * sqrtT
            mu = rd - rf
            v2 = volatility * volatility

            d1 = (lnS0k + mu * tdel + v2 * texp / 2.0) / den
            d2 = (lnS0k + mu * tdel - v2 * texp / 2.0) / den

            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                w = 1
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                w = -1
            else:
                raise FinError("Unknown option type")

            n1 = N(w*d1)
            n2 = N(w*d2)

            pips_spot_delta = w*np.exp(-rf*tdel)*n1
            pips_fwd_delta = w*n1
            pct_spot_delta_prem_adj = w*np.exp(-rd*tdel)*n2*K/S0
            pct_fwd_delta_prem_adj = w*K*n2/F

        return {"pips_spot_delta": pips_spot_delta,
                "pips_fwd_delta": pips_fwd_delta,
                "pct_spot_delta_prem_adj": pct_spot_delta_prem_adj,
                "pct_fwd_delta_prem_adj": pct_fwd_delta_prem_adj}

###############################################################################

    def gamma(self,
              valueDate,
              spotFXRate,  # value of a unit of foreign in domestic currency
              domDiscountCurve,
              forDiscountCurve,
              model):
        ''' This function calculates the FX Option Gamma using the spot delta.
        '''

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("FX Rate must be greater than zero.")

        if model._parentType != FinFXModel:
            raise FinError("Model is not inherited off type FinFXModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        domDf = domDiscountCurve._df(t)
        rd = -np.log(domDf)/t

        forDf = forDiscountCurve._df(t)
        rf = -np.log(forDf)/t

        K = self._strikeFXRate
        S0 = spotFXRate

        if type(model) == FinFXModelBlackScholes:

            volatility = model._volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(S0 / K)
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = rd - rf
            v2 = volatility * volatility
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den
            gamma = np.exp(-rf * t) * nprime(d1)
            gamma = gamma / S0 / den
        else:
            raise FinError("Unknown Model Type")

        return gamma

###############################################################################

    def vega(self,
             valueDate,
             spotFXRate,  # value of a unit of foreign in domestic currency
             domDiscountCurve,
             forDiscountCurve,
             model):
        ''' This function calculates the FX Option Vega using the spot delta.
        '''

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("Spot FX Rate must be greater than zero.")

        if model._parentType != FinFXModel:
            raise FinError("Model is not inherited off type FinEquityModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        domDf = domDiscountCurve._df(t)
        rd = -np.log(domDf)/t

        forDf = forDiscountCurve._df(t)
        rf = -np.log(forDf)/t

        K = self._strikeFXRate
        S0 = spotFXRate

        if type(model) == FinFXModelBlackScholes:

            volatility = model._volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(S0/K)
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = rd - rf
            v2 = volatility * volatility
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den
            vega = S0 * sqrtT * np.exp(-rf * t) * nprime(d1)
        else:
            raise FinError("Unknown Model type")

        return vega

###############################################################################

    def theta(self,
              valueDate,
              spotFXRate,  # value of a unit of foreign in domestic currency
              domDiscountCurve,
              forDiscountCurve,
              model):
        ''' This function calculates the time decay of the FX option. '''

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("Spot FX Rate must be greater than zero.")

        if model._parentType != FinFXModel:
            raise FinError("Model is not inherited off type FinEquityModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        domDf = domDiscountCurve._df(t)
        rd = -np.log(domDf)/t

        forDf = forDiscountCurve._df(t)
        rf = -np.log(forDf)/t

        K = self._strikeFXRate
        S0 = spotFXRate

        if type(model) == FinFXModelBlackScholes:

            vol = model._volatility

            if np.any(vol) < 0.0:
                raise FinError("Volatility should not be negative.")

            vol = np.maximum(vol, 1e-10)

            lnS0k = np.log(S0/K)
            sqrtT = np.sqrt(t)
            den = vol * sqrtT
            mu = rd - rf
            v2 = vol * vol
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den
            d2 = (lnS0k + (mu - v2 / 2.0) * t) / den

            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                v = - S0 * np.exp(-rf * t) * nprime(d1) * vol / 2.0 / sqrtT
                v = v + rf * S0 * np.exp(-rf * t) * N(d1)
                v = v - rd * K * np.exp(-rd * t) * N(d2)
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                v = - S0 * np.exp(-rf * t) * nprime(d1) * vol / 2.0 / sqrtT
                v = v + rd * K * np.exp(-rd * t) * N(-d2)
                v = v - rf * S0 * np.exp(-rf * t) * N(-d1)
            else:
                raise FinError("Unknown option type")

        else:
            raise FinError("Unknown Model Type")

        return v

###############################################################################

    def impliedVolatility(self,
                          valueDate,
                          stockPrice,
                          discountCurve,
                          dividendYield,
                          price):
        ''' This function determines the implied volatility of an FX option
        given a price and the other option details. It uses a one-dimensional
        Newton root search algorith to determine the implied volatility. '''

        argtuple = (self, valueDate, stockPrice,
                    discountCurve, dividendYield, price)

        sigma = optimize.newton(f, x0=0.2, fprime=fvega, args=argtuple,
                                tol=1e-6, maxiter=50, fprime2=None)
        return sigma

###############################################################################

    def valueMC(self,
                valueDate,
                spotFXRate,
                domDiscountCurve,
                forDiscountCurve,
                model,
                numPaths=10000,
                seed=4242):
        ''' Calculate the value of an FX Option using Monte Carlo methods.
        This function can be used to validate the risk measures calculated
        above or used as the starting code for a model exotic FX product that
        cannot be priced analytically. This function uses Numpy vectorisation
        for speed of execution.'''

        if model._parentType == FinFXModel:
            volatility = model._volatility
        else:
            raise FinError("Model Type invalid")

        np.random.seed(seed)
        t = (self._expiryDate - valueDate) / gDaysInYear

        domDF = domDiscountCurve.df(self._expiryDate)
        forDF = forDiscountCurve.df(self._expiryDate)

        rd = -np.log(domDF)/t
        rf = -np.log(forDF)/t

        mu = rd - rf
        v2 = volatility**2
        K = self._strikeFXRate
        sqrtdt = np.sqrt(t)

        # Use Antithetic variables
        g = np.random.normal(0.0, 1.0, size=(1, numPaths))
        s = spotFXRate * np.exp((mu - v2 / 2.0) * t)
        m = np.exp(g * sqrtdt * volatility)
        s_1 = s * m
        s_2 = s / m

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            payoff_a_1 = np.maximum(s_1 - K, 0.0)
            payoff_a_2 = np.maximum(s_2 - K, 0.0)
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            payoff_a_1 = np.maximum(K - s_1, 0.0)
            payoff_a_2 = np.maximum(K - s_2, 0.0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(payoff_a_1) + np.mean(payoff_a_2)
        v = payoff * np.exp(-rd * t) / 2.0
        return v

###############################################################################
