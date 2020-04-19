# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from ...finutils.FinError import FinError

# Non exhaustive list of country codes and currency names

ccyNames = {"AED": ("UNITED ARAB EMIRATES", "DIRHAM"),
            "AUD": ("AUSTRALIA","DOLLAR"),
            "BRL": ("BRAZIL", "REAL"),
            "CAD": ("CANADA", "DOLLAR"),
            "CHF": ("SWITZERLAND","FRANC"),
            "CLP": ("CHILE", "PESO"),
            "CNY": ("CHINA", "RENMIMBI"),
            "COP": ("COLUMBIA", "PESO"),
            "DKK": ("DENMARK", "KRONE"),
            "EUR": ("EUROZONE", "EURO"),
            "GBP": ("UK", "POUND"),
            "HKD": ("HONG KONG", "DOLLAR"),
            "HUF": ("HUNGARY", "FORINT"),
            "IDR": ("INDONESIA", "RUPIAH"),
            "INR": ("INDIA", "RUPEE"),
            "ILS": ("ISRAEL", "SHEKEL"),
            "JPY": ("JAPAN"," YEN"),
            "KRW": ("KOREAN", "WON"),
            "MYR": ("MALYSIA", "RINGIT"),
            "MXN": ("MEXICO", "PESO"),
            "NOK": ("NORWAY", "KRONE"),
            "NZD": ("NEW ZEALAND", "DOLLAR"),
            "PHP": ("PHILIPPINES", "PESO"),
            "PLN": ("POLAND", "ZLOTY"),
            "RON": ("ROMANIA", "LEU"),
            "RUB": ("RUSSIA", "RUBLE"),
            "SAR": ("SAUDI ARABIA", "RIYAL"),
            "SEK": ("SWEDEN", "KRONA"),
            "SGD": ("SINGAPORE", "DOLLAR"),
            "THB": ("THAILAND", "BHAT"),
            "TRY": ("TURKEY", "LIRA"),
            "TWD": ("TAIWAN", "DOLLAR"),
            "USD": ("US","DOLLAR"),
            "ZAR": ("SOUTH AFRICA", "RAND")
            }

ccyQuotes = {"EURUSD": (1, 1, 1),
             "USDJPY": (1, 1, 1),
             "EURJPY": (1, 1, 1),
             "GBPUSD": (1, 1, 1),
             "EURGBP": (1, 1, 1),
             "USDCHF": (1, 1, 1),
             "AUDUSD": (1, 1, 1),
             "NZDUSD": (1, 1, 1),
             "USDCAD": (1, 1, 1),
             "EURNOK": (1, 1, 1),
             "EURSEK": (1, 1, 1),
             "EURDKK": (1, 1, 1),
             "EURHUF": (1, 1, 1),
             "EURPLN": (1, 1, 1),
             "USDTRY": (1, 1, 1),
             "USDZAR": (1, 1, 1),
             "USDMXN": (1, 1, 1),
             "USDBRL": (1, 1, 1),
             "USDSGD": (1, 1, 1)}

###############################################################################


class FinFXRate():

    def __init__(self,
                 ccy1,
                 ccy2,
                 rate):

        if ccy1 in ccyNames:
            self._ccy1 = ccy1
        else:
            raise FinError("Unknown currency code ", ccy1)

        if ccy2 in ccyNames:
            self._ccy2 = ccy2
        else:
            raise FinError("Unknown currency code ", ccy2)

        self._ccy2 = ccy2
        self._rate

###############################################################################
