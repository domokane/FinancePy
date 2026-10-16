#!/usr/bin/env python .

#
# Created on June/2025
# updated on set/2026
# @author: Hernan Lobert
#
# Nota: reinstlar o package a cada vez que altera o codigo fonte conforme abaixo.
#    instala o FinancePy instalado no C:\ com permissao para editar:
#       CD C:\Python\FinancePy
#       uv pip install --system-certs --system --upgrade --editable .
#       uv pip install --system-certs --system --upgrade --editable C:\Python\FinancePy

"""import the packages that we'll use"""

# time
import datetime as dt
import sys

# math libs
import numpy as np
import scipy.optimize as opt

# excel wrapper
import xlwings as xw

# uteis
from . import utils_by_HH as ut

# financial functions
from .market.curves.interpolator import InterpTypes

# CDS functions
from .products.credit.cds import CDS
from .market.curves.cds_curve import CDSCurve

# dates functions
from .utils.date import Date
from .utils.day_count import DayCountTypes
from .utils.frequency import FrequencyTypes, annual_frequency
from .utils.date import from_datetime
from .utils.calendar import CalendarTypes, Calendar

# swap functions
from .products.rates.ibor_swap import SwapTypes

# from .products.rates.ibor_swap import IborSwap, SwapTypes
from .products.rates.ois import OIS
from .products.rates.ois_curve import OISCurve
from .products.rates.ibor_deposit import IborDeposit

# plot
# import matplotlib.pyplot as plt
# Data wrappers
# import pandas as pd
# regressao polinomial
# from numpy.polynomial import Polynomial
# calendario US
# from pandas.tseries.holiday import USFederalHolidayCalendar
# funcoes de DU do pandas
# from pandas.tseries.offsets import CustomBusinessDay
# statistica
# from scipy import stats
# BBG wrapper
# from xbbg import blp


""" global variables """
# curva de juros
curva_OIS = None  # OISCurve
# data e hora da curva de juros
data_curva_juros: Date | None = None  # dt.datetime
# define o modo de debug
verbal_debug = (
    None  # True = mostra todas as informacoes, False = mostra poucas informacoes
)
# matem a linha de debug
linha_debug = 1


def carrega_curva_mercado(linha_debug=1):
    """Carrega a curva de juros do mercado"""
    # dados de verbal
    global verbal_debug
    # nome da aba Debug e inicializa a linha de debug
    linha_debug = ut.print_debug(
        "debug", linha_debug, "### Carregando a curva de juros ###"
    )
    # cria a variavel da planilha
    rc_book = xw.Book.caller()  # xw.Book(planilha)
    # cria a variavel da tab da curva de juros OIS-ISDA
    aba_USD_ISDA_OIS = rc_book.sheets["USD_ISDA_OIS"]
    # checa se eh para mostrar todas as informacoes ou somente as minimas
    if verbal_debug is None:
        verbal_debug = rc_book.sheets["debug"].range("verbal_debug").value
        """if verbal_debug is None:
            verbal_debug = False
        elif not isinstance(
            verbal_debug, bool
        ):  # se nao for um booleano, entao assume False
            verbal_debug = False"""
        # se nao tiver valor ou nao for booleano, entao assume false
        if verbal_debug is None or not isinstance(verbal_debug, bool):
            verbal_debug = False
        # imprime que tipo de debug esta sendo usado
        linha_debug = ut.print_debug(
            "debug", linha_debug, f"   verbal_debug? {verbal_debug}"
        )

    """Calcula a curva de CDS para um unico nome"""
    # data de inicio do swap
    trade_date = from_datetime(
        aba_USD_ISDA_OIS["trade_date"].options(dates=dt.date).value
    )
    settlementDate = from_datetime(
        aba_USD_ISDA_OIS["data_settlement_swap"].options(dates=dt.date).value
    )

    """Deposits"""
    linha_debug = ut.print_debug("debug", linha_debug, "   adicionando os deposits")
    OISdeposits = []  # lista de todos os vertices da curva de deposits
    depoDCCType = DayCountTypes.ACT_360
    calendario = CalendarTypes.UNITED_STATES
    # vertices da curva de swap
    range_data = aba_USD_ISDA_OIS["deposit_maturity"].options(dates=dt.date).value
    range_taxa = aba_USD_ISDA_OIS["deposit_taxa"].value
    # adiciona todos os dados dos deposits
    for i in range(len(range_data)):
        if range_data[i] is None or range_taxa[i] is None:
            # pula se nao tem data ou taxa
            linha_debug = ut.print_debug(
                "debug",
                linha_debug,
                f"   pulando linha {i + 1} - data: {range_data[i]} taxa: {range_taxa[i]}",
                verbal_debug,
            )
        else:
            # info do vertice
            linha_debug = ut.print_debug(
                "debug",
                linha_debug,
                f"   adicionando linha {i + 1} - "
                + f"data: {from_datetime(range_data[i])} taxa: {range_taxa[i]}",
            )
            # adiciona o vertice de deposit
            deposit = IborDeposit(
                trade_date,
                from_datetime(range_data[i]),
                range_taxa[i],
                depoDCCType,
                cal_type=calendario,
            )
            # retorna a variavel do deposit
            linha_debug = ut.print_debug(
                "debug", linha_debug, repr(deposit), verbal_debug
            )
            # adiciona o deposit a lista de deposits
            OISdeposits.append(deposit)

    """OIS Swaps"""
    linha_debug = ut.print_debug("debug", linha_debug, "   adicionanado os swaps")
    OIS_swaps = []  # lista de todos os vertices da curva de swaps
    swapType = SwapTypes.PAY
    dcType = DayCountTypes.ACT_360  # THIRTY_E_360_ISDA
    # dcType = DayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FrequencyTypes.ANNUAL  # SEMI_ANNUAL
    calendario = CalendarTypes.UNITED_STATES
    # vertices da curva de swap
    range_data = aba_USD_ISDA_OIS["swap_maturity"].options(dates=dt.date).value
    range_taxa = aba_USD_ISDA_OIS["swap_taxa"].value
    # verifica se tem dados
    linha_debug = ut.print_debug(
        "debug",
        linha_debug,
        "   Rno de datas: "
        + str(len(range_data))
        + " Rno de taxas: "
        + str(len(range_taxa)),
        verbal_debug,
    )
    # adiciona todos os dados dos swaps
    for i in range(len(range_taxa)):
        if range_data[i] is None or range_taxa[i] is None:
            # pula se nao tem data ou taxa
            linha_debug = ut.print_debug(
                "debug",
                linha_debug,
                f"   pulando linha {i + 1} - data: {range_data[i]} taxa: {range_taxa[i]}",
                verbal_debug,
            )
        else:
            # info do vertice
            linha_debug = ut.print_debug(
                "debug",
                linha_debug,
                f"   adicionando linha {i + 1} - "
                + f"data: {from_datetime(range_data[i])} taxa: {range_taxa[i]}",
                verbal_debug,
            )
            # adiciona o vertice de swap
            swap = OIS(
                settlementDate,
                from_datetime(range_data[i]),
                swapType,
                range_taxa[i],
                fixedFreq,
                dcType,
                1e6,
                2,
                0.0,
                fixedFreq,
                dcType,
                calendario,
            )
            # retorna a variavel do swap
            linha_debug = ut.print_debug("debug", linha_debug, repr(swap), verbal_debug)
            # adiciona o swap a lista de swaps
            OIS_swaps.append(swap)

    """Monta a curva de juros"""
    linha_debug = ut.print_debug("debug", linha_debug, "   montando a curva de juros")
    # cria a curva de juros OIS so com a SOFR e os swaps
    Interest_types = InterpTypes.FLAT_FWD_RATES
    curva_OIS = OISCurve(trade_date, OISdeposits, [], OIS_swaps, Interest_types)

    # mostra no debug a variavel da curva de juros
    # for swap in curva_OIS:
    #    data = swap.maturity_dt  # data
    #    df = curva_OIS.df(data) / curva_OIS.df(settlementDate)  # discount factor
    #    ccZeroRate = curva_OIS.zero_rate(
    #        data, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_360
    #    )  # zero rate (bullet)
    #    linha_debug = ut.print_debug(
    #        "debug",
    #        linha_debug,
    #        f"{data:%d-%b-%Y}, {df:.8f}, {ccZeroRate*100:.6f}%"
    #    )

    # mostra a curva de juros no debug
    linha_debug = ut.print_debug("debug", linha_debug, repr(curva_OIS), verbal_debug)

    # fim
    linha_debug = (
        ut.print_debug("debug", linha_debug + 1, "### curva de juros montada !!! ###")
        + 1
    )
    # marca a data que a curva de juros foi gerada
    global data_curva_juros
    # atualiza a data da curva de juros
    data_curva_juros = from_datetime(dt.datetime.now(tz=dt.timezone.utc))
    # retorna a curva de juros e a linha de debug
    return curva_OIS, data_curva_juros, linha_debug


def cds_single_name():
    """Cria o contrato de CDS e calcula o premio do CDS"""
    # nome da aba Debug e inicializa a linha de debug
    global linha_debug, verbal_debug
    linha_debug = ut.print_debug("debug", linha_debug, "### Calculando o CDS ###")
    # cria a variavel da planilha
    rc_book = xw.Book.caller()  # xw.Book(planilha)
    # checa se eh para mostrar todas as informacoes ou somente as de debug
    verbal_debug = rc_book.sheets["debug"].range("verbal_debug").value
    # se estiver vazio ou nao for um booleano, entao assume False
    if verbal_debug is None or not isinstance(verbal_debug, bool):
        verbal_debug = False
    # imprime que tipo de debug esta sendo usado
    linha_debug = ut.print_debug(
        "debug", linha_debug, f"   verbal_debug? {verbal_debug}"
    )

    """Carrega a curva de juros da ISDA"""
    # forca o uso das variaveis globais
    global curva_OIS, data_curva_juros
    # variaveis da funcao
    coluna = 0
    aba_cds = None
    # verifica se a curva de juros ja foi carregada
    if curva_OIS is None or data_curva_juros is None:
        # se a curva de juros nao foi carregada, entao carrega
        linha_debug = ut.print_debug(
            "debug",
            linha_debug,
            "   Carregando a curva de juros do mercado",
            verbal_debug,
        )
        # chama a funcao que carrega a curva de juros
        curva_OIS, data_curva_juros, linha_debug = carrega_curva_mercado(linha_debug)
    else:
        # se a curva de juros ja foi carregada, entao verifica se a curva de juros esta atualizada
        if data_curva_juros < dt.datetime.now(tz=dt.timezone.utc) - dt.timedelta(
            minutes=15
        ):
            # se a curva de juros nao esta atualizada, entao carrega novamente
            linha_debug = ut.print_debug(
                "debug",
                linha_debug,
                "   A curva de juros nao esta atualizada, carregando novamente",
                verbal_debug,
            )
            # chama a funcao que carrega a curva de juros
            curva_OIS, data_curva_juros, linha_debug = carrega_curva_mercado(
                linha_debug
            )
        else:
            # se a curva de juros esta atualizada, entao usa a curva de juros
            linha_debug = ut.print_debug(
                "debug",
                linha_debug,
                "   Usando a curva de juros da memoria",
                verbal_debug,
            )
            # usa a curva de juros ja carregada

    """Pega a coluna e nome da aba do CDS sendo calculado"""
    coluna = int(rc_book.sheets["debug"].range("coluna_uma_CLN").value)
    if coluna == 0:
        coluna = 2  # se nao tiver definido, assume a coluna 2
    aba_cds = rc_book.sheets["debug"].range("aba_cotacao").value
    if aba_cds is None:
        sys.exit(
            "A aba de cotacao nao foi definida. Favor definir a aba de cotacao na aba debug."
        )
    # cria a variavel da tab da curva de juros OIS-ISDA
    aba_cashflow_CDS = rc_book.sheets[aba_cds]
    linha_debug = ut.print_debug(
        "debug", linha_debug, f"aba_cashflow_CDS: {aba_cashflow_CDS} e coluna: {coluna}"
    )

    """Calcula a curva de CDS para um unico nome"""
    # Dados do CDS
    trade_date = from_datetime(
        aba_cashflow_CDS.range((11, coluna)).options(dates=dt.date).value
    )
    # t+1 do trade date = step-in date
    effective_date = trade_date.add_days(1)
    # data de liquidação do CDS
    # settlement_date = trade_date.add_days(3)

    settlement_date = Calendar(CalendarTypes.UNITED_STATES).add_business_days(
        trade_date, 3
    )

    # recovery rate
    recovery_rate = aba_cashflow_CDS.range((22, coluna)).value
    # CDS Flat Spread negociado com o mercado - spread do CDS em percentual
    cds_flat_spread = aba_cashflow_CDS.range((17, coluna)).value / 10000
    # cria os vertices de CDS
    # cds6m = CDS(effective_date, "6M", cds_flat_spread)  # effective_date
    # cds1y = CDS(effective_date, "1Y", cds_flat_spread)
    # cds2y = CDS(effective_date, "2Y", cds_flat_spread)
    cds3y = CDS(effective_date, "3Y", cds_flat_spread)
    # cds4y = CDS(effective_date, "4Y", cds_flat_spread)
    cds5y = CDS(effective_date, "5Y", cds_flat_spread)
    cds7y = CDS(effective_date, "7Y", cds_flat_spread)
    cds10y = CDS(effective_date, "10Y", cds_flat_spread)
    # Lista dos vertices de CDS
    # cds_list = [cds6m, cds1y, cds2y, cds3y, cds4y, cds5y, cds7y, cds10y]
    cds_list = [cds3y, cds5y, cds7y, cds10y]
    # Curva do CDS
    issuer_curve = CDSCurve(
        trade_date,
        cds_list,
        curva_OIS,
        recovery_rate,
        # DayCountTypes.ACT_360,
        DayCountTypes.ACT_365F,
        InterpTypes.FLAT_FWD_RATES,
    )
    linha_debug = ut.print_debug("debug", linha_debug, repr(issuer_curve), verbal_debug)

    """Cria o contrato de CDS"""
    maturity_date = from_datetime(
        aba_cashflow_CDS.range((19, coluna)).options(dates=dt.date).value
    )
    # spread do CDS em percentual
    running_coupon = aba_cashflow_CDS.range((21, coluna)).value / 10000
    # notional
    notional = aba_cashflow_CDS.range((6, coluna)).value
    # long protection?
    long_protection = True  # True = compra CDS, False = vende CDS
    # coupons' frequency - take the integer part of the value in the cell and convert to FrequencyTypes
    freq_type = FrequencyTypes(int(aba_cashflow_CDS.range((7, coluna)).value))
    # contrato de CDS
    cds_contract = CDS(
        step_in_dt=effective_date,
        maturity_dt_or_tenor=maturity_date,
        running_cpn=running_coupon,
        notional=notional,
        long_protect=long_protection,
        freq_type=freq_type,
    )

    # mostra os dados do contrato de CDS criado acima
    linha_debug = ut.print_debug("debug", linha_debug, repr(cds_contract), verbal_debug)
    # mostra todos os dados da curva de CDS caso o debug seja verbal
    for dados in cds_contract.print_payments_by_HH(
        effective_date, issuer_curve, full_flow=True
    ):
        linha_debug = ut.print_debug("debug", linha_debug, dados, verbal_debug)

    """Calcula o premio do CDS e outros dados"""
    # cds_value = cds_contract.value_by_HH(trade_date, issuer_curve, recovery_rate)
    cds_value = cds_contract.value(trade_date, issuer_curve, recovery_rate)

    # dirty value
    aba_cashflow_CDS.range((27, coluna)).value = cds_value[0]  # "dirty_pv"]
    # clean value
    aba_cashflow_CDS.range((29, coluna)).value = cds_value[1]  # "clean_pv"]
    # accrual: dias e juros
    aba_cashflow_CDS.range((28, coluna + 1)).value = cds_contract.accrued_days(
        trade_date
    )
    aba_cashflow_CDS.range((28, coluna)).value = cds_contract.accrued_interest(
        trade_date
    )
    # premium leg
    # aba_cashflow_CDS.range((30, coluna)).value = cds_contract.premium_leg_pv_by_HH(
    #    trade_date, issuer_curve, recovery_rate
    # )
    aba_cashflow_CDS.range((30, coluna)).value = cds_contract.premium_leg_pv(
        trade_date, issuer_curve, recovery_rate
    )
    # protection leg + Accrued Interest
    # aba_cashflow_CDS.range((31, coluna)).value = cds_contract.prot_leg_pv_by_HH(
    #    trade_date, issuer_curve, recovery_rate
    # )
    aba_cashflow_CDS.range((31, coluna)).value = cds_contract.prot_leg_pv(
        trade_date, issuer_curve, recovery_rate
    )

    # Extras by HH
    aba_cashflow_CDS.range((32, coluna)).value = cds_contract.cash_settlement_amount(
        trade_date, settlement_date, issuer_curve, recovery_rate
    )
    aba_cashflow_CDS.range((33, coluna)).value = cds_contract.upfront(
        trade_date, settlement_date, issuer_curve, recovery_rate
    )
    aba_cashflow_CDS.range((35, coluna)).value = cds_contract.premium_leg_pv_by_HH(
        trade_date, issuer_curve, recovery_rate
    )

    # end
    linha_debug = ut.print_debug(
        "debug",
        linha_debug + 1,
        "### cds_single_name - finalizado ###",
    )
    linha_debug += 1

    # By HH: extras para teste
    temp1 = cds_contract.premium_leg_pv_by_HH(trade_date, issuer_curve, recovery_rate)
    linha_debug = ut.print_debug(
        "debug", linha_debug, f"prm leg by HH: {temp1:.2f}", verbal_debug
    )
    temp1 = cds_contract.prot_leg_pv_plus_accrual(
        trade_date, issuer_curve, recovery_rate
    )
    linha_debug = ut.print_debug(
        "debug", linha_debug, f"prot+accrual leg by HH: {temp1:.2f}", verbal_debug
    )
    temp1 = cds_contract.rpv01(trade_date, issuer_curve, pv01_method=0)
    linha_debug = ut.print_debug(
        "debug",
        linha_debug,
        f"RPV01 dirty: {temp1[0]:.2f}, clean: {temp1[1]:.2f}",
        verbal_debug,
    )
    linha_debug = ut.print_debug(
        "debug",
        linha_debug,
        f"Frequencia: {cds_contract.freq_type} ( {annual_frequency(cds_contract.freq_type):.2f} )",
        verbal_debug,
    )
    temp1 = cds_contract.value_by_HH(trade_date, issuer_curve, recovery_rate)
    linha_debug = ut.print_debug(
        "debug",
        linha_debug,
        f"Value by HH dirty: {temp1[0]:.2f}, clean: {temp1[1]:.2f}",
        verbal_debug,
    )
    temp1 = cds_contract.value(trade_date, issuer_curve, recovery_rate)
    linha_debug = ut.print_debug(
        "debug",
        linha_debug,
        f"Value dirty: {temp1[0]:.2f}, clean: {temp1[1]:.2f}",
        verbal_debug,
    )
    temp1 = cds_contract.get_pcd(trade_date)
    linha_debug = ut.print_debug(
        "debug",
        linha_debug,
        f"pcd: {temp1[0]}, start_index: {temp1[1]}",
        verbal_debug,
    )
    temp1 = cds_contract.accrual_factors
    linha_debug = ut.print_debug(
        "debug",
        linha_debug,
        f"Accrual factors: {temp1}",
        verbal_debug,
    )


def f_premio_cds(
    cds_flat_spread: float,
    trade_date: Date,
    effective_date: Date,
    maturity_date: Date,
    running_coupon: float,
    recovery_rate: float,
    notional: float,
    long_protection: bool,
    curva_OIS,
) -> float:
    """Calcula o premio do CDS da CLN dado o spread

    Args:
        cds_flat_spread (float): spread o CDS
        trade_date (Date): data da negociacao
        effective_date (Date): data de inicio do CDS
        maturity_date (Date): data de vencimento do CDS
        running_coupon (float): coupon do CDS
        recovery_rate (float): taxa de recuperacao
        notional (float): valor nocional
        long_protection (bool): True = compra CDS, False = vende CDS
        curva_OIS (_type_): curva OIS

    Returns:
        float: premio do CDS
    """
    """calcula o premio do CDS embutido com o preco estimado"""
    # verifica se cds_flat_spread eh um array numpy e muda para float
    if isinstance(cds_flat_spread, np.ndarray):
        cds_flat_spread = cds_flat_spread[0].astype(float)
    # cria os vertices de CDS
    cds6m = CDS(effective_date, "6M", cds_flat_spread)
    cds1y = CDS(effective_date, "1Y", cds_flat_spread)
    cds2y = CDS(effective_date, "2Y", cds_flat_spread)
    cds3y = CDS(effective_date, "3Y", cds_flat_spread)
    cds4y = CDS(effective_date, "4Y", cds_flat_spread)
    cds5y = CDS(effective_date, "5Y", cds_flat_spread)
    cds7y = CDS(effective_date, "7Y", cds_flat_spread)
    cds10y = CDS(effective_date, "10Y", cds_flat_spread)
    # Lista dos vertices de CDS
    cds_list = [cds6m, cds1y, cds2y, cds3y, cds4y, cds5y, cds7y, cds10y]
    # Curva do CDS
    issuer_curve = CDSCurve(trade_date, cds_list, curva_OIS, recovery_rate)
    # contrato de CDS
    cds_contract = CDS(
        effective_date, maturity_date, running_coupon, notional, long_protection
    )

    """Calcula o premio do CDS e outros dados"""
    cds_value = cds_contract.value_by_HH(trade_date, issuer_curve, recovery_rate)

    """retorna o premio o CDS - clean value = cash amount (BBG)"""
    return cds_value[1]  # "clean_pv"]


def f_auxiliar_cds_embutido(
    cds_flat_spread: float,
    trade_date: Date,
    effective_date: Date,
    maturity_date: Date,
    running_coupon: float,
    recovery_rate: float,
    notional: float,
    long_protection: bool,
    curva_OIS,
    premio_CDS_mercado: float,
    receita_desejada: float,
) -> float:
    """Calcula o premio do CDS embutido na CLN, soma a receita deseja e desconta o premio do CDs de mercado
    obs.:
        - o premio do CDS de mercado ja vem negativo da planilha
        - a receita deseja vem positiva da planilha
        - o premio do CDS embutido na CLN eh calculado e retorna um valor positivo
    """
    # Calcula o clean value = cash amount (BBG)
    cds_estimado = f_premio_cds(
        cds_flat_spread,
        trade_date,
        effective_date,
        maturity_date,
        running_coupon,
        recovery_rate,
        notional,
        long_protection,
        curva_OIS,
    )
    """retorna o liquido da venda do CDS no mercado + compra do CDS embutido na CLN
    + Receita desejada"""
    return cds_estimado + premio_CDS_mercado - receita_desejada


def procura_preco_cds():
    """Retorna o Spread do CDS de tal forma que a diferenca com o premio do CDS de mercado
    seja igual aa receita desejada"""
    # nome da aba Debug e inicializa a linha de debug
    global linha_debug, verbal_debug
    linha_debug = ut.print_debug(
        "debug", linha_debug, "### Calculando o CDS embutido ###"
    )
    linha_debug += 1  # pula uma linha

    """calcula o premio do CDS de mercado"""
    cds_single_name()

    # cria a variavel da planilha
    rc_book = xw.Book.caller()  # xw.Book(planilha)
    # checa se eh para mostrar todas as informacoes ou somente as de debug
    verbal_debug = rc_book.sheets["debug"].range("verbal_debug").value
    # se estiver vazio ou nao for um booleano, entao assume False
    if verbal_debug is None or not isinstance(verbal_debug, bool):
        verbal_debug = False
    # imprime que tipo de debug esta sendo usado
    linha_debug = ut.print_debug(
        "debug", linha_debug, f"   verbal_debug? {verbal_debug}"
    )

    """Carrega a curva de juros da ISDA"""
    # forca o uso das variaveis globais
    global curva_OIS, data_curva_juros
    # variaveis da funcao
    coluna = 0
    aba_cds = None
    # verifica se a curva de juros ja foi carregada
    if curva_OIS is None or data_curva_juros is None:
        # se a curva de juros nao foi carregada, entao carrega
        linha_debug = ut.print_debug(
            "debug",
            linha_debug,
            "   Carregando a curva de juros do mercado",
            verbal_debug,
        )
        # chama a funcao que carrega a curva de juros
        curva_OIS, data_curva_juros, linha_debug = carrega_curva_mercado(linha_debug)
    else:
        # se a curva de juros ja foi carregada, entao verifica se a curva de juros esta atualizada
        if data_curva_juros < dt.datetime.now(dt.timezone.utc) - dt.timedelta(
            minutes=15
        ):
            # se a curva de juros nao esta atualizada, entao carrega novamente
            linha_debug = ut.print_debug(
                "debug",
                linha_debug,
                "   A curva de juros nao esta atualizada, carregando novamente",
                verbal_debug,
            )
            # chama a funcao que carrega a curva de juros
            curva_OIS, data_curva_juros, linha_debug = carrega_curva_mercado(
                linha_debug
            )
        else:
            # se a curva de juros esta atualizada, entao usa a curva de juros
            linha_debug = ut.print_debug(
                "debug",
                linha_debug,
                "   Usando a curva de juros da memoria",
                verbal_debug,
            )
            # usa a curva de juros ja carregada

    """Pega a coluna e nome da aba do CDS sendo calculado"""
    coluna = int(rc_book.sheets["debug"].range("coluna_uma_CLN").value)
    if coluna == 0:
        coluna = 2  # se nao tiver definido, assume a coluna 2
    aba_cds = rc_book.sheets["debug"].range("aba_cotacao").value
    if aba_cds is None:
        sys.exit(
            "A aba de cotacao nao foi definida. Favor definir a aba de cotacao na aba debug."
        )
    # cria a variavel da tab da curva de juros OIS-ISDA
    aba_cashflow_CDS = rc_book.sheets[aba_cds]
    linha_debug = ut.print_debug(
        "debug", linha_debug, f"aba_cashflow_CDS: {aba_cashflow_CDS} e coluna: {coluna}"
    )

    """Calcula a curva do CDS embutido"""
    # Dados do CDS
    trade_date = from_datetime(
        aba_cashflow_CDS.range((11, coluna)).options(dates=dt.date).value
    )
    # t+1 do trade date
    effective_date = trade_date.add_days(1)
    # CDS Flat Spread negociado com o mercado - spread do CDS em percentual
    # cds_flat_spread = (aba_cashflow_CDS.range((17, coluna)).value / 10000)
    # data de vencimento
    maturity_date = from_datetime(
        aba_cashflow_CDS.range((19, coluna)).options(dates=dt.date).value
    )
    # coupon do CDS - spread do CDS em percentual
    running_coupon = aba_cashflow_CDS.range((21, coluna)).value / 10000
    # notional do CDS
    notional = aba_cashflow_CDS.range((6, coluna)).value
    # Definie sempre como Long o CDS embutido que o Banco compra do cliente via a CLN
    long_protection = False  # True = compra CDS, False = vende CDS
    # recovery rate
    recovery_rate = aba_cashflow_CDS.range((22, coluna)).value

    # premio do CDS de mercado
    premio_CDS_mercado = aba_cashflow_CDS.range((24, coluna)).value
    # mostra no debug o premio do CDS de mercado
    linha_debug = ut.print_debug(
        "debug", linha_debug, f"premio mercado: {premio_CDS_mercado}", verbal_debug
    )

    # receita desejada em moeda
    receita_desejada = aba_cashflow_CDS.range((49, coluna)).value
    # mostra no debug a receita desejada
    linha_debug = ut.print_debug(
        "debug", linha_debug, f"receita desejada: {receita_desejada}", verbal_debug
    )

    # temp = f_auxiliar_cds_embutido(
    #     0.0123,
    #     trade_date,
    #     effective_date,
    #     maturity_date,
    #     running_coupon,
    #     recovery_rate,
    #     notional,
    #     long_protection,
    #     curva_OIS,
    #     premio_CDS_mercado,
    #     receita_desejada,
    # )
    # linha_debug = ut.print_debug("debug", linha_debug, f"teste: {temp}", verbal_debug)

    """faz o goal-seek para achar o valor do CDS que chegue no spread desejado"""
    try:
        # OLD: options = {"disp": False, "xtol": 0.10, "maxfev": 100}
        options = {"disp": False, "maxiter": 100, "xtol": 0.000001, "maxfev": 100}
        # procura a raiz
        resultado_root = opt.root(
            f_auxiliar_cds_embutido,
            x0=0.0123,
            args=(
                trade_date,
                effective_date,
                maturity_date,
                running_coupon,
                recovery_rate,
                notional,
                long_protection,
                curva_OIS,
                premio_CDS_mercado,
                receita_desejada,
            ),
            method="hybr",
            tol=0.000001,  # tolerance of 1 cent
            options=options,
        )
        # coloca o resultado do opt.root
        linha_debug = ut.print_debug(
            "debug", linha_debug, f"resultado_root: {resultado_root}", verbal_debug
        )
        # verifica se o resultado foi encontrado
        if resultado_root.success:
            # retorna o spread do CDS em bps = % * 10000
            aba_cashflow_CDS.range((55, coluna)).value = resultado_root.x[0] * 10000
            # retorna o premio do CDS em moeda
            aba_cashflow_CDS.range((54, coluna + 2)).value = -f_premio_cds(
                resultado_root.x[0],
                trade_date,
                effective_date,
                maturity_date,
                running_coupon,
                recovery_rate,
                notional,
                long_protection,
                curva_OIS,
            )
        # nao achou nenhuma raiz
        else:
            # coloca a mensagem de erro
            aba_cashflow_CDS.range((55, coluna)).value = (
                f"Error: no solution found. {resultado_root.message}"
            )
            # apaga o premio do CDS
            aba_cashflow_CDS.range((55, coluna + 1)).clear_contents()

    # deu algum erro esperado na busca da raiz / nos calculos numericos
    except (TypeError, ValueError, RuntimeError, ZeroDivisionError, OverflowError) as e:
        linha_debug = ut.print_debug("debug", linha_debug, f"error: {e}", verbal_debug)
        aba_cashflow_CDS.range((55, coluna)).value = (
            "Error: no solution found. See debug."
        )
        # apaga o premio do CDS
        aba_cashflow_CDS.range((55, coluna + 1)).clear_contents()

    # end
    linha_debug = ut.print_debug(
        "debug",
        linha_debug + 1,
        "### procura_preco_cds - finalizado ###",
    )
    linha_debug += 1


#############
# Funcoes de teste ou opcionais
#
##############
def retorna_survival_probability():
    """Cria o contrato de CDS e retorna a probabilidade de sobrevivencia do CDS"""
    # variaveis globais
    global curva_OIS, data_curva_juros
    # nome da aba Debug e inicializa a linha de debug
    linha_debug = ut.print_debug("debug", 1, "### Calculando o CDS ###")
    # cria a variavel da planilha
    rc_book = xw.Book.caller()  # xw.Book(planilha)
    # cria a variavel da tab da curva de juros OIS-ISDA
    aba_cashflow_CDS = rc_book.sheets["Cashflow CDS"]
    # checa se eh para mostrar todas as informacoes ou somente as de debug
    verbal_debug = rc_book.sheets["debug"].range("verbal_debug").value
    # se estiver vazio ou nao for um booleano, entao assume False
    if verbal_debug is None or not isinstance(verbal_debug, bool):
        verbal_debug = False
    # imprime que tipo de debug esta sendo usado
    linha_debug = ut.print_debug(
        "debug", linha_debug, f"   verbal_debug? {verbal_debug}"
    )

    """Carrega a curva de juros da ISDA"""
    curva_OIS, data_curva_juros, linha_debug = carrega_curva_mercado(linha_debug)

    """Calcula a curva de CDS para um unico nome"""
    trade_date = from_datetime(
        aba_cashflow_CDS["trade_date_CDS"].options(dates=dt.date).value
    )
    effective_date = from_datetime(
        aba_cashflow_CDS["CDS_effective_date"].options(dates=dt.date).value
    )
    cds_flat_spread = (
        aba_cashflow_CDS["CDS_spread"].value / 10000
    )  # spread do CDS em percentual
    # cria os vertices de CDS
    cds6m = CDS(effective_date, "6M", cds_flat_spread)
    cds1y = CDS(effective_date, "1Y", cds_flat_spread)
    cds2y = CDS(effective_date, "2Y", cds_flat_spread)
    cds3y = CDS(effective_date, "3Y", cds_flat_spread)
    cds4y = CDS(effective_date, "4Y", cds_flat_spread)
    cds5y = CDS(effective_date, "5Y", cds_flat_spread)
    cds7y = CDS(effective_date, "7Y", cds_flat_spread)
    cds10y = CDS(effective_date, "10Y", cds_flat_spread)
    # lista dos vertices de CDS
    cds_list = [cds6m, cds1y, cds2y, cds3y, cds4y, cds5y, cds7y, cds10y]
    # recovery rate
    recovery_rate = aba_cashflow_CDS["CDS_recovery_rate"].value
    # Curva do CDS
    issuer_curve = CDSCurve(trade_date, cds_list, curva_OIS, recovery_rate)
    linha_debug = ut.print_debug("debug", linha_debug, repr(issuer_curve), verbal_debug)

    """Cria o contrato de CDS"""
    maturity_date = from_datetime(
        aba_cashflow_CDS["CDS_maturity"].options(dates=dt.date).value
    )
    running_coupon = (
        aba_cashflow_CDS["CDS_coupon"].value / 10000
    )  # spread do CDS em percentual
    notional = aba_cashflow_CDS["notional"].value
    long_protection = False  # True = compra CDS, False = vende CDS
    # contrato de CDS
    cds_contract = CDS(
        effective_date, maturity_date, running_coupon, notional, long_protection
    )
    linha_debug = ut.print_debug("debug", linha_debug, repr(cds_contract), verbal_debug)
    for dados in cds_contract.print_payments_by_HH(effective_date, issuer_curve):
        linha_debug = ut.print_debug("debug", linha_debug, dados, verbal_debug)

    # end
    linha_debug = (
        ut.print_debug(
            "debug",
            linha_debug + 1,
            "### retorna_survival_probability - finalizado ###",
        )
        + 1
    )


#############
# Main
#############
if __name__ == "__main__":
    # Used for debug
    # xw.Book('RC generic quote Python v0.1.xlsm').set_mock_caller()
    # Used for frozen executable
    carrega_curva_mercado()
