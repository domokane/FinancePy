###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from ...utils.global_vars import g_days_in_year
from ...products.bonds.bond import Bond
from ...utils.date import Date
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.helpers import label_to_string, check_argument_types


# TODO: Examine other exchange conventions.
# TODO: Delivery option model
###############################################################################


class BondFuture:
    """Class for managing futures contracts on government bonds that follows
    CME conventions and related analytics.

    Attributes
    ----------
    ticker_name : str
        Identifier name of the futures contract.
    first_delivery_dt : Date
        The first delivery date of the futures contract (usually IMM date).
    last_delivery_dt : Date
        The last delivery date of the futures contract.
    contract_size : int
        Contract size in currency units (e.g. $100,000).
    cpn : float
        Contract standard coupon rate, used for conversion factor calculations.
    """

    def __init__(
        self,
        ticker_name: str,
        first_delivery_dt: Date,
        last_delivery_dt: Date,
        contract_size: int,
        cpn: float,
    ):
        """
        Initialize BondFuture instance.

        Parameters
        ----------
        ticker_name : str
            The futures contract ticker symbol.
        first_delivery_dt : Date
            First delivery date (IMM date).
        last_delivery_dt : Date
            Last delivery date.
        contract_size : int
            Contract size in currency units.
        cpn : float
            Standard coupon rate of the contract (percent).
        """

        check_argument_types(self.__init__, locals())

        self.ticker_name = ticker_name
        self.first_delivery_dt = first_delivery_dt  # This is the IMM date
        self.last_delivery_dt = last_delivery_dt
        self.contract_size = contract_size
        self.cpn = cpn

    ###########################################################################

    def conversion_factor(self, bond: Bond):
        """
        Compute the CME-style conversion factor for a given bond.

        The conversion factor normalizes bonds with different coupons and
        maturities to a standard contract coupon and maturity schedule.

        Steps:
        - Calculate years to maturity from 1st delivery date to bond maturity.
        - Round this maturity to the nearest 3-month (quarter year) multiple.
        - Construct a hypothetical bond starting on the first delivery date
        with the standard contract coupon and rounded maturity.
        - Calculate the clean price of this hypothetical bond yielding the
        contract coupon.
        - Normalize and round to 4 decimal places.

        Parameters
        ----------
        bond : Bond
            The deliverable bond for which conversion factor is calculated.

        Returns
        -------
        float
            Conversion factor rounded to 4 decimals.
        """

        dy = self.delivery_years(bond)

        roundedt_matInMonths = int(dy * 4) * 3  # DO NOT CHANGE

        new_mat = self.first_delivery_dt.add_months(roundedt_matInMonths)
        ex_div_days = 0

        # Hypothetical issue date for bond is fixed in year 2000 to align
        # with CME convention. This "hack" helps standardize conversion
        # factor calculation
        issue_dt = Date(new_mat.d, new_mat.m, 2000)

        # Construct hypothetical bond with standard contract coupon
        new_bond = Bond(
            issue_dt,
            new_mat,
            bond.cpn,
            bond.freq_type,
            bond.dc_type,
            ex_div_days,
        )

        # Calculate clean price of hypothetical bond given coupon yield
        p = new_bond.clean_price_from_ytm(self.first_delivery_dt, self.cpn)
        p = p / 100.0

        # Convention is to round the conversion factor to 4dp
        p = round(p, 4)
        return p

    ###########################################################################

    def principal_invoice(self, bond: Bond, futures_price: float) -> float:
        """
        Calculate the principal invoice amount for delivering a bond.

        CME defines this as contract size × futures price × conversion factor / 100.

        Parameters
        ----------
        bond : Bond
            The deliverable bond.
        futures_price : float
            The quoted futures price (in percentage points).

        Returns
        -------
        float
            Principal invoice amount in currency units.
        """
        cf = self.conversion_factor(bond)
        pip = self.contract_size * (futures_price * cf) / 100.0
        return pip

    ###########################################################################

    def total_invoice_amount(
        self, settle_dt: Date, bond: Bond, futures_price: float
    ) -> float:
        """
        Calculate total invoice amount paid on delivery, including
        accrued interest.

        Parameters
        ----------
        settle_dt : Date
            Settlement date when the invoice is paid.
        bond : Bond
            The deliverable bond.
        futures_price : float
            Quoted futures price.

        Returns
        -------
        float
            Total invoice amount including accrued interest.
        """

        cf = self.conversion_factor(bond)
        pip = self.contract_size * (futures_price * cf) / 100.0
        accrued = bond.accrued_interest(settle_dt) * self.contract_size / 100.0
        tia = pip + accrued
        return tia

    ###########################################################################

    def delivery_years(self, bond: Bond) -> float:
        """
        Calculate the fractional years from the first delivery date to the
        bond maturity.

        Uses ACT/ACT ISDA day count convention.

        Parameters
        ----------
        bond : Bond
            The deliverable bond.

        Returns
        -------
        float
            Fractional years (e.g., 5.25 years).
        """
        dc = DayCount(DayCountTypes.ACT_ACT_ISDA)

        year_frac, _, _ = dc.year_frac(
            self.first_delivery_dt, bond.maturity_dt
        )

        years = int(year_frac)
        months = int(12 * (year_frac - years))
        del_years = years + months / 12
        return del_years

    ###########################################################################

    def gross_basis(
        self, bond: Bond, clean_price: float, futures_price: float
    ) -> float:
        """
        Compute the gross basis: what is received on delivery vs. bond price.

        gross_basis = clean_price - (conversion_factor × futures_price)

        Parameters
        ----------
        bond : Bond
            The deliverable bond.
        clean_price : float
            The current clean price of the bond.
        futures_price : float
            The quoted futures price.

        Returns
        -------
        float
            Gross basis value.
        """
        cf = self.conversion_factor(bond)
        gross_basis = clean_price - cf * futures_price
        return gross_basis

    ###########################################################################

    def net_basis(
        self,
        bond,
        settle_dt: Date,
        clean_price: float,
        futures_price: float,
        repo_rate: float,
    ):
        """
        Compute net basis, adjusting for financing cost (repo rate).

        net_basis = forward_price - (conversion_factor × futures_price)

        Parameters
        ----------
        bond : Bond
            Deliverable bond.
        settle_dt : Date
            Settlement date.
        clean_price : float
            Current bond clean price.
        futures_price : float
            Futures contract price.
        repo_rate : float
            Repo financing rate (annualized decimal, e.g. 0.02 for 2%).

        Returns
        -------
        float
            Net basis value.

        Raises
        ------
        ValueError
            If any price or repo_rate is negative.
        """

        if clean_price < 0 or futures_price < 0 or repo_rate < 0:
            raise ValueError("Prices and repo rate must be non-negative")

        fwd_date = self.last_delivery_dt
        fwd_price = bond.forward_price(
            settle_dt, fwd_date, clean_price, repo_rate
        )

        cf = self.conversion_factor(bond)
        net_basis = fwd_price - cf * futures_price
        return net_basis

    ###########################################################################

    def implied_repo_rate(
        self,
        bond: Bond,
        settle_dt: Date,
        clean_price: float,
        futures_price: float,
    ):
        """
        Calculate implied repo rate consistent with futures price.

        This is the IRR of financing the bond purchase and delivering at
        futures price, accounting for accrued interest and coupon payments
        between settle and delivery.

        Parameters
        ----------
        bond : Bond
            Deliverable bond.
        settle_dt : Date
            Settlement date.
        clean_price : float
            Current clean price of bond.
        futures_price : float
            Futures price.

        Returns
        -------
        float
            Implied repo rate (annualized decimal).

        Raises
        ------
        ValueError
            If clean_price or futures_price is negative.
        """
        if clean_price < 0 or futures_price < 0:
            raise ValueError("Prices and futures price must be non-negative")

        delivery_dt = self.last_delivery_dt
        YEAR = 360  # AGREES WITH BLOOMBERG

        # Accrued interest on delivery date
        AC_settle = bond.accrued_interest(settle_dt)
        AC_delivery = bond.accrued_interest(delivery_dt)
        full_price = clean_price + AC_settle

        fv_cpns = 0.0
        weighted_days_sum = 0.0
        for dt, amt in zip(bond.cpn_dts[1:], bond.flow_amounts[1:]):
            if settle_dt < dt <= delivery_dt:
                fv_cpns += amt * bond.par
                weighted_days_sum += (dt - settle_dt) * amt * bond.par

        if fv_cpns > 0:
            avg_coupon_days = weighted_days_sum / fv_cpns
        else:
            avg_coupon_days = 0.0

        cf = self.conversion_factor(bond)
        D1 = float(delivery_dt - settle_dt)

        num = futures_price * cf + AC_delivery - full_price + fv_cpns
        denom = full_price * (D1 / YEAR) - fv_cpns * (avg_coupon_days / YEAR)
        irr = num / denom
        return irr

    ###########################################################################

    def ctd(self, bonds: list, bond_clean_prices: list, futures_price: float):
        """
        Determine the Cheapest to Deliver (CTD) bond among candidates.

        The CTD is the bond with the highest gross basis
        (= clean price - cf × futures price).

        Parameters
        ----------
        bonds : list of Bond
            List of deliverable bonds.
        bond_clean_prices : list of float
            Corresponding clean prices of the bonds.
        futures_price : float
            Futures contract price.

        Returns
        -------
        Bond
            The bond identified as Cheapest to Deliver.
        """

        ctd_bond = None
        ctd_net = float("-inf")

        for bond_clean_price, bond in zip(bond_clean_prices, bonds):

            # CHANGE THIS TO NET BASIS ?
            net = self.gross_basis(bond, bond_clean_price, futures_price)

            if net > ctd_net:
                ctd_bond = bond
                ctd_net = net

        return ctd_bond

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("TICKER NAME", self.ticker_name)
        s += label_to_string("FIRST DELIVERY DATE", self.first_delivery_dt)
        s += label_to_string("LAST DELIVERY DATE", self.last_delivery_dt)
        s += label_to_string("CONTRACT SIZE", self.contract_size)
        s += label_to_string("COUPON", self.cpn)
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
