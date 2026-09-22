# ============================================================================
# FINANCEPY EXAMPLES - CDSBasket
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the valuation of nth-to-default CDS baskets.
#
# A CDS basket contains several reference credits. The protection payment is
# triggered by the nth default in the portfolio:
#
#       1st-to-default
#       2nd-to-default
#       ...
#       nth-to-default
#
# Basket valuation therefore depends not only on the individual credit
# curves, but also on the dependence between defaults.
#
# The example considers:
#
#   1. Construction of the interest-rate discount curve
#   2. Construction of homogeneous issuer CDS curves
#   3. CDS index spread measures
#   4. Gaussian-copula Monte Carlo valuation
#   5. One-factor homogeneous Gaussian valuation
#   6. Student-t copula Monte Carlo valuation
#
# The Gaussian and Student-t models allow the effect of default dependence
# and tail dependence on nth-to-default spreads to be examined.
# ============================================================================

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.products.credit.cds_basket import CDSBasket
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
from financepy.utils.date import Date
from financepy.utils.math import corr_matrix_generator
from financepy.utils.format_graphs import set_plot_style

from helpers import build_ibor_curve
from helpers import load_homogeneous_spread_curves

# ============================================================================
# GLOBAL OUTPUT FORMAT
# ============================================================================

LINE = "=" * 110
SUBLINE = "-" * 110
set_plot_style()

# ============================================================================
# SUPPORTING FUNCTION - BUILD IBOR CURVE
# ============================================================================
#
# Construct the interest-rate curve used to discount CDS cash flows.
#
# The original FinancePy example uses five par swaps with maturities from
# one to five years.
# ============================================================================


# ============================================================================
# SUPPORTING FUNCTION - HETEROGENEOUS CREDIT CURVES
# ============================================================================
#
# The original example also provides a loader for CDX.NA.IG Series 7 market
# spreads. It is retained here so that the homogeneous basket can easily be
# replaced by a basket of issuer-specific credit curves.
#
# The main example below uses homogeneous curves, matching the active branch
# of the original FinancePy example.
# ============================================================================


# ============================================================================
# 1. CDS BASKET SETUP
# ============================================================================

print("\n" + LINE)
print("1. CDS BASKET SETUP")
print(LINE)

trade_dt = Date(
    1,
    3,
    2007,
)

step_in_dt = trade_dt.add_days(1)
value_dt = trade_dt

basket_maturity = Date(
    20,
    12,
    2011,
)

num_credits = 5

seed = 1967


print(f"{'Trade Date':<40}: " f"{trade_dt}")

print(f"{'Value Date':<40}: " f"{value_dt}")

print(f"{'Step-In Date':<40}: " f"{step_in_dt}")

print(f"{'Basket Maturity':<40}: " f"{basket_maturity}")

print(f"{'Number of Credits':<40}: " f"{num_credits}")

print(f"{'Monte Carlo Seed':<40}: " f"{seed}")


# ============================================================================
# 2. BUILD INTEREST-RATE CURVE
# ============================================================================

print("\n" + LINE)
print("2. BUILD INTEREST-RATE CURVE")
print(LINE)

libor_curve = build_ibor_curve(trade_dt)

print("IBOR discount curve constructed from " "1Y, 2Y, 3Y, 4Y and 5Y swaps.")


# ============================================================================
# 3. BUILD HOMOGENEOUS CREDIT CURVES
# ============================================================================
#
# All five credits are assigned the same CDS term structure.
#
# The spreads below are decimal rates:
#
#       3Y  = 12 bp
#       5Y  = 25 bp
#       7Y  = 34 bp
#       10Y = 46 bp
#
# Recovery is 40%.
# ============================================================================

print("\n" + LINE)
print("3. BUILD HOMOGENEOUS CREDIT CURVES")
print(LINE)

spd_3yr = 0.0012
spd_5yr = 0.0025
spd_7yr = 0.0034
spd_10yr = 0.0046

issuer_curves = load_homogeneous_spread_curves(
    value_dt,
    libor_curve,
    spd_3yr,
    spd_5yr,
    spd_7yr,
    spd_10yr,
    num_credits,
)


print(f"{'MATURITY':<20}" f"{'CDS SPREAD (bp)':>20}")

print(SUBLINE)

print(f"{'3Y':<20}" f"{spd_3yr * 10000.0:20.6f}")

print(f"{'5Y':<20}" f"{spd_5yr * 10000.0:20.6f}")

print(f"{'7Y':<20}" f"{spd_7yr * 10000.0:20.6f}")

print(f"{'10Y':<20}" f"{spd_10yr * 10000.0:20.6f}")


# ============================================================================
# 4. CDS INDEX PORTFOLIO SPREAD MEASURES
# ============================================================================
#
# Before valuing nth-to-default protection, calculate several spread
# measures for the underlying portfolio:
#
#   Intrinsic spread
#       Portfolio spread obtained from the constituent credit curves.
#
#   Total spread
#       Sum of the constituent spread contributions.
#
#   Minimum spread
#       Lowest constituent-equivalent spread.
#
#   Maximum spread
#       Highest constituent-equivalent spread.
#
# Since this example uses homogeneous issuer curves, the constituent credit
# characteristics are identical.
# ============================================================================

print("\n" + LINE)
print("4. CDS INDEX PORTFOLIO SPREAD MEASURES")
print(LINE)

cds_index = CDSIndexPortfolio()


intrinsic_spread = cds_index.intrinsic_spread(
    value_dt,
    step_in_dt,
    basket_maturity,
    issuer_curves,
)

total_spread = cds_index.total_spread(
    value_dt,
    step_in_dt,
    basket_maturity,
    issuer_curves,
)

minimum_spread = cds_index.min_spread(
    value_dt,
    step_in_dt,
    basket_maturity,
    issuer_curves,
)

maximum_spread = cds_index.max_spread(
    value_dt,
    step_in_dt,
    basket_maturity,
    issuer_curves,
)


print(f"{'SPREAD MEASURE':<40}" f"{'VALUE (bp)':>20}")

print(SUBLINE)

print(f"{'Intrinsic Spread':<40}" f"{intrinsic_spread * 10000.0:20.6f}")

print(f"{'Total Spread':<40}" f"{total_spread * 10000.0:20.6f}")

print(f"{'Minimum Spread':<40}" f"{minimum_spread * 10000.0:20.6f}")

print(f"{'Maximum Spread':<40}" f"{maximum_spread * 10000.0:20.6f}")


# ============================================================================
# 5. GAUSSIAN COPULA
# ============================================================================
#
# Value each nth-to-default contract using:
#
#   1. Gaussian-copula Monte Carlo
#   2. One-factor homogeneous Gaussian model
#
# Two factor loadings are considered:
#
#       beta = 0.0
#       beta = 0.5
#
# In the one-factor construction:
#
#       rho = beta^2
#
# giving asset/default dependence parameters of:
#
#       rho = 0.00
#       rho = 0.25
#
# Comparing the two valuation methods also provides a useful numerical check.
# ============================================================================

print("\n" + LINE)
print("5. GAUSSIAN COPULA")
print(LINE)

basket = CDSBasket(
    value_dt,
    basket_maturity,
)

num_trials = 1000

gaussian_results = []


print(
    f"{'NTD':>5}"
    f"{'RHO':>10}"
    f"{'TRIALS':>12}"
    f"{'MC SPREAD (bp)':>20}"
    f"{'1F SPREAD (bp)':>20}"
    f"{'DIFF (bp)':>16}"
    f"{'TIME (s)':>14}"
)

print(SUBLINE)


for ntd in range(
    1,
    num_credits + 1,
):

    for beta in [
        0.0,
        0.5,
    ]:

        rho = beta * beta

        beta_vector = np.ones(num_credits) * beta

        corr_matrix = corr_matrix_generator(
            rho,
            num_credits,
        )

        start = time.perf_counter()

        gaussian_mc = basket.value_gaussian_mc(
            value_dt,
            ntd,
            issuer_curves,
            corr_matrix,
            libor_curve,
            num_trials,
            seed,
        )

        gaussian_homo = basket.value_1f_gaussian_homo(
            value_dt,
            ntd,
            issuer_curves,
            beta_vector,
            libor_curve,
        )

        elapsed = time.perf_counter() - start

        # The original FinancePy example uses element 2 of the Monte Carlo
        # result and element 3 of the homogeneous one-factor result as the
        # corresponding basket spreads.

        mc_spread = gaussian_mc[2] * 10000.0

        homo_spread = gaussian_homo[3] * 10000.0

        spread_difference = mc_spread - homo_spread

        gaussian_results.append(
            {
                "ntd": ntd,
                "beta": beta,
                "rho": rho,
                "mc_spread": mc_spread,
                "homo_spread": homo_spread,
                "difference": spread_difference,
                "time": elapsed,
            }
        )

        print(
            f"{ntd:5d}"
            f"{rho:10.4f}"
            f"{num_trials:12d}"
            f"{mc_spread:20.6f}"
            f"{homo_spread:20.6f}"
            f"{spread_difference:16.6f}"
            f"{elapsed:14.6f}"
        )


# ============================================================================
# 6. STUDENT-T COPULA
# ============================================================================
#
# The Student-t copula introduces heavier joint tails than the Gaussian
# copula.
#
# The original example compares:
#
#       degrees of freedom = 3
#       degrees of freedom = 4
#
# with the Gaussian-copula result.
#
# Lower degrees of freedom imply heavier tails.
# ============================================================================

print("\n" + LINE)
print("6. STUDENT-T COPULA VERSUS GAUSSIAN COPULA")
print(LINE)

student_results = []


print(f"{'NTD':>5}" f"{'RHO':>10}" f"{'MODEL':>12}" f"{'DOF':>8}" f"{'SPREAD (bp)':>20}" f"{'TIME (s)':>14}")

print(SUBLINE)


for beta in [
    0.0,
    0.5,
]:

    rho = beta**2

    corr_matrix = corr_matrix_generator(
        rho,
        num_credits,
    )

    for ntd in range(
        1,
        num_credits + 1,
    ):

        for degrees_of_freedom in [
            3,
            4,
        ]:

            start = time.perf_counter()

            student_value = basket.value_student_t_mc(
                value_dt,
                ntd,
                issuer_curves,
                corr_matrix,
                degrees_of_freedom,
                libor_curve,
                num_trials,
                seed,
            )

            elapsed = time.perf_counter() - start

            student_spread = student_value[2] * 10000.0

            student_results.append(
                {
                    "ntd": ntd,
                    "rho": rho,
                    "model": "Student-t",
                    "dof": degrees_of_freedom,
                    "spread": student_spread,
                    "time": elapsed,
                }
            )

            print(
                f"{ntd:5d}"
                f"{rho:10.4f}"
                f"{'Student-t':>12}"
                f"{degrees_of_freedom:8d}"
                f"{student_spread:20.6f}"
                f"{elapsed:14.6f}"
            )

        start = time.perf_counter()

        gaussian_value = basket.value_gaussian_mc(
            value_dt,
            ntd,
            issuer_curves,
            corr_matrix,
            libor_curve,
            num_trials,
            seed,
        )

        elapsed = time.perf_counter() - start

        gaussian_spread = gaussian_value[2] * 10000.0

        student_results.append(
            {
                "ntd": ntd,
                "rho": rho,
                "model": "Gaussian",
                "dof": None,
                "spread": gaussian_spread,
                "time": elapsed,
            }
        )

        print(f"{ntd:5d}" f"{rho:10.4f}" f"{'Gaussian':>12}" f"{'-':>8}" f"{gaussian_spread:20.6f}" f"{elapsed:14.6f}")


# ============================================================================
# 7. STUDENT-T COPULA WITH FIVE DEGREES OF FREEDOM
# ============================================================================
#
# Repeat the Student-t calculation with five degrees of freedom, matching
# the final calculation in the original FinancePy example.
# ============================================================================

print("\n" + LINE)
print("7. STUDENT-T COPULA WITH FIVE DEGREES OF FREEDOM")
print(LINE)

degrees_of_freedom = 5

student_5_results = []


print(f"{'NTD':>5}" f"{'RHO':>10}" f"{'TRIALS':>12}" f"{'DOF':>8}" f"{'SPREAD (bp)':>20}" f"{'TIME (s)':>14}")

print(SUBLINE)


for beta in [
    0.0,
    0.5,
]:

    rho = beta**2

    corr_matrix = corr_matrix_generator(
        rho,
        num_credits,
    )

    for ntd in range(
        1,
        num_credits + 1,
    ):

        start = time.perf_counter()

        student_value = basket.value_student_t_mc(
            value_dt,
            ntd,
            issuer_curves,
            corr_matrix,
            degrees_of_freedom,
            libor_curve,
            num_trials,
            seed,
        )

        elapsed = time.perf_counter() - start

        student_spread = student_value[2] * 10000.0

        student_5_results.append(
            {
                "ntd": ntd,
                "rho": rho,
                "spread": student_spread,
                "time": elapsed,
            }
        )

        print(
            f"{ntd:5d}"
            f"{rho:10.4f}"
            f"{num_trials:12d}"
            f"{degrees_of_freedom:8d}"
            f"{student_spread:20.6f}"
            f"{elapsed:14.6f}"
        )


# ============================================================================
# 8. GAUSSIAN NTH-TO-DEFAULT SPREADS
# ============================================================================
#
# Plot the Gaussian Monte Carlo basket spread against default order.
#
# Correlation changes the distribution of joint defaults and therefore
# affects different nth-to-default positions differently.
# ============================================================================

print("\n" + LINE)
print("8. GAUSSIAN NTH-TO-DEFAULT SPREADS")
print(LINE)

for rho in [
    0.0,
    0.25,
]:

    x_values = [result["ntd"] for result in gaussian_results if result["rho"] == rho]

    y_values = [result["mc_spread"] for result in gaussian_results if result["rho"] == rho]

    plt.plot(
        x_values,
        y_values,
        marker="o",
        label=f"rho = {rho:.2f}",
    )

plt.xlabel("Nth Default")

plt.ylabel("Basket Spread (bp)")

plt.title("Gaussian Copula Nth-to-Default Spreads")

plt.xticks(
    range(
        1,
        num_credits + 1,
    )
)

plt.legend()
plt.grid(True)


# ============================================================================
# 9. GAUSSIAN MONTE CARLO VERSUS ONE-FACTOR MODEL
# ============================================================================
#
# Compare the Monte Carlo calculation with the homogeneous one-factor
# Gaussian result.
# ============================================================================

print("\n" + LINE)
print("9. GAUSSIAN MONTE CARLO VERSUS ONE-FACTOR MODEL")
print(LINE)

for rho in [
    0.0,
    0.25,
]:

    x_values = [result["ntd"] for result in gaussian_results if result["rho"] == rho]

    mc_values = [result["mc_spread"] for result in gaussian_results if result["rho"] == rho]

    homo_values = [result["homo_spread"] for result in gaussian_results if result["rho"] == rho]

    plt.figure()

    plt.plot(
        x_values,
        mc_values,
        marker="o",
        label="Gaussian Monte Carlo",
    )

    plt.plot(
        x_values,
        homo_values,
        marker="o",
        label="One-Factor Homogeneous",
    )

    plt.xlabel("Nth Default")

    plt.ylabel("Basket Spread (bp)")

    plt.title(f"Gaussian Basket Valuation: rho = {rho:.2f}")

    plt.xticks(
        range(
            1,
            num_credits + 1,
        )
    )

    plt.legend()
    plt.grid(True)


# ============================================================================
# 10. GAUSSIAN VERSUS STUDENT-T
# ============================================================================
#
# Compare Gaussian and Student-t basket spreads.
#
# The Student-t copula allows stronger joint tail behaviour. The impact can
# be seen by comparing spreads for the same default order and correlation.
# ============================================================================

print("\n" + LINE)
print("10. GAUSSIAN VERSUS STUDENT-T")
print(LINE)

for rho in [
    0.0,
    0.25,
]:

    plt.figure()

    for model_name, dof in [
        ("Gaussian", None),
        ("Student-t", 3),
        ("Student-t", 4),
    ]:

        selected = [
            result
            for result in student_results
            if result["rho"] == rho and result["model"] == model_name and result["dof"] == dof
        ]

        x_values = [result["ntd"] for result in selected]

        y_values = [result["spread"] for result in selected]

        if model_name == "Gaussian":

            label = "Gaussian"

        else:

            label = f"Student-t, dof = {dof}"

        plt.plot(
            x_values,
            y_values,
            marker="o",
            label=label,
        )

    plt.xlabel("Nth Default")

    plt.ylabel("Basket Spread (bp)")

    plt.title(f"Copula Comparison: rho = {rho:.2f}")

    plt.xticks(
        range(
            1,
            num_credits + 1,
        )
    )

    plt.legend()
    plt.grid(True)


# ============================================================================
# 11. SUMMARY
# ============================================================================

print("\n" + LINE)
print("11. SUMMARY")
print(LINE)

print(f"{'Number of Credits':<40}: " f"{num_credits}")

print(f"{'Basket Maturity':<40}: " f"{basket_maturity}")

print(f"{'Monte Carlo Trials':<40}: " f"{num_trials}")

print(f"{'Monte Carlo Seed':<40}: " f"{seed}")

print(f"{'Correlations Examined':<40}: " f"0.00, 0.25")

print(f"{'Student-t Degrees of Freedom':<40}: " f"3, 4, 5")

print(f"{'Nth-to-Default Orders':<40}: " f"1 through {num_credits}")

print("\n" + LINE)
print("END OF CDS BASKET DEMONSTRATION")
print(LINE)


# ============================================================================
# DISPLAY ALL PLOTS
# ============================================================================

plt.show()
