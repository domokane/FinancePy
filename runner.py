#!/usr/bin/env python3
"""runner.py — JSON-in/JSON-out FinancePy bond pricing wrapper.

Role in the FaIR → DicePy → FinancePy pipeline
────────────────────────────────────────────────
DicePy emits a climate risk premium (a decimal yield spread) for one or more
emissions scenarios.  This script prices a fixed-coupon bond using the
FinancePy library, applying that climate risk premium as extra yield above a
baseline YTM.

When run standalone, supply a full bond spec (including ytm or clean_price)
in the input JSON.  When chained after DicePy, only the DicePy output fields
are required — bond parameters fall back to built-in defaults and the YTM is
resolved from adjusted_ytm (preferred) or climate_risk_premium.

TO RUN:
  From inside the repo root:
    python runner.py damage_estimates.json
  Or piped:
    python runner.py < damage_estimates.json

INPUT JSON keys:
  issue_date        list[int]  [year, month, day] of bond issue      (default: [2024,1,1])
  maturity_date     list[int]  [year, month, day] of maturity        (default: [2034,1,1])
  settlement_date   list[int]  [year, month, day] of settlement      (default: [2024,1,3])
  coupon_rate       float      annual coupon as a decimal             (default: 0.04)
  frequency         str        coupon frequency — "annual",           (default: "semi_annual")
                               "semi_annual", "quarterly", "monthly"
  day_count         str        day-count convention — "act_act_isda", (default: "act_act_isda")
                               "thirty_e_360", "act_360", "act_365f"
  ytm               float      yield-to-maturity as a decimal;        (optional)
                               if absent, resolved from DicePy fields
                               below.  Provide either ytm or
                               clean_price (not both).
  clean_price       float      clean price per 100 face value;        (optional)
                               used to back-solve ytm when ytm is
                               absent.  Ignored when ytm is present.
  adjusted_ytm      dict       {scenario: {config: float}} from       (optional)
                               DicePy; collapsed to a mean scalar and
                               used as ytm when ytm is absent.
  climate_risk_premium dict    {scenario: {config: float}} from       (optional)
                               DicePy; used as ytm when both ytm and
                               adjusted_ytm are absent.
  insured_aai_usd   float      Average Annual Insured Loss from        (optional)
                               OasisLMF; combined with
                               total_insured_value_usd (or
                               total_exposed_value_usd) to derive an
                               expected-loss yield spread added to
                               risk_free_rate.
  lambda_annual     float      Modelled annual event rate from         (optional)
                               hurricane-frequency. When supplied
                               alongside event_insured_loss_table,
                               replaces CLIMADA's implicit frequency:
                               AAI is rescaled by lambda_annual /
                               sum(per-event freq).
  event_insured_loss_table list Per-event insured losses from OasisLMF (optional)
                               with per-event frequency. Used for the
                               lambda_annual rescaling above.
  risk_free_rate    float      baseline risk-free rate used when       (optional)
                               deriving YTM from OasisLMF EL.
                               (default: 0.045)

OUTPUT JSON keys:
  clean_price         float   clean price per 100 face value
  dirty_price         float   dirty price (clean + accrued interest)
  accrued_interest    float   accrued coupon since last payment date
  macauley_duration   float   Macaulay duration in years
  modified_duration   float   modified duration (price sensitivity to yield)
  convexity           float   convexity (second-order price/yield sensitivity)
  ytm                 float   yield-to-maturity used (or back-solved) for pricing
"""

import argparse
import json
import sys

try:
    from financepy.utils.date import Date
    from financepy.utils.frequency import FrequencyTypes
    from financepy.utils.day_count import DayCountTypes
    from financepy.products.bonds.bond import Bond
except ModuleNotFoundError as exc:
    missing = getattr(exc, "name", None) or str(exc)
    pyver = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"

    msg_lines = [
        "Failed to import FinancePy dependencies.",
        f"Missing module: {missing}",
        f"Python: {pyver}",
        "",
        "Fix:",
        "- Use Python 3.10–3.12 (FinancePy does not currently support 3.13+).",
        "- Create/activate a venv, then install FinancePy deps (editable install is fine):",
        "    python3.12 -m venv .venv",
        "    source .venv/bin/activate",
        "    python -m pip install -U pip",
        "    python -m pip install -e .",
        "",
        "Then run:",
        "    python runner.py < input.json",
    ]
    print("\n".join(msg_lines), file=sys.stderr)
    raise


# ---------------------------------------------------------------------------
# Default bond parameters used when FinancePy runs as a downstream step
# (e.g., after DicePy in a FaIR → DicePy → FinancePy chain). In that case
# only climate risk fields arrive in the input; bond specifics are absent.
# When running standalone, supply these keys in the input JSON to override.
# ---------------------------------------------------------------------------
DEFAULT_ISSUE_DATE = [2024, 1, 1]
DEFAULT_MATURITY_DATE = [2034, 1, 1]    # 10-year bond
DEFAULT_SETTLEMENT_DATE = [2024, 1, 3]  # T+2 settlement
DEFAULT_COUPON_RATE = 0.04              # 4 % annual coupon
DEFAULT_RISK_FREE_RATE = 0.045          # baseline risk-free rate when deriving YTM from OasisLMF EL


def _mean_nested(d: dict) -> float:
    """Collapse a {scenario: {config: float}} dict to a single mean value.

    DicePy outputs climate_risk_premium and adjusted_ytm in this shape.
    We average across all scenario/config combinations to get one scalar
    yield suitable for a single bond pricing call.
    """
    values = [v for inner in d.values() for v in inner.values()]
    return sum(values) / len(values) if values else 0.0


FREQ_MAP = {
    "annual": FrequencyTypes.ANNUAL,
    "semi_annual": FrequencyTypes.SEMI_ANNUAL,
    "quarterly": FrequencyTypes.QUARTERLY,
    "monthly": FrequencyTypes.MONTHLY,
}

DC_MAP = {
    "act_act_isda": DayCountTypes.ACT_ACT_ISDA,
    "thirty_e_360": DayCountTypes.THIRTY_E_360,
    "act_360": DayCountTypes.ACT_360,
    "act_365f": DayCountTypes.ACT_365F,
}


def run(params: dict) -> dict:
    # Bond parameters — use input values if present, fall back to defaults.
    # Defaults exist so this runner can be chained after DicePy without
    # requiring a full bond spec in the input.
    issue = params.get("issue_date", DEFAULT_ISSUE_DATE)      # [year, month, day]
    mat = params.get("maturity_date", DEFAULT_MATURITY_DATE)
    settle = params.get("settlement_date", DEFAULT_SETTLEMENT_DATE)

    issue_dt = Date(issue[2], issue[1], issue[0])
    maturity_dt = Date(mat[2], mat[1], mat[0])
    settle_dt = Date(settle[2], settle[1], settle[0])

    coupon = params.get("coupon_rate", DEFAULT_COUPON_RATE)
    freq = FREQ_MAP[params.get("frequency", "semi_annual")]
    dc = DC_MAP[params.get("day_count", "act_act_isda")]

    bond = Bond(issue_dt, maturity_dt, coupon, freq, dc)

    # ytm resolution: prefer an explicit scalar in the input (standalone use).
    # When chaining from DicePy, fall back to adjusted_ytm (baseline YTM +
    # climate risk premium) if available, then climate_risk_premium alone.
    # Both are {scenario: {config: float}} dicts — _mean_nested collapses
    # them to a single representative scalar for a single bond pricing call.
    ytm = params.get("ytm")
    if ytm is None:
        for key in ("adjusted_ytm", "climate_risk_premium"):
            val = params.get(key)
            if isinstance(val, dict):
                ytm = _mean_nested(val)
                break
    if ytm is None and "insured_aai_usd" in params:
        # OasisLMF input: derive YTM from expected loss rate + risk-free rate.
        # el_rate = Average Annual Insured Loss / Total Insured Value
        tiv = params.get("total_insured_value_usd") or params.get("total_exposed_value_usd")
        if tiv and float(tiv) > 0:
            aai = float(params["insured_aai_usd"])
            # If hurricane-frequency supplied a modelled lambda alongside
            # OasisLMF's event table, swap CLIMADA's implicit frequency for it:
            # AAI = E[loss|event] * E[events/yr]; we keep severity (encoded in
            # the event table) and replace the frequency. lambda_implicit is
            # CLIMADA's total annual event rate (sum of per-event freq).
            elt = params.get("event_insured_loss_table")
            lam = params.get("lambda_annual")
            if isinstance(elt, list) and elt and isinstance(lam, (int, float)):
                lambda_implicit = sum(float(ev.get("frequency", 0.0)) for ev in elt)
                if lambda_implicit > 0:
                    aai = aai * (float(lam) / lambda_implicit)
            el_rate = aai / float(tiv)
            risk_free = float(params.get("risk_free_rate", DEFAULT_RISK_FREE_RATE))
            ytm = risk_free + el_rate

    clean_price = params.get("clean_price")

    result = {}

    if ytm is not None:
        # Price from yield
        result["clean_price"] = bond.clean_price_from_ytm(settle_dt, ytm)
        result["dirty_price"] = bond.dirty_price_from_ytm(settle_dt, ytm)
        result["accrued_interest"] = bond.accrued_interest(settle_dt)
        result["macauley_duration"] = bond.macauley_duration(settle_dt, ytm)
        result["modified_duration"] = bond.modified_duration(settle_dt, ytm)
        result["convexity"] = bond.convexity_from_ytm(settle_dt, ytm)
        result["ytm"] = ytm
    elif clean_price is not None:
        # Yield from price
        ytm_calc = bond.yield_to_maturity(settle_dt, clean_price)
        result["ytm"] = ytm_calc
        result["clean_price"] = clean_price
        result["dirty_price"] = bond.dirty_price_from_ytm(settle_dt, ytm_calc)
        result["accrued_interest"] = bond.accrued_interest(settle_dt)
        result["macauley_duration"] = bond.macauley_duration(settle_dt, ytm_calc)
        result["modified_duration"] = bond.modified_duration(settle_dt, ytm_calc)
        result["convexity"] = bond.convexity_from_ytm(settle_dt, ytm_calc)
    else:
        raise ValueError("Provide either 'ytm' or 'clean_price'")

    return result


def _load_input_json() -> dict:
    parser = argparse.ArgumentParser(
        description="FinancePy runner: read bond inputs as JSON and emit JSON results."
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        help="Path to primary input JSON file (e.g. OasisLMF damage_estimates), "
             "or '-' to read from stdin (default).",
    )
    parser.add_argument(
        "frequency_input",
        nargs="?",
        default=None,
        help="Optional second input JSON file with lambda_annual "
             "(e.g. hurricane-frequency output). Merged into primary params.",
    )
    args = parser.parse_args()

    if args.input == "-":
        params = json.load(sys.stdin)
    else:
        with open(args.input, "r", encoding="utf-8") as handle:
            params = json.load(handle)

    if args.frequency_input:
        with open(args.frequency_input, "r", encoding="utf-8") as handle:
            freq = json.load(handle)
        # hurricane-frequency keys (lambda_annual, etc.) don't collide with
        # OasisLMF / DicePy keys, so a flat merge is safe.
        params.update(freq)

    return params


if __name__ == "__main__":
    input_data = _load_input_json()
    output_data = run(input_data)
    json.dump(output_data, sys.stdout, indent=2)
    print()  # trailing newline