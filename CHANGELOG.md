## CHANGE LOG
21 August 2026
List of changes
- EquityChooserOption - enforced term structure for q and r
- EquityChooserOption - handled t = t_c
- EquityVolCurve - Amended constructor to hold all internal info
- EquityVolCurve - Added new and improved interpolation schemes
- EquityVolCurve - Added new and improved extrapolation schemes
- black_scholes_analytic - Added fwd, d1, and d2 functions
- black_scholes_analytic - fixed and vectorised American option analytic approximations
- black_scholes_hedging_sim - added new numba BS delta hedging simulator
- black_scholes_mc - moved option_type to last argument in line with black_scholes_analytic
- option_implied_dbn - improved
- EquityForwardStart - new product type
- math.py - vectorisation of M function for bivariates
- Fast solver in IborSingleCurve

10 August 2026
List of changes
- Renaming changes
    - FinCompoundingTypes renamed OISCompoundingTypes
    - FinProcessSimulator renamed ProcessSimulator
    - FinFXATMTypes renamed FXATMTypes
    - FinFXDeltaMethodTypes renamed FXDeltaMethodTypes
    - FinInflationIndexCurve renamed InflationIndexCurve
    - FinInflationSwap renamed InflationSwap
    - FinInflationBond renamed InflationBond
- Migration of types to global_types file
    - InterpTypes moved to global_types.py
    - VolFuncTypes moved to global_types.py
    - FXATMTypes moved to global_types.py
    - FXDeltaMethodTypes moved to global_types.py
    - BlackTypes moved to global_types.py
    - BlackScholesTypes moved to global_types.py
    - CIRNumericalSchemeTypes moved to global_types.py
    - HestonNumericalSchemeTypes moved to global_types
    - LMMModelTypes moved to global_types.py
    - GBMNumericalSchemeTypes moved to global_types.py
    - YTMCalcType moved to global_types.py
    - ProcessTypes moved to global_types.py
    - FXBarrierTypes moved to global_types.py
    - OISCompoundingTypes moved to global_types.py
    - HWEuropeanCalcTypes moved to global_types.py
    - DigitalOptionTypes moved to global_types.py


9 August 2026
- Added more calendars
    - AUSTRALIA_RITS
    - SINGAPORE
    - ZURICH
    - TOKYO
    - US GOVERNMENT SECURITIES
    - US FEDERAL RESERVE
    - TORONTO
    - HONK KONG
- Joint calendars are now possible for following classes
    - CDS
    - CDS_BASKET
    - CDS_TRANCHE
    - CDS_OPTION
    - EQUITY_SWAP
    - EQUITY_SWAP_LEG
    - IBOR_BASIS_SWAP
    - IBOR_FIXED_FLOAT_SWAP
    - IBOR_SWAP
    - IBOR_SWAPTION
    - FinFixedFixedXCCySwap
    - IborIborSwap
    - fixed_fixed_cross_ccy
    - fixed_float_cross_ccy_swap
    - float_float_cross_ccy_swap
    - Schedule
- Upgraded Amount class for potential future use
    - Amount now takes currency
    - Added overloaded math functionality
- Date changes
    - Removed datelite class
    - Moved date_arrays into date class
- Removed singleton class

6 August 2026
RELEASE OF FINANCEPY  V1.1.0
List of changes

- This release renames a wide set of the Discount curves
    - DiscountCurveFlat -> FlatDiscountCurve
    - DiscountCurveNS -> NSDiscountCurve
    - DiscountCurveNSS -> NSSDiscountCurve
    - DiscountCurvePoly -> PolyDiscountCurve
    - DiscountCurvePWFONF -> PWFONFDiscountCurve
    - DiscountCurvePWF -> PWFDiscountCurve
    - DiscountCurveZeros -> ZeroRatesDiscountCurve
- This release moves some product-based curves from product folders to market->curves
    - IborSingleCurve moves
    - IborCurveRiskEngine moves
    - IborSingleCurveParShocker moves
    - IborSingleCurveSmoothingCalibrator moves
    - IborDualCurve moves
    - OISCurve moves
    - CDSCurve moves
- A number of changes in CDS valuation
    - Renamed credit01 to spread dv01 like Bloomberg
    - Renamed rate01 to ir dv01 like Bloomberg
    - Renamed risky_pv01 function name rpv01
    - Removed dicts and define pair returns to be (DIRTY, CLEAN)
    - Improved accuracy of fast approximator
- Tweaked solver1D to be more gentle on convergence
- Added vectorised American pricers to black_scholes_analytic using numerical approximations
- FXBarrierModel code moved to models folder
- Renamed bond_callable.py to bond_embedded_option.py as this is what it contains
- Removed UNITED_KINGDOM from Calendars as UK is England + Wales + Scotland
- Added LONDON to Calendars for England to replace UNITED_KINGDOM
- Repaired intraday date functionality in Date class
- Updated all notebooks to ensure they work
- Updated regression and unit tests

###############################################################################

15 June 2026

- Added bond boostrap discount curve to do exact fit to bond prices
- Added bond parametric discount curve to do parametric best fit to bond prices
- Added bond parametric yield curve to do parametric best fit to yields
- Removed BondFittedZeroCurve (replaced by BondParametricDiscountCurve)
- Removed BondFittedYieldCurve (replaced by BondParametricYieldCurve)
- Removed BondExactZeroCurve (replaced by BondBootstrapDiscountCurve)
- Widespread switch to use times_from_dates to avoid G_DAYS_IN_YEAR to be completed
- Removed LINEAR_FWD_RATE interpolation as incorrect
- Added LINEAR_DISCOUNT interpolation
- Fixed optimiser use in FXVolSurfaces to handle precision
- Bond curve files all moved under market->curves

26 May 2026
- Cleaning up code in models folder - adding validation of inputs, type checking
- Fixed minor bug in Compound Option on tree
- Fixed bug in LSMC - I now regress on ITM paths only
- Tidied up finite_difference module

19 May 2026
- Further tidying up of Discount Curves folder
- Tidying up of vol curve and surfaces
- Extension of BDT tree to cover same cases as BK and HW trees
- Tidying up model code for Bachelier and Black-Scholes

18 May 2026
- Imposed bounds checks on discount factors
- Fixed FXForward value calculations

17 May 2026
- Correct all notebooks except FXBarrierOption which needs to be fixed
- Added methods for continuously compounded zero rates to discount curve

16 May 2026
- Rewrote Interpolator class to be easier to understand
- Ensured it passes test cases
- Reconfigured unit and golden tests. PLEASE SEE READMEs in each folder for new instructions

13 May 2026
- Mosly work on DiscountCurves and all its variants
- Extended inheritance of DiscountCurves to all curves
- All inheriting classes have a bump_parallel method
- Removed alternative calculations of zero, par, swap and forward rates so all such calls use DiscountCurve methods
- Added discount curve level property to determine day count rule for converting dates to times
- Set default value of time_dc_type to ACT_365F - which means we divide by 365
- Fixed bug in swap_float_leg so we use specified accrual basis and not index curve basis
- Set all references to call directly discountCurve._dfs and ._times to avoid copying using accessor methods
- Did same for credit curves to access ._qs survival probabilities
- Added compounding enum to handle rates and separate it from frequency

28 August 2025
- Completed pep8 cleanup
- Ensured caching on as many numba functions as possible
- Getter and Setter functions for discount curve times and discount factors

1 May 2024 version 0.360 released
- Fixed all notebooks to ensure they all work with current version
- Unit tests complete with success
- Gradually removing underscore prefix from class member variable names
- Adjustments to accrued interest calculations for FRN (need a consistent interface)

19 February 2024 version 0.350 released
- A lot of various pep8 fixes - should all be nearly done soon
- Fixed a bug in gauss_approx_tranche_loss

9 December 2023  version 0.34 released
- A lot of various pep8 fixes - should all be nearly done soon
- Some minor bug fixes

10 November 2023  version 0.33 released
- Tidied up key rate code
- Fixed unit tests for pytest
- Fixed vectorisation of barrier options
- Various pep8 fixes

28 August 2023 - version 0.32 released
- Fixed bug in Bond OAS and ASW

24 August 2023 - Version 0.31 released

- Schedule
  - Corrected bug in schedule generation
  - Corrected bug in CDS protection leg integral

- Many Bond Classes have been amended
  - Changed FULL price to DIRTY price in functions UPDATE YOUR CODE PLEASE. APOLS for inconvenience.
  - Removed face amount from bond class - how much you buy is not intrinsic to a bond
  - Made number of ex-dividend days a member of bond class
  - Added adjustment for ex-dividend dates to yield calculations
  - Revised accrued and principal functions to take face amount as input
  - Updated document

29 May 2023 - Version 0.30 released
- Added PrettyPrint to required dependencies

22 Nov 22
Version 0.260 has been released and pushed to PyPI
- Create Date from python datetime
- Zero coupon bond class
- Fixed bug in bond payment date

31-Aug-2022
Version 0.240 has just been released and pushed to PyPI with changes
- Negative terms in date class
- Recovery rates do not default to standard value for CDS curves


