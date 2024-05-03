## CHANGE LOG
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


