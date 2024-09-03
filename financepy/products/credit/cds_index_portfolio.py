##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import pow

from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.error import FinError
from ...products.credit.cds import CDS
from ...products.credit.cds_curve import CDSCurve
from ...utils.helpers import check_argument_types
from ...utils.helpers import label_to_string


###########################################################################
# TODO: Move index spread details into class and then pass in issuer discount
#       to the function when doing the adjustment
###########################################################################


class CDSIndexPortfolio:
    """This class manages the calculations associated with an equally weighted
    portfolio of CDS contracts with the same maturity date."""

    def __init__(
        self,
        freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        day_count_type: DayCountTypes = DayCountTypes.ACT_360,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Create Fincds_indexPortfolio object. Note that all of the inputs
        have a default value which reflects the CDS market standard."""

        check_argument_types(self.__init__, locals())

        self.dc_type = day_count_type
        self.dg_type = dg_type
        self.cal_type = cal_type
        self.freq_type = freq_type
        self.bd_type = bd_type

    ###########################################################################

    def intrinsic_rpv01(
        self, value_dt, step_in_dt, maturity_dt, issuer_curves
    ):
        """Calculation of the risky PV01 of the CDS portfolio by taking the
        average of the risky PV01s of each contract."""

        num_credits = len(issuer_curves)

        cds_contract = CDS(step_in_dt, maturity_dt, 0.0)

        intrinsic_rpv01 = 0.0

        for m in range(0, num_credits):
            ret_value = cds_contract.risky_pv01(value_dt, issuer_curves[m])

            clean_rpv01 = ret_value["clean_rpv01"]

            intrinsic_rpv01 += clean_rpv01

        intrinsic_rpv01 /= num_credits
        return intrinsic_rpv01

    ###########################################################################

    def intrinsic_prot_leg_pv(
        self, value_dt, step_in_dt, maturity_dt, issuer_curves
    ):
        """Calculation of intrinsic protection leg value of the CDS portfolio
        by taking the average sum the protection legs of each contract."""

        num_credits = len(issuer_curves)

        intrinsic_prot_pv = 0.0

        # All contracts have same flows so only need one object
        cds_contract = CDS(step_in_dt, maturity_dt, 0.0, 1.0)

        for m in range(0, num_credits):
            prot_pv = cds_contract.prot_leg_pv(value_dt, issuer_curves[m])

            intrinsic_prot_pv += prot_pv

        intrinsic_prot_pv /= num_credits
        return intrinsic_prot_pv

    ###########################################################################

    def intrinsic_spread(
        self, value_dt, step_in_dt, maturity_dt, issuer_curves
    ):
        """Calculation of the intrinsic spread of the CDS portfolio as the one
        which would make the value of the protection legs equal to the value of
        the premium legs if all premium legs paid the same spread."""

        intrinsic_prot_pv = self.intrinsic_prot_leg_pv(
            value_dt, step_in_dt, maturity_dt, issuer_curves
        )

        intrinsic_rpv01 = self.intrinsic_rpv01(
            value_dt, step_in_dt, maturity_dt, issuer_curves
        )

        intrinsic_spread = intrinsic_prot_pv / intrinsic_rpv01

        return intrinsic_spread

    ###########################################################################

    def average_spread(self, value_dt, step_in_dt, maturity_dt, issuer_curves):
        """Calculates the average par CDS spread of the CDS portfolio."""

        num_credits = len(issuer_curves)

        cds_contract = CDS(step_in_dt, maturity_dt, 0.0)

        average_spread = 0.0

        for m in range(0, num_credits):
            spread = cds_contract.par_spread(value_dt, issuer_curves[m])
            average_spread += spread

        average_spread /= num_credits
        return average_spread

    ###########################################################################

    def total_spread(self, value_dt, step_in_dt, maturity_dt, issuer_curves):
        """Calculates the total CDS spread of the CDS portfolio by summing
        over all of the issuers and adding the spread with no weights."""

        num_credits = len(issuer_curves)

        cds_contract = CDS(step_in_dt, maturity_dt, 0.0)

        total_spd = 0.0

        for m in range(0, num_credits):
            spread = cds_contract.par_spread(value_dt, issuer_curves[m])
            total_spd += spread

        return total_spd

    ###########################################################################

    def min_spread(self, value_dt, step_in_dt, maturity_dt, issuer_curves):
        """Calculates the minimum par CDS spread across all of the issuers in
        the CDS portfolio."""

        num_credits = len(issuer_curves)

        if num_credits < 1:
            raise FinError(
                "Number of credits in index must be > 1 and not"
                + str(num_credits)
            )

        cds_contract = CDS(step_in_dt, maturity_dt, 0.0)

        min_spread = cds_contract.par_spread(value_dt, issuer_curves[0])

        for m in range(1, num_credits):
            spread = cds_contract.par_spread(value_dt, issuer_curves[m])
            if spread < min_spread:
                min_spread = spread

        return min_spread

    ###########################################################################

    def max_spread(self, value_dt, step_in_dt, maturity_dt, issuer_curves):
        """Calculates the maximum par CDS spread across all of the issuers in
        the CDS portfolio."""

        num_credits = len(issuer_curves)

        if num_credits < 1:
            raise FinError(
                "Number of credits in index must be > 1 and not "
                + str(num_credits)
            )

        cds_contract = CDS(step_in_dt, maturity_dt, 0.0)

        max_spread = cds_contract.par_spread(value_dt, issuer_curves[0])

        for m in range(1, num_credits):
            spread = cds_contract.par_spread(value_dt, issuer_curves[m])
            if spread > max_spread:
                max_spread = spread

        return max_spread

    ###########################################################################

    def spread_adjust_intrinsic(
        self,
        value_dt,
        issuer_curves,
        index_cpns,
        index_upfronts,
        index_maturity_dts,
        index_recovery_rate,
        tolerance=1e-6,
    ):
        """Adjust individual CDS discount to reprice CDS index prices.
        This approach uses an iterative scheme but is slow as it has to use a
        CDS curve bootstrap required when each trial spread adjustment is made
        """

        num_credits = len(issuer_curves)

        if num_credits < 1:
            raise FinError(
                "Number of credits in index must be > 1 and not "
                + str(num_credits)
            )

        libor_curve = issuer_curves[0].libor_curve
        num_index_maturity_points = len(index_cpns)

        cds_maturity_dts = []
        for cds in issuer_curves[0].cds_contracts:
            cds_dts = cds.maturity_dt
            cds_maturity_dts.append(cds_dts)

        num_cds_mat_points = len(cds_maturity_dts)

        for issuer_curve in issuer_curves:
            n = len(issuer_curve.cds_contracts)
            if n != len(cds_maturity_dts):
                raise FinError(
                    "All issuer discount must be from same cds maturities"
                )

        cds_spread_multipliers = [1.0] * num_cds_mat_points

        #        spreadDifference = [0.0] * num_cds_mat_points

        adjusted_cds_spreads = [0.0] * num_cds_mat_points

        #######################################################################
        # Set up CDS contracts used to build curve
        #######################################################################

        curve_cds_contracts = []

        for j in range(0, num_cds_mat_points):
            cds_cpn = 1.0

            cds_contract = CDS(value_dt, cds_maturity_dts[j], cds_cpn)

            curve_cds_contracts.append(cds_contract)

        #######################################################################
        # We calibrate the individual CDS discount to fit each index maturity
        #######################################################################

        for i_maturity in range(0, num_index_maturity_points):

            alpha = 0.0
            num_iterations = 0

            while abs(alpha - 1.0) > tolerance:

                num_iterations += 1

                if num_iterations > 20:
                    raise FinError("Num iterations > 20.")

                sum_rpv01 = 0.0
                sum_prot = 0.0

                # This is for the specific index maturity date
                index_maturity_dt = index_maturity_dts[i_maturity]
                cds_index = CDS(value_dt, index_maturity_dt, 0.0, 1.0)

                for i_credit in range(0, num_credits):

                    cds_contracts = issuer_curves[i_credit].cds_contracts
                    recovery_rate = issuer_curves[i_credit].recovery_rate
                    adjusted_cds_contracts = []

                    for j in range(0, num_cds_mat_points):
                        cds_spread = cds_contracts[j].running_cpn
                        adjusted_cds_spreads[j] = (
                            cds_spread * cds_spread_multipliers[j]
                        )
                        curve_cds_contracts[j].running_cpn = (
                            adjusted_cds_spreads[j]
                        )

                    adjusted_issuer_curve = CDSCurve(
                        value_dt,
                        curve_cds_contracts,
                        libor_curve,
                        recovery_rate,
                    )

                    index_prot_pv = cds_index.prot_leg_pv(
                        value_dt, adjusted_issuer_curve, index_recovery_rate
                    )

                    clean_rpv01 = cds_index.risky_pv01(
                        value_dt, adjusted_issuer_curve
                    )["clean_rpv01"]

                    sum_rpv01 += clean_rpv01
                    sum_prot += index_prot_pv

                sum_rpv01 /= num_credits
                sum_prot /= num_credits

                sum_prem = sum_rpv01 * index_cpns[i_maturity]

                numerator = index_upfronts[i_maturity] + sum_prem
                denominator = sum_prot

                alpha = numerator / denominator
                cds_spread_multipliers[i_maturity] *= alpha

        # use spread multipliers to build and store adjusted discount
        adjusted_issuer_curves = []

        for i_credit in range(0, num_credits):

            recovery_rate = issuer_curves[i_credit].recovery_rate

            adjusted_cds_contracts = []
            adjusted_spreads = []

            for j in range(0, num_cds_mat_points):

                unadjusted_spread = (
                    issuer_curves[i_credit].cds_contracts[j].running_cpn
                )

                adjusted_spread = unadjusted_spread * cds_spread_multipliers[j]

                adjusted_cds_contract = CDS(
                    value_dt, cds_maturity_dts[j], adjusted_spread
                )

                adjusted_cds_contracts.append(adjusted_cds_contract)
                adjusted_spreads.append(adjusted_spread)

                adjusted_issuer_curve = CDSCurve(
                    value_dt,
                    adjusted_cds_contracts,
                    libor_curve,
                    recovery_rate,
                )

            adjusted_issuer_curves.append(adjusted_issuer_curve)

        return adjusted_issuer_curves

    ###########################################################################

    def hazard_rate_adjust_intrinsic(
        self,
        value_dt,
        issuer_curves,
        index_cpns,
        index_up_fronts,
        index_maturity_dts,
        index_recovery_rate,
        tolerance=1e-6,
        max_iterations=200,
    ):
        """Adjust individual CDS discount to reprice CDS index prices.
        This approach adjusts the hazard rates and so avoids the slowish
        CDS curve bootstrap required when a spread adjustment is made."""

        if 1 == 0:
            print("=========================================")
            print(value_dt)
            print(index_cpns)
            print(index_up_fronts)
            print(index_maturity_dts)
            print(index_recovery_rate)

        num_credits = len(issuer_curves)

        if num_credits < 1:
            raise FinError("Number of credits must be greater than 1")

        libor_curve = issuer_curves[0].libor_curve
        num_index_maturity_points = len(index_cpns)

        #        hazardRateMultipliers = [1.0] * numIndexMaturityPoints

        adjusted_issuer_curves = []

        # making a copy of the issuer discount
        for issuer_curve in issuer_curves:
            adjusted_issuer_curve = CDSCurve(
                value_dt, [], libor_curve, index_recovery_rate
            )

            adjusted_issuer_curve._times = issuer_curve._times.copy()
            adjusted_issuer_curve._values = issuer_curve._values.copy()
            adjusted_issuer_curves.append(adjusted_issuer_curve)

        # We solve for each maturity point
        for i_maturity in range(0, num_index_maturity_points):

            alpha = 1.0
            ratio = 1.0 + 2.0 * tolerance
            num_iterations = 0

            while abs(ratio - 1.0) > tolerance:

                num_iterations += 1

                if num_iterations > max_iterations:
                    raise FinError("Max Iterations exceeded")

                sum_rpv01 = 0.0
                sum_prot = 0.0

                for i_credit in range(0, num_credits):
                    q1 = adjusted_issuer_curves[i_credit]._values[i_maturity]
                    q2 = adjusted_issuer_curves[i_credit]._values[
                        i_maturity + 1
                    ]
                    q12 = q2 / q1

                    q12_new = pow(q12, ratio)
                    q2_new = q1 * q12_new

                    adjusted_issuer_curves[i_credit]._values[
                        i_maturity + 1
                    ] = q2_new

                    #        if i_maturity == 0 and index_cpns[0] == 0.006:
                    #              print(i_credit, q1, q2_new)

                    index_maturity_dt = index_maturity_dts[i_maturity]

                    # the CDS spreads we extract here
                    # should be to the index maturity dates
                    cds_index = CDS(value_dt, index_maturity_dt, 0, 1.0)

                    index_prot_pv = cds_index.prot_leg_pv(
                        value_dt,
                        adjusted_issuer_curves[i_credit],
                        index_recovery_rate,
                    )

                    rpv01_ret = cds_index.risky_pv01(
                        value_dt, adjusted_issuer_curves[i_credit]
                    )

                    clean_rpv01 = rpv01_ret["clean_rpv01"]

                    sum_rpv01 += clean_rpv01
                    sum_prot += index_prot_pv

                sum_rpv01 /= num_credits
                sum_prot /= num_credits

                spd = sum_prot / sum_rpv01

                sum_prem = sum_rpv01 * index_cpns[i_maturity]

                numerator = index_up_fronts[i_maturity] + sum_prem
                denominator = sum_prot

                ratio = numerator / denominator
                alpha = alpha * ratio

        #  print("Maturity:", i_maturity, "Num:", numerator, "Den:",
        # denominator, "Ratio:", ratio, "Alpha:", alpha)

        return adjusted_issuer_curves

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("DAYCOUNT", self.dc_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUSDAYRULE", self.bd_type)
        s += label_to_string("DATEGENRULE", self.dg_type)
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###########################################################################
