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
    """ This class manages the calculations associated with an equally weighted
    portfolio of CDS contracts with the same maturity date. """

    def __init__(self,
                 freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 day_count_type: DayCountTypes = DayCountTypes.ACT_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create FinCDSIndexPortfolio object. Note that all of the inputs
        have a default value which reflects the CDS market standard. """

        check_argument_types(self.__init__, locals())

        self._day_count_type = day_count_type
        self._date_gen_rule_type = date_gen_rule_type
        self._calendar_type = calendar_type
        self._freq_type = freq_type
        self._bus_day_adjust_type = bus_day_adjust_type

    ###########################################################################

    def intrinsic_rpv01(self,
                        valuation_date,
                        step_in_date,
                        maturity_date,
                        issuer_curves):
        """ Calculation of the risky PV01 of the CDS portfolio by taking the
        average of the risky PV01s of each contract. """

        num_credits = len(issuer_curves)

        cds_contract = CDS(step_in_date,
                           maturity_date,
                           0.0)

        intrinsic_rpv01 = 0.0

        for m in range(0, num_credits):
            retValue = cds_contract.risky_pv01(valuation_date,
                                               issuer_curves[m])

            cleanRPV01 = retValue['clean_rpv01']

            intrinsic_rpv01 += cleanRPV01

        intrinsic_rpv01 /= num_credits
        return (intrinsic_rpv01)

    ###########################################################################

    def intrinsic_protection_leg_pv(self,
                                    valuation_date,
                                    step_in_date,
                                    maturity_date,
                                    issuer_curves):
        """ Calculation of intrinsic protection leg value of the CDS portfolio
        by taking the average sum the protection legs of each contract. """

        num_credits = len(issuer_curves)

        intrinsic_prot_pv = 0.0

        # All contracts have same flows so only need one object
        cds_contract = CDS(step_in_date,
                           maturity_date,
                           0.0,
                           1.0)

        for m in range(0, num_credits):
            protectionPV = cds_contract.protection_leg_pv(valuation_date,
                                                          issuer_curves[m])

            intrinsic_prot_pv += protectionPV

        intrinsic_prot_pv /= num_credits
        return intrinsic_prot_pv

    ###########################################################################

    def intrinsic_spread(self,
                         valuation_date,
                         step_in_date,
                         maturity_date,
                         issuer_curves):
        """ Calculation of the intrinsic spread of the CDS portfolio as the one
        which would make the value of the protection legs equal to the value of
        the premium legs if all premium legs paid the same spread. """

        intrinsic_prot_pv = self.intrinsic_protection_leg_pv(valuation_date,
                                                             step_in_date,
                                                             maturity_date,
                                                             issuer_curves)

        intrinsic_rpv01 = self.intrinsic_rpv01(valuation_date,
                                               step_in_date,
                                               maturity_date,
                                               issuer_curves)

        intrinsic_spread = intrinsic_prot_pv / intrinsic_rpv01

        return (intrinsic_spread)

    ###########################################################################

    def average_spread(self,
                       valuation_date,
                       step_in_date,
                       maturity_date,
                       issuer_curves):
        """ Calculates the average par CDS spread of the CDS portfolio. """

        num_credits = len(issuer_curves)

        cds_contract = CDS(step_in_date,
                           maturity_date,
                           0.0)

        average_spread = 0.0

        for m in range(0, num_credits):
            spread = cds_contract.par_spread(valuation_date, issuer_curves[m])
            average_spread += spread

        average_spread /= num_credits
        return average_spread

    ###########################################################################

    def total_spread(self,
                     valuation_date,
                     step_in_date,
                     maturity_date,
                     issuer_curves):
        """ Calculates the total CDS spread of the CDS portfolio by summing
        over all of the issuers and adding the spread with no weights. """

        num_credits = len(issuer_curves)

        cds_contract = CDS(step_in_date,
                           maturity_date,
                           0.0)

        totalSpread = 0.0

        for m in range(0, num_credits):
            spread = cds_contract.par_spread(valuation_date, issuer_curves[m])
            totalSpread += spread

        return totalSpread

    ###########################################################################

    def min_spread(self,
                   valuation_date,
                   step_in_date,
                   maturity_date,
                   issuer_curves):
        """ Calculates the minimum par CDS spread across all of the issuers in
        the CDS portfolio. """

        num_credits = len(issuer_curves)

        if num_credits < 1:
            raise FinError(
                "Number of credits in index must be > 1 and not" + str(num_credits))

        cds_contract = CDS(step_in_date,
                           maturity_date,
                           0.0)

        min_spread = cds_contract.par_spread(valuation_date, issuer_curves[0])

        for m in range(1, num_credits):
            spread = cds_contract.par_spread(valuation_date, issuer_curves[m])
            if spread < min_spread:
                min_spread = spread

        return min_spread

    ###########################################################################

    def max_spread(self,
                   valuation_date,
                   step_in_date,
                   maturity_date,
                   issuer_curves):
        """ Calculates the maximum par CDS spread across all of the issuers in
        the CDS portfolio. """

        num_credits = len(issuer_curves)

        if num_credits < 1:
            raise FinError(
                "Number of credits in index must be > 1 and not " + str(num_credits))

        cds_contract = CDS(step_in_date,
                           maturity_date,
                           0.0)

        max_spread = cds_contract.par_spread(valuation_date, issuer_curves[0])

        for m in range(1, num_credits):
            spread = cds_contract.par_spread(valuation_date, issuer_curves[m])
            if spread > max_spread:
                max_spread = spread

        return max_spread

    ###########################################################################

    def spread_adjust_intrinsic(self,
                                valuation_date,
                                issuer_curves,
                                index_coupons,
                                indexUpfronts,
                                indexMaturityDates,
                                indexRecoveryRate,
                                tolerance=1e-6):
        """ Adjust individual CDS discount to reprice CDS index prices.
        This approach uses an iterative scheme but is slow as it has to use a
        CDS curve bootstrap required when each trial spread adjustment is made
        """

        num_credits = len(issuer_curves)

        if num_credits < 1:
            raise FinError(
                "Number of credits in index must be > 1 and not " + str(num_credits))

        libor_curve = issuer_curves[0]._libor_curve
        numIndexMaturityPoints = len(index_coupons)

        cdsMaturityDates = []
        for cds in issuer_curves[0]._cds_contracts:
            cdsDates = cds._maturity_date
            cdsMaturityDates.append(cdsDates)

        numCDSMaturityPoints = len(cdsMaturityDates)

        for issuer_curve in issuer_curves:
            n = len(issuer_curve._cds_contracts)
            if n != len(cdsMaturityDates):
                raise FinError(
                    "All issuer discount must be built from same cds maturities")

        cdsSpreadMultipliers = [1.0] * numCDSMaturityPoints

        #        spreadDifference = [0.0] * numCDSMaturityPoints

        adjustedCDSSpreads = [0.0] * numCDSMaturityPoints

        #######################################################################
        # Set up CDS contracts used to build curve
        #######################################################################

        curveCDSContracts = []

        for j in range(0, numCDSMaturityPoints):
            cdsCoupon = 1.0

            cds_contract = CDS(valuation_date,
                               cdsMaturityDates[j],
                               cdsCoupon)

            curveCDSContracts.append(cds_contract)

        #######################################################################

        # We calibrate the individual CDS discount to fit each index maturity
        # point
        for iMaturity in range(0, numIndexMaturityPoints):

            alpha = 0.0
            numIterations = 0

            while abs(alpha - 1.0) > tolerance:

                numIterations += 1

                if numIterations > 20:
                    raise FinError(
                        "Num iterations > 20. Increase limit or reduce tolerance or check inputs.")

                sumRPV01 = 0.0
                sumProt = 0.0

                # This is for the specific index maturity date
                indexMaturityDate = indexMaturityDates[iMaturity]
                cdsIndex = CDS(valuation_date, indexMaturityDate, 0.0, 1.0)

                for iCredit in range(0, num_credits):

                    cds_contracts = issuer_curves[iCredit]._cds_contracts
                    recovery_rate = issuer_curves[iCredit]._recovery_rate
                    adjustedCDSContracts = []

                    for j in range(0, numCDSMaturityPoints):
                        cdsSpread = cds_contracts[j]._running_coupon
                        adjustedCDSSpreads[j] = cdsSpread * \
                            cdsSpreadMultipliers[j]
                        curveCDSContracts[j]._running_coupon = adjustedCDSSpreads[j]

                    adjustedIssuerCurve = CDSCurve(valuation_date,
                                                   curveCDSContracts,
                                                   libor_curve,
                                                   recovery_rate)

                    indexProtectionPV = cdsIndex.protection_leg_pv(valuation_date,
                                                                   adjustedIssuerCurve,
                                                                   indexRecoveryRate)

                    cleanRPV01 = cdsIndex.risky_pv01(valuation_date,
                                                     adjustedIssuerCurve)['clean_rpv01']

                    sumRPV01 += cleanRPV01
                    sumProt += indexProtectionPV

                sumRPV01 /= num_credits
                sumProt /= num_credits

                sumPrem = sumRPV01 * index_coupons[iMaturity]

                numerator = indexUpfronts[iMaturity] + sumPrem
                denominator = sumProt

                alpha = numerator / denominator
                cdsSpreadMultipliers[iMaturity] *= alpha

        # use spread multipliers to build and store adjusted discount
        adjustedIssuerCurves = []

        for iCredit in range(0, num_credits):

            recovery_rate = issuer_curves[iCredit]._recovery_rate

            adjustedCDSContracts = []
            adjustedSpreads = []

            for j in range(0, numCDSMaturityPoints):
                unadjustedSpread = issuer_curves[iCredit]._cds_contracts[j]._running_coupon

                adjustedSpread = unadjustedSpread * cdsSpreadMultipliers[j]

                adjustedcds_contract = CDS(valuation_date,
                                           cdsMaturityDates[j],
                                           adjustedSpread)

                adjustedCDSContracts.append(adjustedcds_contract)
                adjustedSpreads.append(adjustedSpread)

                adjustedIssuerCurve = CDSCurve(valuation_date,
                                               adjustedCDSContracts,
                                               libor_curve,
                                               recovery_rate)

            adjustedIssuerCurves.append(adjustedIssuerCurve)

        return adjustedIssuerCurves

    ###########################################################################

    def hazard_rate_adjust_intrinsic(self,
                                     valuation_date,
                                     issuer_curves,
                                     index_coupons,
                                     indexUpfronts,
                                     indexMaturityDates,
                                     indexRecoveryRate,
                                     tolerance=1e-6,
                                     maxIterations=100):
        """ Adjust individual CDS discount to reprice CDS index prices.
        This approach adjusts the hazard rates and so avoids the slowish
        CDS curve bootstrap required when a spread adjustment is made."""
        num_credits = len(issuer_curves)

        if num_credits < 1:
            raise FinError("Number of credits must be greater than 1")

        libor_curve = issuer_curves[0]._libor_curve
        numIndexMaturityPoints = len(index_coupons)
        #        hazardRateMultipliers = [1.0] * numIndexMaturityPoints
        adjustedIssuerCurves = []

        # making a copy of the issuer discount
        for issuer_curve in issuer_curves:
            adjustedIssuerCurve = CDSCurve(valuation_date,
                                           [],
                                           libor_curve)

            adjustedIssuerCurve._times = issuer_curve._times.copy()
            adjustedIssuerCurve._values = issuer_curve._values.copy()
            adjustedIssuerCurves.append(adjustedIssuerCurve)

        # We solve for each maturity point
        for iMaturity in range(0, numIndexMaturityPoints):

            alpha = 1.0
            ratio = 1.0 + 2.0 * tolerance
            numIterations = 0

            while abs(ratio - 1.0) > tolerance:

                numIterations += 1

                if numIterations == maxIterations:
                    raise FinError("Max Iterations exceeded")

                sumRPV01 = 0.0
                sumProt = 0.0

                for iCredit in range(0, num_credits):
                    q1 = adjustedIssuerCurves[iCredit]._values[iMaturity]
                    q2 = adjustedIssuerCurves[iCredit]._values[iMaturity + 1]
                    q12 = q2 / q1
                    q12NEW = pow(q12, ratio)
                    q2NEW = q1 * q12NEW

                    adjustedIssuerCurves[iCredit]._values[iMaturity + 1] = q2NEW

                    indexMaturityDate = indexMaturityDates[iMaturity]

                    # the CDS spreads we extract here should be the index
                    # maturity dates
                    cdsIndex = CDS(valuation_date, indexMaturityDate, 0, 1.0)

                    indexProtPV = cdsIndex.protection_leg_pv(valuation_date,
                                                             adjustedIssuerCurves[iCredit],
                                                             indexRecoveryRate)

                    rpv01Ret = cdsIndex.risky_pv01(
                        valuation_date, adjustedIssuerCurves[iCredit])

                    cleanRPV01 = rpv01Ret['clean_rpv01']

                    sumRPV01 += cleanRPV01
                    sumProt += indexProtPV

                sumRPV01 /= num_credits
                sumProt /= num_credits
                sumPrem = sumRPV01 * index_coupons[iMaturity]

                numerator = indexUpfronts[iMaturity] + sumPrem
                denominator = sumProt

                ratio = numerator / denominator
                alpha *= ratio

        return adjustedIssuerCurves

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("DAYCOUNT", self._day_count_type)
        s += label_to_string("CALENDAR", self._calendar_type)
        s += label_to_string("BUSDAYRULE", self._bus_day_adjust_type)
        s += label_to_string("DATEGENRULE", self._date_gen_rule_type)
        return s

    ###########################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###########################################################################
