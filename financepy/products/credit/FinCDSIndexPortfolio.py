##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import pow

from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes, FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinError import FinError
from ...products.credit.FinCDS import FinCDS
from ...products.credit.FinCDSCurve import FinCDSCurve
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinHelperFunctions import labelToString

###############################################################################
# TODO: Move index spread details into class and then pass in issuer curves
#       to the function when doing the adjustment
###############################################################################


class FinCDSIndexPortfolio():
    ''' This class manages the calculations associated with an equally weighted
    portfolio of CDS contracts with the same maturity date. '''

    def __init__(self,
                 freqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create FinCDSIndexPortfolio object. Note that all of the inputs
        have a default value which reflects the CDS market standard. '''

        checkArgumentTypes(self.__init__, locals())

        self._dayCountType = dayCountType
        self._dateGenRuleType = dateGenRuleType
        self._calendarType = calendarType
        self._freqType = freqType
        self._busDayAdjustType = busDayAdjustType

###############################################################################

    def intrinsicRPV01(self,
                       valuationDate,
                       stepInDate,
                       maturityDate,
                       issuerCurves):
        ''' Calculation of the risky PV01 of the CDS porfolio by taking the
        average of the risky PV01s of each contract. '''

        numCredits = len(issuerCurves)

        cdsContract = FinCDS(stepInDate,
                             maturityDate,
                             0.0)

        intrinsicRPV01 = 0.0

        for m in range(0, numCredits):

            retValue = cdsContract.riskyPV01(valuationDate,
                                             issuerCurves[m])

            cleanRPV01 = retValue['clean_rpv01']

            intrinsicRPV01 += cleanRPV01

        intrinsicRPV01 /= numCredits
        return(intrinsicRPV01)

###############################################################################

    def intrinsicProtectionLegPV(self,
                                 valuationDate,
                                 stepInDate,
                                 maturityDate,
                                 issuerCurves):
        ''' Calculation of intrinsic protection leg value of the CDS porfolio
        by taking the average sum the protection legs of each contract. '''

        numCredits = len(issuerCurves)

        intrinsicProtPV = 0.0

        # All contracts have same flows so only need one object
        cdsContract = FinCDS(stepInDate,
                             maturityDate,
                             0.0,
                             1.0)

        for m in range(0, numCredits):

            protectionPV = cdsContract.protectionLegPV(valuationDate,
                                                       issuerCurves[m])

            intrinsicProtPV += protectionPV

        intrinsicProtPV /= numCredits
        return intrinsicProtPV

###############################################################################

    def intrinsicSpread(self,
                        valuationDate,
                        stepInDate,
                        maturityDate,
                        issuerCurves):
        ''' Calculation of the intrinsic spread of the CDS portfolio as the one
        which would make the value of the protection legs equal to the value of
        the premium legs if all premium legs paid the same spread. '''

        intrinsicProtPV = self.intrinsicProtectionLegPV(valuationDate,
                                                        stepInDate,
                                                        maturityDate,
                                                        issuerCurves)

        intrinsicRPV01 = self.intrinsicRPV01(valuationDate,
                                             stepInDate,
                                             maturityDate,
                                             issuerCurves)

        intrinsicSpread = intrinsicProtPV / intrinsicRPV01

        return(intrinsicSpread)

###############################################################################

    def averageSpread(self,
                      valuationDate,
                      stepInDate,
                      maturityDate,
                      issuerCurves):
        ''' Calculates the average par CDS spread of the CDS portfolio. '''

        numCredits = len(issuerCurves)

        cdsContract = FinCDS(stepInDate,
                             maturityDate,
                             0.0)

        averageSpread = 0.0

        for m in range(0, numCredits):
            spread = cdsContract.parSpread(valuationDate, issuerCurves[m])
            averageSpread += spread

        averageSpread /= numCredits
        return averageSpread

###############################################################################

    def totalSpread(self,
                    valuationDate,
                    stepInDate,
                    maturityDate,
                    issuerCurves):
        ''' Calculates the total CDS spread of the CDS portfolio by summing
        over all of the issuers and adding the spread with no weights. '''

        numCredits = len(issuerCurves)

        cdsContract = FinCDS(stepInDate,
                             maturityDate,
                             0.0)

        totalSpread = 0.0

        for m in range(0, numCredits):
            spread = cdsContract.parSpread(valuationDate, issuerCurves[m])
            totalSpread += spread

        return totalSpread

###############################################################################

    def minSpread(self,
                  valuationDate,
                  stepInDate,
                  maturityDate,
                  issuerCurves):
        ''' Calculates the minimum par CDS spread across all of the issuers in
        the CDS portfolio. '''

        numCredits = len(issuerCurves)

        if numCredits < 1:
            raise FinError("Number of credits in index must be > 1 and not" + str(numCredits))

        cdsContract = FinCDS(stepInDate,
                             maturityDate,
                             0.0)

        minSpread = cdsContract.parSpread(valuationDate, issuerCurves[0])

        for m in range(1, numCredits):
            spread = cdsContract.parSpread(valuationDate, issuerCurves[m])
            if spread < minSpread:
                minSpread = spread

        return minSpread

###############################################################################

    def maxSpread(self,
                  valuationDate,
                  stepInDate,
                  maturityDate,
                  issuerCurves):
        ''' Calculates the maximum par CDS spread across all of the issuers in
        the CDS portfolio. '''

        numCredits = len(issuerCurves)

        if numCredits < 1:
            raise FinError("Number of credits in index must be > 1 and not " + str(numCredits))

        cdsContract = FinCDS(stepInDate,
                             maturityDate,
                             0.0)

        maxSpread = cdsContract.parSpread(valuationDate, issuerCurves[0])

        for m in range(1, numCredits):
            spread = cdsContract.parSpread(valuationDate, issuerCurves[m])
            if spread > maxSpread:
                maxSpread = spread

        return maxSpread

###############################################################################

    def spreadAdjustIntrinsic(self,
                              valuationDate,
                              issuerCurves,
                              indexCoupons,
                              indexUpfronts,
                              indexMaturityDates,
                              indexRecoveryRate,
                              tolerance=1e-6):
        ''' Adjust individual CDS curves to reprice CDS index prices.
        This approach uses an iterative scheme but is slow as it has to use a
        CDS curve bootstrap required when each trial spread adjustment is made
        '''

        numCredits = len(issuerCurves)

        if numCredits < 1:
            raise FinError("Number of credits in index must be > 1 and not " + str(numCredits))

        liborCurve = issuerCurves[0]._liborCurve
        numIndexMaturityPoints = len(indexCoupons)

        cdsMaturityDates = []
        for cds in issuerCurves[0]._cdsContracts:
            cdsDates = cds._maturityDate
            cdsMaturityDates.append(cdsDates)

        numCDSMaturityPoints = len(cdsMaturityDates)

        for issuerCurve in issuerCurves:
            n = len(issuerCurve._cdsContracts)
            if n != len(cdsMaturityDates):
                raise FinError(
                    "All issuer curves must be built from same cds maturities")

        cdsSpreadMultipliers = [1.0] * numCDSMaturityPoints

#        spreadDifference = [0.0] * numCDSMaturityPoints

        adjustedCDSSpreads = [0.0] * numCDSMaturityPoints

        #######################################################################
        # Set up CDS contracts used to build curve
        #######################################################################

        curveCDSContracts = []

        for j in range(0, numCDSMaturityPoints):

            cdsCoupon = 1.0

            cdsContract = FinCDS(valuationDate,
                                 cdsMaturityDates[j],
                                 cdsCoupon)

            curveCDSContracts.append(cdsContract)

        #######################################################################

        # We calibrate the individual CDS curves to fit each index maturity
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
                cdsIndex = FinCDS(valuationDate, indexMaturityDate, 0.0, 1.0)

                for iCredit in range(0, numCredits):

                    cdsContracts = issuerCurves[iCredit]._cdsContracts
                    recoveryRate = issuerCurves[iCredit]._recoveryRate
                    adjustedCDSContracts = []

                    for j in range(0, numCDSMaturityPoints):

                        cdsSpread = cdsContracts[j]._runningCoupon
                        adjustedCDSSpreads[j] = cdsSpread * \
                            cdsSpreadMultipliers[j]
                        curveCDSContracts[j]._runningCoupon = adjustedCDSSpreads[j]

                    adjustedIssuerCurve = FinCDSCurve(valuationDate,
                                                      curveCDSContracts,
                                                      liborCurve,
                                                      recoveryRate)

                    indexProtectionPV = cdsIndex.protectionLegPV(valuationDate,
                                                                 adjustedIssuerCurve,
                                                                 indexRecoveryRate)

                    cleanRPV01 = cdsIndex.riskyPV01(valuationDate,
                                                    adjustedIssuerCurve)['clean_rpv01']

                    sumRPV01 += cleanRPV01
                    sumProt += indexProtectionPV

                sumRPV01 /= numCredits
                sumProt /= numCredits

                sumPrem = sumRPV01 * indexCoupons[iMaturity]

                numerator = indexUpfronts[iMaturity] + sumPrem
                denominator = sumProt

                alpha = numerator / denominator
                cdsSpreadMultipliers[iMaturity] *= alpha

        # use spread multipliers to build and store adjusted curves
        adjustedIssuerCurves = []

        for iCredit in range(0, numCredits):

            recoveryRate = issuerCurves[iCredit]._recoveryRate

            adjustedCDSContracts = []
            adjustedSpreads = []

            for j in range(0, numCDSMaturityPoints):

                unadjustedSpread = issuerCurves[iCredit]._cdsContracts[j]._runningCoupon

                adjustedSpread = unadjustedSpread * cdsSpreadMultipliers[j]

                adjustedcdsContract = FinCDS(valuationDate,
                                             cdsMaturityDates[j],
                                             adjustedSpread)

                adjustedCDSContracts.append(adjustedcdsContract)
                adjustedSpreads.append(adjustedSpread)

                adjustedIssuerCurve = FinCDSCurve(valuationDate,
                                                  adjustedCDSContracts,
                                                  liborCurve,
                                                  recoveryRate)

            adjustedIssuerCurves.append(adjustedIssuerCurve)

        return adjustedIssuerCurves

###############################################################################

    def hazardRateAdjustIntrinsic(self,
                                  valuationDate,
                                  issuerCurves,
                                  indexCoupons,
                                  indexUpfronts,
                                  indexMaturityDates,
                                  indexRecoveryRate,
                                  tolerance=1e-6,
                                  maxIterations=100):
        ''' Adjust individual CDS curves to reprice CDS index prices.
        This approach adjusts the hazard rates and so avoids the slowish
        CDS curve bootstrap required when a spread adjustment is made.'''
        numCredits = len(issuerCurves)

        if numCredits < 1:
            raise FinError("Number of credits must be greater than 1")

        liborCurve = issuerCurves[0]._liborCurve
        numIndexMaturityPoints = len(indexCoupons)
#        hazardRateMultipliers = [1.0] * numIndexMaturityPoints
        adjustedIssuerCurves = []

        # making a copy of the issuer curves
        for issuerCurve in issuerCurves:

            adjustedIssuerCurve = FinCDSCurve(valuationDate,
                                              [],
                                              liborCurve)

            adjustedIssuerCurve._times = issuerCurve._times.copy()
            adjustedIssuerCurve._values = issuerCurve._values.copy()
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

                for iCredit in range(0, numCredits):

                    q1 = adjustedIssuerCurves[iCredit]._values[iMaturity]
                    q2 = adjustedIssuerCurves[iCredit]._values[iMaturity + 1]
                    q12 = q2 / q1
                    q12NEW = pow(q12, ratio)
                    q2NEW = q1 * q12NEW

                    adjustedIssuerCurves[iCredit]._values[iMaturity + 1] = q2NEW

                    indexMaturityDate = indexMaturityDates[iMaturity]

                    # the CDS spreads we extract here should be the index
                    # maturity dates
                    cdsIndex = FinCDS(valuationDate, indexMaturityDate, 0, 1.0)

                    indexProtPV = cdsIndex.protectionLegPV(valuationDate,
                                                           adjustedIssuerCurves[iCredit],
                                                           indexRecoveryRate)

                    rpv01Ret = cdsIndex.riskyPV01(
                        valuationDate, adjustedIssuerCurves[iCredit])

                    cleanRPV01 = rpv01Ret['clean_rpv01']

                    sumRPV01 += cleanRPV01
                    sumProt += indexProtPV

                sumRPV01 /= numCredits
                sumProt /= numCredits
                sumPrem = sumRPV01 * indexCoupons[iMaturity]

                numerator = indexUpfronts[iMaturity] + sumPrem
                denominator = sumProt

                ratio = numerator / denominator
                alpha *= ratio

        return adjustedIssuerCurves

###############################################################################

    def __repr__(self):

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("FREQUENCY", self._freqType)
        s += labelToString("DAYCOUNT", self._dayCountType)
        s += labelToString("CALENDAR", self._calendarType)
        s += labelToString("BUSDAYRULE", self._busDayAdjustType)
        s += labelToString("DATEGENRULE", self._dateGenRuleType)
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
