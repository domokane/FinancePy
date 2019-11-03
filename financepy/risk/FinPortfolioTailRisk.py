# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 21:57:13 2019

@author: Dominic
"""
'''
public static  portfolioTailRiskAdjBinomial( k1,
                                     k2,
                                     regretLevel,
                                     confidence,
                                     tmat,
                                    numCredits,
                                     weightVector,
                                     hazardRateVector,
                                     recoveryVector,
                                     betaVector,
                                    numPoints)

#if ( DEBUG )
    //         using (ProfilerScoper profiler = new ProfilerScoper(MethodBase.GetCurrentMethod().Name))
    
#endif

        if (numCredits == 0)
            throw new Exception("Num credits is zero")

        if (k2 == 0)
            throw new Exception("k2 should not be zero")

        if (k1 >= k2)
            throw new Exception("K1 must be less than K2")

        HelpfulFunctions.renormaliseWeights(weightVector)

         defaultProbs = np.zeros(numCredits]
         lossDbn = np.zeros(numCredits + 1]

        for (j = 0 j < numCredits j++)
        
            defaultProbs[j] = 1.0 - Math.Exp(-hazardRateVector[j] * tmat)
        

         totalLoss = 0.0

         lossRatio = np.zeros(numCredits]

        for (i = 0 i < numCredits i++)
        
            totalLoss = totalLoss + weightVector[i] * (1.0 - recoveryVector[i])
        

         avgLoss = totalLoss / numCredits

        for (i = 0 i < numCredits i++)
        
            lossRatio[i] = weightVector[i] * (1.0 - recoveryVector[i]) / avgLoss
        

        LossDistributionBuilders.lossDbnHeterogeneousAdjBinomial(numCredits, defaultProbs, lossRatio, betaVector, numPoints, lossDbn)

         expectedShortfall = 0.0
         expectedLoss = 0.0
         expectedRegret = 0.0
         expectedTrancheLoss = 0.0

         tailProbabilityShortfall = 0.0

        numLossUnits = lossDbn.Length
         tailProbability = 0.0

        for (i = 0 i < numLossUnits i++)
        
            tailProbability += lossDbn[i]

             loss = totalLoss * i / (numLossUnits - 1)

             tailLoss = Math.Max(loss - k1, 0)

             trancheLoss = (Math.Max(loss - k1, 0) - Math.Max(loss - k2, 0)) / (k2 - k1)

            expectedLoss += loss * lossDbn[i]

            if (loss < regretLevel)
            
                expectedRegret += tailLoss * lossDbn[i]
            

            if (tailProbability >= confidence)
            
                expectedShortfall += tailLoss * lossDbn[i]
                tailProbabilityShortfall += lossDbn[i]
            

            expectedTrancheLoss += trancheLoss * lossDbn[i]
        

        expectedShortfall /= tailProbability

         output = np.zeros(5]
        output[0] = expectedLoss
        output[1] = expectedShortfall
        output[2] = expectedRegret
        output[3] = expectedTrancheLoss
        output[4] = tailProbabilityShortfall

        return output
'''
