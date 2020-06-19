from financepy.finutils.FinSobol import generateSobol

def test_FinSobol():
    
    numPoints = 1000
    dimensions = 10
    
    points = generateSobol(numPoints, dimensions)

    for d in range(dimensions):
        av = 0.0
        var = 0.0
        
        for point in points[:, d]:
            av += point
            var += point ** 2
            
        av /= numPoints
        var /= numPoints
        
        avError = abs(av - (1/2))
        varError = abs(var - (1/3))
        assert(avError < 0.002)
        assert(varError < 0.002)

test_FinSobol()