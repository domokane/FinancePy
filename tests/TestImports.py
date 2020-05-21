###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


# from financepy.finutils.FinDate import FinDate # Works
# from ..financepy import *  # fails


from financepy.finutils import *
from financepy.products.libor import *
from financepy.market.curves import *

def test_Imports():

    dt = FinDate(4, 7, 2019)

    curve = FinFlatCurve(dt, 0.05, 1)

    print(dt)

###############################################################################

test_Imports()
