########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import sys, os

# Add the *parent of the script's parent directory* to sys.path
sys.path.append(
    os.path.dirname(  # go up two directories
        os.path.abspath(__file__)  # absolute path to current file
    )
)
