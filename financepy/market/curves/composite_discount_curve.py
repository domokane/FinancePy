###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


import numpy as np
from typing import List, Union

###############################################################################

from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.helpers import label_to_string
from ...utils.helpers import check_argument_types
from ...market.curves.discount_curve import DiscountCurve


class CompositeDiscountCurve(DiscountCurve):
    """
    A discount curve that is a sum (in rates) of 'children' discount curves
    """

    ###############################################################################

    def __init__(self, child_curves: List[DiscountCurve]):
        """
        Create a discount curve that is a sum (in rates) of other discount curves
        """

        check_argument_types(self.__init__, locals())
        assert (
            len(child_curves) > 0
        ), "Empty list of child curves is not supported"

        self._children = child_curves

        self.value_dt = self._children[0].value_dt
        assert all(
            c.value_dt == self.value_dt for c in self._children
        ), "Child curves must all have the same vlauation date"

        # Read off the first child
        self.dc_type = self._children[0].dc_type

    ###############################################################################

    def df_t(self, t: Union[float, np.ndarray]):
        """
        Return discount factors given a single or vector of dates.
        ParentRate = Sum of children rates => Parent DF = product of children dfs
        """

        dfs = np.ones_like(np.atleast_1d(t), dtype=float)
        for c in self._children:
            dfc = c.df_t(t)
            dfs *= dfc

        return dfs

    ###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("CHILDREN", (self._children))
        return s

    ###############################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
