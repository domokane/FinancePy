# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


from typing import List, Union

import numpy as np

from ...utils.error import FinError
from ...utils.helpers import label_to_string
from ...utils.helpers import check_argument_types

from ...market.curves.discount_curve import DiscountCurve

###############################################################################


class CompositeDiscountCurve(DiscountCurve):
    """
    A discount curve that is a sum (in rates) of 'children' discount curves
    """

    ###########################################################################

    def __init__(self, child_curves: List[DiscountCurve]) -> None:
        """Create a discount curve that is a sum (in rates) of other
        discount curves.
        """

        check_argument_types(self.__init__, locals())

        if len(child_curves) == 0:
            raise FinError(
                "Empty list of child curves is not supported."
            )

        self._children = child_curves

        # All child curves must have the same anchor date and
        # use the same curve day-count convention.
        self.anchor_dt = self._children[0].anchor_dt
        self.curve_dc_type = self._children[0].curve_dc_type

        for curve in self._children:
            if curve.anchor_dt != self.anchor_dt:
                raise FinError(
                    "Child curves must have the same anchor date."
                )

            if curve.curve_dc_type != self.curve_dc_type:
                raise FinError(
                    "Child curves must have the same curve day count type."
                )

    ###########################################################################

    def df_t(self, t: Union[float, np.ndarray]):
        """
        Return discount factors given a single or vector of dates.
        ParentRate = Sum of children rates => Parent DF = product of
        children dfs
        """

        dfs = np.ones_like(np.atleast_1d(t), dtype=float)
        for c in self._children:
            dfc = c.df_t(t)
            dfs *= dfc

        return dfs

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("CHILDREN", (self._children))
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
