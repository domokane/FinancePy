import numpy as np

from financepy.utils.helpers import grid_index


def test_grid_index_returns_match_without_printing(capsys):
    assert grid_index(0.2, np.array([0.0, 0.2, 0.5])) == 1
    assert capsys.readouterr().out == ""
