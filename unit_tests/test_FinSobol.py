# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.models.sobol import get_uniform_sobol

########################################################################################


def test__fin_sobol():

    num_points = 1000
    dimensions = 3

    points = get_uniform_sobol(num_points, dimensions)

    for d in range(dimensions):
        av = 0.0
        var = 0.0

        for point in points[:, d]:
            av += point
            var += point**2

        av /= num_points
        var /= num_points

        av_error = abs(av - (1 / 2))
        var_error = abs(var - (1 / 3))
        assert av_error < 0.002
        assert var_error < 0.002
