import pytest

from financepy.utils.tenor import Tenor, TenorUnit


@pytest.mark.parametrize('tenor_string,num_periods,units',
                         [
                             (None, 0, TenorUnit.NONE),
                             ('', 0, TenorUnit.NONE),
                             ('on', 1, TenorUnit.DAYS),
                             ('TN', 1, TenorUnit.DAYS),
                             ('5d', 5, TenorUnit.DAYS),
                             ('-7w', -7, TenorUnit.WEEKS),
                             ('12M', 12, TenorUnit.MONTHS),
                             ('-3Y', -3, TenorUnit.YEARS),
                         ])
def test_tenor_create(tenor_string, num_periods, units):
    tenor = Tenor(tenor_string)
    assert tenor._num_periods == num_periods
    assert tenor._units.name == units.name


@pytest.mark.parametrize('tenor_string', ['5d', '-7w', '12M', '-3Y', ])
def test_tenor_to_string(tenor_string):
    tenor = Tenor(tenor_string)
    assert f'{tenor}' == tenor_string.upper()


@pytest.mark.parametrize('tenor,mult,result',
                         [
                             ('5d', 5, '25d'),
                             ('-7w', 3, '-21W'),
                             ('12M', 2, '24M'),
                             ('-3Y', -2, '6Y'),
                         ])
def test_tenor_multiply(tenor, mult, result):
    t1 = Tenor(tenor)
    t2 = t1 * mult
    t3 = mult*t1
    res = Tenor(result)
    assert t2 == res
    assert t3 == res


@pytest.mark.parametrize('tenor1,tenor2,result',
                         [
                             ('5d', '10d', '15d'),
                             ('7w', '3d', '52D'),
                             ('1M', '1D', '29D'),
                             ('1Y', '1W', '49w'),
                             ('12M', '1Y', '24M'),
                             ('-3Y', '-2Y', '-5Y'),
                         ])
def test_tenor_add(tenor1, tenor2, result):
    t1 = Tenor(tenor1)
    t2 = Tenor(tenor2)
    res = Tenor(result)
    assert t1 + t2 == res


# if __name__ == '__main__':
#     # test_tenor_create(None, 0, TenorUnit.NONE)
#     # test_tenor_multiply('5d', 5, '25d')
#     test_tenor_add('5d', '10d', '15d')
