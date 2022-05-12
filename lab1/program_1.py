import numpy as np
import matplotlib.pyplot as plt

actual_points = {
    1920: 105711,
    1930: 123203,
    1940: 131669,
    1950: 150697,
    1960: 179323,
    1970: 203212
}

additional_point = {
    1910: 91772
}


def lagrange_interp(x, inter_points):
    res = 0
    for point in inter_points:
        tmp = 1
        for _point in inter_points:
            if _point is not point:
                tmp = tmp * (x - _point) / (point - _point)
        res = res + tmp * inter_points[point]
    return res


print('L(1910)=', lagrange_interp(1910, actual_points))
print('L(1965)=', lagrange_interp(1965, actual_points))
print('L(2002)=', lagrange_interp(2002, actual_points))

actual_points.popitem()
print('L\'(1965)=', lagrange_interp(1965, {**actual_points, **additional_point}))
print('L\'(2002)=', lagrange_interp(2002, {**actual_points, **additional_point}))
