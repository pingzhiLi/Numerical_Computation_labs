import numpy as np
import matplotlib.pyplot as plt

points = {
    1920: 105711,
    1930: 123203,
    1940: 131669,
    1950: 150697,
    1960: 179323,
    1970: 203212
}


def newton_interp(u, inter_points):
    x, g = [], []
    n = len(inter_points) - 1
    for _x in inter_points:
        x.append(_x)
        g.append(inter_points[_x])
    for k in range(1, n + 1):
        for j in range(n, k - 1, -1):
            g[j] = (g[j] - g[j - 1]) / (x[j] - x[j - k])
    t = 1
    newton = g[0]
    for k in range(1, n + 1):
        t = t * (u - x[k - 1])
        newton = newton + t * g[k]
    return newton


print('N(1965)=', newton_interp(1965, points))
print('N(2012)=', newton_interp(2012, points))
