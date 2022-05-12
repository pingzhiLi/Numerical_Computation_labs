import numpy as np


def y_drv(x, y, z):
    return 0.09 * y * (1 - y / 20) - 0.45 * y * z


def z_drv(x, y, z):
    return 0.06 * z * (1 - z / 15) - 0.001 * y * z


def euler_ode(N, a, b):
    h = (b - a) / N
    x = np.linspace(a, b, N + 1, dtype=np.int16)
    y = np.zeros(N + 1)
    z = np.zeros(N + 1)
    y[0] = 1.6
    z[0] = 1.2
    for i in range(N):
        _y = y[i] + h * y_drv(x[i], y[i], z[i])
        _z = z[i] + h * z_drv(x[i], y[i], z[i])
        y[i + 1] = y[i] + h / 2 * \
            (y_drv(x[i], y[i], z[i]) + y_drv(x[i + 1], _y, _z))
        z[i + 1] = z[i] + h / 2 * \
            (z_drv(x[i], y[i], z[i]) + z_drv(x[i + 1], _y, _z))
    return x, y, z


t, u, v = euler_ode(3, 1, 4)
for i in range(len(t)):
    print("t =", t[i], ", u =", round(u[i], 6), ", v =", round(v[i], 6))
