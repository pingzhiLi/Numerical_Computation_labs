import numpy as np


def y_drv(x, y):
    return np.sin(np.pi * x) * y


def rouge_kutta_2(n, a, b):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = np.zeros(n + 1)
    y[0] = 1
    for i in range(n):
        k1 = y_drv(x[i], y[i])
        k2 = y_drv(x[i] + h, y[i] + h * k1)
        y[i + 1] = y[i] + h * (k1 + k2) / 2
    return x, y


x, y = rouge_kutta_2(10, 0, 1)
for i, j in zip(x, y):
    print("x =", round(i, 2), ", y =", round(j, 5))
