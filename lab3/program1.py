import numpy as np


def f_f(x, y):
    return x**2 + y**2 - 1


def f_g(x, y):
    return x**3 - y


def Jaccobi(x, y):
    return np.array([[2*x, 2*y], [3*x**2, -1]])


def newton_iter(x0, y0, eps):
    x = x0
    y = y0
    k = 1
    while True:
        J = Jaccobi(x, y)
        J_inv = np.linalg.inv(J)
        f = np.array([f_f(x, y), f_g(x, y)])
        delta = -np.dot(J_inv, f)
        x += delta[0]
        y += delta[1]
        print('Step', k, ':', 'x =', x, 'y =', y)
        k += 1
        if np.linalg.norm(delta) < eps:
            break
    return x, y


x, y = newton_iter(0.8, 0.6, 1e-5)
