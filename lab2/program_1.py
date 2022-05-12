from math import log


def sn_calc(a: float, b: float, m: int, f):
    n = 2 * m
    h = (b - a) / n
    sn = f(a) + f(b)
    for i in range(m):
        sn = sn + 4 * f(a + (2 * i + 1) * h)
    for i in range(1, m):
        sn = sn + 2 * f(a + (2 * i) * h)
    sn = sn * h / 3
    return sn


def s2n_from_sn(a: float, b: float, sn: float, n: int, f):
    h = (b - a) / n
    Hn = 0
    H2n = 0
    for i in range(n):
        Hn = Hn + f(a + h * (i + 1 / 2))
    for i in range(2 * n):
        H2n = H2n + f(a + (h / 2) * (i + 1 / 2))
    Hn = Hn * h
    H2n = H2n * h / 2
    s2n = 0.5 * sn + 1 / 6 * (4 * H2n - Hn)
    return s2n


def auto_precision_simpson(a: float, b: float, eps: float, f, initial_m: int):
    m = initial_m
    s2 = sn_calc(a=a, b=b, m=m, f=f)
    s1 = s2 + 1
    iter_times = 0
    while abs(s1 - s2) > eps:
        iter_times += 1
        s1 = s2
        s2 = s2n_from_sn(a=a, b=b, sn=s1, n=2 * m, f=f)
        m = 2 * m
    return s2, iter_times


start = 1
end = 2
e = 10 ** (-4)
for i in range(1, 5):
    integral, i_time = auto_precision_simpson(a=start, b=end, eps=e, f=log, initial_m=i)
    print("Initial n:", 2 * i)
    print("Simpson:", round(sn_calc(a=start, b=end, m=i, f=log), 8))
    print("Auto precision-control Simpson:", round(integral, 8), ".")
    print(i_time, "iterations totally")
