from polynomial import Polynomial
from prettytable import PrettyTable
import numpy as np


def test_func(x, u):
    return x ** 2 + u ** 2


def picard(ns, x, func, pol=False):
    """ Cauchy problem:
        u'(x) = func(x, u) = x^2 + u^2
        u(0) = 0

    y(0) = u(0) """

    y = Polynomial(0, [0.0])
    t = Polynomial(1, [1.0])

    us = []

    for i in range(max(ns)):
        under_int = func(t, y)
        y = under_int.integral_variable_up(0.0)

        if i+1 in ns:
            if pol:
                print("Полином {0} итерации: ".format(i+1), y)
            us.append(round(y.get(x), 7))
    return us


def broken_evident(x, func):
    y = 0
    t = 0

    while t <= x:
        f = func(t, y) * h
        y += f
        t += h

    return y


def runge_kutta_second(x, alpha, h, func):

   # assert(alpha == 0.5 or alpha == 1, "Alpha должна быть равно 1 или 0.5")

    y = 0
    t = 0

    while t <= x:
        f = h * ((1 - alpha) * func(t, y) +
                 alpha * func(t + h / (2 * alpha),
                              y + (h / (2 * alpha)) * func(t, y)))

        y += f
        t += h

    return y


def runge_kutta_fourth(x, h, func):
    y = 0
    t = 0

    while t <= x:
        k1 = func(t, y)
        k2 = func(t + h / 2, y + h / 2 * k1)
        k3 = func(t + h / 2, y + h / 2 * k2)
        k4 = func(t + h, y + h * k3)

        f = h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

        y += f
        t += h

    return y


def calc(ns, a, b, h, alpha, func):
    table_picard = PrettyTable()
    title = ["X"]
    for it in ns:
        title.append("Пикар: {0}-я итерация".format(it))
    title.append("Метод ломанной (явный)")
    title.append("Метод Рунге-Кутты 2-ого порядка (неявный)")
    title.append("Метод Рунге-Кутты 4-ого порядка")

    table_picard.field_names = title

    xs = np.arange(a, b, h)

    for x in xs:
        u = picard(ns, x, func)
        u.insert(0, round(x, 3))
        bevi = broken_evident(x, func)
        u.append(bevi)
        rks = runge_kutta_second(x, alpha, h, func)
        u.append(rks)
        rkf = runge_kutta_fourth(x, h, func)
        u.append(rkf)
        table_picard.add_row(u)

    u = picard(ns, b, func, pol=True)
    u.insert(0, round(b, 3))
    bevi = broken_evident(b, func)
    u.append(bevi)
    rk = runge_kutta_second(b, alpha, h, func)
    u.append(rk)
    rkf = runge_kutta_fourth(b, h, func)
    u.append(rkf)
    table_picard.add_row(u)

    print(table_picard)


if __name__ == "__main__":
    x1 = int(input("Введите X: от "))
    x2 = int(input(" " * 11 + "до "))
    h = float(input(" " * 6 + "с шагом "))
    nl = input("Введите через пробел интересующие итерации для метода  Пикара: ").split(" ")
    a = float(input("Введите α для метода Рунге-Кутты второго порядка: "))
    nl = list(map(lambda x: int(x), nl))
    calc(nl, x1, x2, h, a, test_func)
