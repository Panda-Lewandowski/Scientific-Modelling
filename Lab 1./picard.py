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

    y = Polynomial(0, [0])
    t = Polynomial(1, [1])

    us = []

    for i in range(max(ns)):
        under_int = func(t, y)
        y = under_int.integral_variable_up(0)

        if i+1 in ns:
            if pol:
                print("Полином {0} итерации: ".format(i+1), y)
            us.append(round(y.get(x), 7))
    return us


def calc_picard(ns, a, b, h, func):
    table = PrettyTable()
    title = ["X"]
    for it in ns:
        title.append("Пикар: {0}-я итерация".format(it))

    table.field_names = title

    xs = np.arange(a, b, h)

    for x in xs:
        u = picard(ns, x, func)
        u.insert(0, round(x, 3))
        table.add_row(u)

    u = picard(ns, b, func, pol=True)
    u.insert(0, round(b, 3))
    table.add_row(u)

    print(table)


if __name__ == "__main__":
    print(">>>>> МЕТОД ПИКАРА <<<<<")
    x1 = int(input("Введите X: от "))
    x2 = int(input(" " * 11 + "до "))
    h = float(input(" " * 6 + "с шагом "))
    nl = input("Введите через пробел интересующие итерации: ").split(" ")
    nl = list(map(lambda x: int(x), nl))
    calc_picard(nl, x1, x2, h, test_func)
