from data_getter import Data
from method import thomas_algorithm
from numpy import arange
import matplotlib.pyplot as plt


def left_boundary_conditions():
    k1 = Data.k(Data.x0 + Data.h)
    alpha_half = Data.alpha(Data.x0 + Data.h / 2)

    X_half = (2 * Data.k0 * k1) / (Data.k0 + k1)
    p_half = (2 * alpha_half) / Data.R
    f_half = (2 * alpha_half) / Data.R * Data.Tenv

    p0 = Data.p(Data.x0)
    f0 = Data.f(Data.x0)

    K0 = X_half + Data.h ** 2 / 8 * p_half + Data.h ** 2 / 4 * p0
    M0 = Data.h ** 2 / 8 * p_half - X_half
    P0 = Data.h * Data.F0 + Data.h ** 2 / 4 * (f0 + f_half)

    return K0, M0, P0


def right__boundary_conditions():
    X_half = Data.Xn_minus_half(Data.l)
    pN = Data.p(Data.l)
    pN1 = Data.p(Data.l - Data.h)
    fN = Data.f(Data.l)
    fN1 = (2 * Data.alpha(Data.l - Data.h / 2)) / Data.R * Data.Tenv

    KN = - X_half - Data.h ** 2 / 16 * pN1 - 5 * Data.h ** 2 / 16 * pN
        # Data.a_alpha / Data.Tenv + Data.a_alpha / Data.b_alpha

    MN = X_half - Data.h ** 2 / 16 * pN1 - Data.h ** 2 / 16 * pN
    PN = Data.h * Data.a_alpha - Data.h ** 2 / 4 * (fN - fN1)

    return KN, MN, PN


def calc_coefficients():
    A = []
    B = []
    C = []
    D = []

    for i in arange(Data.x0, Data.l, Data.h):
        An = Data.Xn_minus_half(i) / Data.h
        Cn = Data.Xn_plus_half(i) / Data.h
        Bn = An + Cn + Data.p(i) * Data.h
        Dn = Data.f(i) * Data.h

        A.append(An)
        B.append(Bn)
        C.append(Cn)
        D.append(Dn)

    return A, B, C, D


if __name__ == "__main__":
    a, b, c, d = calc_coefficients()
    # print(a)
    # print(b)
    # print(c)
    # print(d)

    k0, m0, p0 = left_boundary_conditions()
    # print(k0)
    # print(m0)
    # print(p0)

    kN, mN, pN = right__boundary_conditions()
    # print(kN)
    # print(mN)
    # print(pN)

    T = thomas_algorithm(a, b, c, d, k0, m0, p0, kN, mN, pN)
    x = arange(Data.x0, Data.l, Data.h)

    plt.title('Heating the rod')
    plt.grid(True)
    plt.plot(x, T, 'r', linewidth=0.5)
    plt.xlabel("Length (cm)")
    plt.ylabel("Temperature (K)")

    plt.savefig("plot.png")

    plt.show()


