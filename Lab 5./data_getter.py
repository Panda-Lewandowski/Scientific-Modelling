from numpy import arange
from method import two_d_interpolation


class Data:
    x0 = 0
    l = 10           # Length of the rod (cm)
    R = 0.5          # Radius of the rod (cm)
    Tenv = 300       # Ambient temperature (K)
    F0 = 100         # Heat flux density (W / (cm^2 * K))
    alpha0 = 1e-2    # Heat transfer coefficient at the beginning of the rod (W / (cm^2 * K))
    alphaN = 9e-3  # Heat transfer coefficient at the end of the rod (W / (cm^2 * K))
    h = 1e-2
    tau = 1
    eps = 1e-3
    rounding = 4
    b_alpha = (alphaN * l) / (alphaN - alpha0)
    a_alpha = - alpha0 * b_alpha

    T_curr = {round(x, 4): 300 for x in arange(x0, l+h, h)}
    t_curr = 0

    T_past = None

    @staticmethod
    def get_table_of_heat_cap():
        return {300:	2.544,
                400:	2.615,
                500:	2.656,
                600:	2.689,
                700:	2.717,
                800:	2.748,
                900:	2.783,
                1000:	2.817,
                1200:	2.893,
                1400:	2.977,
                1600:	3.070,
                1800:	3.124,
                2000:	3.270,
                2200:	3.381,
                2400:	3.502,
                2600:	3.640,
                2800:	3.792,
                3000:	3.968}

    @staticmethod
    def get_table_of_thermal_cond():
        return {300:    0.163,
                400:	0.156,
                500:	0.146,
                600:	0.137,
                700:	0.130,
                800:	0.124,
                900:	0.120,
                1000:	0.117,
                1200:	0.114,
                1400:	0.111,
                1600:	0.110,
                1800:	0.109,
                2000:   0.108,
                2200:	0.107,
                2400:	0.107,
                2600:	0.107,
                2800:	0.109,
                3000:	0.108}

    @staticmethod
    def alpha(x):
        return Data.a_alpha / (x - Data.b_alpha)

    @staticmethod
    def k(x):
        if x in Data.T_curr.keys():
            T = Data.T_curr[x]
        else:
            T = two_d_interpolation(Data.T_curr, x)
        table = Data.get_table_of_thermal_cond()
        if T in table.keys():
            return table[T]
        else:
            # interpolation
            return two_d_interpolation(table, T)

    @staticmethod
    def C(x):
        if x in Data.T_curr.keys():
            T = Data.T_curr[x]
        else:
            T = two_d_interpolation(Data.T_curr, x)
        table = Data.get_table_of_heat_cap()
        if T in table.keys():
            return table[T]
        else:
            # interpolation
            return two_d_interpolation(table, T)

    @staticmethod
    def Xn_minus_half(x):
        return (Data.k(x - Data.h) + Data.k(x)) / 2

    @staticmethod
    def Xn_plus_half(x):
        return (Data.k(x) + Data.k(x + Data.h)) / 2

    @staticmethod
    def p(x):
        return 2 * Data.alpha(x) / Data.R

    @staticmethod
    def f(x):
        return 2 * Data.alpha(x) / Data.R * Data.Tenv


if __name__ == "__main__":
    print(Data.C(0.49585855))


