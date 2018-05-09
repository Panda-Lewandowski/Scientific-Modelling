class Data:
    x0 = 0
    l = 10           # Length of the rod (cm)
    R = 0.5          # Radius of the rod (cm)
    Tenv = 300       # Ambient temperature (K)
    F0 = 100        # Heat flux density (W / (cm^2 * K))
    k0 = 0.2         # Coefficient of thermal conductivity at the beginning of the rod (W / (cm * K))
    kN = 0.5         # Coefficient of thermal conductivity at the end of the rod (W / (cm * K))
    alpha0 = 1e-2    # Heat transfer coefficient at the beginning of the rod (W / (cm^2 * K))
    alphaN = 9e-3  # Heat transfer coefficient at the end of the rod (W / (cm^2 * K))
    h = 1e-2
    bk = (kN * l) / (kN - k0)
    ak = - k0 * bk
    b_alpha = (alphaN * l) / (alphaN - alpha0)
    a_alpha = - alpha0 * b_alpha


    @staticmethod
    def k(x):
        return Data.ak / (x - Data.bk)

    @staticmethod
    def alpha(x):
        return Data.a_alpha / (x - Data.b_alpha)

    @staticmethod
    def Xn_plus_half(x):
        return (2 * Data.k(x) * Data.k(x + Data.h)) / \
               (Data.k(x) + Data.k(x + Data.h))

    @staticmethod
    def Xn_minus_half(x):
        return (2 * Data.k(x) * Data.k(x - Data.h)) / \
               (Data.k(x) + Data.k(x - Data.h))

    @staticmethod
    def p(x):
        return 2 * Data.alpha(x) / Data.R

    @staticmethod
    def f(x):
        return 2 * Data.alpha(x) / Data.R * Data.Tenv

