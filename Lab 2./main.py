from data_getters import Data
from methods import Runge_Kutta_IV_3var_1step, simpson, dichotomy
from math import pi
from scipy import arange


def get_To_m(I):
    table_I_T0 = Data.get_table_I_T0()
    len_I_T0 = len(table_I_T0)

    In, Ik = 0, 1

    if I > 1200:
        In, Ik = len_I_T0 - 2,  len_I_T0 - 1

    for i in range(len_I_T0 - 1):
        if table_I_T0[i][0] < I < table_I_T0[i + 1][0]:
            In = i
            Ik = i + 1

    X0 = table_I_T0[In][0]
    X1 = table_I_T0[Ik][0]
    Y0 = table_I_T0[In][1]
    Y1 = table_I_T0[Ik][1]
    To = Y0 + (Y1 - Y0) * (I - X0) / (X1 - X0)
    Y0 = table_I_T0[In][2]
    Y1 = table_I_T0[Ik][2]
    m = Y0 + (Y1 - Y0) * (I - X0) / (X1 - X0)
    #print("T: ", To, "\nm:", m, "\n")
    return To, m


def get_T(I, r):
    To, m = get_To_m(I)
    return (Data.T_w - To) * (r / Data.Rad) ** m + To


def get_sigma(I, r, p):
    sigma_table = Data.get_table_sigma()
    len_sigma = len(sigma_table)

    T = get_T(I, r)

    Tn, Tk = 0, 1
    if T > 20000:
        Tn, Tk = len_sigma - 2, len_sigma - 1

    for i in range(len_sigma - 1):
        if sigma_table[i][0] < T < sigma_table[i + 1][0]:
            Tn = i
            Tk = i + 1

    pn, pk = 0, 1
    len_p = len(Data.p)
    if p > 25:
        pn, pk = len_p - 2, len_p - 1

    for i in range(len_p - 1):
        if Data.p[i] < p < Data.p[i + 1]:
            pn = i
            pk = i + 1

    X0 = sigma_table[Tn][0]
    X1 = sigma_table[Tk][0]
    Y0 = sigma_table[Tn][pn]
    Y1 = sigma_table[Tk][pn]
    fy1 = Y0 + (Y1 - Y0) * (T - X0) / (X1 - X0)

    Y0 = sigma_table[Tn][pk]
    Y1 = sigma_table[Tk][pk]
    fy2 = Y0 + (Y1 - Y0) * (T - X0) / (X1 - X0)

    Y0 = Data.p[pn]
    Y1 = Data.p[pk]
    sigma = fy1 + (fy2 - fy1) * (p - Y0) / (Y1 - Y0)
    return sigma


def get_nt(I, r, p):
    table_nt = Data.get_table_nt()
    len_nt = len(table_nt)

    T = get_T(I, r)

    Tn, Tk = 0, 1
    if T > 15000:
        Tn, Tk = len_nt - 2, len_nt - 1

    for i in range(len(table_nt) - 1):
        if table_nt[i][0] < T < table_nt[i + 1][0]:
            Tn = i
            Tk = i + 1

    pn, pk = 0, 1

    len_p = len(Data.p)
    if p > 25:
        pn, pk = len_p - 2, len_p - 1

    for i in range(len_p - 1):
        if Data.p[i] < p < Data.p[i + 1]:
            pn = i
            pk = i + 1

    X0 = table_nt[Tn][0]
    X1 = table_nt[Tk][0]
    Y0 = table_nt[Tn][pn]
    Y1 = table_nt[Tk][pn]
    fy1 = Y0 + (Y1 - Y0) * (T - X0) / (X1 - X0)

    Y0 = table_nt[Tn][pk]
    Y1 = table_nt[Tk][pk]
    fy2 = Y0 + (Y1 - Y0) * (T - X0) / (X1 - X0)

    Y0 = Data.p[pn]
    Y1 = Data.p[pk]
    nt = fy1 + (fy2 - fy1) * (p - Y0) / (Y1 - Y0)
    return nt


def nt_integrand(r, I, p):
    nt = get_nt(I, r, p)
    return nt*r


def p_func(p, I):
    return 2/(Data.Rad * Data.Rad) * simpson(nt_integrand, 0, Data.Rad, 10, *[I, p]) - Data.p_0 * 7242 / Data.T_0


def sigma_integrand(r, I):
    p = dichotomy(p_func, min(Data.p), max(Data.p), Data.eps, *[I])
    sigma = get_sigma(I, r, p)
    # print("R:", r, "\n")
    return sigma * r


def get_Rp(t, I):
    return Data.l_e / (2 * pi * simpson(sigma_integrand, 0, Data.Rad, 10, *[I]))


def dI_dt(t, Uc, I):
    Rp = get_Rp(t, I)
    return (Uc - (Data.R_k + Rp) * I) / Data.L_k


def dUc_dt(t, Uc, I):
    return -1 / Data.C_k * I


if __name__ == '__main__':
    ts = list(arange(Data.t_0, (Data.t_n + Data.tau), Data.tau))
    I, Uc, Rp = [Data.I_0], [Data.U_0], []

    Ii = Data.I_0
    Ui = Data.U_0
    Rpi = get_Rp(Data.t_0, Ii)

    for t in ts[1:]:

        print(t, Ii, Ui, Rpi)
        I_next, Uc_next = Runge_Kutta_IV_3var_1step(dI_dt, dUc_dt, t, Ui, Ii, Data.tau)
        Rpi = get_Rp(t, Ii)

        I.append(I_next)
        Uc.append(Uc_next)
        Rp.append(Rpi)

        Ii = I_next
        Ui = Uc_next



