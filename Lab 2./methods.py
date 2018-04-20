def Runge_Kutta_IV_3var_1step(f, phi, x0, y0, z0, h):
    k1 = f(x0, y0, z0)
    p1 = phi(x0, y0, z0)

    k2 = f(x0 + h/2, y0 + k1 * h / 2, z0 + p1 * h / 2)
    p2 = phi(x0 + h/2, y0 + k1 * h / 2, z0 + p1 * h / 2)

    k3 = f(x0 + h / 2, y0 + k2 * h / 2, z0 + p2 * h / 2)
    p3 = phi(x0 + h / 2, y0 + k2 * h / 2, z0 + p2 * h / 2)

    k4 = f(x0 + h / 2, y0 + k3 * h, z0 + p3 * h)
    p4 = phi(x0 + h / 2, y0 + k3 * h, z0 + p3 * h)

    yn = y0 + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    zn = z0 + h / 6 * (p1 + 2 * p2 + 2 * p3 + p4)

    return yn, zn


def simpson(func, left, right, n, *args):
    h = (right - left) / n
    ans = h / 3
    even = 0.0
    odd = 0.0

    for i in range(1, n):
        if i % 2 == 0:
            even += func(left + h * i, *args)
        else:
            odd += func(left + h * i, *args)

    ans *= (2 * even + 4 * odd + func(left, *args) + func(right, *args))
    return ans


def dichotomy(func, left, right, eps, *args):
    while abs(right - left) > eps:
        c = (left + right) / 2.0

        if func(right, *args) * func(c, *args) <= 0:
            left = c
        else:
            right = c
    return (right + left) / 2.0


def newton_interpolation():
    pass

