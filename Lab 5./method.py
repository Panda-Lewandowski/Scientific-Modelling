def fixed_point_iteration(A, B, C, D, K0, M0, P0, KN, MN, PN):  # Tridiagonal matrix algorithm
    # Initial values
    xi = [None, - M0 / K0]
    eta = [None, P0 / K0]

    # Straight running
    for i in range(1, len(A)):
        x = C[i] / (B[i] - A[i] * xi[i])
        e = (D[i] + A[i] * eta[i]) / (B[i] - A[i] * xi[i])

        xi.append(x)
        eta.append(e)

    # print(xi)
    # print(eta)

    # Reverse running
    y = [(PN - MN * eta[-1]) / (KN + MN * xi[-1])]

    for i in range(len(A) - 1, -1, -1):
        y_i = xi[i + 1] * y[0] + eta[i + 1]

        y.insert(0, y_i)

    return y


def two_d_interpolation(table: dict, x):
    xs = list(table.keys())
    xs.sort()
    if x in xs:
        return x
    elif x > xs[-1]:
        return xs[-1]
    elif x < xs[0]:
        return xs[0]
    else:
        for i in range(len(xs) - 1):
            if xs[i] < x < xs[i + 1]:
                return table[xs[i]] + (table[xs[i + 1]] - table[xs[i]]) / (xs[i + 1] - xs[i]) * (x - xs[i])







