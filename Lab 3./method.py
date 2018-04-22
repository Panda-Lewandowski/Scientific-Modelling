def thomas_algorithm(A, B, C, D, K0, M0, P0, KN, MN, PN):  # Tridiagonal matrix algorithm
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

    for i in range(len(A) - 2, -1, -1):
        y_i = xi[i + 1] *  y[0] + eta[i + 1]

        y.insert(0, y_i)

    return y


