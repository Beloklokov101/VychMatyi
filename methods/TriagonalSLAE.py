import numpy as np
import matplotlib.pyplot as plt


# Dasha's solver
def solveTriagonalSlae(a, b, c, d):
    n = len(b)
    h = d

    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i][j] = b[i]
            if i == j + 1:
                A[i][j] = a[i]
            if i + 1 == j:
                A[i][j] = c[j - 1]
    #print(A, d)
    #Dasha's part
    P = np.zeros((n,1))
    Q = np.zeros((n,1))
    x = np.zeros(n)
    
    P[0] = -A[0, 1] / A[0, 0]
    Q[0] = h[0] / A[0, 0]
    
    #a, b, c, d initialized differently
    for i in range(1,n):
        a, b = A[i,i-1], A[i,i]
        if i != n-1:
            c = A[i,i+1]
            P[i] = c/(-b - a*P[i-1])
        Q[i] = (-h[i]+a*Q[i-1])/(-b - a*P[i-1])
    x[n-1] = Q[n-1]
    
    for i in range(n-2,-1,-1):
        x[i] = P[i]*x[i+1]+Q[i]
    
    return x


# Mine solver
def TriagonalSlae_solver(a, b, c, d, N, M = 1, transpose = False):
    # np.shape(a or b or c or d) == (N, M)

    # a[0] * u[1] + b[0] * u[0] = d[0]
    # a[i] * u[i + 1] + b[i] * u[i] + c[i] * u[i - 1] = d[i], i = range(1, N - 1)
    # b[N - 1] * u[N - 1] + c[N - 1] * u[N - 2] = d[N - 1]

    # (!) There should be c[0] = 0, a[N - 1] = 0
    # This solver can solve M SLAE's simultaneously

    # u[0] = f(m) <==> a[0] = 0, b[0] = 1, d[0] = f(m)
    # u[N - 1] = g(m) <==> c[N - 1] = 0, b[N - 1] = 1, d[N - 1] = g(m)

    dimens = (N, M)
    shape_a = np.shape(a)
    shape_b = np.shape(b)
    shape_c = np.shape(c)
    shape_d = np.shape(d)

    if M != 1:
        if (dimens != shape_a or dimens != shape_b or dimens != shape_c or dimens != shape_d):
            print("\nERROR IN TriagonalSlae_solver")
            print(f"shape_a = {shape_a}, (N, M) = {dimens}\n")
            return False
    else:
        if (dimens[0] != shape_a[0] or dimens[0] != shape_b[0] or dimens[0] != shape_c[0] or dimens[0] != shape_d[0]):
            print("\nERROR IN TriagonalSlae_solver")
            print(f"shape_a = {shape_a}, (N, M) = {dimens}\n")
            return False
    
    if transpose:
        (N, M) = (M, N)
        a, b, c, d = a.T, b.T, c.T, d.T

    # u[i] = alpha[i] * u[i + 1] + beta[i]
    alpha, beta = np.zeros((N - 1, M)), np.zeros((N - 1, M))
    alpha[0] = - a[0] / b[0]
    beta[0] = d[0] / b[0]
    for it in range(1, N - 1):
        divider = b[it] + c[it] * alpha[it - 1]
        alpha[it] = - a[it] / divider
        beta[it] = (d[it] - c[it] * beta[it - 1]) / divider
    
    if M != 1:
        u = np.zeros((N, M))
    else:
        u = np.zeros(N)
    
    u[N - 1] = (d[N - 1] - c[N - 1] * beta[N - 2]) / (b[N - 1] + c[N - 1] * alpha[N - 2])
    for it in range(N - 2, -1, -1):
        u[it] = alpha[it] * u[it + 1] + beta[it]

    if transpose:
        return u.T
    
    return u



# Test SLAE solvers
N = 50
M = 2
a, b, c, d = np.zeros((N, M)), np.zeros((N, M)), np.zeros((N, M)), np.zeros((N, M))
for it in range(N):
    a[it] = np.array([it, it])
    b[it] = np.array([1, 1])
    c[it] = np.array([it ** 2, it ** 2])
    d[it] = np.array([1, 2])

u_mine = TriagonalSlae_solver(a.T, b.T, c.T, d.T, M, N, transpose=True)
u_mine = u_mine.T

u_Dasha_first = solveTriagonalSlae(c[:, 0], b[:, 0], a[:, 0], d[:, 0])
u_Dasha_second = solveTriagonalSlae(c[:, 1], b[:, 1], a[:, 1], d[:, 1])

x_grid = np.linspace(0, 1, N)

plt.figure()
plt.plot(x_grid, u_mine[:, 0] - u_Dasha_first)

plt.figure()
plt.plot(x_grid, u_mine[:, 1] - u_Dasha_second)

plt.show()