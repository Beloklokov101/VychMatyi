import numpy as np
from sympy import *


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
    print(A, d)
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


#   МЕТОД ПРОГОНКИ

v, y, t = symbols("v, y, t")
#a -- u_n-1, b -- u_n, c -- u_n+1
    a_f = lambdify([t, v ,y], 1 + Fdv * h / 2, "numpy")
    b_f = lambdify([t, v ,y], -2 - Fdy * h**2, "numpy")
    c_f = lambdify([t, v ,y], 1 - Fdv * h / 2, "numpy")
    #print(-2 - Fdy * h**2)

    a_n, b_n, c_n, d_n = np.zeros(N + 1), np.zeros(N + 1), np.zeros(N + 1), np.zeros(N + 1)
    a_n[0], c_n[-1] = 0, 0

    b_n[0], c_n[0], d_n[0] = 1, 0, 0
    # DANGER!!!
    a_n[1], b_n[1], c_n[1], d_n[1] = 1, -1, 0, - h #- fdv(T[0], V[0], Y[0]) * h**2 / 2

    for i in range(2, N):
        a_n[i] = a_f(T[i], V[i], Y[i])
        b_n[i] = b_f(T[i], V[i], Y[i])
        c_n[i] = c_f(T[i], V[i], Y[i])
        d_n[i] = 0
    
    a_n[-1], b_n[-1], d_n[-1] = a_f(T[i], V[i], Y[i]), b_f(T[i], V[i], Y[i]), 0

    u_n = solveTriagonalSlae(a_n, b_n, c_n, d_n)
    print(u_n)
