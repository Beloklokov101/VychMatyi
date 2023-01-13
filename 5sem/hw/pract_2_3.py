import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import math as mt


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


#   МЕТОД ПРОГОНКИ

x = symbols("x")

a, b = 0, mt.pi
y_f, y_l = 0, mt.pi**2

N = 1000
h = (b - a) / N
T = np.linspace(a, b, N + 1)

G = x**2 - 3
H = (x**2 - 3)*cos(x)
F = 2 - 6*x + 2*x**3 + (x**2 - 3)*exp(x)*sin(x)*(1 + cos(x)) + cos(x)*(exp(x) + (x**2 - 1) + x**4 - 3*x**2)

#a -- u_n-1, b -- u_n, c -- u_n+1
a_f = lambdify(x, 1 - G * h / 2, "numpy")
b_f = lambdify(x, -2 + H * h**2, "numpy")
c_f = lambdify(x, 1 + G * h / 2, "numpy")
f_f = lambdify(x, F * h**2, "numpy")
#print(-2 - Fdy * h**2)

a_n, b_n, c_n, d_n = np.zeros(N + 1), np.zeros(N + 1), np.zeros(N + 1), np.zeros(N + 1)
a_n[0], c_n[-1] = 0, 0

b_n[0], c_n[0], d_n[0] = 1, 0, y_f
a_n[-1], b_n[-1], d_n[-1] = 0, 1, y_l

for i in range(1, N):
    a_n[i] = a_f(T[i])
    b_n[i] = b_f(T[i])
    c_n[i] = c_f(T[i])
    d_n[i] = f_f(T[i])

Y = solveTriagonalSlae(a_n, b_n, c_n, d_n)

out_dots = np.array([0.5, 1, 1.5, 2, 2.5, 3])
j = 0
for i in range(1, N + 1):
    d = np.sum((out_dots - T[i - 1]) < 0) - np.sum((out_dots - T[i]) < 0)
    if (d != 0):
        num = i - 1 + ((out_dots[j] - T[i - 1]) > (T[i] - out_dots[j]))
        print(f"T[{num}] = {T[num]}")
        print(f"Y[{num}] = {Y[num]}\n")
        j = j + 1

plt.plot(T, Y)
plt.show()