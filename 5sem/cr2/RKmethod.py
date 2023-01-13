import numpy as np
from sympy import *

z = symbols("z")
x, y = symbols('x y', real=True)
"""A = np.array([[1/4, 0], [1/2, 1/4]])
c = np.array([1/4, 3/4])
b = np.array([1/2, 1/2])
"""
A = np.array([[0, 0, 0, 0], [1/2, 0, 0, 0], [0, 1/2, 0, 0], [0, 0, 1, 0]])
c = np.array([0, 1/2, 1/2, 1])
b = np.array([1/6, 1/3, 1/3, 1/6])

n = len(b)

p = 0
print(np.sum(b))
if (1 == np.sum(b)):
    p = 1

if (p == 1):
    sum = 0
    for i in range(n):
        sum = sum + b[i] * np.sum(A, axis=1)[i]
    if (2 * sum == 1):
        p = 2

if (p == 2):
    sum = 0
    for i in range(n):
        for j in range(n):
            sum = sum + b[i] * A[i][j] * np.sum(A, axis=1)[i]
    if (sum * 3 == 1):
        p = 3

print(f"RK approx level p = {p}")

A_s = z * Matrix(A)
Z_s = z * ones(n, 1) * Matrix(b).T
P = (eye(n) - A_s + Z_s).det()
Q = (eye(n) - A_s).det()
print(P, "\n", Q)

#P = 1 + z/3
#Q = 1 - 2*z/3 + z**2 / 6

y = Symbol('y', real=True)
#E = (Abs(Q.subs(z, I*y)))**2 - (Abs(P.subs(z, I*y)))**2
E = (re(Q.subs(z, x + I*y))**2 + im(Q.subs(z, x + I*y))**2) - (re(P.subs(z, x + I*y))**2 + im(P.subs(z, x + I*y))**2)
#E = (P.subs(z, x + I*y)) *  (P.subs(z, x - I*y)) - (Q.subs(z, x + I*y)) *  (Q.subs(z, x - I*y))

#print(Abs(P.subs(z, I*y)))
print(f"E(y) = {E}")