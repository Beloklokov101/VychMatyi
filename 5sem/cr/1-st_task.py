import math
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

a, b = -3, 5
#a, b = 0, 1
N = 3
X = [(a + b) / 2 + (b - a) * math.cos((2 * i + 1) * math.pi / (2 * N)) / 2 for i in range(N)]
X.reverse()
print(X)
x = symbols("x")
F = cos(x)
#F = atan(x)
f = lambdify(x, F, modules="numpy")
P_20 = f(X[0])
P_21 = (f(X[1]) - f(X[0])) * (x - X[0]) / (X[1] - X[0])
f_12 = (f(X[2]) - f(X[1])) / (X[2] - X[1])
f_10 = (f(X[1]) - f(X[0])) / (X[1] - X[0])
P_22 = (x - X[0]) * (x - X[1]) * (f_12 - f_10) / (X[2] - X[0])
P_2 = P_20 + P_21 + P_22
P = simplify(P_2)
print(P)
print("Погрешность = 1/96 = ", 1 / 96)

p = lambdify(x, P)
x_x = np.linspace(-3, 5, 100)
p_x = p(x_x)
#print(x_x, p_x)
plt.plot(x_x, p_x)
plt.plot(x_x, f(x_x))
plt.show()