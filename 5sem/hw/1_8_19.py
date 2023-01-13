import math
from sympy import *
import numpy as np


x = symbols("x")
F = cos(x)

e = 1e-3
n = 0
a, b = 0, 1

max_end = max(abs(a), abs(b))
F_diff = diff(F, x, 1)
X = np.linspace(a, b, 1000)
f_diff = lambdify(x, F_diff, modules="numpy")
f_diff_max = np.max(abs(f_diff(X)))
R_n = max_end * f_diff_max
n = 0
print(R_n)

while(R_n > e):
    n = n + 1
    F_diff = diff(F, x, n + 1)
    f_diff = lambdify(x, F_diff, modules="numpy")
    f_diff_max = abs(np.max(f_diff(X)))
    R_n = max_end ** (n + 1) * f_diff_max / math.factorial(n + 1)

print(n)