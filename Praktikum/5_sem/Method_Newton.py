import numpy as np
import math
from sympy import *


#Segment with root
a, b = 7.6, 7.8
#Initial approximation
x_prev = a
#Accuracy of root search
e = 1e-6
#Number of iterations
n = 0

x = symbols("x")
F = sin(x) - exp(1 / (1 - x ** 2))
f = lambdify(x, F, modules="numpy")

X = np.linspace(a, b, int(1e6))
flag = False
if (np.max(f(X)) * np.min(f(X)) <= 0):  flag = True

if(flag):
    F_diff = diff(F, x)
    f_diff = lambdify(x, F_diff)

    x = x_prev - f(x_prev) / f_diff(x_prev)
    n = 1
    while(abs(x - x_prev) > e):
        n = n + 1
        x_prev = x
        x = x_prev - f(x_prev) / f_diff(x_prev)
        #print(abs(x - x_prev))

print(x, n)