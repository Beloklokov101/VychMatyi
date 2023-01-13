import math
from sympy import *
import numpy as np


t = symbols("t")

F = sin(t) * exp((- t ** 2) / 1e4) * cos(t / 50)


func = lambdify(t, F, modules="numpy")
grid = 1001
A, B = 0, 300
count_of_seg = 3 * (grid - 1)
h = (B - A) / count_of_seg

knots = np.zeros(count_of_seg + 1)
for i in range(count_of_seg + 1):
    knots[i] = A + i * h
f = func(knots)
I = 0
#print((3 * h) * (f[0] + 3 * f[1] + 3 * f[2] + f[3]) / 8)
for i in range(grid - 1):
    I_i = (3 * h) * (f[3 * i] + 3 * f[3 * i + 1] + 3 * f[3 * i + 2] + f[3 * i + 3]) / 8
    I = I + I_i

print("Integral equals: ", I / 100)