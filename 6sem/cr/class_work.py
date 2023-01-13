import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import math as mt


x = symbols("x")
A = mt.pi
T = 10

num_x = 100 * 3
num_t = 2000
dx = A / num_x
dt = T / num_t
#print(dx, dt)

 #Начальные условия
u_ = sin(x)
v_ = cos(x)
u = lambdify(x, u_, modules="numpy")
v = lambdify(x, v_, modules="numpy")

p = 0

 #Определение пирамидального массива для будущих значений U и V
arr = np.array([np.zeros(num_x + 1 + 2 * num_t - 2 * j) for j in range(num_t + 1)])
#print(arr[0])
U, V = arr, arr
#print(U[0].shape[0])

#print(u(np.linspace(start=0, stop=float(A), num=int(U[0].shape[0]), dtype=float)).shape)

U[0] = u(np.linspace(0, A, U[0].shape[0]))
V[0] = v(np.linspace(0, A, V[0].shape[0]))

#for t in range(1, num_t + 1):
