import numpy as np
from sympy import *
import matplotlib.pyplot as plt


# k = 4
def ode_4(f, g, t0, v0, y0, h):
    #print(t0, v0, y0)
    k1 = f(t0, v0, y0)
    m1 = g(t0, v0, y0)
    #print(k1, m1)
    
    k2 = f(t0 + h/2, v0 + k1*h/2, y0 + m1*h/2)
    m2 = g(t0 + h/2, v0 + k1*h/2, y0 + m1*h/2)
    #print(k2, m2)
    
    k3 = f(t0 + h/2, v0 + k2*h/2, y0 + m2*h/2)
    m3 = g(t0 + h/2, v0 + k2*h/2, y0 + m2*h/2)
    #print(k3, m3)
    
    k4 = f(t0 + h, v0 + k3*h, y0 + m3*h)
    m4 = g(t0 + h, v0 + k3*h, y0 + m3*h)
    #print(k4, m4, "\n")
    
    #print("\n", k1, k2, k3, k4, y0, "\n")
    v1 = v0 + h*(k1 + 2*k2 + 2*k3 + k4)/6
    y1 = y0 + h*(m1 + 2*m2 + 2*m3 + m4)/6
    t1 = t0 + h

    return (t1, v1, y1)


def ode2(f, g, t0, v0, y0, h):
    k1 = f(t0, v0, y0)
    m1 = g(t0, v0, y0)

    k2 = f(t0 + h/2, v0 + k1*h/2, y0 + m1*h/2)
    m2 = g(t0 + h/2, v0 + k1*h/2, y0 + m1*h/2)

    v1 = v0 + h*k2
    y1 = y0 + h*m2
    t1 = t0 + h

    return (t1, v1, y1)


def ode4sys(f, g, a, b, v, y0, N):
    h = (b - a) / N
    T, V, Y = np.zeros(N + 1), np.zeros(N + 1), np.zeros(N + 1)
    T[0], V[0], Y[0] = a, v, y0
    for i in range(1, N + 1):
        T[i], V[i], Y[i] = ode_4(f, g, T[i-1], V[i-1], Y[i-1], h)
    
    return (T, V, Y)


v, y, t = symbols("v, y, t")

mu = 1000
F = mu * (1 - v**2) * v - y
#F = -2*v - y
G = v
f = lambdify([t, v, y], F, "numpy")
g = lambdify([t, v, y], G, "numpy")

a, b = 0, 1
y0, dy0 = 0, 0.001

N = 100000
h = (b - a) / N

for count in range(1):
    T, V, Y = ode4sys(f, g, a, b, dy0, y0, N)
    plt.plot(T, Y)
    plt.plot(T, V)
    
"""h = (b - a) / N
T, V, Y = np.zeros(N + 1), np.zeros(N + 1), np.zeros(N + 1)
T[0], V[0], Y[0] = a, dy0, y0
for i in range(1, N + 1):
    T[i], V[i], Y[i] = ode2(f, g, T[i-1], V[i-1], Y[i-1], h)
"""
plt.plot(T, Y)
plt.plot(T, V)
plt.show()