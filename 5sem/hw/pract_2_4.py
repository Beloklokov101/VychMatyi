import numpy as np
from sympy import *
import matplotlib.pyplot as plt

# k = 4
def ode_4(f, g, t0, v0, y0, h):
    k1 = f(t0, v0, y0)
    m1 = g(t0, v0, y0)
    #print(t0, v0, y0, h, k1, m1)
    
    k2 = f(t0 + h/2, v0 + k1*h/2, y0 + m1*h/2)
    m2 = g(t0 + h/2, v0 + k1*h/2, y0 + m1*h/2)
    
    k3 = f(t0 + h/2, v0 + k2*h/2, y0 + m2*h/2)
    m3 = g(t0 + h/2, v0 + k2*h/2, y0 + m2*h/2)
    
    k4 = f(t0 + h, v0 + k3*h, y0 + m3*h)
    m4 = g(t0 + h, v0 + k3*h, y0 + m3*h)
    
    #print("\n", k1, k2, k3, k4, v0, "\n")
    #print(m1, m2, m3, m4, y0, "\n")
    v1 = v0 + h*(k1 + 2*k2 + 2*k3 + k4)/6
    y1 = y0 + h*(m1 + 2*m2 + 2*m3 + m4)/6
    t1 = t0 + h

    return (t1, v1, y1)


def ode4sys(f, g, a, v, y_f, N, h):
    T, V, Y = np.zeros(N + 1), np.zeros(N + 1), np.zeros(N + 1)
    T[0], V[0], Y[0] = a, v, y_f
    for i in range(1, N + 1):
        T[i], V[i], Y[i] = ode_4(f, g, T[i-1], V[i-1], Y[i-1], h)
    
    return (T, V, Y)


v, y, t = symbols("v, y, t")

F = sqrt(1/t**2 + exp(1)*y**2 / log(t) - exp(v) * y)
G = v
f = lambdify([t, v, y], F, "numpy")
g = lambdify([t, v, y], G, "numpy")

a, b = float(E), float(E**2) 
y_f, y_l = float(E), float(2*E**2)

alpha = 1

N = 10000
h = float((b - a) / N)


start_alpha = alpha + 1e3

while(abs(start_alpha - alpha) > 1e-3):
    start_alpha = alpha
    T, V, Y = ode4sys(f, g, a, alpha, y_f, N, h)
    #plt.plot(T, Y)
    #print(Y)
    #plt.plot(T, V)
    #print(f"alpha = {alpha}")
    r0 = Y[-1] - y_l
    #print(f"r0 = {r0}")

    da = alpha / 1e3
    new_alpha = alpha + da
    
    T, V, Y = ode4sys(f, g, a, new_alpha, y_f, N, h)
    #plt.plot(T, Y)
    #plt.plot(T, V)

    r1 = Y[-1] - y_l
    print(f"r1 = {r1}")

    dr = (r1 - r0) / da
    print(f"dr = {dr}")
    alpha = alpha - r0 / dr

    print(f"new alpha = {alpha}", "\n")

plt.plot(T, Y)
plt.show()