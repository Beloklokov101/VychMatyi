import numpy as np
from sympy import *
import matplotlib.pyplot as plt


# k = 4
def ode_4(f, g, t0, v0, y0, h):
    k1 = f(t0, v0, y0)
    m1 = g(t0, v0, y0)
    
    k2 = f(t0 + h/2, v0 + k1*h/2, y0 + m1*h/2)
    m2 = g(t0 + h/2, v0 + k1*h/2, y0 + m1*h/2)
    
    k3 = f(t0 + h/2, v0 + k2*h/2, y0 + m2*h/2)
    m3 = g(t0 + h/2, v0 + k2*h/2, y0 + m2*h/2)
    
    k4 = f(t0 + h, v0 + k3*h, y0 + m3*h)
    m4 = g(t0 + h, v0 + k3*h, y0 + m3*h)
    
    #print("\n", k1, k2, k3, k4, y0, "\n")
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


y1, y2, x = symbols("y1, y2, x")

# A, B, D = -5, 2, 100
A, B, D = 0, 5, 10

 #Theoretical answer

C1 = A/10 + B/4
C2 = A/10 - B/4

y1th = C1 * 5 * exp(x) + C2 * 5 * exp(-199 * x)
y2th = C1 * 2 * exp(x) - C2 * 2 * exp(-199 * x)
Y1th = lambdify(x, y1th, "numpy")
Y2th = lambdify(x, y2th, "numpy")

 #Plotting

F = - 99 * y1 + 250 * y2
G = 40 * y1 - 99 * y2
f = lambdify([x, y1, y2], F, "numpy")
g = lambdify([x, y1, y2], G, "numpy")

y1f, y2f = A, B
a, b = 0, D
length = b - a
err = 1e-6
# N = 160000
N = 16000

h = length / N

X, Y1, Y2 = ode4sys(f, g, a, y1f, y2f, N, h)
err_step = err * 2
X0, Y10, Y20 = X, Y1, Y2
"""
while(err_step >= err):
    N = 2 * N
    h = length / N
    X, Y1, Y2 = ode4sys(f, g, a, y1f, y2f, N, h)

    n_step = int(N / 10 / 2)
    
    err_set1 = [abs(Y10[j * n_step] - Y1[2 * j * n_step]) for j in range(11)]
    err_set2 = [abs(Y20[j * n_step] - Y2[2 * j * n_step]) for j in range(11)]
    #err_set1 = [abs(Y1[2 * j * n_step] - Y1th(D * j / 10)) for j in range(11)]
    #err_set2 = [abs(Y2[2 * j * n_step] - Y2th(D * j / 10)) for j in range(11)]
    
    print(f"Errors in {N}: {max(err_set1)}, {max(err_set2)}")
    err_step = max(max(err_set1), max(err_set2))
    X0, Y10, Y20 = X, Y1, Y2
"""
n_step = N / 10 / 2


plt.plot(X, Y1, label="Y1exp")
plt.plot(X, Y2, label="Y2exp")

plt.plot(X, Y1th(np.linspace(0, D, N + 1)), label="Y1th")
plt.plot(X, Y2th(np.linspace(0, D, N + 1)), label="Y2th")


f = open("C:/My_Progs/VychMatyi/Praktikum/6_sem/1-st_lab/1-st_lab.txt", "w")
f.write("Laboratory work number: 1\n")
f.write("Variant number: 16\n\n")
f.write(f"Number of segments N = {N}\n\n")

for i in range(11):
    f.write(f"X = {X[int(2 * i * n_step)]:.3f}\n")
    f.write("{:>2}Y1_th = {}\n".format("", Y1th(D * i / 10)))
    f.write("{:>2}Y1_ex = {}\n".format("", Y1[int(2 * i * n_step)]))
    f.write("{:>2}Diff1 = {}\n\n".format("", abs(Y1[int(2 * i * n_step)] - Y1th(D * i / 10))))

    f.write("{:>2}Y2_th = {}\n".format("", Y2th(D * i / 10)))
    f.write("{:>2}Y2_ex = {}\n".format("", Y2[int(2 * i * n_step)]))
    f.write("{:>2}Diff2 = {}\n\n".format("", abs(Y2[int(2 * i * n_step)] - Y2th(D * i / 10))))
        

f.close()
plt.legend()
plt.show()

