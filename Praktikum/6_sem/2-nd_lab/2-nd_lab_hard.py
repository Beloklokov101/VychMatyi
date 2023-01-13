import matplotlib.pyplot as plt
import numpy as np
from sympy import *

def therm_eq(x0, L, k_a, k_b, q_a, q_b, f_a, f_b, u0, u1):
    h = 1 / L
    x_grid = np.linspace(0, 1, L + 1)

    l_a = int(x0 / h)
    l_b = l_a + 1

        #a, b, c, d
    a = [0] + [k_a(x_grid[i] + h/2) for i in range(1, l_a)] + [0, 0] + [k_b(x_grid[i] + h/2) for i in range(l_b + 1, L)] + [0]
    b = [0] + [- (k_a(x_grid[i] + h/2) + k_a(x_grid[i] - h/2) + q_a(x_grid[i]) * h**2) for i in range(1, l_a)] 
    b = b + [0, 0] + [- (k_b(x_grid[i] + h/2) + k_b(x_grid[i] - h/2) + q_b(x_grid[i]) * h**2) for i in range(l_b + 1, L)] + [0]
    c = [0] + [k_a(x_grid[i] - h/2) for i in range(1, l_a)] + [0, 0] + [k_b(x_grid[i] - h/2) for i in range(l_b + 1, L)] + [0]
    d = [0] + [- f_a(x_grid[i]) * h**2 for i in range(1, l_a)] + [0, 0] + [- f_b(x_grid[i]) * h**2 for i in range(l_b + 1, L)] + [0]

        #Alpha and beta
    alpha = [0 for i in range(L + 1)]
    beta = [u0] + [0 for i in range(1, L)] + [u1]

    for l in range(1, l_a):
        alpha[l] = - a[l] / (b[l] + c[l] * alpha[l - 1])
        beta[l] = (d[l] - c[l] * beta[l - 1]) / (b[l] + c[l] * alpha[l - 1])

    for l in range(L - 1, l_b, -1):
        alpha[l] = - c[l] / (b[l] + a[l] * alpha[l + 1])
        beta[l] = (d[l] - a[l] * beta[l + 1]) / (b[l] + a[l] * alpha[l + 1])

    #print(alpha, beta)

        #Calculating u[]
    u = [u0] + [0 for i in range(1, L)] + [u1]

    u[l_a] = (k_a(x_grid[l_a]) * beta[l_a - 1] + k_b(x_grid[l_b]) * beta[l_b + 1]) / (k_a(x_grid[l_a]) * (1 - alpha[l_a - 1]) + k_b(x_grid[l_b]) * (1 - alpha[l_b + 1]))
    u[l_b] = u[l_a]

    for l in range(l_a - 1, 0, -1):
        u[l] = alpha[l] * u[l + 1] + beta[l]

    for l in range(l_b + 1, L):
        u[l] = alpha[l] * u[l - 1] + beta[l]
    
    return u



x = symbols("x")

    #initial datas
x0 = 0.525
u0, u1 = 0, 1
k_a = lambdify(x, x)
k_b = lambdify(x, x**2 + 1)
q_a = lambdify(x, exp(-x))
q_b = lambdify(x, exp(-x))
f_a = lambdify(x, x**3)
f_b = lambdify(x, 1)

err = 1e-4

N = 1e2
L = int(10 * N)
L = 163841 * 2 - 1

u = therm_eq(x0, L, k_a, k_b, q_a, q_b, f_a, f_b, u0, u1)
err_step = err * 2

#while(err_step >= err and L < 10000000):
for t in range(1):
    u_prev = u
    L = 2 * L
    u = therm_eq(x0, L, k_a, k_b, q_a, q_b, f_a, f_b, u0, u1)
    err_set = [0 for i in range(11)]
    for i in np.linspace(0, L, 11):
        err_set[int(i * 10 / L)] = abs(u[int(i)] - u_prev[int(i/2)])
    err_step = max(err_set)
    print(err_step, L)

x_grid = np.linspace(0, 1, L + 1)

plt.figure()
plt.plot(x_grid, u)


    #Output in text

f = open("C:/My_Progs/VychMatyi/Praktikum/6_sem/2-nd_lab/2-nd_lab_hard.txt", "w")
f.write("Laboratory work number: 2\n")
f.write("Variant number: 14\n\n")
f.write(f"Number of segments L = {L}\n\n")

diff = [0 for i in range(11)]

for i in np.linspace(0, L, 11):
    f.write(f"X[{int(i)}] = {x_grid[int(i)]}\n")
    f.write(f"u[{int(i)}] = {u[int(i)]}\n")
    f.write(f"u_prev[{int(i/2)}] = {u_prev[int(i/2)]}\n")
    f.write(f"diff = {abs(u[int(i)] - u_prev[int(i/2)])}\n\n")
    diff[int(i * 10 / L)] = abs(u[int(i)] - u_prev[int(i/2)])
    
print(f"Max error = {max(diff)}\n")

plt.show()