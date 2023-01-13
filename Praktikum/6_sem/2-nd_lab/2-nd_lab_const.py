import matplotlib.pyplot as plt
import numpy as np
from sympy import *

def therm_eq_const(x0, L, k_a, k_b, q_a, q_b, f_a, f_b, u0, u1):
    h = 1 / L

    l_a = int(x0 / h)
    l_b = l_a + 1

        #a, b, c, d
    a = [0] + [k_a for i in range(1, l_a)] + [0, 0] + [k_b for i in range(l_b + 1, L)] + [0]
    b = [0] + [- 2 * k_a - q_a * h**2 for i in range(1, l_a)] + [0, 0] + [- 2 * k_b - q_b * h**2 for i in range(l_b + 1, L)] + [0]
    c = [0] + [k_a for i in range(1, l_a)] + [0, 0] + [k_b for i in range(l_b + 1, L)] + [0]
    d = [0] + [- f_a * h**2 for i in range(1, l_a)] + [0, 0] + [- f_b * h**2 for i in range(l_b + 1, L)] + [0]

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
    u[l_a] = (k_a * beta[l_a - 1] + k_b * beta[l_b + 1]) / (k_a * (1 - alpha[l_a - 1]) + k_b * (1 - alpha[l_b + 1]))
    u[l_b] = u[l_a]

    for l in range(l_a - 1, 0, -1):
        u[l] = alpha[l] * u[l + 1] + beta[l]

    for l in range(l_b + 1, L):
        u[l] = alpha[l] * u[l - 1] + beta[l]

    return u



x = symbols("x")

    #Initial datas
x0 = 0.525
u0, u1 = 0, 1
k_a, k_b = x0, x0**2 + 1
q_a, q_b = np.exp(-x0), np.exp(-x0)
f_a, f_b = x0**3, 1

    #Theoretical solution

    #Lambda and mu
lambda_a, lambda_b = np.sqrt(q_a / k_a), np.sqrt(q_b / k_b)
m_a, m_b = f_a / q_a, f_b / q_b

    #A and B
a11 = np.exp(- lambda_a * x0) - np.exp(lambda_a * x0)
a12 = np.exp(lambda_b * (2 - x0)) - np.exp(lambda_b * x0)
a21 = k_a * lambda_a * (np.exp(lambda_a * x0) + np.exp(- lambda_a * x0))
a22 = k_b * lambda_b * (np.exp(lambda_b * (2 - x0)) + np.exp(lambda_b * x0))

b1 = m_b - m_a + (m_a - u0) * np.exp(lambda_a * x0) - (m_b - u1) * np.exp(lambda_b * (1 - x0))
b2 = k_a * lambda_a * (u0 - m_a) * np.exp(lambda_a * x0) + k_b * lambda_b * (u1 - m_b) * np.exp(lambda_b *(1 - x0))

    #C
c1 = (((u0 - m_a) * a11 - b1) * a22 - ((u0 - m_a) * a21 - b2) * a12) / (a11 * a22 - a12 * a21)
c2 = (b1 * a22 - b2 * a12) / (a11 * a22 - a12 * a21)
c3 = (b2 * a11 - b1 * a21) / (a11 * a22 - a12 * a21)
c4 = (u1 - m_b) * np.exp(lambda_b) - c3 * np.exp(2 * lambda_b)
#print(c1, c2, c3, c4)

    #Theor. solution

x_set_first = np.linspace(0, x0, 100)
x_set_second = np.linspace(x0, 1, 100)

plt.figure()
plt.plot(x_set_first, c1 * np.exp(lambda_a * x_set_first) + c2 * np.exp(- lambda_a * x_set_first) + m_a)
plt.plot(x_set_second, c3 * np.exp(lambda_b * x_set_second) + c4 * np.exp(- lambda_b * x_set_second) + m_b)

u_theor_a = lambdify(x, c1 * exp(lambda_a * x) + c2 * exp(- lambda_a * x) + m_a)
u_theor_b = lambdify(x, c3 * exp(lambda_b * x) + c4 * exp(- lambda_b * x) + m_b)

    #Numerical solution

N = 1e3
L = int(10 * N)
L = 81921
err = 1e-5

u = therm_eq_const(x0, L, k_a, k_b, q_a, q_b, f_a, f_b, u0, u1)
err_step = err * 2

# for t in range(1):
#     u_prev = u
#     L = 2 * L
#     u = therm_eq_const(x0, L, k_a, k_b, q_a, q_b, f_a, f_b, u0, u1)
#     err_set = [0 for i in range(11)]
#     for i in np.linspace(0, L, 11):
#         err_set[int(i * 10 / L)] = abs(u[int(i)] - u_prev[int(i/2)])
#     err_step = max(err_set)
#     print(err_step, L)

x_grid = np.linspace(0, 1, L + 1)

    #Theory values
u_theor = [0 for i in range(L + 1)]
for i in range(len(x_grid)):
    if x_grid[i] <= x0:
        u_theor[i] = u_theor_a(x_grid[i])
    else:
        u_theor[i] = u_theor_b(x_grid[i])


f = open("C:/My_Progs/VychMatyi/Praktikum/6_sem/2-nd_lab/2-nd_lab_const.txt", "w")
f.write("Laboratory work number: 2\n")
f.write("Variant number: 14\n\n")
f.write(f"Number of segments L = {L}\n\n")

diff = [0 for i in range(11)]

for i in np.linspace(0, L, 11):
    f.write(f"X[{int(i)}] = {x_grid[int(i)]}\n")
    f.write(f"u_theor = {u_theor[int(i)]}\n")
    f.write(f"u[{int(i)}] = {u[int(i)]}\n")
    f.write(f"diff = {abs(u[int(i)] - u_theor[int(i)])}\n\n")
    diff[int(i * 10 / L)] = abs(u[int(i)] - u_theor[int(i)])
    #f.write(f"u_prev[{int(i/2)}] = {u_prev[int(i/2)]}\n")
    #f.write(f"diff = {abs(u[int(i)] - u_prev[int(i/2)])}\n\n")

print(f"Max error = {max(diff)}\n")
#plt.figure()
plt.plot(x_grid, u)

plt.show()
f.close()