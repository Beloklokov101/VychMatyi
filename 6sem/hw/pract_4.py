import numpy as np
import matplotlib.pyplot as plt
import math as mt

def TriagonalSlae_solver(a, b, c, d, N, M = 1, transpose = False):
    # np.shape(a or b or c or d) == (N, M)

    # a[0] * u[1] + b[0] * u[0] = d[0]
    # a[i] * u[i + 1] + b[i] * u[i] + c[i] * u[i - 1] = d[i], i = range(1, N - 1)
    # b[N - 1] * u[N - 1] + c[N - 1] * u[N - 2] = d[N - 1]

    # (!) There should be c[0] = 0, a[N - 1] = 0
    # This solver can solve M SLAE's simultaneously

    # u[0] = f(m) <==> a[0] = 0, b[0] = 1, d[0] = f(m)
    # u[N - 1] = g(m) <==> c[N - 1] = 0, b[N - 1] = 1, d[N - 1] = g(m)

    dimens = (N, M)
    shape_a = np.shape(a)
    shape_b = np.shape(b)
    shape_c = np.shape(c)
    shape_d = np.shape(d)

    if M != 1:
        if (dimens != shape_a or dimens != shape_b or dimens != shape_c or dimens != shape_d):
            print("\nERROR IN TriagonalSlae_solver")
            print(f"shape_a = {shape_a}, (N, M) = {dimens}\n")
            return False
    else:
        if (dimens[0] != shape_a[0] or dimens[0] != shape_b[0] or dimens[0] != shape_c[0] or dimens[0] != shape_d[0]):
            print("\nERROR IN TriagonalSlae_solver")
            print(f"shape_a = {shape_a}, (N, M) = {dimens}\n")
            return False
    
    if transpose:
        (N, M) = (M, N)
        a, b, c, d = a.T, b.T, c.T, d.T

    # u[i] = alpha[i] * u[i + 1] + beta[i]
    alpha, beta = np.zeros((N - 1, M)), np.zeros((N - 1, M))
    alpha[0] = - a[0] / b[0]
    beta[0] = d[0] / b[0]
    for it in range(1, N - 1):
        divider = b[it] + c[it] * alpha[it - 1]
        alpha[it] = - a[it] / divider
        beta[it] = (d[it] - c[it] * beta[it - 1]) / divider
    
    if M != 1:
        u = np.zeros((N, M))
    else:
        u = np.zeros(N)
    
    u[N - 1] = (d[N - 1] - c[N - 1] * beta[N - 2]) / (b[N - 1] + c[N - 1] * alpha[N - 2])
    for it in range(N - 2, -1, -1):
        u[it] = alpha[it] * u[it + 1] + beta[it]

    if transpose:
        return u.T
    
    return u


L_x = 50
L_y = 50
N = 50
h_x = 1 / L_x
h_y = 1 / L_y
tau = 1 / N
time_out = np.array([0.2, 0.4, 0.5, 0.6, 0.8, 1.0])
time_out_size = time_out.shape[0]

x_grid, y_grid = np.linspace(0, 1, L_x + 1), np.linspace(0, 1, L_y + 1)
t_grid = np.linspace(0, 1, N + 1)

lamb_const = 1e-4
sigma_x = lamb_const * tau / h_x**2
sigma_y = lamb_const * tau / h_y**2
a, c = 25 * sigma_x, 25 * sigma_x
e, d = sigma_y, sigma_y
b = - 25 * 2 * sigma_x - 2 * sigma_y - 1


    #for j = const
a_i, c_i = a * np.ones(L_x + 1), c * np.ones(L_x + 1)
e_i, d_i = e * np.ones(L_x + 1), d * np.ones(L_x + 1)
b_i = b * np.ones(L_x + 1)

    # for i = const
a_j, c_j = a * np.ones(L_y + 1), c * np.ones(L_y + 1)
e_j, d_j = e * np.ones(L_y + 1), d * np.ones(L_y + 1)
b_j = b * np.ones(L_y + 1)

b_i[0], b_i[-1] = 1, 1
a_i[0], a_i[-1], c_i[0], c_i[-1] = 0, 0, 0, 0
e_i[0], e_i[-1], d_i[0], d_i[-1] = 0, 0, 0, 0

b_j[0], b_j[-1] = 1, 1
a_j[0], a_j[-1], c_j[0], c_j[-1] = 0, 0, 0, 0
e_j[0], e_j[-1], d_j[0], d_j[-1] = 0, 0, 0, 0


phi_0 = np.outer(np.cos(mt.pi * x_grid), np.sin(5 * mt.pi * y_grid))
phi_next = phi_0.copy()

phi_time_out_errors = np.zeros(time_out_size)
time_count = 0
for n in range(N):
    phi_open = phi_next.copy()
    phi_next = np.zeros((L_x + 1, L_y + 1))

    # a * phi[i - 1, j] + b * phi[i, j] + c * phi[i + 1, j] + d * phi[i, j - 1] + e * phi[i, j + 1]

    f = - phi_open.copy()

    phi_next_k = phi_open.copy()
    for k in range(1):
        phi_open_k = phi_next_k.copy()
        phi_next_k = np.zeros((L_x + 1, L_y + 1))

        # j = const, ..._i
        phi_next_k[:, 0], phi_next_k[:, -1] = 0, 0
        for j in range(1, L_y):
            d_triag = np.zeros(L_x + 1)
            d_triag[0] = np.sin(5 * mt.pi * y_grid[j]) * np.exp(- 50 * mt.pi**2 * lamb_const * t_grid[n + 1])
            d_triag[-1] = - d_triag[0]
            d_triag[1 : L_x] = f[1 : L_x, j] - d * phi_next_k[1 : L_x, j - 1] - e * phi_open_k[1 : L_x, j + 1]
            phi_next_k[:, j] = TriagonalSlae_solver(c_i, b_i, a_i, d_triag, N = L_x + 1)

        phi_open_k = phi_next_k.copy()
        phi_next_k = np.zeros((L_x + 1, L_y + 1))

        # i = const, ...j
        phi_next_k[0] = np.sin(5 * mt.pi * y_grid) * np.exp(- 50 * mt.pi**2 * lamb_const * t_grid[n + 1])
        phi_next_k[-1] = - phi_next_k[0]
        for i in range(1, L_x):
            d_triag = np.zeros(L_y + 1)
            d_triag[0], d_triag[-1] = 0, 0
            d_triag[1 : L_y] = f[i, 1 : L_y] - a * phi_next_k[i - 1, 1 : L_y] - c * phi_open_k[i + 1, 1 : L_y]
            phi_next_k[i] = TriagonalSlae_solver(e_j, b_j, d_j, d_triag, N = L_y + 1)

        phi_open_k = phi_next_k.copy()
        phi_next_k = np.zeros((L_x + 1, L_y + 1))

        # BACK j = const, ..._i
        phi_next_k[:, 0], phi_next_k[:, -1] = 0, 0
        for j in range(L_y - 1, 0, -1):
            d_triag = np.zeros(L_x + 1)
            d_triag[0] = np.sin(5 * mt.pi * y_grid[j]) * np.exp(- 50 * mt.pi**2 * lamb_const * t_grid[n + 1])
            d_triag[-1] = - d_triag[0]
            d_triag[1 : L_x] = f[1 : L_x, j] - d * phi_open_k[1 : L_x, j - 1] - e * phi_next_k[1 : L_x, j + 1]
            phi_next_k[:, j] = TriagonalSlae_solver(c_i, b_i, a_i, d_triag, N = L_x + 1)

        phi_open_k = phi_next_k.copy()
        phi_next_k = np.zeros((L_x + 1, L_y + 1))

        # BACK i = const, ...j
        phi_next_k[0] = np.sin(5 * mt.pi * y_grid) * np.exp(- 50 * mt.pi**2 * lamb_const * t_grid[n + 1])
        phi_next_k[-1] = - phi_next_k[0]
        for i in range(L_x - 1, 0, -1):
            d_triag = np.zeros(L_y + 1)
            d_triag[0], d_triag[-1] = 0, 0
            d_triag[1 : L_y] = f[i, 1 : L_y] - a * phi_open_k[i - 1, 1 : L_y] - c * phi_next_k[i + 1, 1 : L_y]
            phi_next_k[i] = TriagonalSlae_solver(e_j, b_j, d_j, d_triag, N = L_y + 1)
    
    phi_next = phi_next_k.copy()

    if t_grid[n + 1] >= time_out[time_count]:
        phi_theory_n = phi_0 * np.exp(- 50 * mt.pi**2 * lamb_const * t_grid[n + 1])
        phi_err_max_n = np.max(np.abs(phi_next - phi_theory_n))
        phi_time_out_errors[time_count] = phi_err_max_n
        time_count += 1

    if n % 5 == 0: print(f"n = {n}")

if time_count != time_out_size:
    phi_theory_n = phi_0 * np.exp(- 50 * mt.pi**2 * lamb_const * t_grid[n + 1])
    phi_err_max_n = np.max(np.abs(phi_next - phi_theory_n))
    phi_time_out_errors[time_count] = phi_err_max_n
    time_count += 1

phi_final = phi_next.copy()

    # Saving errors
# f = open("C:/My_Progs/VychMatyi/6sem/hw/pract_4_errors.txt", "a")
with open("C:/My_Progs/VychMatyi/6sem/hw/pract_4_errors.txt", "a") as f:
    f.write(f"L = {L_x}\n")
    for it in range(time_out_size):
        f.write(f"{phi_time_out_errors[it]} ")
        print(f"time = {time_out[it]}, err = {phi_time_out_errors[it]}")
    f.write("\n")

    # Theory
phi_theory_final = phi_0 * np.exp(- 50 * mt.pi**2 * lamb_const)

    # Plotting
XMeshgrid, YMeshgrid = np.meshgrid(x_grid, y_grid)

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.contour3D(XMeshgrid, YMeshgrid, phi_final, 75, cmap="viridis")
ax.set_title("phi_exp")

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.contour3D(XMeshgrid, YMeshgrid, phi_theory_final, 75, cmap="viridis")
ax.set_title("phi_theory")

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.contour3D(XMeshgrid, YMeshgrid, abs(phi_final - phi_theory_final), 75, cmap="viridis")
ax.set_title("phi_err")

# f.close()
plt.show()