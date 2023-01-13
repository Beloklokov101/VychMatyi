import numpy as np
import math as mt
import sympy as sym
import matplotlib.pyplot as plt

def convection_equation_solver(u_0, psi_set, psi_dt_set, psi_2dt_set, psi_3dt_set, L, N, h, tau, t_grid_0):
    a_n = - (2 * t_grid_0 + 3)
    a_dt = - 2

    u_dx = - psi_dt_set / a_n
    u_2dx = (psi_2dt_set + a_dt * u_dx) / (a_n ** 2)
    u_3dx = - psi_3dt_set / (a_n ** 3) + 3 * a_dt * u_2dx / (a_n ** 2)

    u_L = psi_set
    u_L_1 = u_L - h * u_dx + h**2 * u_2dx / 2 - h**3 * u_3dx / 6
    u_L_2 = u_L - 2 * h * u_dx + 2 * h**2 * u_2dx - 4 * h**3 * u_3dx / 3

    u_next = u_0.copy()
    for n in range(1, N + 1):
        # print(n)
        u_open = u_next.copy()

        u_next = np.zeros(L + 1)
        u_next[L] = u_L[n]
        u_next[L - 1] = u_L_1[n]
        u_next[L - 2] = u_L_2[n]

        for l in range(L - 2):
            u_l_3 = u_open[l + 3]
            u_l_2 = u_open[l + 2]
            u_l_1 = u_open[l + 1]
            u_l = u_open[l]
            first_sum = (2 * n * tau + 3 + tau) * (2 * u_l_3 - 9 * u_l_2 + 18 * u_l_1 - 11 * u_l) * tau / (6 * h)
            second_sum = (2 * n * tau + 3) * (2 * n * tau + 3 + 2 * tau) * (- u_l_3 + 4 * u_l_2 - 5 * u_l_1 + 2 * u_l) * tau**2 / (2 * h**2)
            third_sum = (2 * n * tau + 3)**3 * (u_l_3 - 3 * u_l_2 + 3 * u_l_1 - u_l) * tau**3 / (6 * h**3)
            u_next[l] = u_l + first_sum + second_sum + third_sum
    
    return u_next


# Initial datas
N_L = 4     # N / L
tau_h = 1 / N_L     #tau / h

L_start = 20
L_end = 40
L_step = 20
count_L = int((L_end - L_start) / L_step) + 1

N_start = N_L * L_start
N_end = N_L * L_end
N_step = N_L * L_step
# N = N_L * L

# tau = 1 / N
# h = 1 / L

# x_grid_0 = np.linspace(0, 1, L + 1)
# t_grid_0 = np.linspace(0, 1, N + 1)

x, t = sym.symbols("x t")
phi_fun = sym.ln(1 + x ** 2)
psi_fun = sym.ln(1 + (1 + 3 * t + t ** 2) ** 2)

phi = sym.lambdify(x, phi_fun, modules="numpy")

psi = sym.lambdify(t, psi_fun, modules="numpy")
psi_dt = sym.lambdify(t, sym.diff(psi_fun, t, 1), modules="numpy")
psi_2dt = sym.lambdify(t, sym.diff(psi_fun, t, 2), modules="numpy")
psi_3dt = sym.lambdify(t, sym.diff(psi_fun, t, 3), modules="numpy")

# plt.plot(t_grid_0, psi_set)
# plt.show()


f = open("C:/My_Progs/VychMatyi/Praktikum/6_sem/4-th_lab/4-th_lab.txt", "w")
f.write("Laboratory work number: 4\n")
f.write("Variant number: 5\n\n")
f.write(f"Current number tau / h = {tau_h}\n")
f.write(f"Number of distance segments L is in {L_start} : {L_end} : {L_step}\n")
f.write(f"Number of time segments N is in {N_start} : {N_end} : {N_step}\n\n")

all_errors = np.zeros(count_L)
j = 0
for L in range(L_start, L_end + L_step, L_step):
    print(f"L = {L}")
    f.write(f"L = {L}\n\n")
    N = N_L * L
    tau = 1 / N
    h = 1 / L   

    x_grid_0 = np.linspace(0, 1, L + 1)
    t_grid_0 = np.linspace(0, 1, N + 1)

    u_0 = phi(x_grid_0)

    psi_set = psi(t_grid_0)
    psi_dt_set = psi_dt(t_grid_0)
    psi_2dt_set = psi_2dt(t_grid_0)
    psi_3dt_set = psi_3dt(t_grid_0)

    u_last_t1 = convection_equation_solver(u_0, psi_set, psi_dt_set, psi_2dt_set, psi_3dt_set, L, N, h, tau, t_grid_0)
    answer_theory = np.log(1 + (x_grid_0 + 4) ** 2)

    err_stepL_set = np.zeros(11)
    i = 0
    for it in np.linspace(0, L, 11):
        diff = abs(u_last_t1[int(it)] - answer_theory[int(it)])
        err_stepL_set[i] = diff
        # print(f"err[{int(it)}] = {diff}")

        f.write(f"x_grid[{int(it)}] = {x_grid_0[int(it)]}\n")
        f.write(f"u_exp[{int(it)}] = {u_last_t1[int(it)]}\n")
        f.write(f"u_theory[{int(it)}] = {answer_theory[int(it)]}\n")
        f.write(f"err[{int(it)}] = {diff}\n\n")

        i += 1

    err_max = np.max(err_stepL_set)
    all_errors[j] = err_max
    print(f"Max error = {err_max}\n")
    f.write(f"Max error = {err_max}\n\n")

    j += 1

plt.figure()
plt.plot(x_grid_0, u_last_t1, "-r", label="exp")
plt.plot(x_grid_0, answer_theory, ".-b", label="theory")

plt.figure()
err_x_grid = np.linspace(L_start, L_end, count_L)
plt.plot(err_x_grid, all_errors, "*-r", label=f"tau/h = {tau_h}")

plt.legend()
plt.show()
f.close()