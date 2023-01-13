import numpy as np
import pandas as pd
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
    
    if (dimens != shape_a or dimens != shape_b or dimens != shape_c or dimens != shape_d):
        print("\nERROR IN TriagonalSlae_solver\n")
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
    
    u = np.zeros((N, M))
    u[N - 1] = (d[N - 1] - c[N - 1] * beta[N - 2]) / (b[N - 1] + c[N - 1] * alpha[N - 2])
    for it in range(N - 2, -1, -1):
        u[it] = alpha[it] * u[it + 1] + beta[it]

    if transpose:
        return u.T
    
    return u


L, M = 5, 5
N = 6
r_max, teta_max, time_max = 1, mt.pi / 2, 1
h_r, h_teta = r_max / L, teta_max / M
tau = time_max / N

# My case (N.15)
mu = 1 / 2
C_t = 14

f_error = open("C:/My_Progs/VychMatyi/Praktikum/6_sem/6-th_lab/thermal_cond_eq_active_errors.txt", "w")
f_results = open("C:/My_Progs/VychMatyi/Praktikum/6_sem/6-th_lab/thermal_cond_eq_results.txt", "w")

f_error.write("Laboratory number: 6\n")
f_error.write("Variant number: 2.15\n\n")
f_error.write(f"L = {L}, M = {M}\nN = {N}\n\n")
f_error.write(f"mu = {mu}, C_t = {C_t}\n")

f_results.write("Laboratory number: 6\n")
f_results.write("Variant number: 15\n\n")
f_results.write(f"L = {L}, M = {M}\nN = {N}\n\n")
f_results.write(f"mu = {mu}, C_t = {C_t}\n")

print(f"L = {L}, M = {M}\nN = {N}")

r_grid = np.linspace(0, r_max, L + 1)
teta_grid = np.linspace(0, teta_max, M + 1)
time_grid = np.linspace(0, time_max, N + 1)
sin_grid = np.sin(teta_grid)

u_0 = np.outer(r_grid, sin_grid) ** (2 / mu) * C_t ** (- 1 / mu)
u_next = u_0.copy()

for n in range(N):
    f_error.write(f"n = {n}:\n")
    u_open = u_next.copy()
    u_next_k = u_next.copy()
    for k in range(3):
        f_error.write(f" k = {k}:\n")

        u_open_k = u_next_k.copy()

        # Boundary conditions:
        u_Lm = (sin_grid ** (2 / mu)) * ((C_t - 4 * (mu + 1) * time_grid[n + 1] / mu) ** (- 1 / mu))
        u_lM = (r_grid ** (2 / mu)) * ((C_t - 4 * (mu + 1) * time_grid[n + 1] / mu) ** (- 1 / mu))

        #  First half-step, [inter]:
        u_inter = np.zeros((L + 1, M + 1))

        # a, b, c, d .shape == l = 0, L; m = 1, M - 1
        a_inter, b_inter = np.zeros((L + 1, M - 1)), np.zeros((L + 1, M - 1))
        c_inter, d_inter = np.zeros((L + 1, M - 1)), np.zeros((L + 1, M - 1))

        # c_inter[0] and a_inter[L] already equal 0
        # I have initial and final condition's, so: 
        # (a_inter[0] and c_inter[L] already equal 0)
        b_inter[0], b_inter[L] = np.ones(M - 1), np.ones(M - 1)
        # u_inter[0] = 0 <-- already done (d_inter[0] = 0)
        d_inter[L] = u_Lm[1 : M]
        
        for l in range(1, L):
            multip_inter = tau / (2 * (r_grid[l] * h_r) ** 2)
            multip_inter_a = u_open_k[l + 1, 1 : M] ** mu + u_open_k[l, 1 : M] ** mu
            multip_inter_c = u_open_k[l, 1 : M] ** mu + u_open_k[l - 1, 1 : M] ** mu
            a_inter[l] = - (r_grid[l] + h_r / 2) ** 2 * multip_inter_a * multip_inter
            c_inter[l] = - (r_grid[l] - h_r / 2) ** 2 * multip_inter_c * multip_inter
            b_inter[l] = 1 - a_inter[l] - c_inter[l]
            d_inter[l] = u_open[l, 1 : M]
        
        # u_inter[:, 0] and [:, M] will not be used
        u_inter[:, 1 : M] = TriagonalSlae_solver(a_inter, b_inter, c_inter, d_inter, L + 1, M - 1)

        #  Second hald-step, [next_k]:
        u_next_k = np.zeros((L + 1, M + 1))

        # a, b, c, d .shape == l = 1, L - 1; m = 0, M
        a_next_k, b_next_k = np.zeros((L - 1, M + 1)), np.zeros((L - 1, M + 1))
        c_next_k, d_next_k = np.zeros((L - 1, M + 1)), np.zeros((L - 1, M + 1))
        
        # c_next_k[:, 0] and a_next_k[:, M] already equal 0
        # I have initial and final condition's, so: 
        # (a_next_k[:, 0] and c_next_k[:, M] already equal 0)
        b_next_k[:, 0], b_next_k[:, M] = np.ones(L - 1), np.ones(L - 1)
        # u_next_k[:, 0] = 0 <-- already done (d_next_k[:, 0] = 0)
        d_next_k[:, M] = u_lM[1 : L]

        for m in range(1, M):
            multip_next_k = tau / (2 * sin_grid[m] * (r_grid[1 : L] * h_teta) ** 2)
            multip_next_k_a = u_open_k[1 : L, m + 1] ** mu + u_open_k[1 : L, m] ** mu
            multip_next_k_c = u_open_k[1 : L, m] ** mu + u_open_k[1 : L, m - 1] ** mu
            a_next_k[:, m] = - np.sin(teta_grid[m] + h_teta / 2) * multip_next_k_a * multip_next_k
            c_next_k[:, m] = - np.sin(teta_grid[m] - h_teta / 2) * multip_next_k_c * multip_next_k
            b_next_k[:, m] = 1 - a_next_k[:, m] - c_next_k[:, m]
            d_next_k[:, m] = u_inter[1 : L, m]
        
        u_next_k[1 : L, :] = TriagonalSlae_solver(a_next_k, b_next_k, c_next_k, d_next_k, L - 1, M + 1, transpose = True)

        # u_next_k[0] = 0 <-- already done
        u_next_k[L] = u_Lm

        #  Error at k-th step:
        rel_err_k = (u_next_k[1 : L, 1 : M] - u_open_k[1 : L, 1 : M]) / u_next_k[1 : L, 1 : M]
        max_rel_err_k = np.max(rel_err_k)
        f_error.write(f" rel_err_k = {max_rel_err_k}\n")

    u_next = u_next_k.copy()

    #  Error at n-th step:
    rel_err_n = (u_next[1 : L, 1 : M] - u_open[1 : L, 1 : M]) / u_next[1 : L, 1 : M]
    max_rel_err_n = np.max(rel_err_n)
    f_error.write(f"\n rel_err_n = {max_rel_err_n}\n\n")

u_final = u_next.copy()
u_theory = np.outer(r_grid, sin_grid) ** (2 / mu) * (C_t - 4 * (mu + 1) * time_max / mu) ** (- 1 / mu)

    # Saving results
r_indexes = np.arange(0, L + 1, int(L / 5))
teta_indexes = np.arange(0, M + 1, int(M / 5))
# print(u_theory[r_indexes, teta_indexes])

    # (!) INT only for this case
columns_indexes = np.degrees(teta_grid[teta_indexes])
df_theory = pd.DataFrame(u_theory[r_indexes, :][:, teta_indexes], index=r_grid[r_indexes], columns=np.degrees(teta_grid[teta_indexes]))
df_exp = pd.DataFrame(u_final[r_indexes, :][:, teta_indexes], index=r_grid[r_indexes], columns=np.degrees(teta_grid[teta_indexes]))

# f_results.write(f"Theory, time = {time_max}:\n")
# f_results.write(df_theory)
# f_results.write(f"\nExperimental, time = {time_max}:\n")
# f_results.write(df_exp)
print("Theory, time = {time_max}:")
print(df_theory)
print("Experimental, time = {time_max}:")
print(df_exp)
df_theory.to_csv("C:/My_Progs/VychMatyi/Praktikum/6_sem/6-th_lab/therm_cond_eq_theory.txt", sep="\t")
df_exp.to_csv("C:/My_Progs/VychMatyi/Praktikum/6_sem/6-th_lab/therm_cond_eq_experiment.txt", sep="\t")

    # Plotting:
RMeshgrid, TetaMeshgrid = np.meshgrid(r_grid, teta_grid)

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.contour3D(RMeshgrid, TetaMeshgrid, u_final, 75, cmap="viridis")
ax.contour3D(RMeshgrid, TetaMeshgrid, u_theory, 75, cmap="binary")
ax.set_title("u_exp")


    # ERRORS
u_error = np.absolute(u_theory - u_final)

df_error = pd.DataFrame(u_error[r_indexes, :][:, teta_indexes], index=r_grid[r_indexes], columns=np.degrees(teta_grid[teta_indexes]))
df_exp.to_csv("C:/My_Progs/VychMatyi/Praktikum/6_sem/6-th_lab/therm_cond_eq_errors.txt", sep="\t")

# print(u_error, "\n")
# print(u_theory, "\n")
u_rel_error = u_error[1 : L, 1 : M] / u_theory[1 : L, 1 : M]
# print(u_rel_error)

u_error_max = np.max(u_error)
u_rel_error_max = np.max(u_rel_error)
rel_err_max_index = np.unravel_index(np.argmax(u_rel_error), u_rel_error.shape)
r_rel_err_max = r_grid[rel_err_max_index[0] + 1]
teta_rel_err_max = teta_grid[rel_err_max_index[1] + 1]

f_error.write(f"Max error WHOLE GRID = {u_error_max}\n")
print(f"Max error WHOLE GRID = {u_error_max}")
print(f"Max error 5x5 = {np.max(u_error[r_indexes, :][:, teta_indexes])}")
f_error.write(f"Max relative error = {u_rel_error_max}\n")
print(f"Max relative error = {u_rel_error_max}")
f_error.write(f"Relative error index: r = {r_rel_err_max}, teta = {teta_rel_err_max}\n")
print(f"Relative error index: r = {r_rel_err_max}, teta = {teta_rel_err_max}")

print(f"u_error_max_rel = {u_error[rel_err_max_index[0] + 1, rel_err_max_index[1] + 1]}")
print(f"u_theory_max_rel = {u_theory[rel_err_max_index[0] + 1, rel_err_max_index[1] + 1]}")

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.contour3D(RMeshgrid, TetaMeshgrid, u_error, 75, cmap="viridis")
ax.set_title("u_err")

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.contour3D(RMeshgrid[1 : M, 1 : L], TetaMeshgrid[1 : M, 1 : L], u_rel_error, 75, cmap="viridis")
ax.set_title("u_rel_err")

"""    # Printing results in point
degree_out = 54
r_out = 0.8
exp_out = df_exp[degree_out][r_out]
theory_out = df_theory[degree_out][r_out]
print(f"exp = {exp_out}")
print(f"theory = {theory_out}")
print(f"Error = {np.abs(exp_out - theory_out)}")"""

f_error.close()
f_results.close()
plt.show()