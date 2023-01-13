import numpy as np
import math as mt
import matplotlib.pyplot as plt

def get_mu_min_max(L, M, h_x, h_y):
    mu_min = 4 * ((mt.sin(mt.pi / (2 * L)) / h_x)**2 + (mt.sin(mt.pi / (2 * M)) / h_y)**2)
    mu_max = 4 * ((mt.sin(mt.pi * (L - 1) / (2 * L)) / h_x)**2 + (mt.sin(mt.pi * (M - 1) / (2 * M)) / h_y)**2)
    return mu_min, mu_max


def get_tau_set(mu_min, mu_max, N, cheb_seq):
    tau_arr = 2 / ((mu_min + mu_max) + (mu_max - mu_min) * np.cos(mt.pi * cheb_seq / (2 * N)))
    return np.array(tau_arr)


def chebyshev_seq_even(teta_N, N, overline):
    teta_2N = np.zeros(2 * N)
    for i in range(1, 2 * N + 1):
        if i % 2 != 0:
            teta_2N[i - 1] = teta_N[int((i + 1) / 2 - 1)]
        else:
            teta_2N[i - 1] = 4 * N + 2 * overline - teta_2N[i - 2]
    return teta_2N


def chebyshev_seq_not_even(teta_2N, N2, overline):
    N = int(N2 / 2)
    teta_2N_1 = np.zeros(2 * N + 1)
    teta_2N_1[ : 2 * N] = teta_2N[ : 2 * N]
    teta_2N_1[2 * N] = (4 - 2 * overline) * N + 1
    return teta_2N_1


def chebyshev_sequence(N):
    N_it = int(N)
    flag = True
    N_seq = []
    while(flag):
        if N_it % 2 == 0:
            if N_it & (N_it - 1) == 0:
                degree_2 = int(mt.log(N_it, 2))
                flag = False
            else:
                N_seq.append(0)
                N_it = int(N_it / 2)
        else:
            N_seq.append(1)
            N_it -= 1
    # print(N_seq, degree_2)

    N_seq = N_seq[::-1]
    # print(f"N_seq = {N_seq}, {degree_2}")

    teta = np.array([1])
    N_check = 1
    if len(N_seq) < 2:
        for it in range(degree_2 - 1):
            teta = chebyshev_seq_even(teta, N_check, overline=False)
            N_check *= 2

        if len(N_seq) == 1:
            teta = chebyshev_seq_even(teta, N_check, overline=True)
            N_check *= 2
            teta = chebyshev_seq_not_even(teta, N_check, overline=True)
            N_check += 1    
        else:
            teta = chebyshev_seq_even(teta, N_check, overline=False)
            N_check *= 2
    else:
        for it in range(degree_2):
            teta = chebyshev_seq_even(teta, N_check, overline=False)
            N_check *= 2
        
        # print("degree", teta, N_check)
        flag = False
        teta = chebyshev_seq_not_even(teta, N_check, overline=False)
        N_check += 1
        it = 1

        while it < len(N_seq) - 2:
            # print(f"it = {it}")
            if (N_seq[it] == 0 and N_seq[it - 1] == 1 and N_seq[it + 1] == 1):
                # print("YES", teta, N_check)
                flag = True
                teta = chebyshev_seq_even(teta, N_check, overline=True)
                N_check *= 2
                teta = chebyshev_seq_not_even(teta, N_check, overline=True)
                N_check += 1  
                it += 2
            else:
                # print("NO", teta, N_check)
                flag = False
                if (N_seq[it] == 0):
                    teta = chebyshev_seq_even(teta, N_check, overline=False)
                    N_check *= 2  
                else:
                    teta = chebyshev_seq_not_even(teta, N_check, overline=False)
                    N_check += 1
                it += 1
        
        # print("AFTER", teta, N_check)
        # print("it = ", it, flag)
        # it = len(N_seq) - 2
        if ((flag == False or it - flag == len(N_seq) - 3) and N_seq[it] == 0 and N_seq[it + 1] == 1):
            # print("First case")
            teta = chebyshev_seq_even(teta, N_check, overline=True)
            N_check *= 2
            teta = chebyshev_seq_not_even(teta, N_check, overline=True)
            N_check += 1  
            it += 2
        else:
            if (flag and it == len(N_seq) - 1):
                teta = chebyshev_seq_even(teta, N_check, overline=False)
                N_check *= 2
                it += 1
            else:
                for it in range(len(N_seq) - 2, len(N_seq)):
                    # print(f"it = {it}")
                    # print(teta, N_check)
                    if (N_seq[it] == 0):
                        teta = chebyshev_seq_even(teta, N_check, overline=False)
                        N_check *= 2  
                        it += 1
                    else:
                        teta = chebyshev_seq_not_even(teta, N_check, overline=False)
                        N_check += 1
                        it += 1
    
    if N_check == N: print(f"OK")
    return teta


L, M = 20, 20
x_max, y_max = 1, 1
h_x = x_max / L
h_y = y_max / M
err = 1e-3

mu_min, mu_max = get_mu_min_max(L, M, h_x, h_y)

    # min N
mmin_2, mmax_2 = mt.sqrt(mu_min), mt.sqrt(mu_max)
fract = (mmax_2 + mmin_2) / (mmax_2 - mmin_2)
N_min = mt.log(2 / err) / mt.log(fract)
# print(N_min)

N = int(N_min) + 1
print(f"N_theory = {N}")
N = 512
# N = int(N_min) + 50

    # tau set
cheb_seq = chebyshev_sequence(N)
tau_set = get_tau_set(mu_min, mu_max, N, cheb_seq)
print(f"N = {N}")
# print(f"tau_set = {tau_set}")


    # Solver
u_0 = np.zeros((L + 1, M + 1))
# u_open = u_0.copy()
u_next = u_0.copy()

it = 0
for tau in tau_set:
    u_open = u_next.copy()
    u_next = np.zeros((L + 1, M + 1))
    for l in range(1, L):
        for m in range(1, M):
            first_add = (u_open[l - 1, m] - 2 * u_open[l, m] + u_open[l + 1, m]) / h_x**2
            second_add = (u_open[l, m - 1] - 2 * u_open[l, m] + u_open[l, m + 1]) / h_y**2
            third_add = 2 * (l * h_x * (1 - l * h_x) + m * h_y * (1 - m * h_y))
            u_next[l, m] = u_open[l, m] + tau * (first_add + second_add + third_add)
    if it % 50 == 0: print(f"it = {it}")
    it += 1



# f = open("C:/My_Progs/VychMatyi/Praktikum/6_sem/5-th_lab/puasson.txt", "w")
# f.write("Laboratory work number: 5\n")
# f.write("Variant number: 1\n\n")

u_final = u_next.copy()
np.savetxt('C:/My_Progs/VychMatyi/Praktikum/6_sem/5-th_lab/puasson_exp.txt', u_final)
    # format: u[x, y]

# f.write(u_final)


x_grid = np.linspace(0, x_max, L + 1)
y_grid = np.linspace(0, y_max, M + 1)

# plt.subplot(2, 1, 1)
# plt.plot(x_grid, u_final[:, int(M / 2)])
# plt.title(f"y = {y_max / 2}")

# plt.subplot(2, 1, 2)
# plt.plot(y_grid, u_final[int(L / 2), :])
# plt.title(f"x = {x_max / 2}")

X_meshgrid, Y_meshgrid = np.meshgrid(x_grid, y_grid)

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.contour3D(X_meshgrid, Y_meshgrid, u_final, 50, cmap="viridis")
ax.set_title("u_exp")

u_theory = X_meshgrid * (1 - X_meshgrid) * Y_meshgrid * (1 - Y_meshgrid)
np.savetxt('C:/My_Progs/VychMatyi/Praktikum/6_sem/5-th_lab/puasson_theory.txt', u_theory)
err_field = abs(u_theory - u_final)
np.savetxt('C:/My_Progs/VychMatyi/Praktikum/6_sem/5-th_lab/puasson_err.txt', err_field)
print("Max error index = ", np.unravel_index(np.argmax(err_field), err_field.shape))
print("Max error = ", np.max(err_field))


fig = plt.figure()
ax = plt.axes(projection="3d")
ax.contour3D(X_meshgrid, Y_meshgrid, err_field, 50, cmap="viridis")
ax.set_title("err")


err_relative = err_field[1 : -1, 1 : -1] / u_theory[1 : -1, 1 : -1]
print("Max rel_error index = ", np.unravel_index(np.argmax(err_relative), err_relative.shape))
print("Max rel_error = ", np.max(err_relative))

# f.close()
plt.show()