import matplotlib.pyplot as plt
import numpy as np
import cmath as mt
import matplotlib.animation as animation


    #Initial datas
L = 10
T = 0.02

gamma = 5/3
gamma_1 = gamma - 1
u_left_0 = 0
ro_left_0 = 13
P_left_0 = 10 * 101325
u_right_0 = 0
ro_right_0 = 1.3
P_right_0 = 1 * 101325

e_left_0 = P_left_0 / ((gamma_1) * ro_left_0)
e_right_0 = P_right_0 / ((gamma_1) * ro_right_0)

    #Time and x parameters
tau = 1e-5
result_out_time = 0.015
NX = 100
x_grid = np.linspace(-L, L, NX + 1)
h = 2 * L / NX
sigma = tau / h
dt_print = 0.0005

    #Making values set
v_prev = np.zeros((NX, 3))
for i in range(NX):
    if x_grid[i] < 0:
        v_prev[i] = np.array([ro_left_0, 
                            u_left_0, 
                            e_left_0])
    else:
        v_prev[i] = np.array([ro_right_0, 
                            u_right_0, 
                            e_right_0])

w_prev = np.zeros((NX, 3))
for i in range(NX):
    if x_grid[i] < 0:
        w_prev[i] = np.array([ro_left_0, 
                            ro_left_0 * u_left_0, 
                            ro_left_0 * e_left_0])
    else:
        w_prev[i] = np.array([ro_right_0, 
                            ro_right_0 * u_right_0, 
                            ro_right_0 * e_right_0])

t_moment = 0
CFL_max = 0
t_print = 0
while t_moment <= T:
    w_open = np.zeros((NX, 3))
    v_open = np.zeros((NX, 3))
    CFL_moment = np.zeros(NX)

    for i in range(1, NX - 1):
        e_i = v_prev[i][2]
        u_i = v_prev[i][1]
        c__2 = gamma * (gamma_1) * e_i
        c = np.sqrt(c__2)

        omega_T = np.array([[-u_i * c, c, gamma_1], 
                        [-c__2, 0, gamma_1], 
                        [u_i * c, - c, gamma_1]])
        omega_inv_T = np.array([[1 / (2 * c__2), - 1 / c__2, 1 / (2 * c__2)],
                            [(u_i + c) / (2 * c__2), - u_i / c__2, (u_i - c) / (2 * c__2)],
                            [1 / (2 * (gamma_1)), 0, 1 / (2 * (gamma_1))]])
        lambda_abs = np.array([[abs(u_i + c), 0, 0],
                            [0, abs(u_i), 0],
                            [0, 0, abs(u_i - c)]])
        A = np.array([[0, 1, 0],
                    [ - u_i ** 2, 2 * u_i, gamma_1],
                     [- gamma * u_i * e_i, gamma * e_i, u_i]])
        
        first_summand = sigma * A @ (w_prev[i + 1] - w_prev[i - 1]) / 2
        second_summand = sigma * omega_inv_T @ lambda_abs @ omega_T @ (w_prev[i + 1] - 2 * w_prev[i] + w_prev[i - 1]) / 2
        w_open[i] = w_prev[i] - first_summand + second_summand
        v_open[i] = np.array([w_open[i][0],
                            w_open[i][1] / w_open[i][0],
                            w_open[i][2] / w_open[i][0]])
        
        CFL_moment[i] = sigma * np.max(np.diag(lambda_abs))

    if np.max(CFL_moment) > CFL_max:
        CFL_max = np.max(CFL_moment)
    
    w_open[0] = w_open[1]
    w_open[NX - 1] = w_open[NX - 2]
    v_open[0] = v_open[1]
    v_open[NX - 1] = v_open[NX - 2]

    if t_moment >= result_out_time - tau / 2 and t_moment <= result_out_time + tau / 2:
        v_out = v_open.copy()
        print("SUCCESS!!!", t_moment)
    
    w_prev = w_open.copy()
    v_prev = v_open.copy()

    if t_moment >= t_print:
        print(t_moment)
        print(CFL_max)
        t_print += dt_print

    t_moment += tau

    #Printing result
x_set = np.linspace(- L + h/2, L - h/2, NX)

fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True)

for i in range(3):
    axes[int(i // 2)][int(i % 2)].plot(x_set, v_out[:, i])
    axes[int(i // 2)][int(i % 2)].set_title(f"v_out[{i}]")

axes[1][1].plot(x_set, gamma_1 * v_out[:, 0] * v_out[:, 2])
axes[1][1].set_title("P")


plt.show()