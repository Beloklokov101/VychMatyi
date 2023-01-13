import matplotlib.pyplot as plt
import numpy as np
import cmath as mt
import matplotlib.animation as animation

    #Initial datas
a = 1   #Velocity
T = 18
L = 20
h = 0.5
CFL = 1.01
tau = h * CFL
t_grid_num = int(T / tau)
x_set_num = int(L / h)
dt_out_results = 1
results_out_num = int(T / dt_out_results) + 1 + int(np.sign(T % dt_out_results))
x_grid = np.linspace(0, L, x_set_num + 1)
u_0 = np.sin(4 * mt.pi * x_grid / L)

    #Scheme corner

# plt.figure()
# plt.plot(x_grid, u_0)

u_open = u_0.copy()
u_next = u_0.copy()
U_answer_corn = np.zeros((results_out_num, x_set_num + 1))
U_next_it = 0

time_next_output = 0
t_moment = 0
while t_moment <= T + tau / 2:
    u_open = u_next.copy()
    for i in range(1, x_set_num + 1):
        u_next[i] = u_open[i] - a * CFL * (u_open[i] - u_open[i - 1])
    u_next[0] = u_next[x_set_num]

    if t_moment >= time_next_output or t_moment > T - tau / 2:
        U_answer_corn[U_next_it] = u_open.copy()
        U_next_it += 1
        time_next_output += dt_out_results
        print(t_moment)
        # plt.plot(x_grid, u_0)
    t_moment += tau
    
# plt.figure()
# for i in range(results_out_num):
#     plt.plot(x_grid, U_answer_corn[i])

#   Scheme Lax-Wendroff

u_open = u_0.copy()
u_next = u_0.copy()
U_answer_LW = np.zeros((results_out_num, x_set_num + 1))
U_next_it = 0

time_next_output = 0
t_moment = 0
while t_moment <= T + tau / 2:
    u_open = u_next.copy()
    for i in range(x_set_num):
        u_next[i] = u_open[i] - a * CFL * (u_open[i + 1] - u_open[i - 1]) / 2 + a**2 * CFL**2 * (u_open[i + 1] - 2 * u_open[i] + u_open[i - 1]) / 2
    u_next[x_set_num] = u_open[x_set_num] - a * CFL * (u_open[0] - u_open[x_set_num - 1]) / 2 + a**2 * CFL**2 * (u_open[0] - 2 * u_open[x_set_num] + u_open[x_set_num - 1]) / 2

    if t_moment >= time_next_output or t_moment > T - tau / 2:
        U_answer_LW[U_next_it] = u_open.copy()
        U_next_it += 1
        time_next_output += dt_out_results
        print(t_moment)
    t_moment += tau


    #Plotting

fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)

def animate_corner(i):
    ax[0].clear()
    ax[0].set_title(f"Corner scheme, CFL = {CFL}")
    line = ax[0].plot(x_grid, U_answer_corn[i])
    return line

def animate_LW(i):
    ax[1].clear()
    ax[1].set_title(f"LW scheme, CFL = {CFL}")
    line = ax[1].plot(x_grid, U_answer_LW[i])
    return line

interval_animation = 50
repeat_animation = True
corner_animation = animation.FuncAnimation(fig, 
                                      animate_corner, 
                                      frames = np.linspace(0, results_out_num - 1, results_out_num, dtype=int),
                                      interval = interval_animation,
                                      repeat = repeat_animation)

LW_animation = animation.FuncAnimation(fig, 
                                      animate_LW, 
                                      frames = np.linspace(0, results_out_num - 1, results_out_num, dtype=int),
                                      interval = interval_animation,
                                      repeat = repeat_animation)


# plt.figure()
# for i in range(results_out_num):
#     plt.plot(x_grid, U_answer_LW[i])
# plt.plot(x_grid, U_answer_LW[0])
# # plt.plot(x_grid, u_0)

plt.show()