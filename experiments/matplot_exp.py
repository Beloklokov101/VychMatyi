import numpy as np
import math as mt
import matplotlib.pyplot as plt
import matplotlib.animation as animation


a = 1   #Velocity
T = 18
L = 20
h = 0.5
CFL = 1.0
tau = h * CFL
t_grid_num = int(T / tau)
x_set_num = int(L / h)
dt_out_results = 1
results_out_num = int(T / dt_out_results) + 1 + int(np.sign(T % dt_out_results))
#print(results_out_num)
x_grid = np.linspace(0, L, x_set_num + 1)
u_0 = np.sin(4 * mt.pi * x_grid / L)

    #Scheme corner

# plt.figure()
# print(x_grid, u_0)
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
        U_answer_corn[U_next_it] = u_open
        U_next_it += 1
        time_next_output += dt_out_results
        #print(t_moment)
        # plt.plot(x_grid, u_0)
    t_moment += tau

fig = plt.figure(figsize=(5, 5), facecolor="pink")
ax = plt.subplot(frameon=False)

def animate(i):
    ax.clear()
    ax.text(1, 1, "DOLBAEB", color="b", fontsize=16)
    
    line = ax.plot(x_grid, U_answer_corn[i], color="w")
    return line


fun_animation = animation.FuncAnimation(fig, 
                                      animate, 
                                      frames = np.linspace(0, results_out_num - 1, results_out_num, dtype=int),
                                      interval = 20,
                                      repeat = True)

# fun_animation.save('animation.gif',
#                  writer='Pillow', 
#                  fps=30)
# print(U_answer_corn[[-1]])

plt.show()
