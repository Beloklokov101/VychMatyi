import numpy as np
import matplotlib.pyplot as plt


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



def get_ro(P):
    ro_0 = 1000
    c_f = 1e-4  # [atm]^-1
    p_0 = 120
    ro = ro_0 * (1 + c_f * (P / 101325 - p_0))
    return ro



L = 500
T = 10 * 24 * 60 * 60
count_x = 100
count_t = 10 * 24
h = L / count_x
tau = T / count_t
timestamps_out_days = np.array([0.1, 0.25, 0.5, 1, 1.5, 2, 3, 10])
timestamps_out = timestamps_out_days * 24 * 60 * 60
timestamps_out_size = timestamps_out.shape[0]
x_grid = np.linspace(0, L, count_x + 1)
t_grid = np.linspace(0, T, count_t + 1)


k = 1e-14
mu = 1e-3
phi = 0.2
c_f = 1e-4 / 101325
ro_0 = 1000
P_left = 150 * 101325
P_right = 50 * 101325

P_0_number = 100 * 101325
P_0 = P_0_number * np.ones(count_x + 1)
P_0[0] = P_left
P_0[-1] = P_right

P_next = P_0.copy()
P_result_timestamps = np.zeros((timestamps_out_size, count_x + 1))
timestamps_count = 0
for n in range(count_t):
    P_open = P_next.copy()
    P_next = np.zeros(count_x + 1)
    
    a, b = np.zeros(count_x + 1), np.zeros(count_x + 1)
    c, d = np.zeros(count_x + 1), np.zeros(count_x + 1)

    # a[i] * P[i + 1] + b[i] * P[i] + c[i] * P[i - 1] = d[i], i = range(1, count_x - 1)
    # a[count_x] = c[0] = 0
    # a[0] = c[count_x] = 0
    b[0], b[count_x] = 1, 1
    d[0], d[count_x] = P_left, P_right
    
    for it in range(1, count_x):
        p_it = P_open[it]
        p_it_next = P_open[it + 1]
        p_it_prev = P_open[it - 1]

        if p_it >= p_it_next:
            ro_half_plus = get_ro(p_it)
        else:
            ro_half_plus = get_ro(p_it_next)
        
        if p_it_prev >= p_it:
            ro_half_minus = get_ro(p_it_prev)
        else:
            ro_half_minus = get_ro(p_it)
        
        multip_1 = k / (mu * h**2)
        multip_2 = phi * c_f * ro_0 / tau
        a[it] = multip_1 * ro_half_plus
        c[it] = multip_1 * ro_half_minus
        b[it] = (- c[it] - a[it] - multip_2)
        d[it] = - multip_2 * p_it

    P_next = TriagonalSlae_solver(a, b, c, d, count_x + 1)

    if t_grid[n] >= timestamps_out[timestamps_count]:
        P_result_timestamps[timestamps_count] = P_next.copy()
        timestamps_count += 1
        print(f"timestamps_count = {timestamps_count}")

    if n % 10 == 0: print(f"n = {n}")

if timestamps_count != timestamps_out_size:
    P_result_timestamps[timestamps_count] = P_next.copy()
    timestamps_count += 1
    print(f"timestamps_count = {timestamps_count}")

fig, ax = plt.subplots(4, 2, sharex="all", sharey="all")

for it in range(8):
    row = it // 2
    col = it % 2
    ax[row, col].plot(x_grid, P_result_timestamps[it], label=f"P for t = {timestamps_out_days[it]} days")
    ax[row, col].legend()

plt.show()