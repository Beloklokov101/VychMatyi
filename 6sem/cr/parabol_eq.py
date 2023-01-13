import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
import math as mt

def plotting(U, x_grid):
    fig, ax = plt.subplots(1, 2, sharex=True)

    ax[0].plot(x_grid, U[:, 0])
    ax[0].set_title("u")

    ax[1].plot(x_grid, U[:, 1])
    ax[1].set_title("v")


T = 1e-2

N = 1    # Time
L = 100    # x
tau = T / N
h = 1 / L
x_grid = np.linspace(0, 1, L + 1)

m_Lamb = np.array([[1, 2], [2, 1]])

    # Initial datas on u and v
u_0 = 0
u_left = 1
u_right = 0

v_0 = 0
v_left = 0
v_right = 1

    # A, B, C
sigma = tau / h**2
m_A = - sigma * m_Lamb
m_B = 2 * sigma * m_Lamb
m_C = np.eye(2) - sigma * m_Lamb

    # L, U
# L_c = np.zeros((L - 1, 2, 2))
# L_c[1 : L - 1] = m_A
L_w = np.zeros((L - 1, 2, 2))
L_w[0] = m_B
U_q = np.zeros((L - 1, 2, 2))

for it in range(1, L - 1):
    U_q[it - 1] = LA.inv(L_w[it - 1]) @ m_C
    L_w[it] = m_B - m_A @ U_q[it - 1]

det_min = 1000
for w in L_w:
    det_w = abs(LA.det(w))
    if det_w < det_min: det_min = det_w

    # U_start   
U_0 = np.zeros((L + 1, 2))
U_0[:, 0] = u_0
U_0[:, 1] = v_0

U_left = np.array([u_left, v_left])
U_right = np.array([u_right, v_right])

# U_open = U_0.copy()
U_next = U_0.copy()

# plotting(U_0, x_grid)

for t in range(1, N + 1):
    U_open = U_next.copy()
    U_next = np.zeros((L + 1, 2))

    D = np.zeros((L - 1, 2))
    D[0] = U_open[1] - m_A @ U_left
    D[1 : L - 2] = U_open[2 : L - 1]
    D[L - 2] = U_open[L - 1] - m_C @ U_right

    g = np.zeros((L - 1, 2))
    g[0] = LA.inv(L_w[0]) @ D[0]

    for it in range(1, L - 1):
        g[it] = LA.inv(L_w[it]) @ (D[it] - m_A @ g[it - 1])
    
    U_next[0] = U_left
    U_next[L] = U_right
    U_next[L - 1] = g[L - 2]

    for it in range(L - 2, 0, -1):
        U_next[it] = g[it - 1] - U_q[it - 1] @ U_next[it + 1]

    if t % 50 == 0: print("time = ", t)

U_final = U_next.copy()

# fig = plt.figure()
plotting(U_final, x_grid)

plt.show()