import numpy as np
import math as mt
import sympy as sym

def get_a(tau, n):
    return - (2 * tau * n + 3)


def get_psi(tau, n):
    t_n = tau * n
    t = sym.symbols("t")
    psi = sym.lambdify(t, sym.ln(1 + (1 + 3 * t + t ** 2) ** 2))

def convection_eq_next_step(u_prev, L, N, n):
    u_next = np.zeros(L + 1)
    tau = 1 / N
    h = 1 / L
    u_dt = - 2
    a_n = get_a(tau, n)

    u_dx = (- u_dt) / a_n


    return u_next


    #Initial datas
L = 1000
N = 10 * L
tau = 1 / N
h = 1 / L

x_grid_0 = np.linspace(0, 1, L + 1)
u_0 = np.ln(1 + x_grid_0 ** 2)
