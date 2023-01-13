import numpy as np
import math
from sympy import *

def sign_change(z):
    z_new = np.zeros(len(z))
    for i in range(len(z) - 1):
        z_new[i] = z[i + 1] - z[i]
    z_new = (z_new / 2 ) ** 2
    return np.sum(z_new)


def Shturm_method(F, z_min, z_max):
    F1 = diff(F, x)
    Functions = [F, F1]

    flag, i = True, 0
    while(flag):
        q, r = div(Functions[i], Functions[i + 1], domain='RR')
        Functions.append(- r)
        if(diff(r, x) == 0): flag = False
        i = i + 1

    functions = [lambdify(x, i, "numpy") for i in Functions]

    Z = [[z_min, z_max]]
    next_Z = []
    localized_roots = []
    while(Z != []):
        for seg in Z:
            seg_min, seg_max, length = seg[0], seg[1], seg[1] - seg[0]
            seg_value = [[f(z) for f in functions] for z in seg]
            seg_val_mod = np.sign(seg_value)
            k = int(abs(sign_change(seg_val_mod[0]) - sign_change(seg_val_mod[1])))
            if k != 0 and k != 1:
                for i in range(k):
                    next_Z.append([seg_min + i * length / k, seg_min + (i + 1) * length / k])
            if k == 1:
                localized_roots.append(seg)
        Z = next_Z
        next_Z = []

    return localized_roots


def half_div_method(F, localized_roots, eps_0):
    func = lambdify(x, F, modules="numpy")
    roots = []
    for seg in localized_roots:
        eps = (seg[1] - seg[0]) / 2
        while(eps > eps_0):
            seg_min, seg_max, length = seg[0], seg[1], seg[1] - seg[0]
            eps = (seg_max - seg_min) / 2
            medium = (seg_max + seg_min) / 2
            if func(seg_min) * func(medium) < 0:
                seg = [seg_min, medium]
            else:
                seg = [medium, seg_max]
        roots.append(seg)
    return roots


#First data
g_0 = 5/3
g_3 = 7/5
ro_0 = 1.694e-4
P_0 = 1.013e6
#P_3 = 1.6768e6
P_3 = 1e8
C_3 = 3.6537e4
U_0, U_3 = 0, 0

#Second data
ro_3 = g_3 * P_3 / (C_3)**2
alpha_0 = (g_0 + 1) / (g_0 - 1)
n = 2 * g_3 / (g_3 - 1)
#Only in my case mu = 0
mu = 0
nu = 2 / (g_3 - 1) * math.sqrt(g_3 * (g_0 - 1) * P_3 * ro_0 / (2 * P_0 * ro_3))
X = P_3 / P_0

#Equation coefficients
K_deg_eq = int(2 * n)
eq_coeff = [0 for i in range(K_deg_eq + 1)]
eq_coeff[0] = X**2
eq_coeff[K_deg_eq - int(n + 2)] = - alpha_0 * nu**2 * X
eq_coeff[K_deg_eq - int(n + 1)] = 2 * alpha_0 * nu * (mu + nu) * X
eq_coeff[K_deg_eq - int(n)] = - (2 + (mu + nu)**2 * alpha_0) * X
eq_coeff[K_deg_eq - 2] = - nu**2
eq_coeff[K_deg_eq - 1] = 2 * nu * (mu + nu)
eq_coeff[K_deg_eq] = 1 - (mu + nu)**2

#Circle of roots
B = max([eq_coeff[i] for i in range(K_deg_eq)])
A = max([eq_coeff[i] for i in range(1, K_deg_eq + 1)])
z_min = abs(eq_coeff[K_deg_eq]) / (abs(eq_coeff[K_deg_eq]) + B)
z_max = 1 + A / abs(eq_coeff[0])

Eq_coeff = np.array(eq_coeff)

#Sympy function
x = symbols("x")
F = 0
for i in range(K_deg_eq + 1):
    F = F + Eq_coeff[i] * x ** (K_deg_eq - i)


localized_roots = Shturm_method(F, z_min, z_max)

loc_prec_roots = half_div_method(F, localized_roots, 1e-8)
print(loc_prec_roots)

roots = [(loc_prec_roots[i][0] + loc_prec_roots[i][1]) / 2 for i in range(len(loc_prec_roots))]
Z = np.array(roots)

P_1 = P_3 * Z ** n
U_1 = U_3 + 2 * C_3 * (1 - Z) / (g_3 - 1)
#print(P_1, U_1)

b = np.array(1 - U_1)
c = np.array((P_1 - P_0) / ro_0)
#print(b, c)

Dis = b ** 2 - 4 * c
#print(Dis)
good_bad = [d >= 0 for d in Dis]
#print(good_bad)
D_0 = []
for i in range(len(good_bad)):
    if (good_bad[i]):
        D_0.append([(-b[i] + (-1) ** n * math.sqrt(Dis[i])) / 2 for n in range(2)])

print(D_0)