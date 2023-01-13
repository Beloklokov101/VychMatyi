import sympy as sym
import matplotlib.pyplot as plt
import numpy as np


def div_diff_set(X, f):
    F = []
    n = len(X)
    F_curr = [(f[i + 1] - f[i]) / (X[i + 1] - X[i]) for i in range(n - 1)]
    F.append(F_curr)
    for i in range(n - 2):
        F_curr = [(F[i][j + 1] - F[i][j]) / (X[(j + 1) + (i + 1)] - X[j]) for j in range(len(F[i]) - 1)]
        F.append(F_curr)
    return F


def splyne_val(x_grid, X, S):
    f_grid = np.zeros(len(x_grid))
    j = 0
    for x in x_grid:
        i = 1
        while(x >= X[i] and i < N): 
            #print(x, X[i], i)
            i = i + 1
        f_grid[j] = S[i - 1](x)
        j = j + 1
    return f_grid


x = sym.symbols("x")
#F = x
#f_func = lambdify(x, F, modules="numpy")

"""
X = [0.52360,
    0.87267,
    1.22173,
    1.57080,
    1.91986,
    2.26893,
    2.61799]

F_set = [4e-5,
        0.00068,
        0.00518,
        0.02554,
        0.09624,
        0.30046,
        0.81548]
"""

X = [-3,
    -2, 
    -1,
    0,
    1,
    2,
    3]

F_set = [-120,
        -30,
        -4,
        0,
        0,
        -10,
        -60]

#Building Newton's interpolant
F = div_diff_set(X, F_set)
F_0 = [F[i][0] for i in range(len(F))]
#print(F_0)

P_n = F_set[0]
P_curr = 1
for i in range(len(F_0)):
    P_curr = 1
    for j in range(i + 1):
        P_curr = P_curr * (x - X[j])
    P_curr = P_curr * F_0[i]
    P_n = P_n + P_curr

#print(P_n)

"""
for i in range(10):
    print(P(X[i]) - F_set[i])

plt.plot(X, F_set)

x_x = np.linspace(X[0], X[-1], 100)
p_p = P(x_x)
plt.plot(x_x, p_p)

plt.show()
"""

P = sym.lambdify(x, P_n, modules="numpy")
print(P(1.95))
P_diff = sym.lambdify(x, sym.diff(P_n), modules="numpy")
N = len(X) - 1
#print(N)
a = np.zeros((4, int(N)), dtype=float)
S = [0 for k in range(N)]
#print(a)

#Calculating coeff "a" and building splyne S
for i in range(N):
    xf, xs = X[i], X[i + 1]
    dx = xs - xf
    ff, fs = F_set[i], F_set[i + 1]
    df = fs - ff
    pdf = P_diff(xf)
    pds = P_diff(xs)

    a[3][i] = (pds*dx - 2*df + pdf*dx) / dx**3
    a[2][i] = (-pds*dx*(xs + 2*xf) + 3*df*(xs + xf) - pdf*dx*(xf + 2*xs)) / dx**3
    a[1][i] = (pds*xf*(2*xs + xf)*dx - 6*df*xf*xs + pdf*xs*(xs + 2*xf)*dx) / dx**3
    a[0][i] = (-pds*(xf**2)*xs*dx + fs*(xf**2)*(3*xs - xf) + ff*(xs**2)*(xs - 3*xf) - pdf*xf*(xs**2)*dx) / dx**3

    S_f = a[0][i] + a[1][i]*x + a[2][i]*x**2 + a[3][i]*x**3
    S[i] = sym.lambdify(x, S_f, modules="numpy")

#print(a)

#Sending coeff "a" to the file
f = open("C:/My_Progs/VychMatyi/Praktikum/2-nd_lab/coeff_splyne.txt", "w")
f.write("My variant is I.6\n\n")
for i in range(N):
    f.write("S[{}]: ".format(i) + '\n')
    for j in range(4):
        f.write("{:>5}: ".format("a[{}]".format(j)) + str(a[j][i]) + '\n')
    f.write('\n')
  

#Plotting all graphs
grid_shape = 100
x_grid = np.linspace(X[0], X[-1], num=grid_shape)
#x_grid = np.linspace(X[0], 5, num=grid_shape)
f_grid = splyne_val(x_grid, X, S)
print(splyne_val([1.95], X, S))


plt.plot(X, F_set, 'b.')

p_p = P(x_grid)
plt.plot(x_grid, p_p, 'g-.')

plt.plot(x_grid, f_grid, 'r-')

P_n_sympl = sym.simplify(P_n)
print(P_n_sympl)
all_coeff = sym.Poly(P_n_sympl, x).all_coeffs()[::-1]
f.write("Coeff of interpolation polynomial are (from x**0 to x**{}):\n{}\n\n".format(len(all_coeff) - 1, all_coeff))   

#print(np.max(abs(f_grid - p_p)))
#x_val = np.array([1, 1, 1.5, 2, 2.5])
#print("Value of splyne in point(s) x = {} is S(x) = {}".format(x_val, splyne_val(x_val, X, S)))
#f.write("Spline values in some points \nx = {} are \nS(x) = {}".format(x_val, splyne_val(x_val, X, S)))

f.close() 
plt.show()