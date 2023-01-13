import math as mt
from sympy import *
import numpy as np

x = symbols("x")
F = sin(1000 * x) / (1 + x)
a, b = 0, mt.pi / 2
f = lambdify(x, F, "numpy")

    #4-th derrivative

F_d4 = diff(F, x, 4)
#print(F_d4)
f_d4 = lambdify(x, F_d4, "numpy")

x = np.linspace(a, b, 10000)
M4 = np.max(f_d4(x))

    #Newton-Cortes 3/8

#h_nc = root(684 / M4, 5) / 10
#print(h_nc)
#print(M4 * (0.001**5) / 6480)
h_nc = 0.001
seg = mt.ceil((b - a) / h_nc)
print(seg)
h_nc = (b - a) / seg
print(f"new h_nc = {h_nc}")
print("e_nc = ", M4 * (h_nc**5) / 6480)
dots = 3 * seg + 1
#print(dots)
x_nc = np.linspace(a, b, dots)
f_nc = f(x_nc)

I_nc = 0
for i in range(seg):
    I_i = h_nc * (f_nc[3 * i] + 3 * f_nc[3 * i + 1] + 3 * f_nc[3 * i + 2] + f_nc[3 * i + 3]) / 8
    I_nc = I_nc + I_i

print("diff = ", abs(0.000611014 - I_nc))
print(f"Newton-Cortes 3/8: I_nc = {I_nc}\n")


    #Gauss with 2 points

#h_nc = root(432 / M4, 5) / 10
#print(h_nc)
h_g = 0.001
seg = mt.ceil((b - a) / h_g)

seg = 1

print(seg)
dots = 2 * seg
h_g = (b - a) / seg
print(f"new h_g = {h_g}")
print("e_g = ", M4 * (h_g**5) / 4320)

I_g = 0
for i in range(seg):
    a_i, b_i = a + i*h_g, a + (i + 1)*h_g
    dif = (b_i - a_i)/2
    med = (b_i + a_i)/2
    x_abs_kn = sqrt(mt.pi**2 - 3/(2*500**2))
    c_abs = mt.pi / (2 * 500 * x_abs_kn)
    x1_i = med - x_abs_kn * dif / mt.pi
    x2_i = med + x_abs_kn * dif / mt.pi
    print(x1_i, " ", x2_i)
    #x1_i = med + dif/mt.sqrt(3) 
    #x2_i = med - dif/mt.sqrt(3) 
    #print(x1_i)
    #ans = f(x1_i)
    I_i = (c_abs / (1 + float(x1_i)) - c_abs / (1 + float(x2_i)))
    I_g = I_g + I_i

print("diff = ", abs(0.000611014 - I_g))
print(f"Gauss: I_g = {I_g}")

f = open("C:/My_Progs/VychMatyi/Praktikum/5_sem/integral/integral.txt", "w")
f.write("Laboratory work number: 5\n")
f.write("Variant number: 7\n\n")

f.write("Were used two methods to calculate integral, Newton-Cotes 3/8 and Gauss with two nodes\n")
f.write(f"Step length are h = {h_nc} (NC), {h_g} (G)\n")
f.write(f"Inaccuracy e = 1e-6\n\n")

f.write(" Newton-Cotes\n")
f.write(f"N-C theoretical error e_nc = {M4 * (h_nc**5) / 6480}\n")
f.write(f"I_nc = {I_nc}\n\n")

#f.write(f"Gauss theoretical error e_g = {M4 * (h_g**5) / 4320}\n")
f.write(" GAUSS\n")
f.write("Function :: sin(1000x), n = 500\n")
f.write("Nodes coordinate are: \n")
f.write("{:>2}x1 = {}\n".format("", x1_i))
f.write("{:>2}x2 = {}\n".format("", x2_i))
f.write("Nodes wheights are: \n")
f.write("{:>2}c1,2 = +-{}\n".format("", c_abs))
f.write(f"I_g = {I_g}\n\n")

f.close()