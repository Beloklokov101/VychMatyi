from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import math as mt

# k = 4
def ode_4(f, x0, y0, h):
    f1 = f(x0, y0)
    f2 = f(x0 + h/2, y0 + f1*h/2)
    f3 = f(x0 + h/2, y0 + h*f2/2)
    f4 = f(x0 + h, y0 + h*f3)
    y1 = y0 + h*(f1 + 2*f2 + 2*f3 + f4)/6
    x1 = x0 + h
    return (x1, y1)

# k = 3
def Heun_ns(f, x0, y0, h):
    f1 = f(x0, y0)
    f2 = f(x0 + h/3, y0 + f1*h/3)
    f3 = f(x0 + h*2/3, y0 + h*f2*2/3)
    y1 = y0 + h*(f1 + 3*f3)/4
    x1 = x0 + h
    return (x1, y1)


x = symbols("x")
y = symbols("y")

#Problem statement
expr = (2*x**3 + x**2 - y**2)/(2*x**2*y)
#expr = (y**2 - 5*x)/(2*x*y)
#expr = (y - x*y**2) / x
#x0, y0 = 1, 1
x0, y0 = 1, 0.02
a, b = 1, 2
length = 0.2
err = 1e-4
k = 4
dots_num = 11

f = lambdify([x, y], expr, "numpy")

#First run
n = dots_num - 1
h = (length)/n
#print(n, " ", h)
X1, Y1 = np.zeros(n + 1), np.zeros(n + 1)
X1[0], Y1[0] = x0, y0
for i in range(1, n + 1):
    X1[i], Y1[i] = ode_4(f, X1[i - 1], Y1[i - 1], h)
err_step = err * 2
X0, Y0 = X1, Y1

#Running until we get enough accuracy
while(err_step >= err):
#for j in range(1):
    X, Y = X1, Y1

    n = 2 * n
    #n = 10000
    h = (length)/n
    #print(n, " ", h)
    X1, Y1 = np.zeros(n + 1), np.zeros(n + 1)
    X1[0], Y1[0] = x0, y0
    for i in range(1, n + 1):
        X1[i], Y1[i] = ode_4(f, X1[i - 1], Y1[i - 1], h)
    step_lev = int(n / (dots_num - 1))
    err_set = [abs(Y[int(i*step_lev/2)] - Y1[int(i*step_lev)])/(1 - 2**(-k)) for i in range(dots_num)]
    #err_set = [abs(Y[int(i)] - Y1[int(i*step_lev)])/(1 - 2**(-k)) for i in range(dots_num)]
    
    err_step = max(err_set)

#Write information to the file
f = open("C:/My_Progs/VychMatyi/Praktikum/3-d_lab/R-K_methods.txt", "w")
f.write("Laboratory work number: 6\n")
f.write("Variant number: 1\n\n")

f.write(f"Was used method with accuracy degree k = {k}\n")
f.write("For {} segments vectors of values are X0, Y0\n".format(dots_num - 1))
f.write("For {} segments vectors of values are X_, Y_\n".format(int(n / 2)))
f.write("For {} segments vectors of values are X1, Y1\n\n".format(n))

f.write(f"{dots_num} equidistant dots with their values:\n\n")
for i in range(dots_num):
    f.write(f"{i}\n")
    f.write("{:>2}X0 = {:.3f}\n".format("", X0[i]))
    f.write("{:>2}Y0 = {}\n".format("", Y0[i]))
    #f.write("{:>2}Y_ = {}\n".format("", Y[int(i*step_lev/2)]))
    f.write("{:>2}Y1 = {}\n".format("", Y1[int(i*step_lev)]))
    f.write("{:>2}|| || = {}\n\n".format("", abs(Y0[i] - Y1[int(i*step_lev)])))
    
    """
    f.write("{:>2}X1[{}], Y1[{}] = {:.3f}, {:.6f}\n".format("", int(i*step_lev), int(i*step_lev), X1[int(i*step_lev)], Y1[int(i*step_lev)]))
    #f.write("{:>2}X[{}] , Y[{}]  = {:.3f}, {:.6f}\n".format("", int(i*step_lev/2), int(i*step_lev/2), X[int(i*step_lev/2)], Y[int(i*step_lev/2)]))
    f.write("{:>2}X[{}] , Y[{}]  = {:.3f}, {:.6f}\n".format("", int(i), int(i), X[int(i)], Y[int(i)]))
    
    #f.write("{:>2}diff = {:.2e}\n\n".format("", abs(Y[int(i*step_lev/2)] - Y1[int(i*step_lev)])))
    f.write("{:>2}diff = {:.2e}\n\n".format("", abs(Y[int(i)] - Y1[int(i*step_lev)])))
    """


#Plotting
#Result
plt.plot(X1, Y1, 'b.')
print("Grid shape: ", n)
print("Error: ", err_step)


#Theory
x_th = np.linspace(1, 1.2)
#res_th = x
res_th = sqrt(x**2 - 0.367732*mt.e**(1/x))
#res_th = sqrt(x*(1 - 5*log(x)))
#res_th = 2/x

f_th = lambdify(x, res_th, "numpy")
y_th = f_th(x_th)
plt.plot(x_th, y_th, 'r-.')

plt.show()
f.close()