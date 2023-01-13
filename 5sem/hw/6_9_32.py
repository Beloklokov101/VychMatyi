import matplotlib.pyplot as plt
import sympy as sym
import numpy as np


def solveTriagonalSlae(A, h, n):
    P = np.zeros((n,1))
    Q = np.zeros((n,1))
    x = np.zeros((n,1))
    
    P[0] = -A[0, 1] / A[0, 0]
    Q[0] = h[0] / A[0, 0]
    
    for i in range(1,n):
        a, b = A[i,i-1], A[i,i]
        if i != n-1:
            c = A[i,i+1]
            P[i] = c/(-b - a*P[i-1])
        Q[i] = (-h[i]+a*Q[i-1])/(-b - a*P[i-1])
    x[n-1] = Q[n-1]
    
    for i in range(n-2,-1,-1):
        x[i] = P[i]*x[i+1]+Q[i]
    
    return x


class CubicSpline():
    def __init__(self,h,p,s):
        self.x = h
        self.a = p
        self.n = np.size(p)
        self.s = s
        
    def interpolate(self):
        n = self.n
        self.h = np.diff(self.x)
        self.A = np.zeros((self.n-2,self.n-2))
        self.f = np.zeros(self.n-2) 
        for i in range(1,self.n-1):
            a2 = self.a[i+1]
            h0, h1 = self.h[i], self.h[i-1]
            a1 = self.a[i]
            a0 = self.a[i-1]
            
            self.f[i-1] = ((a2 - a1)/h0 - (a1 - a0)/h1)/(h0+h1)*6
            for j in range(1,n-1):               
                if i == j:
                    self.A[i-1,j-1] = 2
                if j == i - 1:
                    self.A[i-1,j-1] = h0/(h1 + h0)
                if j == i + 1:
                    self.A[i-1,j-1] = h1/(h1 + h0)
        
        self.c = solveTriagonalSlae(self.A, self.f, self.n-2)
        self.c2 = np.zeros(n-1)
        self.c2[0:n-2] = self.c.transpose()[0]
        self.c = self.c2
        self.b, self.d = np.zeros(n-1), np.zeros(n-1)
        for i in range(1,n):
            c,h,c1 = self.c[i-1],self.h[i-1], self.c[i-2]
            a1 = self.a[i-1]
            a0 = self.a[i]
            if i == 1:
                self.b[0] = (2*c*h + 6*(a0 - a1)/h)/6
                self.d[0] = self.c[0]/h
            else:
                self.b[i-1] = (2*c*h + c1*h + 6*(a0 - a1)/h)/6
                self.d[i-1] = (c - c1)/h
        
        self.X = np.linspace(self.x[0], self.x[self.n-1], int((-self.x[0] + self.x[self.n-1])/self.s)+1)
        N = np.size(self.X)
        self.F = np.zeros(N)
        j = 0
        for i in range(n-1):
            while j < N and (self.X[j] <= self.x[i+1]) :
                
                k = self.X[j] - self.x[i+1]
                self.F[j] = self.a[i+1] + self.b[i]*k + 1/2*self.c[i]*k**2 + 1/6*self.d[i]*k**3
                j = j+1
            


def div_diff_func(X, f):
    F = []
    n = len(X)
    F_curr = [(f(X[i + 1]) - f(X[i])) / (X[i + 1] - X[i]) for i in range(n - 1)]
    F.append(F_curr)
    for i in range(n - 2):
        F_curr = [(F[i][j + 1] - F[i][j]) / (X[(j + 1) + (i + 1)] - X[j]) for j in range(len(F[i]) - 1)]
        F.append(F_curr)
    return F


def div_diff_set(X, f):
    F = []
    n = len(X)
    F_curr = [(f[i + 1] - f[i]) / (X[i + 1] - X[i]) for i in range(n - 1)]
    F.append(F_curr)
    for i in range(n - 2):
        F_curr = [(F[i][j + 1] - F[i][j]) / (X[(j + 1) + (i + 1)] - X[j]) for j in range(len(F[i]) - 1)]
        F.append(F_curr)
    return F


x = sym.symbols("x")
#F = x
#f_func = lambdify(x, F, modules="numpy")

X = [1910 + 10 * i for i in range(10)]
F_set = [92228496,
        106021537,
        123202624,
        132164569,
        151325798,
        179323175,
        203211926,
        226545805,
        248709873,
        281421906]

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


P = sym.lambdify(x, P_n, modules="numpy")
"""
for i in range(10):
    print(P(X[i]) - F_set[i])
"""
plt.plot(X, F_set)

x_x = np.linspace(X[0], 2010, 100)
p_p = P(x_x)
plt.plot(x_x, p_p)

print("By Newton, population of USA in 2010 =", int(P(2010)))

X_np = np.array(X)
F_np = np.array(F_set)
C = CubicSpline(X_np, F_np, 0.1)
C.interpolate()
#print(C.F)
plt.plot(C.X, C.F)

x_ext = np.linspace(X[8], 2020, 100)
a = C.a[-1]
b = C.b[-1]
c = C.c[-1]
d = C.d[-1]

x_for_ext = x_ext - X[9]
f_ext = sym.lambdify(x, a + b * x + (c * x ** 2) / 2 + (d * x ** 3) / 6)
y_ext = f_ext(x_for_ext)
plt.plot(x_ext, y_ext)

print("With spline, population of USA in 2010 =", f_ext(2010 - 2000))

plt.show()