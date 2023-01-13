#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
from sympy import *
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate


# In[8]:


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


# In[9]:


data = pd.read_csv("Atmosphere_1.csv", skiprows = 1, header = None)
h, p = data[0].to_numpy(), data[1].to_numpy()
data2 = pd.read_csv("Atmosphere_2.csv", skiprows = 1, header = None)
h2, p2 = data2[0].to_numpy(), data2[1].to_numpy()


# In[7]:


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
        for i in range(1,self.n-1): #здесь все ок
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
            


# In[10]:


#график интерполяции,плотность от высоты
C = CubicSpline(h,p,250)
C.interpolate()
plt.plot(C.X, C.F)


# In[14]:


#первые четыре значения. рыжая линия соединяет узлы, синяя - сплайн
C = CubicSpline(h[0:4],p[0:4],25)
C.interpolate()
plt.plot(C.X, C.F,h[0:4],p[0:4])


# In[25]:


#значение высоты, где плотность 2*10^(-8)
C = CubicSpline(h,p,20)
C.interpolate()
print(C.X[2102])


# In[228]:


#ошибка как отклонение от значений из Atmosphere_2
dp = abs(p2 - C.F)
plt.plot(h2[0:100],dp[0:100])

