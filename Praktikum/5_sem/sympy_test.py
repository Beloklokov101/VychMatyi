from sympy import *

x = symbols("x")

f = x**3 - 7 * x + 7
F = [f, diff(f, x)]
print(F)

flag, i = True, 0
while(flag):
#for i in range(K_deg_eq - 1):
    q, r = div(F[i], F[i + 1], domain='RR')
    F.append(- r)
    #print(r, "\n")
    if(diff(r, x) == 0): flag = False
    i = i + 1

print(*F, sep="\n")