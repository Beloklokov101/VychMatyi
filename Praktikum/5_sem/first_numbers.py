import math
import numpy as np


def next_element(previous_element, element):
    pi = float(format(math.pi, '.15f'))
    e = float(format(math.e, '.15f'))
    """
    print("element: ", element)
    print("previous_element", previous_element)
    print("12*element: ", 12*element)
    print("5*previous_element: ", 5 * previous_element)
    print("pi - ...: ", (pi - (12 * element) - (5 * previous_element)))
    print("(pi - ... )/ e: ", (pi - (12 * element) - (5 * previous_element)) / e, "\n\n")
    """
    return ((pi - (12 * element) - (5 * previous_element)) / e)


N = int(input())
pi = float(format(math.pi, '.15f'))
e = float(format(math.e, '.15f'))

y = []
y0 = pi / (17 + e)
y.append(y0)
y.append(y0)

for i in range(N):
    y.append(next_element(y[i], y[i+1]))

print()
print(*y, sep='\n')

#for i in y:
#    print('%.32f' % i, sep='\n')
#for i in y:
#    format(i, '.24')

a, b, c = e, 12, 5
A = np.zeros((N + 2, N + 2))
A[0, 0], A[1, 1] = 1, 1

for i in range(2, N + 2):
    A[i, i - 2] = c
    A[i, i - 1] = b
    A[i, i] = a

A_1 = np.linalg.inv(A)
A_1_norm = np.max(np.sum(np.abs(A_1), axis=1))
#E = A@A_1

#print(A_1)
print("\nA^-1 norm: ", A_1_norm)

y.clear()
#dy = y0 / 10**6
dy = 10**-17
y.append(y0 + dy)
y.append(y0 + dy)

for i in range(N):
    y.append(next_element(y[i], y[i+1]))

y_th = (-0.420861 * (-3.94873)**(N + 1) + 1.42086 * (-0.46582)**(N + 1)) * dy + y0

print()
print(*y, sep='\n')
#print("\nTheoretical last {0}: ".format(N + 1), y_th)