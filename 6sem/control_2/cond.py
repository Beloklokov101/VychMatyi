import numpy as np
import numpy.linalg as lng


A = np.array([[5, 8, 13],
            [8, 13, 5],
            [13, 5, 8]])

print(f"A^-1:\n{lng.inv(A)}")
A_cond_1 = lng.cond(A, p=np.inf)
A_cond_2 = lng.cond(A, p=1)
A_cond_3 = lng.cond(A, p=2)

print(f"cond:\n1 = {A_cond_1}\n2 = {A_cond_2}\n3 = {A_cond_3}")