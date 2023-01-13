import numpy as np
import numpy.linalg as lng

A = np.array([[18, -6, -7],
            [-6, 6, 0],
            [-7, 0, 6]])
f = np.array([-13, 6, 6])
N_dim = A.shape[0]

L, U, D = np.zeros(A.shape), np.zeros(A.shape), np.zeros(A.shape)

for i in range(N_dim):
    for j in range(N_dim):
        if i > j:
            L[i, j] = A[i, j]
        else:
            if i < j:
                U[i, j] = A[i, j]
            else:
                D[i, j] = A[i, j]

print(f"L:\n{L}\nU:\n{U}\nD:\n{D}")

    #Jacoby
P = - lng.inv(D) @ (L + U)
print(f"Jacoby_P:\n{P}")
print(f"Jacoby_f:\n{lng.inv(D) @ f}")
Jacoby_lambda = lng.eigvals(P)
print(f"Jacoby_lambda = {Jacoby_lambda}")

    #Zeidel
P = - lng.inv(L + D) @ U
print(f"Zeidel_P:\n{P}")
print(f"Zeidel_f:\n{lng.inv(L + D) @ f}")
Zeidel_lambda = lng.eigvals(P)
print(f"Zeidel_lambda = {Zeidel_lambda}")