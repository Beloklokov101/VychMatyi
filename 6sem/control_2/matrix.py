import numpy as np
import numpy.linalg as lng

A = np.array([[3, 5],
              [2, 1]])
f = np.array([4, 1])
x_first = np.array([0, 0])
err = 1e-3

N_dim = A.shape[0]

    # Symmetrizing 
if np.array_equal(A, A.T):
    A_new = A
    f_new = f
else:
    A_new = A.T @ A
    f_new = A.T @ f

print(f"A_new = {A_new}")

lambda_A_new = lng.eigvals(A_new)
print(f"lambda_A_new = {lambda_A_new}")

print(f"l_max + l_min = {np.max(lambda_A_new) + np.min(lambda_A_new)}")
tau_opt = 2 / (np.max(lambda_A_new) + np.min(lambda_A_new))

B = np.eye(N_dim) - tau_opt * A_new
print(f"B = {B}")

lambda_B = lng.eigvals(B)
norm_B = lng.norm(B, ord=2)
print(f"lambda_B = {lambda_B}")

norm_A_new = lng.norm(A_new, ord=2)
r_0 = A_new @ x_first - f_new
norm_r_0 = lng.norm(r_0, ord=2)
print(f"norm_A_new(3) = {norm_A_new}")

N = np.abs(np.log(err * norm_A_new / norm_r_0) / np.log(norm_B))
print(f"N_opt = {N}")
