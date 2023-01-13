import numpy as np
import math as mt


def sht_liu_model_solver(lamb, h, alpha, L):
    A = np.zeros((L - 2, lamb.shape[0]))
    B = np.zeros((L - 2, lamb.shape[0]))
    sigma = lamb * h**2 / alpha
    
    divider_0 = 2 * sigma - 7
    A[0] = (sigma - 4) / divider_0
    B[0] = 1 / divider_0
    divider_1 = 2 * sigma - 6 - A[0] * (sigma - 4)
    A[1] = (sigma - 4 + B[0] * (sigma - 4)) / divider_1
    B[1] = 1 / divider_1

    for l in range(2, L - 3):
        divider_i = A[l - 1] * (4 - A[l - 2]) - 6 - B[l - 2] - (A[l - 1] - 2) * sigma
        A[l] = (B[l - 1] * (A[l - 2] - 4) - 4 + (1 + B[l - 1]) * sigma) / divider_i
        B[l] = 1 / divider_i
    
    A[L - 3] = (5 - 2 * sigma + B[L - 4]) / (4 - A[L - 4] - sigma)
    last_coeff = sigma - 4 + A[L - 3] * (6 + B[L - 5] - 2 * sigma) + (A[L - 4] * A[L - 3] + B[L - 4]) * (A[L - 5] - 4 + sigma)
    return last_coeff


def sht_liu_model(alpha, L, h, prec, lamb):
        #First step
    last_coeff = sht_liu_model_solver(lamb, h, alpha, L)
    if last_coeff[1] * last_coeff[0] <= 0:
            lamb[2] = lamb[1]
            lamb[1] = (lamb[0] + lamb[2]) / 2
            dist = np.sqrt(lamb / alpha)[2] - np.sqrt(lamb / alpha)[0]
    else:
        if last_coeff[1] * last_coeff[2] <= 0:
            lamb[0] = lamb[1]
            lamb[1] = (lamb[0] + lamb[2]) / 2
            dist = np.sqrt(lamb / alpha)[2] - np.sqrt(lamb / alpha)[0]
        else:
            dist = -1

    lamb_middle = np.array([lamb[1]])

        #All next steps
    while(dist > prec):
        last_coeff[1] = sht_liu_model_solver(lamb_middle, h, alpha, L)

        if last_coeff[1] * last_coeff[0] <= 0:
            lamb[2] = lamb[1]
            lamb[1] = (lamb[0] + lamb[2]) / 2
            dist = np.sqrt(lamb / alpha)[2] - np.sqrt(lamb / alpha)[0]
        else:
            if last_coeff[1] * last_coeff[2] <= 0:
                lamb[0] = lamb[1]
                lamb[1] = (lamb[0] + lamb[2]) / 2
                dist = np.sqrt(lamb / alpha)[2] - np.sqrt(lamb / alpha)[0]
            else:
                dist = -1
        
        lamb_middle = np.array([lamb[1]])
    
    if dist >= 0:
        return True, lamb
    else:
        return False, 0


    #initial datas
alpha = 1
L = 1000
h = 1 / L
prec = 1e-5

k_left_border = mt.pi * 2 / 2
k_right_border = mt.pi * 3   / 2
lamb_left_0 = alpha * k_left_border ** 2
lamb_right_0 = alpha * k_right_border ** 2
lamb = np.array([lamb_left_0, (lamb_left_0 + lamb_right_0) / 2, lamb_right_0])
lamb_step = 5
lambda_result_num = 4

file = open("C:/My_Progs/VychMatyi/Praktikum/6_sem/3-d_lab/3-d_lab_model.txt", "w")
file.write("Laboratory work number: 2\n")
file.write("Variant number: 2\n\n")
file.write(f"L = {L}  prec = {prec}\n\n")

print(f"L = {L}  prec = {prec}\n")

root_number = 0
while root_number < lambda_result_num:
# for it in range(1):
    flag, lamb_root = sht_liu_model(alpha, L, h, prec, lamb.copy())

    if flag:
        k_range = np.sqrt(lamb_root / alpha)
        print(f"Lambda[{root_number}] is located in\n{lamb_root[0]} | {lamb_root[2]}")
        print(f"K[{root_number}] is located in:\n{k_range[0]} | {k_range[2]}")
        print(f"Dist = {k_range[2] - k_range[0]}\n")

        file.write(f"Lambda[{root_number}] is located in\n{lamb_root[0]} | {lamb_root[2]}\n")
        file.write(f"K[{root_number}] is located in:\n{k_range[0]} | {k_range[2]}\n")
        file.write(f"Dist = {k_range[2] - k_range[0]}\n\n")

        root_number += 1
    # else:
        # print("There are no roots\n")

    
    lamb[0] = lamb[2]
    lamb[2] = lamb[0] + lamb_step
    lamb[1] = (lamb[0] + lamb[2]) / 2
