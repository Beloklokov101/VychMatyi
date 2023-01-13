import numpy as np
import math as mt


def k_fun(x):
    return 2 - x
    # return 1


def coeff_abc(k_fun, lamb, h, l):
    sigma = lamb * h**2
    k_l_next = k_fun(h * (l + 1/2))
    k_l_prev = k_fun(h * (l - 1/2)) 
    k_l_open = k_fun(h * l)
    a = - 2 * k_l_next - 2 * k_l_open + sigma
    b = k_l_next + 4 * k_l_open + k_l_prev - 2 * sigma
    c = - 2 * k_l_open - 2 * k_l_prev + sigma
    return a, b, c


def sht_liu_hard_solver(lamb, h, k_fun, L):
    A = np.zeros((L - 2, lamb.shape[0]))
    B = np.zeros((L - 2, lamb.shape[0]))

    a1, b1, c1 = coeff_abc(k_fun, lamb, h, 1)
    divider_0 = - b1 - k_fun(h * 1/2)
    A[0] = a1 / divider_0
    B[0] = k_fun(h * 3 / 2) / divider_0

    a2, b2, c2 = coeff_abc(k_fun, lamb, h, 2)
    divider_1 = - b2 - c2 * A[0]
    A[1] = (a2 + c2 * B[0]) / divider_1
    B[1] = k_fun(h * 5 / 2) / divider_1

    for l in range(2, L - 3):
        al, bl, cl = coeff_abc(k_fun, lamb, h, l)
        divider_i = - bl - A[l - 1] * cl - k_fun(h * (l - 1/2)) * (A[l - 2] * A[l - 1] + B[l - 2])
        A[l] = (al + B[l - 1] * cl + k_fun(h * (l - 1/2)) * A[l - 2] * B[l - 1]) / divider_i
        B[l] = k_fun(h * (l + 1/2)) / divider_i
    
    aL_1, bL_1, cL_1= coeff_abc(k_fun, lamb, h, L - 1)
    A[L - 3] = - (bL_1 - k_fun(h * (L - 1/2)) + k_fun(h * (L - 3/2)) * B[L - 4]) / (cL_1 + k_fun(h * (L - 3/2)) * A[L - 4])

    aL_2, bL_2, cL_2= coeff_abc(k_fun, lamb, h, L - 2)
    last_coeff = aL_2 + A[L - 3] * bL_2 + B[L - 5] * A[L - 3] * k_fun(h * (L - 5/2)) + (A[L - 4] * A[L - 3] + B[L - 4]) * (cL_2 + A[L - 5] * k_fun(h * (L - 5/2)))
    return last_coeff


def sht_liu_hard(k_fun, L, h, prec, lamb):
        #First step
    last_coeff = sht_liu_hard_solver(lamb, h, k_fun, L)
    # print(lamb, "\n", last_coeff, "\n")
    if last_coeff[1] * last_coeff[0] <= 0:
            lamb[2] = lamb[1]
            lamb[1] = (lamb[0] + lamb[2]) / 2
            dist = lamb[2] - lamb[0]
    else:
        if last_coeff[1] * last_coeff[2] <= 0:
            lamb[0] = lamb[1]
            lamb[1] = (lamb[0] + lamb[2]) / 2
            dist = lamb[2] - lamb[0]
        else:
            dist = -1

    lamb_middle = np.array([lamb[1]])

        #All next steps
    while(dist > prec):
        last_coeff[1] = sht_liu_hard_solver(lamb_middle, h, k_fun, L)
        # print(lamb, "\n", last_coeff, "\n")

        if last_coeff[1] * last_coeff[0] <= 0:
            lamb[2] = lamb[1]
            lamb[1] = (lamb[0] + lamb[2]) / 2
            dist = lamb[2] - lamb[0]
        else:
            if last_coeff[1] * last_coeff[2] <= 0:
                lamb[0] = lamb[1]
                lamb[1] = (lamb[0] + lamb[2]) / 2
                dist = lamb[2] - lamb[0]
            else:
                dist = -1
        
        lamb_middle = np.array([lamb[1]])

    if dist >= 0:
        return True, lamb
    else:
        return False, 0


    #Initial datas
L = 10000
h = 1 / L
prec = 1e-4


lamb_left_0 = 20
lamb_right_0 = 25
lamb = np.array([lamb_left_0, (lamb_left_0 + lamb_right_0) / 2, lamb_right_0])
lamb_step = 3
lambda_result_num = 1

file = open("C:/My_Progs/VychMatyi/Praktikum/6_sem/3-d_lab/3-d_lab_hard.txt", "w")
file.write("Laboratory work number: 2\n")
file.write("Variant number: 2\n\n")
file.write(f"L = {L}  prec = {prec}\n\n")

print(f"L = {L}  prec = {prec}\n")

root_number = 0
while root_number < lambda_result_num:
# for it in range(1):
    flag, lamb_root = sht_liu_hard(k_fun, L, h, prec, lamb.copy())

    if flag:
        print(f"Lambda[{root_number}] is located in:\n{lamb_root[0]} | {lamb_root[2]}")
        print(f"Dist = {lamb_root[2] - lamb_root[0]}\n")

        file.write(f"Lambda[{root_number}] is located in:\n{lamb_root[0]} | {lamb_root[2]}\n")
        file.write(f"Dist = {lamb_root[2] - lamb_root[0]}\n\n")
        root_number += 1
    # else:
    #     print("There are no roots")
    
    lamb[0] = lamb[2]
    lamb[2] = lamb[0] + lamb_step
    lamb[1] = (lamb[0] + lamb[2]) / 2

file.close()