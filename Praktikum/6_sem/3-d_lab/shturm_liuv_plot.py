import numpy as np
import math as mt
import matplotlib.pyplot as plt


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
L = 50000
lamb_set_num = 500
h = 1 / L

lamb_left_0 = 25
lamb_right_0 = 35
lamb = np.linspace(lamb_left_0, lamb_right_0, lamb_set_num + 1)

last_coeff = sht_liu_hard_solver(lamb.copy(), h, k_fun, L)

plt.plot(lamb, last_coeff, color="b")
less_then_0 = last_coeff < 0
plt.plot(lamb[less_then_0], last_coeff[less_then_0], color="r")
plt.show()