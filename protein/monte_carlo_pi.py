import numpy as np
import math as mt


N = 10000
M = 10000
pi_set = np.zeros(M)
for i in range(M):
    rand_set = np.random.uniform(-1, 1, size=(2, N))

    dist = np.sqrt(rand_set[0]**2 + rand_set[1]**2)

    dots_in_circle = np.sum(dist <= 1)

    pi_exp = 4 * dots_in_circle / N
    pi_set[i] = pi_exp

pi_average = np.mean(pi_set)
# print(pi_set)

f = open("C:/My_Progs/VychMatyi/protein/monte_carlo_datas.txt", "a")
f.write(f"N = {N}, M = {M}:\npi_average = {pi_average}\n\n")
print(f"pi_average = {pi_average}")
print(f"difference = {abs(mt.pi - pi_average)}")