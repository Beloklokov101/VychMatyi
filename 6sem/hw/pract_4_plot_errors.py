import numpy as np
import matplotlib.pyplot as plt


# f = open("C:/My_Progs/VychMatyi/6sem/hw/pract_4_errors.txt", "r")
with open("C:/My_Progs/VychMatyi/6sem/hw/pract_4_errors.txt", "r") as f:
    file_data = f.readlines()

time_out = np.array([0.2, 0.4, 0.5, 0.6, 0.8, 1.0])
time_size = time_out.shape[0]

length = int(len(file_data) / 2)

L_set = np.zeros(length)
L_count = 0
errors_set = np.zeros((length, time_size))
errors_count = 0

for it in range(len(file_data)):
    line = file_data[it].split(" ")

    if it % 2 == 0:
        L_set[L_count] = float(line[-1])
        L_count += 1
    else:
        errors_set[errors_count] = np.array(line[:-1])
        errors_count += 1

fig, ax = plt.subplots(3, 2, sharex="all")
for it in range(time_size):
    row = it // 2
    col = it % 2
    ax[row, col].plot(L_set, - np.log(errors_set[:, it]), ".-")
    ax[row, col].set_title(f"time = {time_out[it]}")

plt.show()