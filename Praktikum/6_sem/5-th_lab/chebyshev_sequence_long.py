import numpy as np
import math as mt

def chebyshev_seq_even(teta_N, N, overline):
    teta_2N = np.zeros(2 * N)
    for i in range(1, 2 * N + 1):
        if i % 2 != 0:
            teta_2N[i - 1] = teta_N[int((i + 1) / 2 - 1)]
        else:
            teta_2N[i - 1] = 4 * N + 2 * overline - teta_2N[i - 2]
    return teta_2N


def chebyshev_seq_not_even(teta_2N, N2, overline):
    N = int(N2 / 2)
    teta_2N_1 = np.zeros(2 * N + 1)
    teta_2N_1[ : 2 * N] = teta_2N[ : 2 * N]
    teta_2N_1[2 * N] = (4 - 2 * overline) * N + 1
    return teta_2N_1


def chebyshev_sequence(N):
    N_it = int(N)
    flag = True
    N_seq = []
    while(flag):
        if N_it % 2 == 0:
            if N_it & (N_it - 1) == 0:
                degree_2 = int(mt.log(N_it, 2))
                flag = False
            else:
                N_seq.append(0)
                N_it = int(N_it / 2)
        else:
            N_seq.append(1)
            N_it -= 1
    # print(N_seq, degree_2)

    N_seq = N_seq[::-1]
    print(f"N_seq = {N_seq}, {degree_2}")

    teta = np.array([1])
    N_check = 1
    if len(N_seq) < 2:
        for it in range(degree_2 - 1):
            teta = chebyshev_seq_even(teta, N_check, overline=False)
            N_check *= 2

        if len(N_seq) == 1:
            teta = chebyshev_seq_even(teta, N_check, overline=True)
            N_check *= 2
            teta = chebyshev_seq_not_even(teta, N_check, overline=True)
            N_check += 1    
        else:
            teta = chebyshev_seq_even(teta, N_check, overline=False)
            N_check *= 2
    else:
        for it in range(degree_2):
            teta = chebyshev_seq_even(teta, N_check, overline=False)
            N_check *= 2
        
        # print("degree", teta, N_check)
        flag = False
        teta = chebyshev_seq_not_even(teta, N_check, overline=False)
        N_check += 1
        it = 1

        while it < len(N_seq) - 2:
            print(f"it = {it}")
            if (N_seq[it] == 0 and N_seq[it - 1] == 1 and N_seq[it + 1] == 1):
                print("YES", teta, N_check)
                flag = True
                teta = chebyshev_seq_even(teta, N_check, overline=True)
                N_check *= 2
                teta = chebyshev_seq_not_even(teta, N_check, overline=True)
                N_check += 1  
                it += 2
            else:
                print("NO", teta, N_check)
                flag = False
                if (N_seq[it] == 0):
                    teta = chebyshev_seq_even(teta, N_check, overline=False)
                    N_check *= 2  
                else:
                    teta = chebyshev_seq_not_even(teta, N_check, overline=False)
                    N_check += 1
                it += 1
        
        print("AFTER", teta, N_check)
        print("it = ", it, flag)
        # it = len(N_seq) - 2
        if ((flag == False or it - flag == len(N_seq) - 3) and N_seq[it] == 0 and N_seq[it + 1] == 1):
            print("First case")
            teta = chebyshev_seq_even(teta, N_check, overline=True)
            N_check *= 2
            teta = chebyshev_seq_not_even(teta, N_check, overline=True)
            N_check += 1  
            it += 2
        else:
            if (flag and it == len(N_seq) - 1):
                teta = chebyshev_seq_even(teta, N_check, overline=False)
                N_check *= 2
                it += 1
            else:
                for it in range(len(N_seq) - 2, len(N_seq)):
                    print(f"it = {it}")
                    print(teta, N_check)
                    if (N_seq[it] == 0):
                        teta = chebyshev_seq_even(teta, N_check, overline=False)
                        N_check *= 2  
                        it += 1
                    else:
                        teta = chebyshev_seq_not_even(teta, N_check, overline=False)
                        N_check += 1
                        it += 1

    # print(N_check)
    if N_check == N: print(f"OK")
    return teta


print(chebyshev_sequence(45))