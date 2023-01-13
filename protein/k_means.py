import numpy as np
import matplotlib.pyplot as plt


def get_normal_distr(N, center, distance, size):
    X, Y = np.zeros(np.sum(size)), np.zeros(np.sum(size))
    j = 0
    for it in range(N):
        center_it = center[it]
        distance_it = distance[it]
        size_it = size[it]

        # X[j : j + size_it] = np.random.uniform(center_it[0] - distance_it / 2, center_it[0] + distance_it / 2, size=size_it)
        X[j : j + size_it] = np.random.normal(center_it[0], distance_it, size=size_it)
        Y[j : j + size_it] = np.random.normal(center_it[1], distance_it, size=size_it)
        j += size_it

    return X, Y


N = 3

centers = np.random.uniform(- 100, 100, size=(N, 2))
# center1 = (50, 60)
# center2 = (100, 20)
size_set = np.random.randint(40, 60, size=N)
size_sum = np.sum(size_set)
distance = np.random.uniform(0, 10, size=N)

X, Y = get_normal_distr(N, centers, distance, size_set)

# print(X, Y)
# plt.scatter(X, Y)
# plt.show()

X_min, X_max = np.min(X), np.max(X)
Y_min, Y_max = np.min(Y), np.max(Y)

clust_centers_x = np.random.uniform(X_min, X_max, size=N)
clust_centers_y = np.random.uniform(Y_min, Y_max, size=N)

clust_centers = np.column_stack((clust_centers_x, clust_centers_y))

clust_disp_max = 1
while clust_disp_max > 0.1:
    dist_all = np.zeros((size_sum, N))
    for n in range(N):
        X_n = X - clust_centers[n, 0]
        Y_n = Y - clust_centers[n, 1]
        dist_n = np.sqrt(X_n**2 + Y_n**2)
        dist_all[:, n] = dist_n
    
    dist_min_indeces = np.argmin(dist_all, axis=1)

    clust_displace = np.zeros(N)
    for n in range(N):
        dots_bool_n = dist_min_indeces == n
        if np.sum(dots_bool_n) != 0:            
            x_mass_center = np.sum(X[dots_bool_n]) / np.sum(dots_bool_n)
            y_mass_center = np.sum(Y[dots_bool_n]) / np.sum(dots_bool_n)

            clust_displace[n] = np.sqrt((clust_centers[n, 0] - x_mass_center)**2 + (clust_centers[n, 1] - y_mass_center)**2)            
            clust_centers[n] = np.array([x_mass_center, y_mass_center])
        else:
            clust_displace[n] = 0
    clust_disp_max = np.max(clust_displace)
    print(clust_disp_max)


plt.scatter(X, Y)
print(clust_centers)
plt.scatter(clust_centers[:, 0], clust_centers[:, 1], c="r")
plt.show()