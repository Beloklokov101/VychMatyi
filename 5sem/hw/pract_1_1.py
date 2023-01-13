import math

def next_step(X_prev, n):
    X = [0, 0]
    det = 2 * (X_prev[0] + X_prev[1] / (math.cos(X_prev[0]) ** 2))
    #print(det)
    X[0] = X_prev[0] - (X_prev[0] ** 2 + X_prev[1] ** 2 - 1 - 2 * X_prev[1] * ( X_prev[1] - math.tan(X_prev[0]))) / det
    X[1] = X_prev[1] - ((X_prev[0] ** 2 + X_prev[1] ** 2 - 1) / (math.cos(X_prev[0]) ** 2) + 2 * X_prev[0] * ( X_prev[1] - math.tan(X_prev[0]))) / det
    return X, n + 1


X = []
n = 0
X_prev = [0.6, 0.7]
X, n = next_step(X_prev, n)
while(abs(X[0] - X_prev[0]) > 1e-6 or abs(X[1] - X_prev[1]) > 1e-6):
    X_prev = X
    X, n = next_step(X_prev, n)

print("Root 1 is: ", X)

n = 0
X_prev = [- 0.6, - 0.7]
X, n = next_step(X_prev, n)
while(abs(X[0] - X_prev[0]) > 1e-6 or abs(X[1] - X_prev[1]) > 1e-6):
    X_prev = X
    X, n = next_step(X_prev, n)

print("Root 2 is: ", X)
