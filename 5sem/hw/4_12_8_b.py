import math

x1_prev = 0
x1 = math.exp(x1_prev ** 2) / (2 * math.sqrt(2 * math.e))
n = 1
while(abs(x1 - x1_prev) > 1e-3 / 2):
    x1_prev = x1
    x1 = math.exp(x1_prev ** 2) / (2 * math.sqrt(2 * math.e))
    n = n + 1

x1_med = (x1 + x1_prev) / 2
print("First root: ", x1_med, n)

x2_prev = 1
x2 = math.sqrt(math.log(x2_prev * 2 * math.sqrt(2 * math.e)))
n = 1
while(abs(x2 - x2_prev) > 1e-3 / 2):
    x2_prev = x2
    x2 = math.sqrt(math.log(x2_prev * 2 * math.sqrt(2 * math.e)))
    n = n + 1

x2_med = (x2 + x2_prev) / 2
print("Second root: ", x2_med, n)

print("Length is: ", x2_med - x1_med)
