file_data = open("VychMatyi/Sasha_task/res.csv")
# print(file_data)
lines = file_data.readlines()
# print(len(lines))

sum_forward = 0
for row in lines:
    # print(row)
    sum_forward += float(row)

sum_reverse = 0
for row in lines[::-1]:
    sum_reverse += float(row)

print(f"Forward sum = {sum_forward}")
print(f"Reverse sum = {sum_reverse}")
file_data.close()