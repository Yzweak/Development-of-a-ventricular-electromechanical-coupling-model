import numpy as np
from math import sqrt
def load_data():
    f0 = np.zeros((19944))
    file1 = open("./fiber_direction_before.txt",'r')

    data1 = file1.readlines()

    for i in range(len(data1)):
        tmp_list = data1[i].split('\n')
        f0[i] = float(tmp_list[0])

    f0 = np.reshape(f0,(1662,12))

    res = []
    for line in f0:
        x,y,z = 0,0,0
        for i in range(4):
            x = x + line[3*i]
            y = y + line[3*i + 1]
            z = y + line[3*i + 2]
        norm = sqrt(x**2 + y**2 + z**2)
        res.append((x/norm, y/norm, z/norm))
    return res



if __name__ == '__main__':
    f0 = load_data()
    with open('fiber_direction.txt','w') as f:
        for item in f0:
            f.write(str(item[0]) + ",")
            f.write(str(item[1]) + ",")
            f.write(str(item[1]) + "\n")

    print("Done")