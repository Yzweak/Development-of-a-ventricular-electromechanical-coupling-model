import numpy as np
from math import fabs,sqrt
def openreadtxt(file_name,mode):
    data = []
    file = open(file_name, 'r')  # 打开文件
    file_data = file.readlines()  # 读取所有行
    for row in file_data:
        tmp_list = row.split(',')  # 按‘，’切分每行的数据
        tmp_list[-1] = tmp_list[-1].replace('\n','') #去掉换行符
        if mode == 1:
            tmp_list = list(map(float, tmp_list))
        else:
            tmp_list = list(map(int, tmp_list))
        data.append(tmp_list)
    return data



if __name__ == '__main__':
    center = openreadtxt("center_location.txt",1)
    fiber = openreadtxt("fiber_direction.txt",1)
    matrix = openreadtxt("matrix.txt", 2)
    num_cell = len(center)
    dot_matrix = np.zeros((num_cell, 4))
    for i in range(num_cell):
        index = 0
        for j in matrix[i]:
            if j != -1:
                location_vec = (center[j][0]-center[i][0], center[j][1]-center[i][1], center[j][2]-center[i][2])
                fiber_vec = (fiber[i][0], fiber[i][1], fiber[i][2])
                norm_location = sqrt(location_vec[0]**2 + location_vec[1]**2 + location_vec[2]**2)
                norm_fiber = sqrt(fiber_vec[0]**2 + fiber_vec[1]**2 + fiber_vec[2]**2)
                dot_res = (fiber_vec[0]*location_vec[0] + fiber_vec[1]*location_vec[1] + fiber_vec[2]*location_vec[2]) / (norm_fiber*norm_location)
                dot_matrix[i][index] = fabs(dot_res)
                index = index + 1

    with open("dot_matrix.txt",'w') as f:
        for line in dot_matrix:
            for i in range(4):
                if i < 4-1:
                    f.write(str(line[i])+',')
                else:
                    f.write(str(line[i]) + '\n')

    print("Done!")