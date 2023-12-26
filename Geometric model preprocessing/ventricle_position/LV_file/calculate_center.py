import numpy as np
from math import fabs
def openreadtxt(file_name,mode):
    data = []
    file = open(file_name, 'r')  # 打开文件
    file_data = file.readlines()  # 读取所有行
    for row in file_data:
        tmp_list = row.split('\t')  # 按‘，’切分每行的数据
        tmp_list[-1] = tmp_list[-1].replace('\n','') #去掉换行符
        if mode == 1:
            tmp_list = list(map(float, tmp_list))
        else:
            tmp_list = list(map(int, tmp_list))
        data.append(tmp_list)
    return data



if __name__ == '__main__':
    points = openreadtxt('..\\fiber\\Points.txt',1)
    cells = openreadtxt('..\\fiber\\Cells.txt',2)
    num_cells = len(cells)
    num_points = len(points)
    center = []



    # 计算相邻矩阵
    for i in range(num_cells):
        p1, p2, p3, p4 = cells[i][0], cells[i][1], cells[i][2], cells[i][3]
        center_x = (points[p1][0] + points[p2][0] + points[p3][0] + points[p4][0]) / 4.0
        center_y = (points[p1][1] + points[p2][1] + points[p3][1] + points[p4][1]) / 4.0
        center_z = (points[p1][2] + points[p2][2] + points[p3][2] + points[p4][2]) / 4.0
        center.append((center_x,center_y,center_z))


    with open("center_location.txt",'w') as f:
        for line in center:
            f.write(str(line[0])+',')
            f.write(str(line[1])+',')
            f.write(str(line[2])+'\n')

