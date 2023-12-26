import numpy as np
import functools
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



def compare(a, b):
    if a[0][2] != b[0][2]:
        return -1 if a[0][2] < b[0][2] else 1
    else:
        if a[0][1] != b[0][1]:
            return -1 if a[0][1] < b[0][1] else 1
        else:
            return -1 if a[0][0] < b[0][0] else 1


if __name__ == '__main__':
    stim_cells = []
    pos_order = []
    res = []
    endo = openreadtxt('endo.txt',2)
    myo = openreadtxt('myo.txt',2)
    points = openreadtxt('Points.txt', 1)
    cells = openreadtxt('Cells.txt', 2)
    cell_num = len(cells)
    endo_num = len(endo)

    x_left,x_right = 8.40012,19.6031
    y_back,y_front = 10.2868,18.5685
    z_low,z_high = -5.02939,0.0

    for endo_cell in endo:
        index = endo_cell[0]
        p1,p2,p3,p4 = points[cells[index][0]],points[cells[index][1]],points[cells[index][2]],points[cells[index][3]]
        center = ((p1[0]+p2[0]+p3[0]+p4[0])/4.0, (p1[1]+p2[1]+p3[1]+p4[1])/4.0, (p1[2]+p2[2]+p3[2]+p4[2])/4.0)
        if center[2] <= 1.0*(z_high+z_low)/4.0:
            stim_cells.append(index)

    for index in range(cell_num):
        p1, p2, p3, p4 = points[cells[index][0]], points[cells[index][1]], points[cells[index][2]], points[cells[index][3]]
        center = ((p1[0] + p2[0] + p3[0] + p4[0]) / 4.0, (p1[1] + p2[1] + p3[1] + p4[1]) / 4.0,(p1[2] + p2[2] + p3[2] + p4[2]) / 4.0)
        res.append((center,index))

    s = sorted(res, key=functools.cmp_to_key(compare))

    for i in range(cell_num):
        pos_order.append(s[i][1])

    with open("stim_cells.txt",'w') as f:
        for index in stim_cells:
            f.write(str(index)+'\n')

    with open("order.txt",'w') as f:
        for index in pos_order:
            f.write(str(index)+'\n')