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

def calculate_near(cell_neighbor,cells,i,j):
    vertex_i = cells[i]
    vertex_j = cells[j]
    same_point = 0
    for v_i in vertex_i:
        for v_j in vertex_j:
            if v_i == v_j:
                same_point = same_point + 1
                break
    if same_point == 3 :
        cell_neighbor.append(j)



if __name__ == '__main__':
    points = openreadtxt('Points.txt',1)
    cells = openreadtxt('Cells.txt',2)
    num_cells = len(cells)
    num_points = len(points)
    matrix = []
    epi,myo,endo = [],[],[]
    outside_cell,inside = [],[]
    epi_pt,endo_pt = [],[]

    x_left,x_right = 8.40012,19.6031
    y_back,y_front = 10.2868,18.5685
    z_low,z_high = -5.02939,0.0

    # x_left,x_right = -0.392433,0.396316
    # y_back,y_front = -0.401817,0.405352
    # z_low,z_high = -0.782033,0.0

    # 计算相邻矩阵
    for i in range(num_cells):
        cell_neighbor = []
        for j in range(num_cells):
            if i != j:
                calculate_near(cell_neighbor,cells,i,j)
        matrix.append(cell_neighbor)
        p1, p2, p3, p4 = cells[i][0], cells[i][1], cells[i][2], cells[i][3]
        z_cell_1 = 1 if fabs(points[p1][2]) < 0.00001 else 0
        z_cell_2 = 1 if fabs(points[p2][2]) < 0.00001 else 0
        z_cell_3 = 1 if fabs(points[p3][2]) < 0.00001 else 0
        z_cell_4 = 1 if fabs(points[p4][2]) < 0.00001 else 0
        # 四个内部共享面或者基面
        if len(cell_neighbor) == 4 or z_cell_1+z_cell_2+z_cell_3+z_cell_4 == 3:
            myo.append(i)
            inside.append(p1)
            inside.append(p2)
            inside.append(p3)
            inside.append(p4)
        elif len(cell_neighbor) == 3:
            outside_cell.append(i)
    # myo得到内点，遍历outer layer对其中的内点进行距离判断
    for i in outside_cell:
        p1, p2, p3, p4 = cells[i][0], cells[i][1], cells[i][2], cells[i][3]
        if p1 in inside:
            center_tri = ((points[p2][0] + points[p3][0] + points[p4][0]) / 3,
                          (points[p2][1] + points[p3][1] + points[p4][1]) / 3,
                          points[p1][2])
            global_center = ((x_left + x_right) / 2.0, (y_front + y_back) / 2.0, center_tri[2])
            dist_1 = (center_tri[0] - global_center[0]) * (center_tri[0] - global_center[0]) + \
                     (center_tri[1] - global_center[1]) * (center_tri[1] - global_center[1])
            dist_2 = (points[p1][0] - global_center[0]) * (points[p1][0] - global_center[0]) + \
                     (points[p1][1] - global_center[1]) * (points[p1][1] - global_center[1])
            if dist_1 < dist_2:
                endo.append(i)
                endo_pt.append([p2, p3, p4])
            else:
                epi.append(i)
                epi_pt.append([p2, p3, p4])

        elif p2 in inside:
            center_tri = ((points[p1][0] + points[p3][0] + points[p4][0]) / 3,
                          (points[p1][1] + points[p3][1] + points[p4][1]) / 3,
                          points[p2][2])
            global_center = ((x_left + x_right) / 2.0, (y_front + y_back) / 2.0, center_tri[2])
            dist_1 = (center_tri[0] - global_center[0]) * (center_tri[0] - global_center[0]) + \
                     (center_tri[1] - global_center[1]) * (center_tri[1] - global_center[1])
            dist_2 = (points[p2][0] - global_center[0]) * (points[p2][0] - global_center[0]) + \
                     (points[p2][1] - global_center[1]) * (points[p2][1] - global_center[1])
            if dist_1 < dist_2:
                endo.append(i)
                endo_pt.append([p1, p3, p4])
            else:
                epi.append(i)
                epi_pt.append([p1, p3, p4])

        elif p3 in inside:
            center_tri = ((points[p1][0] + points[p2][0] + points[p4][0]) / 3,
                          (points[p1][1] + points[p2][1] + points[p4][1]) / 3,
                          points[p3][2])
            global_center = ((x_left + x_right) / 2.0, (y_front + y_back) / 2.0, center_tri[2])
            dist_1 = (center_tri[0] - global_center[0]) * (center_tri[0] - global_center[0]) + \
                     (center_tri[1] - global_center[1]) * (center_tri[1] - global_center[1])
            dist_2 = (points[p3][0] - global_center[0]) * (points[p3][0] - global_center[0]) +\
                     (points[p3][1] - global_center[1]) * (points[p3][1] - global_center[1])
            if dist_1 < dist_2:
                endo.append(i)
                endo_pt.append([p1, p2, p4])
            else:
                epi.append(i)
                epi_pt.append([p1, p2, p4])

        else:
            center_tri = ((points[p1][0] + points[p2][0] + points[p3][0]) / 3,
                          (points[p1][1] + points[p2][1] + points[p3][1]) / 3,
                           points[p4][2])
            global_center = ((x_left + x_right) / 2.0, (y_front + y_back) / 2.0, center_tri[2])
            dist_1 = (center_tri[0] - global_center[0]) * (center_tri[0] - global_center[0]) + \
                     (center_tri[1] - global_center[1]) * (center_tri[1] - global_center[1])
            dist_2 = (points[p4][0] - global_center[0]) * (points[p4][0] - global_center[0]) + \
                     (points[p4][1] - global_center[1]) * (points[p4][1] - global_center[1])
            if dist_1 < dist_2:
                endo.append(i)
                endo_pt.append([p1, p2, p3])
            else:
                epi.append(i)
                epi_pt.append([p1, p2, p3])

    # removeList = []
    # # 过滤myo层
    # for i in range(len(myo)):
    #     inEndo = False
    #     p1, p2, p3, p4 = cells[myo[i]][0], cells[myo[i]][1], cells[myo[i]][2], cells[myo[i]][3]
    #     if p1 in endo_pt or p2 in endo_pt or p3 in endo_pt or p4 in endo_pt:
    #         endo.append(myo[i])
    #         removeList.append(myo[i])
    #         inEndo = True
    #
    #     if inEndo == False:
    #         if p1 in epi_pt or p2 in epi_pt or p3 in epi_pt or p4 in epi_pt:
    #             epi.append(myo[i])
    #             removeList.append(myo[i])

    # for index in removeList:
    #     myo.remove(index)

    with open("matrix.txt",'w') as f:
        for near_cells in matrix:
            num = len(near_cells)
            for i in range(len(near_cells)):
                if i<3:
                    f.write(str(near_cells[i])+',')
                elif i==3:
                    f.write(str(near_cells[i]) + '\n')
            for i in range(0,4-num):
                if i<3-num:
                    f.write("-1,")
                else:
                    f.write("-1\n")

    with open("endo.txt",'w') as f:
        for index in endo:
            f.write(str(index)+'\n')

    with open("myo.txt", 'w') as f:
        for index in myo:
            f.write(str(index) + '\n')

    with open("epi.txt", 'w') as f:
        for index in epi:
            f.write(str(index) + '\n')

