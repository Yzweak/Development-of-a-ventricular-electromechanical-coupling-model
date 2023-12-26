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
    Faces = openreadtxt('Faces.txt',2)
    Flags = openreadtxt('Face_flag.txt',2)

    num_cells = len(cells)
    num_points = len(points)
    matrix = []
    epi,myo,endo = [],[],[]
    epi_pt,endo_pt = [],[]

    for i in range(len(Flags)):
        p_1,p_2,p_3 = Faces[i][0], Faces[i][1], Faces[i][2]
        # EPI
        if Flags[i][0] == 1:
            epi_pt.append(p_1)
            epi_pt.append(p_2)
            epi_pt.append(p_3)
        # ENDO
        elif Flags[i][0] == 2:
            endo_pt.append(p_1)
            endo_pt.append(p_2)
            endo_pt.append(p_3)




    # 计算相邻矩阵
    for i in range(num_cells):
        cell_neighbor = []
        for j in range(num_cells):
            if i != j:
                calculate_near(cell_neighbor,cells,i,j)
        matrix.append(cell_neighbor)
        p1, p2, p3, p4 = cells[i][0], cells[i][1], cells[i][2], cells[i][3]

        if p1 in epi_pt and p2 and epi_pt and p3 in epi_pt or\
           p1 in epi_pt and p2 and epi_pt and p4 in epi_pt or\
           p1 in epi_pt and p3 and epi_pt and p4 in epi_pt or\
           p2 in epi_pt and p3 and epi_pt and p4 in epi_pt:
            epi.append(i)
        elif p1 in endo_pt and p2 and endo_pt and p3 in endo_pt or\
           p1 in endo_pt and p2 and endo_pt and p4 in endo_pt or\
           p1 in endo_pt and p3 and endo_pt and p4 in endo_pt or\
           p2 in endo_pt and p3 and endo_pt and p4 in endo_pt:
            endo.append(i)
        # if p1 in epi_pt or p2 in epi_pt or p3 in epi_pt or p4 in epi_pt:
        #     epi.append(i)
        # elif p1 in endo_pt or p2 in endo_pt or p3 in endo_pt or p4 in endo_pt:
        #     endo.append(i)
        else:
            myo.append(i)


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

