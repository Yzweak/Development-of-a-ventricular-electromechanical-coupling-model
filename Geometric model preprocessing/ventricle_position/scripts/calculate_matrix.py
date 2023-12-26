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

    if same_point == 3:
        cell_neighbor.append(j)

if __name__ == '__main__':
    points = openreadtxt('Points.txt',1)
    cells = openreadtxt('Cells.txt',2)
    edges = openreadtxt('Edges.txt',2)
    num_cells = len(cells)
    num_points = len(points)
    num_edges = len(edges)
    matrix = []
    epi,myo,endo = [],[],[]
    inside,outside = [],[]
    epi_p,endo_p,base = [],[],[]
    pts = dict()

    # 计算内外边界顶点
    for i in range(num_edges):
        first = edges[i][0]
        second = edges[i][1]
        if pts.get(first) == None:
            pts[first] = 1
        else:
            pts[first] = pts[first] + 1
        if pts.get(second) == None:
            pts[second] = 1
        else:
            pts[second] = pts[second] + 1
    z_min = 0
    apex = 1
    for i in range(num_points):
        # 计算base面上的顶点
        if points[i][2] == 0:
            base.append(i)
        if pts[i] <= 6:
            outside.append(i)
            if points[i][2] < z_min:
                zmin = points[i][2]
                apex = i
        else:
            inside.append(i)
    # 从外边界中删除base顶点，只剩下epi和endo顶点
    for b_p in base:
        outside.remove(b_p)

    # 计算单元临界矩阵
    for i in range(num_cells):
        cell_neighbor = []
        for j in range(num_cells):
            if i != j:
                calculate_near(cell_neighbor,cells,i,j)
        matrix.append(cell_neighbor)
        # 四个内部共享面
        if len(cell_neighbor) == 4:
            myo.append(i)


    for i in range(num_cells):
        p1, p2, p3, p4 = cells[i][0], cells[i][1], cells[i][2], cells[i][3]
        z_cell_1 = 1 if fabs(points[p1][2]) < 0.00001 else 0
        z_cell_2 = 1 if fabs(points[p2][2]) < 0.00001 else 0
        z_cell_3 = 1 if fabs(points[p3][2]) < 0.00001 else 0
        z_cell_4 = 1 if fabs(points[p4][2]) < 0.00001 else 0
        # 基面
        if z_cell_1+z_cell_2+z_cell_3+z_cell_4 == 3:
            myo.append(i)
        # 三面共享
        elif len(matrix[i]) == 3:
            if p1 in inside:
                center_tri = ((points[p2][0]+points[p3][0]+points[p4][0])/3,(points[p2][1]+points[p3][1]+points[p4][1])/3,(points[p2][2]+points[p3][2]+points[p4][2])/3)
                global_center = ((2.00346+2.817)/2.0, (2.13965+2.94121)/2.0, center_tri[2])
                dist_1 = (center_tri[0]-global_center[0])*(center_tri[0]-global_center[0]) + (center_tri[1]-global_center[1])*(center_tri[1]-global_center[1]) + (center_tri[2]-global_center[2])*(center_tri[2]-global_center[2])
                dist_2 = (points[p1][0]-global_center[0])*(points[p1][0]-global_center[0]) + (points[p1][1]-global_center[1])*(points[p1][1]-global_center[1]) + (points[p1][2]-global_center[2])*(points[p1][2]-global_center[2])
                if dist_1 < dist_2:
                    endo.append(i)
                else:
                    epi.append(i)
            elif p2 in inside:
                center_tri = ((points[p1][0]+points[p3][0]+points[p4][0])/3,(points[p1][1]+points[p3][1]+points[p4][1])/3,(points[p1][2]+points[p3][2]+points[p4][2])/3)
                global_center = ((2.00346+2.817)/2.0, (2.13965+2.94121)/2.0, center_tri[2])
                dist_1 = (center_tri[0]-global_center[0])*(center_tri[0]-global_center[0]) + (center_tri[1]-global_center[1])*(center_tri[1]-global_center[1]) + (center_tri[2]-global_center[2])*(center_tri[2]-global_center[2])
                dist_2 = (points[p2][0]-global_center[0])*(points[p2][0]-global_center[0]) + (points[p2][1]-global_center[1])*(points[p2][1]-global_center[1]) + (points[p2][2]-global_center[2])*(points[p2][2]-global_center[2])
                if dist_1 < dist_2:
                    endo.append(i)
                else:
                    epi.append(i)
            elif p3 in inside:
                center_tri = ((points[p1][0]+points[p2][0]+points[p4][0])/3,(points[p1][1]+points[p2][1]+points[p4][1])/3,(points[p1][2]+points[p2][2]+points[p4][2])/3)
                global_center = ((2.00346+2.817)/2.0, (2.13965+2.94121)/2.0, center_tri[2])
                dist_1 = (center_tri[0]-global_center[0])*(center_tri[0]-global_center[0]) + (center_tri[1]-global_center[1])*(center_tri[1]-global_center[1]) + (center_tri[2]-global_center[2])*(center_tri[2]-global_center[2])
                dist_2 = (points[p3][0]-global_center[0])*(points[p3][0]-global_center[0]) + (points[p3][1]-global_center[1])*(points[p3][1]-global_center[1]) + (points[p3][2]-global_center[2])*(points[p3][2]-global_center[2])
                if dist_1 < dist_2:
                    endo.append(i)
                else:
                    epi.append(i)
            else:
                center_tri = ((points[p1][0]+points[p2][0]+points[p3][0])/3,(points[p1][1]+points[p2][1]+points[p3][1])/3,(points[p1][2]+points[p2][2]+points[p3][2])/3)
                global_center = ((2.00346+2.817)/2.0, (2.13965+2.94121)/2.0, center_tri[2])
                dist_1 = (center_tri[0]-global_center[0])*(center_tri[0]-global_center[0]) + (center_tri[1]-global_center[1])*(center_tri[1]-global_center[1]) + (center_tri[2]-global_center[2])*(center_tri[2]-global_center[2])
                dist_2 = (points[p4][0]-global_center[0])*(points[p4][0]-global_center[0]) + (points[p4][1]-global_center[1])*(points[p4][1]-global_center[1]) + (points[p4][2]-global_center[2])*(points[p4][2]-global_center[2])
                if dist_1 < dist_2:
                    endo.append(i)
                else:
                    epi.append(i)

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

