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




if __name__ == '__main__':
    cells = openreadtxt('Cells.txt',2)
    points = openreadtxt('Points.txt', 1)
    order = openreadtxt('order.txt', 2)
    normal,ischemia_1a,ischemia_1b,MI_short,MI_long = [],[],[],[],[]
    cell_num = len(cells)
    point_num = len(points)

    x_left,x_right = -0.392433,0.396316
    y_back,y_front = -0.401817,0.405352
    z_low,z_high = -0.782033,0.0

    R = 2.0 * (z_high - z_low) / 16.0
    R1 = 1.0 * R / 4.0
    R2 = 2.0 * R / 4.0
    R3 = 3.0 * R / 4.0
    # circle_center = (0.0, 0.0, -R)
    circle_center = (0.0, 0.0, -4*R)
    for item in order:


        index = item[0]
        p1,p2,p3,p4 = cells[index][0],cells[index][1],cells[index][2],cells[index][3]
        center_x = (points[p1][0] + points[p2][0] + points[p3][0] + points[p4][0]) / 4.0
        center_y = (points[p1][1] + points[p2][1] + points[p3][1] + points[p4][1]) / 4.0
        center_z = (points[p1][2] + points[p2][2] + points[p3][2] + points[p4][2]) / 4.0
        # if center_y < 0:
        #     normal.append(index)
        # else:
        pt_project = (center_x,0,center_z)
        distance = (pt_project[0]-circle_center[0])**2 + (pt_project[1]-circle_center[1])**2 +(pt_project[2]-circle_center[2])**2
        #     # if distance <= R1**2:
        #     #     MI_long.append(index)
        #     # elif distance <= R2**2:
        #     #     MI_short.append(index)
        #     # elif distance <= R3**2:
        #     #     ischemia_1b.append(index)
        #     # elif distance <= R**2:
        #     #     ischemia_1a.append(index)
        #     # else:
        #     #     normal.append(index)
        # if(distance >R**2):
        #     normal.append(index)
        # elif pt_project[2]-circle_center[2] > R2 :
        #     MI_long.append(index)
        # elif (pt_project[2]-circle_center[2] <= R2 and pt_project[2]-circle_center[2]>0):
        #     MI_short.append(index)
        # elif (pt_project[2]-circle_center[2] > -R2 and pt_project[2]-circle_center[2]<=0) :
        #     ischemia_1b.append(index)
        # elif pt_project[2]-circle_center[2] <-R2 :
        #     ischemia_1a.append(index)
        pt_project = (center_x, 0, center_z)
        distance = (pt_project[0] - circle_center[0]) ** 2 + (pt_project[1] - circle_center[1]) ** 2 + (
                    pt_project[2] - circle_center[2]) ** 2
        if(abs(pt_project[2]-circle_center[2])>R ):
            normal.append(index)
        elif pt_project[2]-circle_center[2] > R2 :
            ischemia_1a.append(index)
        elif (pt_project[2]-circle_center[2] <= R2 and pt_project[2]-circle_center[2]>0):
            ischemia_1b.append(index)
        elif (pt_project[2]-circle_center[2] > -R2 and pt_project[2]-circle_center[2]<=0) :
            MI_short.append(index)
        elif pt_project[2]-circle_center[2] <-R2 :
            MI_long.append(index)

    with open("normal_cell.txt",'w') as f:
        for index in normal:
            f.write(str(index)+'\n')

    with open("ischemia-1a.txt",'w') as f:
        for index in ischemia_1a:
            f.write(str(index)+'\n')

    with open("ischemia-1b.txt",'w') as f:
        for index in ischemia_1b:
            f.write(str(index)+'\n')

    with open("MI-short.txt",'w') as f:
        for index in MI_short:
            f.write(str(index)+'\n')

    with open("MI-long.txt",'w') as f:
        for index in MI_long:
            f.write(str(index)+'\n')