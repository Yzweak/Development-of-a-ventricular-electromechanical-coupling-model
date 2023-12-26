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

    for item in order:
        index = item[0]
        p1,p2,p3,p4 = cells[index][0],cells[index][1],cells[index][2],cells[index][3]
        min_x = (points[p1][0] + points[p2][0] + points[p3][0] + points[p4][0]) / 4.0
        if min_x >=4.2 and min_x <=5.8:
            MI_long.append(index)
        elif (min_x >= 3.4 and min_x <4.2) or (min_x <= 6.6 and min_x > 5.8):
            MI_short.append(index)
        elif (min_x >= 2.6 and min_x <3.4) or (min_x <= 7.4 and min_x > 6.6):
            ischemia_1b.append(index)
        elif (min_x >= 1.8 and min_x <2.6) or (min_x <= 8.2 and min_x > 7.4):
            ischemia_1a.append(index)
        elif (min_x >= 1.0 and min_x <1.8) or (min_x <= 9.0 and min_x > 8.2):
            normal.append(index)
        else:
            normal.append(index)


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