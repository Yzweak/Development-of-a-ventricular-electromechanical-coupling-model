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
    return -1 if a[0] < b[0] else 1


if __name__ == '__main__':
    pos_order = []
    res = []

    points = openreadtxt('Points.txt', 1)
    cells = openreadtxt('Cells.txt', 2)
    cell_num = len(cells)



    for index in range(cell_num):
        p1, p2, p3, p4 = points[cells[index][0]], points[cells[index][1]], points[cells[index][2]], points[cells[index][3]]
        x_center = (p1[0] + p2[0] + p3[0] + p4[0]) / 4.0
        res.append((x_center,index))

    s = sorted(res, key=functools.cmp_to_key(compare))

    for i in range(cell_num):
        pos_order.append(s[i][1])

    with open("order.txt",'w') as f:
        for index in pos_order:
            f.write(str(index)+'\n')