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

def calculate_relate(face,cells):
    for i in range(len(cells)):
        if face[0] in cells[i] and face[1] in cells[i] and face[2] in cells[i]:
            return i

if __name__ == '__main__':
    points = openreadtxt('Points.txt',1)
    cells = openreadtxt('Cells.txt',2)
    facets = openreadtxt('Facets.txt',2)
    res = []
    for face in facets:
        res.append(calculate_relate(face,cells))


    with open("faces.txt",'w') as f:
        for index in res:
            f.write(str(index)+'\n')


