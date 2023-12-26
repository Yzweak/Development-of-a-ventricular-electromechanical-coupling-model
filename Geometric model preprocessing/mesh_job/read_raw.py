import numpy as np


def read_raw(file: str, shape: tuple, dtype):
    '''
    读取raw图
    :param file: 文件名
    :param shape: 读取的数据排列，(row,col,channel)
    :param dtype: raw文件类型
    :return: 读取的数据
    '''
    # 从raw文件中读取数据
    data = np.fromfile(file, dtype=dtype)
    # 将读取到的数据重新排列
    data = np.reshape(data, newshape=shape)
    # 返回数据
    return data


def write_raw(file: str, data: np.ndarray):
    '''
    保存raw图
    :param file: 文件名
    :param data: 保存的数据
    :return: 无返回值
    '''
    data.tofile(file)  # 保存数据data到文件file中



if __name__ == '__main__':
    data = read_raw("3D_island.raw",(325,325,425),np.uint8)
    for i in range(325):
        for j in range(325):
            for k in range(425):
                if data[i][j][k] != 20 and data[i][j][k] != 22:
                    data[i][j][k] = 0

    write_raw("left_ventricle.raw",data)
    print("Done!")