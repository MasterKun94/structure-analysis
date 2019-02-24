"""
除了 mat_functions 其他一些可能需要的函数
或者是后续分析中会用到的函数
------
"""
import sys
import time
import warnings

import numpy as np
import scipy.io
import xlsxwriter


class Monster:

    desk = None

    def prepare(self, new_food=None, weight=1, food=None):
        if food is None:
            food = FoodForMonster(new_food, weight)
        if self.desk is None:
            self.desk = [food]
        elif self.desk[-1].weight == food.weight:
            self.magic_desk(food)
        else:
            self.desk.append(food)

    def magic_desk(self, other):
        self.desk[-1].add_food(other)
        try:
            if self.desk[-2].weight == self.desk[-1].weight:
                self.prepare(food=self.desk.pop())
        except IndexError:
            pass

    def eat_all(self):
        if self.desk is None:
            return None
        elif len(self.desk) == 1:
            self.desk = self.desk[0]
        else:
            food = self.desk.pop(-1)
            self.prepare(food.food, food.weight * 2)
            self.eat_all()


class FoodForMonster:
    def __init__(self, food, weight):
        self.weight = weight
        self.food = food

    def add_food(self, other):
        self.food += other.food
        self.weight += other.weight


def time_counter(func):

    def inner(*args, **kwargs):

        start_time = time.perf_counter()
        name = func.__name__
        print('Running <' + name + '>')

        # done = False
        # # here is the animation
        # def animate():
        #     for c in itertools.cycle(['|', '/', '-', '\\']):
        #         if done:
        #             break
        #         sys.stdout.write('\rloading ' + c)
        #         sys.stdout.flush()
        #         time.sleep(0.1)
        #     sys.stdout.write('\rDone!     ')
        # t = threading.Thread(target=animate)
        # t.start()

        res = func(*args, **kwargs)
        # done = True
        print('\nFunction <' + name + '> time cost:', time.perf_counter() - start_time,
              '\n --------------------------------------------------')

        return res
    return inner


def param_factory(func):

    def gen_int_float_str(param):
        while True:
            yield param

    def gen_list(param):
        for ii in param:
            for jj in ii:
                yield jj

    def gen_ele(ele):
        while True:
            try:
                newparam = {}
                for ikey, ivalue in ele.items():
                    newparam[ikey] = next(ivalue)
                yield newparam
            except StopIteration:
                break

    def inner(self, lat_name, lat_param, ele_name, ele_param, ele_size, **kwargs):
        for key, value in ele_param.items():
            if isinstance(value, (int, float, str)):
                value = gen_int_float_str(value)
            elif isinstance(value, list):
                value = gen_list(value)
            ele_param[key] = value
        ele_param = gen_ele(ele_param)
        goods = func(self, lat_name, lat_param, ele_name, ele_param, ele_size, **kwargs)
        return goods
    return inner


def save_exl(data, path, sheet='sheet1', magnitude=1e-10):
    """
    将数据（矩阵）保存到 excel 中
    ------
    :param data: 数据
    :param path: 路径
    :param sheet:
    :param magnitude: 某个数值如果绝对值小于magnitude 就可以看成 0, 用来略掉浮点运算产生的误差，不然用excel看太丑了
    :return: 一个 xlsx 文件
    """
    f = xlsxwriter.Workbook(path)  # 创建工作簿
    sheet1 = f.add_worksheet(sheet)  # 创建sheet
    [h, l] = data.shape  # h为行数，l为列数
    for ii in range(0, h):
        for jj in range(0, l):
            if abs(data[ii, jj]) > magnitude:
                sheet1.write(ii, jj, data[ii, jj])
            else:
                sheet1.write(ii, jj, 0)
    f.close()


def save_mat(value, value_name, path):
    """
    将数据（矩阵）保存成 matlab 的 mat 文件
    注意：这个文件是一个字典，参数名称value_name是键，参数的值value是指
    ------
    :param value: 数据，字典的值
    :param value_name: 参数名，字典的键
    :param path: 路径
    :return: 一个 mat 文件
    """
    scipy.io.savemat(path, {value_name: value})


def screw_solve_fy(distance, pvalue, qvalue=0):
    py = np.arctan((pvalue + qvalue) / distance)
    return py


def build_wrench(f, r, q, reciprocal='yes'):
    if reciprocal is 'no':
        wrench = np.vstack((f.T, np.cross(r, f).T + q * f.T))
    elif reciprocal is 'yes':
        wrench = np.vstack((np.cross(r, f).T + q * f.T, f.T))
    else:
        warnings.warn('函数 build_wrench 有问题')
        sys.exit()
    return wrench


def build_twist(w, c, p, reciprocal='no'):
    """
    构造一条旋量
    ------
    :param w: 螺旋轴线
        类型：1 * 3 的矩阵
    :param c: 轴线上一点的位置
        类型：1 * 3 的矩阵
    :param p: 螺距
    :param reciprocal: 是否交换主部和副部
    :return: 旋量，一个 6 * 1 的矩阵
    """
    if reciprocal is 'no':
        twist = np.vstack((w.T, np.cross(c, w).T + p * w.T))
    elif reciprocal is 'yes':
        twist = np.vstack((np.cross(c, w).T + p * w.T, w.T))
    else:
        warnings.warn('函数 build_twist 有问题')
        sys.exit()
    return twist


if __name__ == '__main__':
    monster = Monster()
    print(monster.desk)
    for i in range(10):
        monster.prepare(np.array([[i, i, i], [1, 2, 3]]))
        # monster.prepare(i)
        print(i, " ", monster.desk[-1].food, "len", len(monster.desk))
    monster.eat_all()
    print(monster.desk[-1].food, "len", len(monster.desk))
