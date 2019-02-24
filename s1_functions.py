"""
structure_class.py 中的类 Strcture_1() 的相关函数
"""
import sys
import warnings

import numpy as np

import process_functions as pf


def build_pelb_sp(row, nx, ny, nz, direction, elementclass, pelb=None):

    def node_identifier_pe(i, j, k):
        if i <= 2 * nx and j <= 2 * ny and k <= 2 * nz and (i % 2 + j % 2 + k % 2) == 1:
            if k % 2 == 0:
                return (3 * nx * ny + 2 * nx + 2 * ny + 1) * int(k / 2) + ny * int((i + 1) / 2) + (
                        ny + 1) * int(i / 2) + int(j / 2)
            else:
                return ((nx + 1) * ny + (ny + 1) * nx) * int((k + 1) / 2) + (nx + 1) * (ny + 1) * int(
                        (k - 1) / 2) + (ny + 1) * int(i / 2) + int(j / 2)
        else:
            raise ValueError('nd', i, j, k, ' 编号不对')

    def ele_identifier_pe():

        if direction == 0:
            elb['rotation'] = np.array([0, np.pi / 2, 0])
            for line in range(0, nz):
                for col in range(0, ny):
                    yield [[
                        node_identifier_pe(2 * nx - 2 * row, 2 * col + 1, 2 * line),
                        node_identifier_pe(2 * nx - 2 * row, 2 * col + 2, 2 * line + 1),
                        node_identifier_pe(2 * nx - 2 * row, 2 * col + 1, 2 * line + 2),
                        node_identifier_pe(2 * nx - 2 * row, 2 * col, 2 * line + 1)
                    ]]
        elif direction == 1:
            elb['rotation'] = np.array([-np.pi / 2, 0, np.pi / 2])
            for line in range(0, nz):
                for col in range(0, nx):
                    yield [[
                        node_identifier_pe(2 * nx - 2 * col - 1, 2 * ny - 2 * row, 2 * line),
                        node_identifier_pe(2 * nx - 2 * col - 2, 2 * ny - 2 * row, 2 * line + 1),
                        node_identifier_pe(2 * nx - 2 * col - 1, 2 * ny - 2 * row, 2 * line + 2),
                        node_identifier_pe(2 * nx - 2 * col, 2 * ny - 2 * row, 2 * line + 1)
                    ]]
        elif direction == 2:
            elb['rotation'] = np.array([0, 0, 0])
            for line in range(0, nx):
                for col in range(0, ny):
                    yield [[
                        node_identifier_pe(2 * line, 2 * col + 1, 2 * row),
                        node_identifier_pe(2 * line + 1, 2 * col + 2, 2 * row),
                        node_identifier_pe(2 * line + 2, 2 * col + 1, 2 * row),
                        node_identifier_pe(2 * line + 1, 2 * col, 2 * row)
                    ]]
        else:
            warnings.warn('参数 ', str(direction), ' 有问题')
            sys.exit()

    if pelb is None:
        pelb = []
    elb, elb_list = {'name': direction, 'class': elementclass}, []
    for ele in ele_identifier_pe():
        pelb = pelb + ele
    elb['elb'] = np.array(pelb)
    elb_list.append(elb)
    return elb_list


def build_xy_s1(ele_d, nx, ny, nz):
    xy, num= np.mat([[0, 0, ele_d / 2, 0]]), 0
    for k in range(0, 2 * nz + 1):
        for i in range(0, 2 * nx + 1):
            for j in range(0, 2 * ny + 1):
                if (i % 2 + j % 2 + k % 2) == 1:
                    if num != 0:
                        xy = np.row_stack((xy, np.mat([[num, ele_d / 2 * i, ele_d / 2 * j, -ele_d / 2 * k]])))
                        num += 1
                    elif num == 0:
                        num = 1
    return xy


def choose_sita(nx, ny, nz, ele_d, p):
    sitalist = [
        pf.screw_solve_fy((nx / 2 - np.arange(0, nx + 1)) * ele_d, p),
        pf.screw_solve_fy((ny / 2 - np.arange(0, ny + 1)) * ele_d, p),
        pf.screw_solve_fy((nz / 2 - np.arange(0, nz + 1)) * ele_d, p),
    ]
    return sitalist
