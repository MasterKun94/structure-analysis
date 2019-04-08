"""
所有的矩阵相关的函数都写在这里
------
下面的函数会用到几个重要的参数：
------
elb: 每个单元的旋转角 以及结点在结构中的编号
        结构的每个单元和每个节点都有编号，每个单元自身的结点也有自己的编号，通过 elb 这个列表可以
        将单元的结点和结构的结点一一对应起来，这样就可以构造刚度矩阵
    格式：
        elb = list([[elementenumber, rx, ry, rz, id0, id1, id2,  ..., idn], ... ])
        elementnumber 为单元的编号，从0 开始，rx, ry, rz 为单元的旋转角， id为每个单元的结点在结构中的编号
    例子：
        elb = [
            [0, math.pi / 2, 0, 0, 1, 2, 3],
            [1, 0, math.pi / 2, 0, 2, 3, 4],
            [2, 0, 0, math.pi / 2, 3, 4, 1],
            [3, 0, 0, 0, 3, 1, 2],
        ]
------
xy: 结构的所有结点的位置
        另方面，xy为变形前的结点坐标，求解得到结点位移后又可以知道变形后的结点坐标，可以实现变形的可视化
    格式：
        xy = numpy.mat([[nodenumber, tx, ty, tz], ... ], ... )
        nodenumber 为结点的编号，从 0 开始，tx, ty, tz 表示坐标位置
    例子：
        elb = np.mat([
            [0, 0， 0， 0],
            [1, 0， 10， 0],
            [2, 0， 10， 10],
            [3, 10， 10， 10],
        ])
    几个注意事项：
        1、 0 编号单元上应有 0 编号结点
        2、 相邻单元在 elb 中也应该相邻， 最好按编号顺序排列
------
load_nod：受载荷结点的编号以及载荷信息
    格式：
        load_nod = numpy.mat([[nodenumber, mx, my, mz, fx, fy, fz]],  ... )
        nodenumber 为结点的编号，mx, my, mz 表示力矩，fx, fy, fz 表示力
    例子：
        elb = numpy.mat([
            [0, 0， 0， 0， 0， 0， 0],
            [3, 0， 12， 0, 0， 17， 0],
            [7, 0， 15， 10, 0， 30， 0],
            [8, 10， 10， 10, 0， 10， 0],
        ])
------
constraind_nod: 所有受约束结点的编号组成的 列表
------
rigid_nod: 与刚性单元连接的结点的 列表
------
redundant_nod: 冗余结点的 列表
------
注意：
    1) 所有带有编号信息的列表，数组或矩阵，请按编号从小到大排列，不然有些情况下会出问题，或者增加运算时间
    2) 所有编号都从 0 开始一次往下。虽然一般都习惯性的从1开始编号，但是由于python 本身的特点，这里从 0 开始编号会更方便
    3) 尽可能使用向量化编程
------
"""

import copy
import warnings

import numpy as np


def adjoint(tx, ty, tz, rx, ry, rz):
    """
    坐标变换矩阵
    ------
    :param tx: x 轴方向 位移
    :param ty: y 轴方向 位移
    :param tz: z 轴方向 位移
    :param rx: x 轴方向 旋转
    :param ry: y 轴方向 旋转
    :param rz: z 轴方向 旋转
    :return: 矩阵 Ad
    """
    r = rotation_xyz(rx, ry, rz)
    t = np.mat([
        [0, -tz, ty],
        [tz, 0, -tx],
        [-ty, tx, 0]
    ])

    ad = np.vstack([[r, np.zeros((3, 3))],
                    [t * r, r]])
    return ad


def adjoint_k(tx, ty, tz, rx, ry, rz):
    """
    坐标变换矩阵，上个函数 adjoint 的逆矩阵，
    适合用于刚度矩阵的坐标变换
    ------
    :param tx: x 轴方向 位移
    :param ty: y 轴方向 位移
    :param tz: z 轴方向 位移
    :param rx: x 轴方向 旋转
    :param ry: y 轴方向 旋转
    :param rz: z 轴方向 旋转
    :return: 矩阵 Ad 的逆矩阵
    """
    r = rotation_xyz(rx, ry, rz).T
    t = np.mat([
        [0, -tz, ty],
        [tz, 0, -tx],
        [-ty, tx, 0]
    ])

    ad_k = np.bmat([[r, np.zeros((3, 3))],
                    [(-r) * t, r]])
    return ad_k


def rotation_xyz(rx, ry, rz):
    """
    旋转矩阵
    ------
    :param rx: x 轴方向 旋转
    :param ry: y 轴方向 旋转
    :param rz: y 轴方向 旋转
    :return:
    """
    cx, cy, cz = np.cos(rx), np.cos(ry), np.cos(rz)
    sx, sy, sz = np.sin(rx), np.sin(ry), np.sin(rz)
    a = np.mat([
        [1, 0, 0],
        [0, cx, -sx],
        [0, sx, cx]
    ])
    b = np.mat([
        [cy, 0, sy],
        [0, 1, 0],
        [-sy, 0, cy]
    ])
    c = np.mat([
        [cz, -sz, 0],
        [sz, cz, 0],
        [0, 0, 1]
    ])
    r = a * b * c

    return r


def build_kmat(elb, kmat=None, shape=None):

    """
    构造有限元刚度矩阵
    构造刚度矩阵需要将每个单元的刚度矩阵输入到结构刚度矩阵合适的位置中，通过 elb 可以做到，xy 确定了结点的数量也就是刚度矩阵的维数。
    最简单的情况：
        如果用来构造的单元只有一种，那么只需要一个 elb 列表，kmat = None 使用一次 该函数 就可以得到结构的刚度矩阵
    很多情况下：
        结构由多个不同的单元组成，这时候 需要多个 elb 列表，每个elb 表示一种单元，并且多次使用 该函数 来构造刚度矩阵。
        第一次使用该函数时，参数 kmat = None，后面使用该函数时，参数 kmat 为上一次返回的值，这样就可以不断叠加 最后返回的
        kmat 为最终的刚度矩阵。
        具体例子可见 structure_class.py 的 line (15 - 26)
    ------
    :param elb: 每个单元的旋转角 以及结点在结构中的编号
        类型：列表
        格式：具体请看 开头说明
        结构的每个单元和每个节点都有编号，每个单元自身的结点也有自己的编号，通过 elb 这个列表可以
        将单元的结点和结构的结点一一对应起来，这样就可以构造刚度矩阵

    :param shape:
    :param kmat: 构造前的结构刚度矩阵
        类型：矩阵
    :return: 结构刚度矩阵
    """
    if shape is None:
        if kmat is None:
            if len(elb) == 1:
                maxnum = elb[0]['elb'].max()
            else:
                # warnings.warn('请定义 shape ，不然可能影响运算速度')
                print('请定义 shape ，不然可能影响运算速度')
                maxnum = 0
                for i in range(len(elb)):
                    maxnum0 = elb[i]['elb'].max()
                    if maxnum < maxnum0:
                        maxnum = maxnum0
            kmat = np.zeros((maxnum * 6 + 6, maxnum * 6 + 6))
    else:
        if kmat is None:
            kmat = np.zeros(shape)

    for ielb in elb:
        elb_ele = ielb['elb']
        elb_rot = ielb['rotation']
        kele = ielb['class'].kmat()
        hele, lele = elb_ele.shape
        hrot = elb_rot.shape

        if hrot == (3, ) or hrot[0] == 1:
            ad = adjoint_k(0, 0, 0, elb_rot[0], elb_rot[1], elb_rot[2])
            for ele in range(0, hele):
                for nod1 in range(0, lele):
                    ele_nod1 = elb_ele[ele][nod1]
                    for nod2 in range(0, lele):
                        ele_nod2 = elb_ele[ele][nod2]
                        kmat[6 * ele_nod1:6 * ele_nod1 + 6, 6 * ele_nod2:6 * ele_nod2 + 6] += \
                            ad.T * kele[6 * nod1:6 * nod1 + 6, 6 * nod2: 6 * nod2 + 6] * ad

        elif hrot[0] == hele:
            for ele in range(0, hele):
                ad = adjoint_k(0, 0, 0, elb_rot[ele][0], elb_rot[ele][1], elb_rot[ele][2])
                for nod1 in range(0, lele):
                    ele_nod1 = elb_ele[ele][nod1]
                    for nod2 in range(0, lele):
                        ele_nod2 = elb_ele[ele][nod2]
                        kmat[6 * ele_nod1:6 * ele_nod1 + 6, 6 * ele_nod2:6 * ele_nod2 + 6] += \
                            ad.T * kele[6 * nod1:6 * nod1 + 6, 6 * nod2: 6 * nod2 + 6] * ad
    return kmat


def solve_kmat(constraind_nod, load_nod, kmat):  # ？？？？？？？？？？？
    """

    知道结构的约束 constrand_nod 后，通过删去刚度矩阵中对应的部分，就可以求逆得到柔度矩阵
    知道载荷 load_node 后就可以得到 结点的位移
    -------
    :param constraind_nod: 约束节点所组成的列表
        类型：列表
    :param load_nod: 每个受载荷结点的受力信息
        类型： 矩阵， 第一列为结点编号， 后面 6 列为载荷 >>> [[nodenumber, mx, my, mz, fx, fy, fz]]
    :param kmat: 刚度矩阵
        类型：矩阵
    :return: 位移矩阵
    """
    num = kmat.shape[0]
    load_nod_list = load_nod[:, 0].T.tolist()[0]
    for nod in constraind_nod:
        if nod in load_nod_list:
            warnings.warn("受载荷结点 " + str(nod) + " 与约束结点重合")
        if nod >= num/6:
            warnings.warn("受约束结点 " + str(nod) + " 编号不对")
    for nod in load_nod_list:
        if nod >= num/6:
            warnings.warn("受载荷结点 " + str(nod) + " 编号不对")

    cons = np.array(constraind_nod) * 6
    cons = sorted(np.hstack((cons, cons + 1, cons + 2, cons + 3, cons + 4, cons + 5)))
    kmat = np.mat(np.delete(np.delete(kmat, cons, axis=0), cons, axis=1))
    # pf.save_exl(kmat, 'C:\work\python\structure1\Try.xlsx')
    wrench_mat = np.zeros((num, 1))
    nl = np.arange(0, load_nod.shape[0])
    for n in [0, 1, 2, 3, 4, 5]:
        wrench_mat[load_nod[nl, 0] * 6 + n, 0] = load_nod[nl, n + 1]
    wrench_mat = np.delete(wrench_mat, cons, axis=0)
    # pf.save_exl(kmat, 'C:\work\python\structure1\Try.xlsx')
    twist_mat = kmat.I * wrench_mat
    return twist_mat


# def build_xy(elb, xy_ele):  # ？？？？？？？？？？
#     """
#     构造 xy 结点坐标
#     ------
#     :param elb: 单元的结点在结构中的编号   列表
#         几个注意事项：
#             1、 0 编号单元上应有 0 编号结点
#             2、 相邻单元在 elb 中也应该相邻， 最好按编号顺序排列
#     :param xy_ele: 单元结点的坐标    矩阵
#     :return: 结构的所有节点坐标
#     """
#     xy_ele = xy_ele[np.lexsort(np.array(xy_ele)[:, ::-1].T)][:, 1:]
#     xy_13, xy_0, ele_nod_num = np.mat([[0, 0, 0]]), [0], xy_ele.shape[0]
#
#     for en in range(0, len(elb)):
#         elb_en = elb[en]
#         if en != 0:                                                         # 编号除 0 以外的单元
#             i = list(set(elb_en[4:]) & set(elb[en - 1][4:]))
#             if i:                                                           # 如果两个单元有相同结点
#                 xy_ind = xy_13[xy_0.index(i[0]), :]                         # 根据 xy_ele 和 elb 就可以知道结点坐标
#                 xy_per = (xy_ele - xy_ele[elb_en[4:].index(i[0]), :]) * rotation_xyz(
#                     -elb_en[1], -elb_en[2], -elb_en[3])
#                 for enn in range(0, ele_nod_num):
#                     if elb_en[enn + 4] not in xy_0:
#                         xy_13 = np.row_stack((xy_13, xy_per[enn, :] + xy_ind))
#                         xy_0.append(elb_en[enn + 4])
#             else:                                                           # 如果没有相同结点，那么 elb 有问题
#                 warnings.warn("相邻编号的单元应相连，也就是说有相同的节点")
#
#         elif en in elb_en[4:]:                                               # 编号 0 的单元
#             xy_per = (xy_ele - xy_ele[elb_en[4:].index(0), :]) * \
#                      rotation_xyz(-elb_en[1], -elb_en[2], -elb_en[3])
#             for enn in range(0, ele_nod_num):
#                 if elb_en[enn + 4] not in xy_0:
#                     xy_13 = np.row_stack((xy_13, xy_per[enn, :]))
#                     xy_0.append(elb_en[enn + 4])
#                 # else:
#         else:
#             warnings.warn("编号 0 的单元应有编号 0 的结点")
#     xy = np.array(np.hstack((np.mat(xy_0).T, xy_13)))           # 将结点编号xy_0和对应坐标xy_13组合
#     xy = xy[np.lexsort(xy[:, ::-1].T)]                          # 按第一列顺序排序
#     return xy


def new_xy(xy, twist_mat, constraind_nod):
    """
    所有结点在变形以后的坐标， 以及结点位移
    ------
    :param xy: 变形前的结点坐标   矩阵
    :param twist_mat: 矩阵运算得到的未约束结点位移矩阵
    :param constraind_nod: 受约束结点
    :return: 变形后的xy
            结点的位移
    """
    noncons_list = [i for i in range(0, np.shape(xy)[0]) if i not in constraind_nod]
    twist_mat = np.reshape(twist_mat, (len(noncons_list), 6))
    twist_mat_all = np.zeros(np.shape(xy)[0], 6)
    twist_mat_all[noncons_list] = twist_mat
    xy[:, 1:4] = twist_mat_all[:, 3:6]

    return xy,  twist_mat_all                       # 返回结点坐标和位移列表，放在一个列表里


def rigidelement_k(xy, featurexy, kmat, rigid_node):

    """
    刚性单元的刚度矩阵变换过程
    如果使用该函数后需要重新返回原有维数的柔度矩阵，或是需要重新观测连接刚性单元的结点位移，
    强烈介意将这些结点的编号设置为最前，也就是 connectnode = [0, 1, 2, 3, ......]
    如果只是作为刚性连接，则不用在意
    ------
    :param xy: 结点的坐标    矩阵
    :param featurexy: 特征结点的坐标    矩阵
    :param kmat: 变换前的刚度矩阵
    :param rigid_node: 刚性单元连接的结点的编号  列表
    :return: 一个列表，第一个为消除了刚性连接结点的 xy ，特征结点的编号为 connectnode 中 第一个结点 的编号
             第二个为变换后的刚度矩阵
    """

    mapto = xy[:, 1:4] - featurexy
    adbig, cnd = [], []
    for i in rigid_node:
        cnd = cnd + [6 * i, 6 * i + 1, 6 * i + 2, 6 * i + 3, 6 * i + 4, 6 * i + 5]
        adbig = adbig + [adjoint_k(mapto[i, 0], mapto[i, 1], mapto[i, 2], 0, 0, 0)]
    cnd = sorted(cnd)
    adbig = np.vstack(adbig)
    kmat = copy.deepcopy(kmat)  # http://www.cnblogs.com/monkey-moon/p/9347505.html
    kmat[:, cnd[0:6]] = kmat[:, cnd] * adbig
    kmat[cnd[0:6], :] = adbig.T * kmat[cnd, :]
    kmat = np.delete(np.delete(kmat, cnd[6:], axis=1), cnd[6:], axis=0)   # 删除多余的行和列
    xy = np.delete(xy, rigid_node[1:], axis=0)
    xy[rigid_node[0], 1:4] = featurexy

    return xy, kmat


def redundantnode_simplify(xy, kmat, redundantnode):
    """
    冗余结点单元的刚度矩阵变换过程
    ------
    :param xy: 结点的坐标    矩阵
    :param kmat: 变换前的刚度矩阵
    :param redundantnode: 冗余结点的编号  列表
    :return: 一个列表，第一个为消除了冗余结点的 xy
             第二个为变换后的刚度矩阵
    """
    rnd, nnd, kmat = [], list(range(0, 6 * xy.shape[0])), np.mat(kmat)
    for i in redundantnode:
        rnd = rnd + [6 * i, 6 * i + 1, 6 * i + 2, 6 * i + 3, 6 * i + 4, 6 * i + 5]
    rnd = sorted(rnd)
    for i in rnd:
        nnd.remove(i)
    kmat = kmat[nnd, :][:, nnd] - kmat[nnd, :][:, rnd] * kmat[rnd, :][:, rnd].I * kmat[rnd, :][:, nnd]
    xy = np.delete(xy, redundantnode, axis=0)

    return xy, kmat


def build_constriand(kmat, constraindnode):
    cons = np.array(constraindnode) * 6
    cons = sorted(np.hstack((cons, cons + 1, cons + 2, cons + 3, cons + 4, cons + 5)))
    kmat = np.mat(np.delete(np.delete(kmat, cons, axis=0), cons, axis=1))
    return kmat
