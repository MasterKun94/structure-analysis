"""
这里是所有基本单元的类，输入尺寸参数来建立基本单元的各种信息
------
一般都是两结点的单元
"""
import math

import numpy as np
import scipy.sparse as sp

import mat_functions as mf


class StrbeamCC:
    """
    圆截面直梁
    梁长方向为 x-轴 方向，y-z 平面与截面平行

    假设：x-轴向右，y-轴向上，z-轴垂直纸面向上，则以下：
    ------
    compliance_mat(self):
        若梁的左端接地，求右端的柔度矩阵
    ------
    xy(self):
        梁的两个结点坐标， 右端结点编号为 0 ，坐标为原点； 左端编号为 1
    ------
    kmat_ele(self):
        单元刚度矩阵，由于有两个结点，所以是12 * 12 的矩阵
    ------
    """
    def __init__(self, rad, lens, ela, v):
        """
        初始化
        :param rad: 圆截面半径   单位：mm
        :param lens: 梁长度     单位：mm       %%%%%%%% 坐标系的x轴方向 %%%%%%%%%%%%
        :param ela: 弹性模量     单位：Mpa
        :param v: 泊松比
        """
        self.rad = rad
        self.lens = lens
        self.ela = ela
        self.v = v

    def compliance_mat(self):

        ix = math.pi * self.rad ** 4 / 4
        j = 2 * ix
        a = math.pi * self.rad ** 2
        g = self.ela / (2 * (1 + self.v))
        c = np.mat([
            [self.lens / (g * j), 0, 0, 0, 0, 0],
            [0, self.lens / (self.ela * ix), 0, 0, 0, 0],
            [0, 0, self.lens / (self.ela * ix), 0, 0, 0],
            [0, 0, 0, self.lens / (self.ela * a), 0, 0],
            [0, 0, 0, 0, self.lens ** 3 / (12 * self.ela * ix), 0],
            [0, 0, 0, 0, 0, self.lens ** 3 / (12 * self.ela * ix)]
        ])
        ad = mf.adjoint(-self.lens / 2, 0, 0, 0, 0, 0)
        c = ad * c * ad.T

        return c

    def xy(self):
        xy = np.mat([
            [0, 0,          0, 0],
            [1, -self.lens, 0, 0],
        ])
        return xy

    def kmat(self):
        xy_ele = self.xy()
        ad12 = mf.adjoint_k(
            xy_ele[0, 1] - xy_ele[1, 1],
            xy_ele[0, 2] - xy_ele[1, 2],
            xy_ele[0, 3] - xy_ele[1, 3],
            0, 0, 0)
        k11 = self.compliance_mat().I
        k12 = k11 * ad12
        k22 = ad12.T * k12
        k21 = k22 * ad12.I
        kele = np.vstack((np.hstack((k11, -k12)),
                          np.hstack((-k21, k22))))

        return kele

    def kmat_sparse(self):
        return sp.bsr_matrix(self.kmat(), blocksize=(6, 6))


class StrbeamSQ:
    """
    矩形截面直梁
    梁长方向为 x-轴 方向，梁宽为 y-轴 方向，梁厚为 z-轴方向

    假设：x-轴向右，y-轴向上，z-轴垂直纸面向上，则以下：
    ------
    compliance_mat(self):
        若梁的左端接地，求右端的柔度矩阵
    ------
    xy(self):
        梁的两个结点坐标， 右端结点编号为 0 ，坐标为原点； 左端编号为 1
    ------
    kmat_ele(self):
        单元刚度矩阵，由于有两个结点，所以是12 * 12 的矩阵
    ------
    """

    def __init__(self, thc, wid, lens, ela, v):
        """
        初始化
        :param thc: 圆截面厚    单位：mm       %%%%%%%% 坐标系的 z-轴方向 %%%%%%%%%%%%
        :param wid: 圆截面宽    单位：mm       %%%%%%%% 坐标系的 y-轴方向 %%%%%%%%%%%%
        :param lens: 梁长度     单位：mm       %%%%%%%% 坐标系的 x-轴方向 %%%%%%%%%%%%
        :param ela: 弹性模量    单位：Mpa
        :param v: 泊松比
        """
        self.thc = thc
        self.wid = wid
        self.lens = lens
        self.ela = ela
        self.v = v

    def compliance_mat(self):
        eix = self.ela * self.thc * self.wid ** 3 / 12
        eiy = self.ela * self.wid * self.thc ** 3 / 12
        area = self.thc * self.wid
        x = 1 / (2 * (1 + self.v))
        y = 12 * (1 / 3 - 0.21 * self.thc / self.wid * (1 - 1 / 12 * (self.thc / self.wid) ** 4))
        c = np.mat([
            [self.lens / (x * y * eiy), 0, 0, 0, 0, 0],
            [0, self.lens / eiy, 0, 0, 0, -(self.lens ** 2 / (2 * eiy))],
            [0, 0, self.lens / eix, 0, self.lens ** 2 / (2 * eix), 0],
            [0, 0, 0, self.lens / (self.ela * area), 0, 0],
            [0, 0, self.lens ** 2 / (2 * eix), 0, self.lens ** 3 / (3 * eix), 0],
            [0, -(self.lens ** 2 / (2 * eiy)), 0, 0, 0, self.lens ** 3 / (3 * eiy)]
        ])

        return c

    def xy(self):
        xy = np.mat([
            [0, 0,          0, 0],
            [1, -self.lens, 0, 0],
        ])
        return xy

    def kmat(self):  # 构造单元的刚度矩阵
        xy_ele = self.xy()
        ad12 = mf.adjoint_k(
            xy_ele[0, 1] - xy_ele[1, 1],
            xy_ele[0, 2] - xy_ele[1, 2],
            xy_ele[0, 3] - xy_ele[1, 3],
            0, 0, 0)

        k11 = self.compliance_mat().I
        k12 = k11 * ad12
        k22 = ad12.T * k12
        k21 = k22 * ad12.I
        kele = np.vstack((np.hstack((k11, -k12)),
                          np.hstack((-k21, k22))))

        return kele

    def kmat_sparse(self):
        return sp.bsr_matrix(self.kmat(), blocksize=(6, 6))


class CurbeamSQ:
    """
    矩形截面曲梁
    若按顺时针方向画圆弧线为粱的中心线，柔度矩阵的坐标系在弧线的止点
    弧线在止点的最终方向为 x 轴，梁宽为 y-轴 方向，梁厚为 z-轴方向

    假设：x-轴向右，y-轴向上，z-轴垂直纸面向上，则以下：
    ------
    compliance_mat(self):
        弧线的起点接地，求止点的柔度矩阵
    ------
    xy(self):
        梁的两个结点坐标， 止点结点编号为 0 ，坐标为原点； 起点编号为 1
    ------
    kmat_ele(self):
        单元刚度矩阵，由于有两个结点，所以是12 * 12 的矩阵
    ------
    """
    def __init__(self, radl, thc, wid, sita, ela, v):
        """
        :param radl: 粱的弧线的半径        单位：mm
        :param thc: 横截面的厚            单位：mm         ###坐标系 z 轴方向
        :param wid: 横截面的宽            单位：mm         ###坐标系 y 轴方向
        :param sita: 弧线的旋转角度
        :param ela: 弹性模量              单位：Mpa
        :param v: 泊松比
        """
        self.radl = radl
        self.thc = thc
        self.wid = wid
        self.sita = sita
        self.ela = ela
        self.v = v

    def compliance_mat(self):
        g = self.ela / (2 * (1 + self.v))
        iy = 1 / 12 * self.wid * self.thc ** 3
        iz = 1 / 12 * self.thc * self.wid ** 3
        jx = iy + iz
        ax = self.thc * self.wid
        ay = ax * 5 / 6
        az = ay

        eax = self.ela * ax
        gay = g * ay
        gaz = g * az
        eiy = self.ela * iy
        eiz = self.ela * iz
        gjx = g * jx
        sin = math.sin(self.sita)
        cos = math.cos(self.sita)
        eps1 = (2 * self.sita + math.sin(2 * self.sita)) / 4
        eps2 = (2 * self.sita - math.sin(2 * self.sita)) / 4
        eps3 = self.sita - sin

        c = np.zeros((6, 6))
        c[3, 3] = self.radl * eps1 / eax + self.radl * eps2 / gay + self.radl ** 3 * (2 * eps3 - eps2) / eiz
        c[3, 4] = self.radl * sin ** 2 / (2 * eax) - self.radl * sin ** 2 / (2 * gay) - self.radl ** 3 * (
                2 - 2 * cos - sin ** 2) / (2 * eiz)
        c[4, 3] = c[3, 4]
        c[3, 2] = -self.radl ** 2 * eps3 / eiz
        c[2, 3] = c[3, 2]
        c[4, 4] = self.radl * eps2 / eax + self.radl * eps1 / gay + self.radl ** 3 * eps2 / eiz
        c[4, 2] = self.radl ** 2 * (1 - cos) / eiz
        c[2, 4] = c[4, 2]
        c[5, 5] = self.radl * self.sita / gaz + self.radl ** 3 * (2 * eps3 - eps2) / gjx + self.radl ** 3 * eps2 / eiy
        c[5, 0] = -self.radl ** 2 * (sin - eps1) / gjx + self.radl ** 2 * eps2 / eiy
        c[0, 5] = c[5, 0]
        c[5, 1] = -self.radl ** 2 * (2 - 2 * cos - sin ** 2) / (2 * gjx) - self.radl ** 2 * sin ** 2 / (2 * eiy)
        c[1, 5] = c[5, 1]
        c[0, 0] = self.radl * eps1 / gjx + self.radl * eps2 / eiy
        c[0, 1] = self.radl * sin ** 2 / (2 * gjx) - self.radl * sin ** 2 / (2 * eiy)
        c[1, 0] = c[0, 1]
        c[1, 1] = self.radl * eps1 / eiy + self.radl * eps2 / gjx
        c[2, 2] = self.radl * self.sita / eiz

        return np.mat(c)

    def xy(self):
        xy = np.mat([
            [0, 0, 0, 0],
            [1, 0, -self.radl - self.radl * math.cos(self.sita), -self.radl * math.sin(self.sita)],
        ])
        return xy

    def kmat(self):
        xy_ele = self.xy()
        ad12 = mf.adjoint_k(
            xy_ele[0, 1] - xy_ele[1, 1],
            xy_ele[0, 2] - xy_ele[1, 2],
            xy_ele[0, 3] - xy_ele[1, 3],
            0, 0, 0)

        k11 = self.compliance_mat().I
        k12 = k11 * ad12
        k22 = ad12.T * k12
        k21 = k22 * ad12.I
        kele = np.vstack((np.hstack((k11, -k12)),
                          np.hstack((-k21, k22))))

        return kele

    def kmat_sparse(self):
        return sp.bsr_matrix(self.kmat(), blocksize=(6, 6))
