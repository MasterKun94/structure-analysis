import math

import numpy as np

import mat_functions as mf


def c_strbeam_cc(r, l, e, v):  # 圆形直梁的柔度矩阵
    # r ： 横截面的半径
    # l ： 梁的长度  ######坐标系的z轴方向#####
    # e ： 弹性模量
    # v ： 泊松比

    ix = math.pi * r ** 4 / 4
    j = 2 * ix
    a = math.pi * r ** 2
    g = e / (2 * (1 + v))
    c = np.mat([
        [l / (g * j),  0,            0,           0,                      0,                      0],
        [0,            l / (e * ix), 0,           0,                      0,                      0],
        [0,            0,            l / (e * ix),0,                      0,                      0],
        [0,            0,            0,           l / (e * a), 0,                      0],
        [0,            0,            0,           0,                      l ** 3 / (12 * e * ix), 0],
        [0,            0,            0,           0,                      0,            l ** 3 / (12 * e * ix)]
    ])
    ad = mf.adjoint(0, 0, -l/2, 0, 0, 0)
    c = ad * c * ad.T

    return c


def xy__strbeam_cc(l):
    xy = np.mat([
        [0, 0, 0, 0],
        [1, 0, 0, -l],
    ])
    return xy


def kele_strbeam_cc(r, l, e, v):
    c = c_strbeam_cc(r, l, e, v)
    ad = mf.adjoint(0, 0, l, 0, 0, 0)


def c_strbeam_sq(t, w, l, e, v):  # 矩形直梁的柔度矩阵
    # t ： 横截面的厚  ###坐标系 x 轴方向
    # w ： 横截面的宽  ###坐标系 y 轴方向
    # l ： 梁的长度   ###坐标系 z 轴方向
    # E ： 弹性模量
    # v ： 泊松比

    eix = e * t * w ** 3 / 12
    eiy = e * w * t ** 3 / 12
    a = t * w
    x = 1 / (2 * (1 + v))
    y = 12 * (1 / 3 - 0.21 * t / w * (1 - 1 / 12 * (t / w) ** 4))
    c = np.mat([
        [l / eix,               0,                  0,                 0,                  -(l ** 2 / (2 * eix)), 0],
        [0,                     l / eiy,            0,                 l ** 2 / (2 * eiy), 0,                     0],
        [0,                     0,                  l / (x * y * eiy), 0,                  0,                     0],
        [0,                     l ** 2 / (2 * eiy), 0,                 l**3 / (3 * eiy),   0,                     0],
        [-(l ** 2 / (2 * eix)), 0,                  0,                 0,                  l**3 / (3 * eix),      0],
        [0,                     0,                  0,                 0,                  0,           l / (e * a)]
    ])

    return c


def xy__strbeam_sq(l):
    xy = np.mat([
        [0, 0, 0, 0],
        [1, 0, 0, -l],
    ])
    return xy


def c_curbeam_sq(r, t, w, sita, e, v):  # 弯曲矩形粱的柔度矩阵
    # 若按顺时针方向画圆弧线为粱的中心线，柔度矩阵的坐标系在弧线的止点
    # 弧线在止点的最终方向为 z 轴
    # r ： 粱中心线的半径
    # t ： 横截面的厚  ###坐标系 x 轴方向
    # w ： 横截面的宽  ###坐标系 y 轴方向
    # E ： 弹性模量
    # v ： 泊松比

    g = e / (2 * (1 + v))
    iy = 1 / 12 * w * t ** 3
    iz = 1 / 12 * t * w ** 3
    jx = iy + iz
    ax = t * w
    ay = ax * 5 / 6
    az = ay

    eax = e * ax
    gay = g * ay
    gaz = g * az
    eiy = e * iy
    eiz = e * iz
    gjx = g * jx
    sin = math.sin(sita)
    cos = math.cos(sita)
    eps1 = (2 * sita + math.sin(2 * sita)) / 4
    eps2 = (2 * sita - math.sin(2 * sita)) / 4
    eps3 = sita - sin

    c = np.zeros((6, 6))
    c[3, 3] = r * eps1 / eax + r * eps2 / gay + r ** 3 * (2 * eps3 - eps2) / eiz
    c[3, 4] = r * sin ** 2 / (2 * eax) - r * sin ** 2 / (2 * gay) - r ** 3 * (2 - 2 * cos - sin ** 2) / (2 * eiz)
    c[4, 3] = c[3, 4]
    c[3, 2] = -r ** 2 * eps3 / eiz
    c[2, 3] = c[3, 2]
    c[4, 4] = r * eps2 / eax + r * eps1 / gay + r ** 3 * eps2 / eiz
    c[4, 2] = r ** 2 * (1 - cos) / eiz
    c[2, 4] = c[4, 2]
    c[5, 5] = r * sita / gaz + r ** 3 * (2 * eps3 - eps2) / gjx + r ** 3 * eps2 / eiy
    c[5, 0] = -r ** 2 * (sin - eps1) / gjx + r ** 2 * eps2 / eiy
    c[0, 5] = c[5, 0]
    c[5, 1] = -r ** 2 * (2 - 2 * cos - sin ** 2) / (2 * gjx) - r ** 2 * sin ** 2 / (2 * eiy)
    c[1, 5] = c[5, 1]
    c[0, 0] = r * eps1 / gjx + r * eps2 / eiy
    c[0, 1] = r * sin ** 2 / (2 * gjx) - r * sin ** 2 / (2 * eiy)
    c[1, 0] = c[0, 1]
    c[1, 1] = r * eps1 / eiy + r * eps2 / gjx
    c[2, 2] = r * sita / eiz

    return c


def xy__curbeam_sq(r, sita):
    xy = np.mat([
        [0, 0, 0, 0],
        [1, 0, r - r * math.cos(sita), -r * math.sin(sita)],
    ])