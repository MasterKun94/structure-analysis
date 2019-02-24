"""
这里是所有高级单元组成的类，高级单元由基本单元或者其他高级单元组成
------
"""

import numpy as np

import element_class as ec
import process_functions as pf
import structurebuilder as sb


class PlanarElement(sb.BMStructureBuilder):
    """
    平面单元1
    具体见：
    ？？？？
    """
    def __init__(self, thc, wid, ela, v, ele_sita, ele_zita, ele_d):
        """

        :param thc: 横截面的厚
        :param wid: 横截面的宽
        :param ela: 弹性模量
        :param v: 泊松比
        :param ele_sita:
        :param ele_zita:
        :param ele_d:
        """

        ssita, csita = np.sin(ele_sita), np.cos(ele_sita)
        # print(np.tan(ele_zita))
        l0 = ele_d / 2 * (csita + ssita / np.tan(ele_zita))
        # print(l0)
        ls, lc = l0 * ssita, l0 * csita
        straightbeamsq = ec.StrbeamSQ(thc, wid, l0, ela, v)
        elb = [{
            'name': 0,
            'class': straightbeamsq,
            'rotation': np.array([
                [0, 0, np.pi - ele_sita],
                [0, 0, np.pi / 2 - ele_sita],
                [0, 0, -ele_sita],
                [0, 0, -np.pi / 2 - ele_sita]
            ]),
            'elb': np.array([
                [0, 4],
                [1, 5],
                [2, 6],
                [3, 7]
            ])
        }]

        xy = np.array([
            [0, 0, 0, 0],
            [1, ele_d / 2, ele_d / 2, 0],
            [2, 0, ele_d, 0],
            [3, ele_d / 2, -ele_d / 2, 0],
            [4, lc, -ls, 0],
            [5, ele_d / 2 - ls, ele_d / 2 - lc, 0],
            [6, ele_d - lc, ls, 0],
            [7, ele_d / 2 + ls, -ele_d / 2 + lc, 0],
        ])
        super(PlanarElement, self).__init__(elb, rigid_connect=([4, 5, 6, 7], np.array([0, 0, 0])), redundant_node=[4],
                                            xy=xy, shape=(48, 48))


if __name__ == '__main__':
    pele = PlanarElement(1, 1, 2000, 0.3, np.pi/6, np.pi/2, ele_d=20)
    pele.kmat_sparse()
    # print(pele.elb)
    pf.save_exl(pele.kmat_sparse().todense(), 'C:\work\python\structure1\Try2.xlsx', magnitude=1e-10)
