"""
这里是所有结构组成的类，由一个或者多个基本单元或者高级单元组成

        elb = [{'name': 0,
                'class': straightbeamsq,
                'rotation': np.array([[0, 0, math.pi - sita],
                                      [0, 0, math.pi / 2 - sita],
                                      [0, 0, -sita],
                                      [0, 0, -math.pi / 2 - sita]]),
                'elb': np.array([[0, 4],
                                 [1, 5],
                                 [2, 6],
                                 [3, 7]])
                }]

        elb = [{'name': 0, 'class': straightbeamsq, 'rotation': [0, 0, math.pi - sita], 'elb':[[0, 4]]},
               {'name': 1, 'class': straightbeamsq, 'rotation': [0, 0, math.pi/2 - sita], 'elb':[[1, 5]]},
               {'name': 2, 'class': straightbeamsq, 'rotation': [0, 0, -sita], 'elb':[[2, 6]]},
               {'name': 3, 'class': straightbeamsq, 'rotation': [0, 0, -math.pi/2 - sita], 'elb':[[3, 7]]}]

------
"""

import numpy as np

import element_highclass as ehc
# import s1_functions as cf
# import sparse_functions as sf
# from process_functions import timecounter
# from process_functions import paramfactory
# from itertools import count
import structurebuilder as sb


# import latticebuilder as lb


class Structure3D(sb.LatStructure):
    def __init__(self, nx, ny, nz, thc, wid, ela, v, ele_d, sitalist):
        nxyz = {'nx': nx, 'ny': ny, 'nz': nz}
        eleparam = {'thc': thc, 'wid': wid, 'ela': ela, 'v': v, 'ele_d': ele_d, 'ele_sita': sitalist,
                    'ele_zita': np.pi / 2}
        super(Structure3D, self).__init__(
            'LatticeWith2DDiamondEle', nxyz, 'PlanarElement', eleparam, ele_d, per_plane=True)
        self.add_rigid(self.lattice.top_node(), feature_node=None)
        self.add_constraint(self.lattice.ground_node())


class PlaneStructure(sb.SPStructureBuilder):
    def __init__(self, nx, ny, thc, wid, ela, v, ele_sita, ele_zita, ele_d):
        elblist = []
        for h in range(nx):
            nyh = (2 * ny + 1) * h + ny
            for l in range(ny):
                nyhl = nyh + l
                elblist += [[
                    nyhl - ny,
                    nyhl + 1,
                    nyhl + ny + 1,
                    nyhl
                ]]

        elb = [{'name': 'PlanarElement',
                'class': ehc.PlanarElement(thc, wid, ela, v, ele_sita, ele_zita, ele_d),
                'rotation': np.array([0, 0, 0]),
                'elb': np.array(elblist)}]
        # print(elb[0]['class'].kmat_sparse())
        self.xy, num, d = [], 0, ele_d / 2
        for h in range(2 * nx + 1):
            for l in range(2 * ny + 1):
                if (h + l) % 2 == 1:
                    self.xy += [[
                        num,
                        h * d,
                        l * d,
                        0
                    ]]
                    num += 1
        rigid_node = list(range(ny))
        feature_xy = np.array([0, 10 * ny, 0])
        super(PlaneStructure, self).__init__(
            elb, rigid_connect=(rigid_node, feature_xy), shape=(num * 6, num * 6), xy=np.array(self.xy))
        self.add_constraint(list(range(len(self.xy) - ny, len(self.xy))))


if __name__ == '__main__':
    # inx, iny, inz, ithc, iwid, iela, iv, iele_sita, iele_zita, iele_d = \
    #     1, 1, 1, 1, 1, 2000, 0.3, np.pi / 6, np.pi / 2, 20
    # #
    # # struct = PlaneStructure(nx, ny, thc, wid, ela, v, ele_sita, ele_zita, ele_d)
    # # kmat = struct.solve_nodecmat(0)
    # # print(kmat.shape)
    # # pf.save_exl(kmat, 'C:\work\python\structure1\kmat.xlsx', magnitude=0)
    #
    # struct = Strcture3(inx, iny, inz, ithc, iwid, iela, iv, iele_d, iele_sita)
    # print(struct)
    pass
