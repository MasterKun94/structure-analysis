from itertools import count

import numpy as np


class LatticeBuilder:

        pass


class LatticeWith2DDiamondEle(LatticeBuilder):
    def __init__(self, nxyz):
        self.nx = nxyz['nx']
        self.ny = nxyz['ny']
        self.nz = nxyz['nz']

    def build_xy(self, element_zize):
        ixy, inum, eld = [], count(0, 1), element_zize / 2
        for k in range(0, 2 * self.nz + 1):
            for i in range(0, 2 * self.nx + 1):
                for j in range(0, 2 * self.ny + 1):
                    if (i % 2 + j % 2 + k % 2) == 1:
                        ixy += [[next(inum), eld * i, eld * j, -eld * k]]
        return np.array(ixy)

    def generater_elb(self, axis, per_plane=False, per_direction=False):
        def map_node(i, j, k):
            if k % 2 == 0:
                return (3 * self.nx * self.ny + 2 * self.nx + 2 * self.ny + 1) * (k // 2) + \
                       self.ny * ((i + 1) // 2) + (self.ny + 1) * (i // 2) + j // 2
            else:
                return ((self.nx + 1) * self.ny + (self.ny + 1) * self.nx) * ((k + 1) // 2) + \
                       (self.nx + 1) * (self.ny + 1) * ((k - 1) // 2) + (self.ny + 1) * (i // 2) + j // 2

        def gen_elelstx():
            ii = 2 * self.nx - 2 * row
            for line in range(self.nz):
                kk = 2 * line
                for col in range(self.ny):
                    jj = 2 * col
                    yield [map_node(ii, jj + 1, kk),
                           map_node(ii, jj + 2, kk + 1),
                           map_node(ii, jj + 1, kk + 2),
                           map_node(ii, jj, kk + 1)]

        def gen_elelsty():
            jj = 2 * self.ny - 2 * row
            for line in range(self.nz):
                kk = 2 * line
                for col in range(self.nx):
                    ii = 2 * self.nx - 2 * col
                    yield [map_node(ii - 1, jj, kk),
                           map_node(ii - 2, jj, kk + 1),
                           map_node(ii - 1, jj, kk + 2),
                           map_node(ii, jj, kk + 1)]

        def gen_elelstz():
            kk = 2 * row
            for line in range(self.nx):
                ii = 2 * line
                for col in range(self.ny):
                    jj = 2 * col
                    yield [map_node(ii, jj + 1, kk),
                           map_node(ii + 1, jj + 2, kk),
                           map_node(ii + 2, jj + 1, kk),
                           map_node(ii + 1, jj, kk)]

        if per_direction is True:
            ele_list = []
            if axis == 0:
                for row in range(self.nx + 1):
                    ele_list += list(gen_elelstx())
                yield ele_list
            elif axis == 1:
                for row in range(self.ny + 1):
                    ele_list += list(gen_elelsty())
                yield ele_list
            elif axis == 2:
                for row in range(self.nz + 1):
                    ele_list += list(gen_elelstz())
                yield ele_list

        elif per_plane is True:
            if axis == 0:
                for row in range(self.nx + 1):
                    yield list(gen_elelstx())
            elif axis == 1:
                for row in range(self.ny + 1):
                    yield list(gen_elelsty())
            elif axis == 2:
                for row in range(self.nz + 1):
                    yield list(gen_elelstz())

        else:
            if axis == 0:
                for row in range(self.nx + 1):
                    for iele in gen_elelstx():
                        yield [iele]
            elif axis == 1:
                for row in range(self.ny + 1):
                    for iele in gen_elelsty():
                        yield [iele]
            elif axis == 2:
                for row in range(self.nz + 1):
                    for iele in gen_elelstz():
                        yield [iele]

    def generater_elementdic(self, axis=None, name=None, per_plane=False, per_direction=False):
        element = yield None
        if axis:
            for ielb in self.generater_elb(axis, per_plane, per_direction):
                element = yield {'name': name,
                                 'class': element,
                                 'rotation': [[0, np.pi / 2, 0], [-np.pi / 2, 0, np.pi / 2], [0, 0, 0]][int(axis)],
                                 'elb': ielb}

        else:
            for i in range(0, 3):
                for ielb in self.generater_elb(i, per_plane, per_direction):
                    element = yield {'name': name,
                                     'class': element,
                                     'rotation': [[0, np.pi / 2, 0], [-np.pi / 2, 0, np.pi / 2], [0, 0, 0]][int(i)],
                                     'elb': ielb}

    def build_elementlist(self):
        pass

    def element_number(self):
        x, y, z = self.nx, self.ny, self.nz
        return 3 * (x * y * z) + x * y + y * z + x * z

    def node_number(self):
        x, y, z = self.nx, self.ny, self.nz
        return 3 * x * y * z + 2 * (x * y + x * z + y * z) + x + y + z

    def top_node(self):
        b = (self.ny + 1) * self.nx + (self.nx + 1) * self.ny
        return list(range(b))

    def ground_node(self):
        a = (3 * self.nx * self.ny + 2 * self.nx + 2 * self.ny + 1) * self.nz
        b = (self.ny + 1) * self.nx + (self.nx + 1) * self.ny
        return list(range(a - b + 1, a + 1))


if __name__ == '__main__':

    pass
