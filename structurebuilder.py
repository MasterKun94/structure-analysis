"""

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

    rigidconnect = ([0, 1, 3, 5, 7], np.array([0.0, 5.2, 1.8]))
        rigidnode = [0, 1, 3, 5, 7]
        featurexy = np.array([0.0, 5.2, 1.8]

    redundantnode = [8, 9, 12]
"""

# import element_class
from scipy.sparse import csr_matrix

import element_highclass as ehc
import latticebuilder as lb
import mat_functions as mf
import sparse_functions as sf
from process_functions import time_counter, param_factory


class SuperBuilder:

    def __init__(self, ele_list, xy, kmat, method):
        if ele_list is None:
            self.ele_list = []
        else:
            self.ele_list = ele_list
        if xy is None:
            # self.xy = 'build_xy(elb) working on'
            pass
        else:
            self.xy = xy
        self.kmat = kmat
        self.method = method

    def __str__(self):

        return (f'StructureBuilder - {self.method!r} matrix method\n'
                f'  Structure name : {self.__class__.__name__}\n'
                f'     Node number : {self.get_node_number()!r} \n'
                f'  Element number : {self.get_element_number()!r} \n'
                f'      Kmat shape : {self.kmat.shape!r} \n'
                )

    def __repr__(self):
        return (f'StructureBuilder - {self.method} matrix method\n'
                f'  Structure name : {self.__class__.__name__}\n')

    def get_element_number(self):
        try:
            return self.element_number
        except AttributeError:
            number = 0
            for ele in self.ele_list:
                number += len(ele['elb'])
            return number

    def get_node_number(self):
        try:
            return self.node_number
        except AttributeError:
            return self.xy.shape[0]

# *****************************


class BMStructureBuilder(SuperBuilder):

    def __init__(self, ele_list, rigid_connect=None, redundant_node=None, shape=None, xy=None):
        super(BMStructureBuilder, self).__init__(ele_list, xy, mf.build_kmat(ele_list, shape=shape), 'BASIC')

        if rigid_connect:
            rigid_node, feature_xy = rigid_connect
            self.xy, self.kmat = mf.rigidelement_k(self.xy, feature_xy, self.kmat, rigid_node)
        if redundant_node:
            self.xy, self.kmat = mf.redundantnode_simplify(self.xy, self.kmat, redundant_node)

    def add_rigid(self, new_rigid_connect, out=False):
        rigid_node, feature_xy = new_rigid_connect
        if out is True:
            return mf.rigidelement_k(self.xy, feature_xy, self.kmat, rigid_node)
        elif out is False:
            self.xy, self.kmat = mf.rigidelement_k(self.xy, feature_xy, self.kmat, rigid_node)

    def add_redundant(self, new_redundant_node, out=False):
        if out is True:
            return mf.redundantnode_simplify(self.xy, self.kmat, new_redundant_node)
        elif out is False:
            self.xy, self.kmat = mf.redundantnode_simplify(self.xy, self.kmat, new_redundant_node)

    def add_element(self, new_ele_list):
        self.kmat = mf.build_kmat(new_ele_list, self.kmat, shape=None)
        # self.xy =
        self.ele_list += new_ele_list

    def add_constraint(self, constraint_node):
        return mf.build_constriand(self.kmat, constraint_node)

    def solve_deformation(self):
        pass

    def solve_node_cmat(self):
        pass

    def kmat_sparse(self):
        # print('using BM')
        return csr_matrix(self.kmat)


class SPStructureBuilder(SuperBuilder):

    def __init__(self, ele_list, rigid_connect=None, redundant_node=None, shape=None, xy=None):
        super(SPStructureBuilder, self).__init__(ele_list, xy, sf.build_kmat(ele_list, shape=shape), 'SPARSE')
        if rigid_connect:
            rigid_node, feature_xy = rigid_connect
            self.kmat = sf.rigid(self.xy, feature_xy, self.kmat, rigid_node)
            self.xy = sf.rigid_xy(self.xy, feature_xy, rigid_node)
        if redundant_node:
            self.xy, self.kmat = sf.redundant(self.xy, self.kmat, redundant_node)

    def add_rigid(self, new_rigid_connect, out=False, reduce_xy=False):
        rigid_node, feature_xy = new_rigid_connect
        if feature_xy is None:
            feature_xy = sum(self.xy[rigid_node])[1:4] / len(rigid_node)
        if out is True:
            if reduce_xy is True:
                return sf.rigid_xy(
                    self.xy, feature_xy, rigid_node), sf.rigid(self.xy, feature_xy, self.kmat, rigid_node)
            else:
                return sf.rigid(self.xy, feature_xy, self.kmat, rigid_node)
        elif out is False:
            self.kmat = sf.rigid(self.xy, feature_xy, self.kmat, rigid_node)
            if reduce_xy is True:
                self.xy = sf.rigid_xy(self.xy, feature_xy, rigid_node)

    def add_redundant(self, new_redundant_node: list, out=False):
        if out is True:
            return sf.redundant(self.xy, self.kmat, new_redundant_node)
        elif out is False:
            self.xy, self.kmat = sf.redundant(self.xy, self.kmat, new_redundant_node)

    def add_element(self, new_ele_list: "list or generator", shape=None):
            self.kmat = sf.build_kmat(new_ele_list, self.kmat, shape=shape)
            self.ele_list += new_ele_list

    def add_constraint(self, constraint_node: list, out=False):
        if out is True:
            return sf.build_constraint(self.kmat, constraint_node)
        elif out is False:
            self.kmat = sf.build_constraint(self.kmat, constraint_node)

    @time_counter
    def solve_deformation(self):
        print('working on')
        return None

    @time_counter
    def solve_node_cmat(self, request_node, constraint_node=None):
        if constraint_node is None:
            # try:
            return sf.sp_inv(self.kmat, request_node)
            # except:
            #     print(self.__str__(), '未施加约束')
        else:
            return sf.sp_inv(self.add_constraint(constraint_node, out=True), request_node)

    def kmat_sparse(self):
        # print('using SP')
        return self.kmat.tocsr()


class LatStructure(SPStructureBuilder):
    @time_counter
    @param_factory
    def __init__(self, lat_name: str, lat_param: dict, ele_name: str, ele_param: dict,
                 ele_size: float, per_plane=False, per_direction=False):
        def _gen_ele():
            dic = self.lattice.generater_elementdic(per_plane=per_plane, per_direction=per_direction)
            ele_class = getattr(ehc, ele_name)
            next(dic)
            for param in ele_param:
                try:
                    yield dic.send(ele_class(**param))
                except StopIteration:
                    break

        self.lattice = getattr(lb, lat_name)(lat_param)
        self.node_number = self.lattice.node_number()
        super(LatStructure, self).__init__(ele_list=[], shape=(self.node_number * 6, self.node_number * 6))
        self.add_element(_gen_ele())
        self.xy = self.lattice.build_xy(ele_size)
        self.element_number = self.lattice.element_number()
        self.method = "SPARSE LATTICE"

# ----------------------------------------------
#
#
# class ElementList:
#     def __init__(self, ele_list=None, **kwargs):
#         kwargs = self.check(kwargs)
#         if ele_list is None or []:
#             self.ele_list = kwargs
#         else:
#             if isinstance(ele_list, list):
#                 self.ele_list = ele_list + kwargs
#             elif isinstance(ele_list, dict):
#                 self.ele_list = [ele_list] + kwargs
#
#     def __repr__(self):
#         return 'ElementList\n', self.ele_list
#
#     def __len__(self):
#         return len(self.ele_list)
#
#     def __getitem__(self, item):
#         try:
#             row, key = item
#             if isinstance(row, int):
#                 return self.ele_list[row][key]
#             elif isinstance(row, slice):
#                 elb = []
#                 for small_list in self.ele_list[row]:
#                     elb += small_list[key]
#                 return elb
#         except TypeError:
#             cls = type(self)
#             if isinstance(item, int) or isinstance(item, slice):
#                 return cls(self.ele_list[item])
#             elif isinstance(item, str):
#                 warnings.warn('still working on')
#                 return None
#
#     def __add__(self, other):
#         cls = type(self)
#         if isinstance(other, list):
#             return cls(self.ele_list + other)
#         elif isinstance(other, dict):
#             return cls(self.ele_list + [other])
#         else:
#             try:
#                 return cls(self.ele_list + other.ele_list)
#
#             except TypeError:
#                 raise TypeError('参数必须是字典，列表或者 ElementList 类')
#
#     def __radd__(self, other):
#         return self + other
#
#     def __iadd__(self, other):
#         return self + other
#
#     def __iter__(self):
#         return iter(self.ele_list)
#
#     def max_node(self):
#         maxnode = 0
#         for ele in self:
#             maxnode = max(np.max(np.array(ele['elb'])), maxnode)
#         return maxnode
#
#     def check(self, dic):
#         if dic is not None:
#             if 'class' in dic and 'rotation' in dic and 'elb' in dic:
#                 return [dic]
#             else:
#                 warnings.warn('参数信息不完整')
#                 return None
#         else:
#             return None
#
#     def append(self, ele_list=None, **kwargs):
#         kwargs = self.check(kwargs)
#         if ele_list is None or []:
#             self.ele_list += kwargs
#         else:
#             if isinstance(ele_list, list):
#                 self.ele_list += ele_list + kwargs
#             elif isinstance(ele_list, dict):
#                 self.ele_list += [ele_list] + kwargs


if __name__ == '__main__':

    pass
