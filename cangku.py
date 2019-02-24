# class PlanarElement1:
#     """
#     平面单元1
#     具体见：
#     ？？？？
#     """
#     def __init__(self, thc, wid, ela, v, ele_sita, ele_d):
#         """
#
#         :param thc: 横截面的厚
#         :param wid: 横截面的宽
#         :param ela: 弹性模量
#         :param v: 泊松比
#         :param ele_sita:
#         :param ele_d:
#         """
#         ssita, csita = np.sin(ele_sita), np.cos(ele_sita)
#         l0 = ele_d * csita / 2
#         ls, lc = l0 * ssita, l0 * csita
#         straightbeamsq = ec.StrbeamSQ(thc, wid, l0, ela, v)
#         elb = [
#             [0, 0, 0, np.pi - ele_sita, 0, 4],
#             [1, 0, 0, np.pi/2 - ele_sita, 1, 5],
#             [2, 0, 0, -ele_sita, 2, 6],
#             [3, 0, 0, -np.pi/2 - ele_sita, 3, 7],
#         ]
#         xy = np.array([
#             [0, 0, 0, 0],
#             [1, ele_d / 2, ele_d / 2, 0],
#             [2, 0, ele_d, 0],
#             [3, ele_d / 2, -ele_d / 2, 0],
#             [4, lc, -ls, 0],
#             [5, ele_d / 2 - ls, ele_d / 2 - lc, 0],
#             [6, ele_d - lc, ls, 0],
#             [7, ele_d / 2 + ls, -ele_d / 2 + lc, 0],
#         ])
#         kele = mf.build_kmat(elb, xy, straightbeamsq.kmat_ele())
#         kele = mf.rigidelement_k(xy, np.mat([[0, -ele_d/2, 0]]), kele, [4, 5, 6, 7])
#         self.kele = mf.redundantnode_simplify(kele[0], kele[1], [4])
#
#     def kmat_kele(self):
#
#         return self.kele[1]
#
#     def xy(self):
#
#         return self.kele[0]
#
#     def kmat_sparse(self):
#
#         return sparse.csr_matrix(self.kele[1])
#
#
# class PlanarElement1Update:
#     """
#     平面单元1
#     具体见：
#     ？？？？
#     """
#     def __init__(self, thc, wid, ela, v, ele_sita, ele_d):
#         """
#
#         :param thc: 横截面的厚
#         :param wid: 横截面的宽
#         :param ela: 弹性模量
#         :param v: 泊松比
#         :param ele_sita:
#         :param ele_d:
#         """
#
#         self.thc = thc
#         self.wid = wid
#         self.ela = ela
#         self.v = v
#         self.ele_sita = ele_sita
#         self.ele_d = ele_d
#         # self.kele = self.ele()
#
#     def ele(self):
#         ssita, csita = np.sin(self.ele_sita), np.cos(self.ele_sita)
#         l0 = self.ele_d * csita / 2
#         ls, lc = l0 * ssita, l0 * csita
#         straightbeamsq = ec.StrbeamSQ(self.thc, self.wid, l0, self.ela, self.v)
#         elb = [{'name': 0,
#                 'class': straightbeamsq,
#                 'rotation': np.array([
#                     [0, 0, np.pi - self.ele_sita],
#                     [0, 0, np.pi / 2 - self.ele_sita],
#                     [0, 0, -self.ele_sita],
#                     [0, 0, -np.pi / 2 - self.ele_sita]
#                 ]),
#                 'elb': np.array([
#                     [0, 4],
#                     [1, 5],
#                     [2, 6],
#                     [3, 7]
#                 ])
#                 }]
#         # elb = [
#         #     {'name': 0, 'class': straightbeamsq, 'rotation': [0, 0, math.pi - self.ele_sita], 'elb':[[0, 4]]},
#         #     {'name': 1, 'class': straightbeamsq, 'rotation': [0, 0, math.pi/2 - self.ele_sita], 'elb':[[1, 5]]},
#         #     {'name': 2, 'class': straightbeamsq, 'rotation': [0, 0, -self.ele_sita], 'elb':[[2, 6]]},
#         #     {'name': 3, 'class': straightbeamsq, 'rotation': [0, 0, -math.pi/2 - self.ele_sita], 'elb':[[3, 7]]}
#         # ]
#         xy = np.array([
#             [0, 0, 0, 0],
#             [1, self.ele_d / 2, self.ele_d / 2, 0],
#             [2, 0, self.ele_d, 0],
#             [3, self.ele_d / 2, -self.ele_d / 2, 0],
#             [4, lc, -ls, 0],
#             [5, self.ele_d / 2 - ls, self.ele_d / 2 - lc, 0],
#             [6, self.ele_d - lc, ls, 0],
#             [7, self.ele_d / 2 + ls, -self.ele_d / 2 + lc, 0],
#         ])
#         kele = mf.build_kmat(elb, xy, straightbeamsq.kmat_ele())
#         kele = mf.rigidelement_k(xy, np.mat([[0, -self.ele_d/2, 0]]), kele, [4, 5, 6, 7])
#         kele = mf.redundantnode_simplify(kele[0], kele[1], [4])
#         return kele
#
#     def kmat_kele(self):
#         kele = self.ele()[1]
#         return kele
#
#     def xy(self):
#         kele = self.ele()[0]
#         return kele
#
#     def kmat_sparse(self):
#         kele = self.ele()[1]
#         return sparse.csr_matrix(kele)
#
#
# class PlanarElement2:
#     """
#     平面单元1
#     具体见：
#     ？？？？
#     """
#     def __init__(self, thc, wid, ela, v, ele_sita, ele_zita, ele_d):
#         """
#
#         :param thc: 横截面的厚
#         :param wid: 横截面的宽
#         :param ela: 弹性模量
#         :param v: 泊松比
#         :param ele_sita:
#         :param ele_zita:
#         :param ele_d:
#         """
#
#         self.thc = thc
#         self.wid = wid
#         self.ela = ela
#         self.v = v
#         self.ele_sita = ele_sita
#         self.ele_zita = ele_zita
#         self.ele_d = ele_d
#         # self.kele = self.ele()
#
#     def ele(self):
#         ssita, csita = np.sin(self.ele_sita), np.cos(self.ele_sita)
#         l0 = self.ele_d / 2 * (csita + ssita / np.tan(self.ele_zita))
#         ls, lc = l0 * ssita, l0 * csita
#         straightbeamsq = ec.StrbeamSQ(self.thc, self.wid, l0, self.ela, self.v)
#         elb = [{'name': 0,
#                 'class': straightbeamsq,
#                 'rotation': np.array([
#                     [0, 0, np.pi - self.ele_sita],
#                     [0, 0, np.pi / 2 - self.ele_sita],
#                     [0, 0, -self.ele_sita],
#                     [0, 0, -np.pi / 2 - self.ele_sita]
#                 ]),
#                 'elb': np.array([
#                     [0, 4],
#                     [1, 5],
#                     [2, 6],
#                     [3, 7]
#                 ])
#                 }]
#         # elb = [
#         #     {'name': 0, 'class': straightbeamsq, 'rotation': [0, 0, math.pi - self.ele_sita], 'elb':[[0, 4]]},
#         #     {'name': 1, 'class': straightbeamsq, 'rotation': [0, 0, math.pi/2 - self.ele_sita], 'elb':[[1, 5]]},
#         #     {'name': 2, 'class': straightbeamsq, 'rotation': [0, 0, -self.ele_sita], 'elb':[[2, 6]]},
#         #     {'name': 3, 'class': straightbeamsq, 'rotation': [0, 0, -math.pi/2 - self.ele_sita], 'elb':[[3, 7]]}
#         # ]
#         xy = np.array([
#             [0, 0, 0, 0],
#             [1, self.ele_d / 2, self.ele_d / 2, 0],
#             [2, 0, self.ele_d, 0],
#             [3, self.ele_d / 2, -self.ele_d / 2, 0],
#             [4, lc, -ls, 0],
#             [5, self.ele_d / 2 - ls, self.ele_d / 2 - lc, 0],
#             [6, self.ele_d - lc, ls, 0],
#             [7, self.ele_d / 2 + ls, -self.ele_d / 2 + lc, 0],
#         ])
#         kele = mf.build_kmat(elb, xy, straightbeamsq.kmat_ele())
#         kele = mf.rigidelement_k(xy, np.mat([[0, -self.ele_d/2, 0]]), kele, [4, 5, 6, 7])
#         kele = mf.redundantnode_simplify(kele[0], kele[1], [4])
#         return kele
#
#     def kmat_kele(self):
#         kele = self.ele()[1]
#         return kele
#
#     def xy(self):
#         kele = self.ele()[0]
#         return kele
#
#     def kmat_sparse(self):
#         kele = self.ele()[1]
#         return sparse.csr_matrix(kele)
#
#
# class PlanarElement2Update:
#     """
#     平面单元1
#     具体见：
#     ？？？？
#     """
#     def __init__(self, thc, wid, ela, v, ele_sita, ele_zita, ele_d):
#         """
#
#         :param thc: 横截面的厚
#         :param wid: 横截面的宽
#         :param ela: 弹性模量
#         :param v: 泊松比
#         :param ele_sita:
#         :param ele_zita:
#         :param ele_d:
#         """
#
#         # self.thc = thc
#         # self.wid = wid
#         # self.ela = ela
#         # self.v = v
#         # self.ele_sita = ele_sita
#         # self.ele_zita = ele_zita
#         # self.ele_d = ele_d
#         # # self.kele = self.ele()
#
#     # def ele(self):
#         ssita, csita = np.sin(ele_sita), np.cos(ele_sita)
#         l0 = ele_d / 2 * (csita + ssita / np.tan(ele_zita))
#         ls, lc = l0 * ssita, l0 * csita
#         straightbeamsq = ec.StrbeamSQ(thc, wid, l0, ela, v)
#         self.elb = [{
#             'name': 0,
#             'class': straightbeamsq,
#             'rotation': np.array([
#                 [0, 0, np.pi - ele_sita],
#                 [0, 0, np.pi / 2 - ele_sita],
#                 [0, 0, -ele_sita],
#                 [0, 0, -np.pi / 2 - ele_sita]
#             ]),
#             'elb': np.array([
#                 [0, 4],
#                 [1, 5],
#                 [2, 6],
#                 [3, 7]
#             ])
#         }]
#
#         self.xy = np.array([
#             [0, 0, 0, 0],
#             [1, ele_d / 2, ele_d / 2, 0],
#             [2, 0, ele_d, 0],
#             [3, ele_d / 2, -ele_d / 2, 0],
#             [4, lc, -ls, 0],
#             [5, ele_d / 2 - ls, ele_d / 2 - lc, 0],
#             [6, ele_d - lc, ls, 0],
#             [7, ele_d / 2 + ls, -ele_d / 2 + lc, 0],
#         ])
#
#         self.structure = sb.BMStructureBuilder(
#             self.elb, ([4, 5, 6, 7], np.array([0, 0, 0])), [4], xy=self.xy, shape=(48, 48))
#
#     def init_structure(self):
#
#         self.structure = sb.BMStructureBuilder(
#             self.elb, xy=self.xy, shape=(48, 48))
#
#     def kmat_sparse(self):
#         try:
#             return sparse.csr_matrix(self.structure.kmat)
#         except AttributeError:
#             self.init_structure()
#             return sparse.csr_matrix(self.structure.kmat)
#
#     def add_rigid(self, newrigid):
#         try:
#             self.structure.add_rigid(newrigid)
#         except AttributeError:
#             self.init_structure()
#             self.structure.add_rigid(newrigid)
#
#     def add_redundant(self, newredundant):
#         try:
#             self.structure.add_redundant(newredundant)
#         except AttributeError:
#             self.init_structure()
#             self.structure.add_redundant(newredundant)
#
#
# # ------------------------------------------------
#
#
# class Strcture1:
#     def __init__(self, nx, ny, nz, thc, wid, ela, v, ele_d, sitalist):
#
#         self.xy = cf.build_xy_s1(ele_d, nx, ny, nz)
#         self.featurexy = np.mat([[ele_d * nx / 2, ele_d * ny / 2, 0]])
#         self.nx = nx
#         self.ny = ny
#         self.nz = nz
#         self.thc = thc
#         self.wid = wid
#         self.ela = ela
#         self.v = v
#         self.ele_d = ele_d
#         self.sitalist = sitalist
#
#     def kmat_beforerigid_sparse(self):
#         time_start = time.time()
#         print('runnung srtuctureclass.kmat_beforerigid_sparse()')
#
#         num = self.xy.shape[0] * 6
#         kmat = None
#
#         for sitax in range(0, len(self.sitalist[0])):
#             planar_ele = ehc.PlanarElement1Update(
#                 self.thc, self.wid, self.ela, self.v, self.sitalist[0][sitax], self.ele_d)
#             pelb = cf.build_pelb_sp(sitax, self.nx, self.ny, self.nz, 'x', planar_ele)
#             kmat = sf.build_kmat(pelb, kmat, shape=(num, num))
#             # pf.save_exl(kmat.todense(), 'C:\work\python\structure1\kkk.xlsx')
#         for sitay in range(0, len(self.sitalist[1])):
#             planar_ele = ehc.PlanarElement1Update(
#                 self.thc, self.wid, self.ela, self.v, self.sitalist[1][sitay], self.ele_d)
#             pelb = cf.build_pelb_sp(sitay, self.nx, self.ny, self.nz, 'y', planar_ele)
#             kmat = sf.build_kmat(pelb, kmat, shape=(num, num))
#         for sitaz in range(0, len(self.sitalist[2])):
#             planar_ele = ehc.PlanarElement1Update(
#                 self.thc, self.wid, self.ela, self.v, self.sitalist[2][sitaz], self.ele_d)
#             pelb = cf.build_pelb_sp(sitaz, self.nx, self.ny, self.nz, 'z', planar_ele)
#             kmat = sf.build_kmat(pelb, kmat, shape=(num, num))
#         time_end = time.time()
#         print('Structure1 inittime cost: ', time_end - time_start, '\n -------------------------')
#         return kmat
#
#     def kmat_beforerigid(self):
#         kmat = self.kmat_beforerigid_sparse()
#         return kmat.todense()
#
#     def kmat_afterrigid_sparse(self):
#         kmat = self.kmat_beforerigid_sparse()
#         xy = self.xy
#         featurexy = self.featurexy
#         rignd = list(range(0, (self.ny + 1) * self.nx + (self.ny + 1) * self.nx))
#         kmat = sf.rigid(xy, featurexy, kmat, rignd)
#         return kmat.todense()
#
#     def kmat_afterrigid(self):
#         kmat = self.kmat_afterrigid_sparse()
#         return kmat.todense()
#
#     def cmat_featurenode(self):
#         a = (3 * self.nx * self.ny + 2 * self.nx + 2 * self.ny + 1) * self.nz
#         b = (self.ny + 1) * self.nx + (self.nx + 1) * self.ny
#         rignd = list(range(0, b))
#
#         kmat = sf.sp_inv(
#             sf.build_constriand(
#                 sf.rigid(
#                     self.xy, self.featurexy, self.kmat_beforerigid_sparse(), rignd),
#                 np.array(list(range(a - len(rignd) + 1, a + b - len(rignd) + 1)))),
#             0)
#
#         # kmat = self.kmat_beforerigid_sparse()
#         # kmat = sf.rigid(self.xy, self.featurexy, kmat, rignd)
#         # kmat = sf.build_constriand(kmat, np.array(list(range(a - len(rignd) + 1, a + b - len(rignd) + 1))))
#         # kmat = sf.sp_inv(kmat, 0)
#         return kmat