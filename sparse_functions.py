"""

"""

import copy
import sys
import warnings

import numpy as np
from scipy import sparse as spa
from scipy.sparse import linalg as sla

import mat_functions as matf
from process_functions import Monster


def build_kmat(elb, kmat=None, shape=None):

    def adele_for_build_kmat(adj):
        adbig = spa.block_diag([adj] * l_ele)
        kele00 = spa.bsr_matrix((adbig.transpose().dot(kele0)).dot(adbig), blocksize=(6, 6))
        return kele00.indptr, kele00.indices, kele00.data

    def builder_for_build_kmat():
        indptr, indices, data = ind
        indices2 = ele_elb[i, indices]
        kele = spa.bsr_matrix(
            (data, indices2, indptr), shape=(indptr.shape[0] * 6 - 6, shape[1]), blocksize=(6, 6)).transpose()
        kkmat = spa.bsr_matrix(
            (kele.data, indices2, kele.indptr), shape=shape, blocksize=(6, 6)).tocsr()
        return kkmat

    if shape is None:
        if kmat is None:
            if isinstance(elb, list):
                if len(elb) == 1:
                    maxnum = elb[0]['elb'].max()
                else:
                    # warnings.warn('请定义 shape ，不然可能影响运算速度')
                    warnings.warn('请定义 shape ，不然可能影响运算速度')
                    maxnum = 0
                    for i in range(len(elb)):
                        maxnum = max(np.array(elb[i]['elb']).max(), maxnum)
                shape = (maxnum * 6 + 6, maxnum * 6 + 6)
            else:
                raise ValueError("elb 如果为生成器类型那么需要定义参数 shape")
        else:
            shape = kmat.shape
    if kmat is None:
        kmat = spa.csr_matrix(shape)
    if elb is None or []:
        return kmat
    monster = Monster()
    for i_elb in elb:
        kele0 = i_elb['class'].kmat_sparse()
        ele_elb = np.array(i_elb['elb'])
        h_ele, l_ele = ele_elb.shape
        rot = np.array(i_elb['rotation'])
        h_rot = rot.shape
        if h_rot == (3, ) or (1, 3):
            ind = adele_for_build_kmat(matf.adjoint_k(0, 0, 0, rot[0], rot[1], rot[2]))
            for i in range(h_ele):
                kmat_food = builder_for_build_kmat()
                monster.prepare(kmat_food)

        elif h_rot[0] == h_ele:
            for i in range(h_ele):
                ind = adele_for_build_kmat(matf.adjoint_k(0, 0, 0, rot[i, 0], rot[i, 1], rot[i, 2]))
                kmat_food = builder_for_build_kmat()
                monster.prepare(kmat_food)
        else:
            raise TypeError("字典 elb['rotation'] 的值存在问题")
    monster.eat_all()
    if monster.desk is None:
        return kmat
    else:
        kmat += monster.desk.food
        return kmat

# --------------------------------------------------------


def rigid(xy, feature_xy, kmat, rigid_node):
    # time_start = time.time()
    # print('runnung sparse_functions.rigid()')

    if not isinstance(rigid_node, list):
        rigid_node = list(rigid_node)

    lens = len(rigid_node) * 6
    shape_bf = kmat.shape[0]
    shape_af = shape_bf - lens
    kmat = spa.bsr_matrix(kmat, blocksize=(6, 6))
    maper = xy[:, 1:4] - feature_xy
    ad_big = spa.csr_matrix(
        np.hstack([matf.adjoint_k(maper[i, 0], maper[i, 1], maper[i, 2], 0, 0, 0).T for i in rigid_node]))

    kmat00, kmat11 = row_block_bsr(kmat, rigid_node, returnall=True)
    kmat = spa.bsr_matrix(
        spa.vstack((ad_big.dot(kmat00), kmat11)), blocksize=(6, 6), shape=(shape_af + 6, shape_bf)).transpose()
    kmat00, kmat11 = row_block_bsr(kmat, rigid_node, returnall=True)
    kmat = spa.bsr_matrix(
        spa.vstack((ad_big.dot(kmat00), kmat11)), blocksize=(6, 6), shape=(shape_af + 6, shape_af + 6))
    # print('rigid time cost: ', time.time() - time_start, '\n -------------------------')
    return kmat


def rigid_xy(xy, feature_xy, rigid_node):
    xy = np.delete(xy, rigid_node[1:], axis=0)
    xy[rigid_node[0], 1:4] = feature_xy
    return xy


def redundant(xy, kmat, redundant_node):
    """
    冗余结点单元的刚度矩阵变换过程
    ------
    :param xy: 结点的坐标    矩阵
    :param kmat: 变换前的刚度矩阵
    :param redundant_node: 冗余结点的编号  列表
    :return: 一个列表，第一个为消除了冗余结点的 xy
             第二个为变换后的刚度矩阵
    """

    kmat_bb, kmat_aa = row_block_bsr(kmat, redundant_node, returnall=True)
    kmat_ba, kmat_aa = row_block_bsr(kmat_aa.transpose(), redundant_node, returnall=True)
    kmat_bb, kmat_ab = row_block_bsr(kmat_bb.transpose(), redundant_node, returnall=True)
    xy = np.delete(xy, redundant_node, axis=0)
    return xy, kmat_aa - (kmat_ab.dot(sla.inv(kmat_bb.tocsc()))).dot(kmat_ba)


def build_constraint(kmat, constraint_nod):
    return row_delete_bsr(row_delete_bsr(kmat, constraint_nod).transpose(), constraint_nod)


def sp_inv(kmat, node):

    if not spa.isspmatrix_csr(kmat):
        kmat = spa.csc_matrix(kmat)

    position = list(range(6 * node, 6 * node + 6))
    ii = spa.csc_matrix(
        ([1, 1, 1, 1, 1, 1], (position, [0, 1, 2, 3, 4, 5])), shape=(kmat.shape[0], 6))
    # Ainv = sla.splu(kmat).solve(ii.toarray())
    ainv = sla.spsolve(kmat, ii)[position]
    return ainv


def row_delete_bsr(kmat, ptr, blocksize=(6, 6)):
    kmat = spa.bsr_matrix(kmat, blocksize=blocksize)
    ptr = np.array(ptr)

    indptr, indices, data = kmat.indptr, kmat.indices, kmat.data
    indptr2 = copy.copy(indptr)
    nod, nnn = np.array(range(0, indptr[-1])), []
    for i in ptr:
        nnn += list(range(indptr[i], indptr[i + 1]))
        indptr2[i + 1:] -= indptr[i + 1] - indptr[i]
    nod = np.delete(nod, nnn)
    indices = indices[nod]
    data = data[nod, :, :]
    indptr2 = np.delete(indptr2, ptr + 1)
    kmat = spa.bsr_matrix((data, indices, indptr2), blocksize=(6, 6))

    return kmat


def row_block_bsr(kmat, ptr, blocksize=(6, 6), returnall=False, shape=None):
    if shape is None:
        shape = kmat.shape
    kmat = spa.bsr_matrix(kmat, blocksize=blocksize, shape=shape)
    ptr = np.array(ptr)
    indptr, indices, data = kmat.indptr, kmat.indices, kmat.data
    indptr_b, indptr_a = copy.copy(indptr), np.zeros(len(indptr), dtype=int)
    nod, nnn = np.array(range(0, indptr[-1])), []
    for i in ptr:
        nnn += list(range(indptr[i], indptr[i + 1]))
        indptr_b[i + 1:] -= indptr[i + 1] - indptr[i]
        indptr_a[i + 1:] += indptr[i + 1] - indptr[i]
    indices_a = indices[nnn]
    data_a = data[nnn, :, :]
    indptr_a = indptr_a[np.hstack((0, ptr + 1))]
    kmat_a = spa.bsr_matrix((data_a, indices_a, indptr_a), blocksize=(6, 6), shape=(len(ptr) * 6, shape[1]))
    if returnall is False:
        return kmat_a
    else:
        nod = np.delete(nod, nnn)
        indices_b = indices[nod]
        data_b = data[nod, :, :]
        indptr_b = np.delete(indptr_b, ptr + 1)
        kmat_b = spa.bsr_matrix(
            (data_b, indices_b, indptr_b), blocksize=(6, 6), shape=(shape[0] - len(ptr) * 6, shape[1]))
        return kmat_a, kmat_b


def build_xy():

    print('Function: build_xy()  still working on')


def bsrtoboo(indices, indptr, shape=None):
    row, col = [], indices
    for i in range(indptr.shape[0] - 1):
        row += [i] * (indptr[i + 1] - indptr[i])
    if shape is None:
        return np.array(row), col
    else:
        return np.array(row), col, shape


def bootobsr(row, col, shape=None):

    #  ?????
    indices = col
    a0, indptr = 0, np.array([0])
    for i in range(0, row.shape[0]):
        a1 = row[i]
        aa = a1 - a0
        if aa > 0:
            indptr = np.append(indptr, np.array([i] * aa))
            a0 = a1
        if aa < 0:
            warnings.warn('row 要从小到大排列')
            sys.exit()

    if shape is None:
        indptr = np.append(indptr, [col.shape[0]])
    else:
        indptr = np.append(indptr, [col.shape[0]] * (int(shape[0] / 6) - len(indptr) + 1))
    return indices, indptr
