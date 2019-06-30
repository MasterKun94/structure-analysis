import os
import time

import xlsxwriter

import s1_functions as cf
import structure_class as sc

nx, ny, nz = 1, 1, 1
thc, wid, ele_d = 1, 1, 20
ela, v = 2000, 0.3
p = 0.5
path = os.path.abspath('./data/data.xlsx')
file = xlsxwriter.Workbook(path)
sheet1, sheet2, sheet3 = file.add_worksheet('param'), file.add_worksheet('sitalist'), file.add_worksheet('cmat')

for i in range(4):
    print('*****************' + ' round << ' + str(i) + ' >> *******************')
    time_start = time.perf_counter()
    sitalist = cf.choose_sita(nx, ny, nz, ele_d, p)
    # print(sitalist)
    structure1 = sc.Structure3D(nx, ny, nz, thc, wid, ela, v, ele_d, sitalist)
    cmat = structure1.solve_node_cmat(0)
    print(structure1)

    k, [h, l] = 0, cmat.shape
    sheet2.write(4 * i, 0, i)
    sheet3.write(7 * i, 0, i)
    for j in [i, nx, nx, nz, thc, wid, ele_d, ela, v, p]:
        sheet1.write(i + 1, k, j)
        k += 1
    for j in [0, 1, 2]:
        for k in range(0, len(sitalist[j])):
            sheet2.write(4 * i + j, k + 1, sitalist[j][k])
    for j in range(0, h):
        for k in range(0, l):
            if abs(cmat[j, k]) > 1e-10:
                sheet3.write(j + 7 * i, k + 1, cmat[j, k])
            else:
                sheet3.write(j + 7 * i, k + 1, 0)

    thc, wid, ele_d = thc * nx / (nx + 2), wid * nx / (nx + 2), ele_d * nx / (nx + 2)
    nx += 2
    ny += 2
    nz += 2
    print('  total time cost: ', time.perf_counter() - time_start,
          '\n---------------------------------------------------')

k = 1
for i in ['nx', 'ny', 'nz', 'thc', 'wid', 'ele_d', 'ela', 'v', 'p']:
    sheet1.write(0, k, i)
    k += 1

file.close()
print('File saved at:', path)

# if __name__ == '__main__':
#     slist = [[1 , 2, 2, 1],
#              [1, 2, 2, 1],
#              [1, 2, 2, 1] ]
#     structure = sc.Structure3D(3, 3, 3, 1, 4, 500, 0.5, 20, slist)
#     c = structure.solve_node_cmat(0)
#     print(c.toarray()[])