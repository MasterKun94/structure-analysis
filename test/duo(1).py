import math

import matplotlib.pyplot as plt
import numpy as np

import structure_class as sc


class PSO(object):
    def __init__(self, population_size, max_steps):
        self.w = 0.4  # 惯性权重
        self.c1 = self.c2 = 2
        self.population_size = population_size  # 粒子群数量
        self.dim = 2 # 搜索空间的维度
        self.max_steps = max_steps  # 迭代次数
        self.x_bound = [0, math.pi / 2]  # 解空间范围
        self.x = np.random.uniform(self.x_bound[0], self.x_bound[1],
                                   (self.population_size, self.dim))  # 初始化粒子群位置
        self.x = np.zeros((population_size, self.dim))
        self.x[:, 0] = 0.5
        self.x[:, 1] = 1.0

        self.x += np.random.rand(population_size, self.dim) * 0.1

        self.v = np.random.rand(self.population_size, self.dim)  # 初始化粒子群速度
        fitness = self.calculate_fitness(self.x)
        self.p = self.x  # 个体的最佳位置
        self.pg = self.x[np.argmin(fitness)]  # 全局最佳位置
        self.individual_best_fitness = fitness  # 个体的最优适应度
        self.global_best_fitness = np.max(fitness)  # 全局最佳适应度

    # def calculate_fitness(self, x):
    #     finalV = 0
    #     for element in x:
    #         #开头的代码{import structure_class as sc}表示导入struce_class这个包并改名为{sc}这里调用用的却是原来的名字，所以会有报错，
    #         #只要将这里改成sc，或者把{import structure_class as sc} 改成{import structure_class}就可以了
    #         #structure = structure_class.Structure3D(3, 3, 3, 1, 4, 500, 0.5, 20, slist)
    #
    #
    #         structure = sc.Structure3D(3, 3, 3, 1, 4, 500, 0.5, 20, slist) #像这样
    #
    #         c = structure.solve_node_cmat(0)
    #         valueA = (c[2, 2] / c[5, 2] * c[2,5] / c[5, 5])
    #         valueB = (c[2, 2]-c[5, 5]+ math.squrt((c[5, 5]-c[2, 2])**2+4*c[5, 2]))/2*c[5, 2]
    #         #对获得的两个值平方相加，值 1 和 0.3 分别为想要迭代得到的理想值，也就是说，值 finalV越小越好
    #         finalV += (valueA - 1) ** 2 + (valueB - 0.3) ** 2
    #     return  finalV

    def calculate_fitness(self, x):
        finalV = np.zeros(x.shape[0])
        i = 0
        for element in x:
            print(element)
            print(self.x_bound)
            booa = element[0] > self.x_bound[1]
            boob = element[0] < self.x_bound[0]
            booc = element[1] > self.x_bound[1]
            bood = element[1] < self.x_bound[0]
            if booa or boob or booc or bood:
                finalV[i] = self.individual_best_fitness[i]
            # 开头的代码{import structure_class as sc}表示导入struce_class这个包并改名为{sc}这里调用用的却是原来的名字，所以会有报错，
            # 只要将这里改成sc，或者把{import structure_class as sc} 改成{import structure_class}就可以了
            # structure = structure_class.Structure3D(3, 3, 3, 1, 4, 500, 0.5, 20, x)

            # 这里的参数element是一个一维列表，即[a, b],但是需要的格式是[[a, b, b, a], [a, b, b, a], [a, b, b, a]] 需要做一个转换
            a = np.append(element, -element[::-1])
            slist = [a, a, a]
            print(slist)
            structure = sc.Structure3D(3, 3, 3, 1, 4, 500, 0.5, 20, slist)

            c = structure.solve_node_cmat(0)
            valueA = (c[2, 2] / c[5, 2] * c[5, 5]) / c[2, 5]
            valueB = -c[2, 5] / c[5, 5]
            # valueB = (c[2, 2] - c[5, 5] + math.sqrt((c[5, 5] - c[2, 2]) ** 2 + 4 * c[5, 2])) / 2 * c[5, 2]
            # 对获得的两个值平方相加，值 1 和 0.3 分别为想要迭代得到的理想值，也就是说，值 finalV越小越好
            finalV[i] = (valueA - 1) ** 2 + (valueB - 0.3) ** 2
            i += 1
        return finalV

    def evolve(self):
        fig = plt.figure()
        for step in range(self.max_steps):
            print("step: " + str(step))
            r1 = np.random.rand(self.population_size, self.dim)
            r2 = np.random.rand(self.population_size, self.dim)
            # 更新速度和权重
            self.v = self.w * (self.v + self.c1 * r1 * (self.p - self.x) + self.c2 * r2 * (self.pg - self.x))
            self.x = self.v + self.x
            plt.clf()
            plt.scatter(self.x[:, 0], self.x[:, 1], s=30, color='k')
            # plt.xlim(-1, 1)
            # plt.ylim(-1, 1)
            plt.xlim(self.x_bound[0], self.x_bound[1])
            plt.ylim(self.x_bound[0], self.x_bound[1])
            plt.pause(0.01)
            fitness = self.calculate_fitness(self.x)
            print(fitness)
            # 需要更新的个体
            update_id = np.greater(self.individual_best_fitness, fitness)
            print(update_id)
            print("update id")
            self.p[update_id] = self.x[update_id]
            self.individual_best_fitness[update_id] = fitness[update_id]
            # 新一代出现了更小的fitness，所以更新全局最优fitness和位置
            if np.min(fitness) < self.global_best_fitness:
                self.pg = self.x[np.argmin(fitness)]
                print(self.pg)
                self.global_best_fitness = np.min(fitness)
            if self.global_best_fitness <0.000000001:
                print("global best fitness")
                break
            print('best fitness: %.5f, mean fitness: %.5f' % (self.global_best_fitness, np.mean(fitness)))

if __name__ == '__main__':

    pso = PSO(50, 50)
    pso.evolve()
    plt.show()

    x = [0.8,1]

    slist = [[0.8, 1, -1, -0.8],
             [0.8, 1, -1, -0.8]]
