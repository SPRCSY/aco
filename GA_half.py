# !/usr/bin/env python
# -*- coding:utf-8 -*-
# Author: wsw
#遗传算法求解y=f(x1,x2,x3,x4)=1/(x1^2+x2^2+x3^2+x4^2+1)其中-5<=x1,x2,x3,x4<=5, y的最大值
import random
import math
import numpy as np
from numpy import random

population_size=5 #种群数量
generations=1 #迭代次数
pc=0.88 #交配概率
pm=0.1 #变异概率
population =[] #种群对应的十进制数值，并标准化范围到（-5,5）
fitness =[] #适应度
fitness_sum=[]
optimum_solution=[] #每次迭代所获得的最优解
population_proportion=[]#每个染色体适应度总和的比

#生成初始种群
def init_population():
    N = np.random.uniform(-5, 5, size=(5, 5))
    for i in range(5):
        N[i][0] = i + 1
    population.extend(N)
    #print(population)
    return population
print(population)
#计算每个染色体的适应度
def calculate_fitness():
    for i in range(population_size):
        function_value=1/(population[i][1]**2+population[i][2]**2+population[i][3]**2+population[i][4]**2+1)
        fitness.append(function_value)
    return population[i][1],population[i][2],population[i][3],population[i][4]



#获取最大适应度的个体
def best_value():
    max_fitness=fitness[0]
    for i in range(population_size):
        if fitness[i]>max_fitness:
            max_fitness=fitness[i]
    return max_fitness


#采用轮盘赌算法进行选择过程
def selection():
    fitness_sum = 0
    for i in range(population_size):
        fitness_sum += fitness[i]
        # 计算生存率
    for i in range(population_size):
        population_proportion.append(fitness[i] / fitness_sum)
    pie_fitness = []
    cumsum = 0.0
    for i in range(population_size):
        pie_fitness.append(cumsum + population_proportion[i])
        cumsum += population_proportion[i]
    # 生成随机数在轮盘上选点[0, 1)
    random_selection = []
    for i in range(population_size):
        random_selection.append(random.random())
    # 选择新种群
    new_population = []
    random_selection_id = 0
    global population
    for i in range(population_size):
        while random_selection_id < population_size and random_selection[random_selection_id] < pie_fitness[i]:
            new_population.append(population[i])
            random_selection_id += 1
    population = new_population


#进行交配
def crossover():
    for i in range(0,population_size-1,2):
        if random.random()<pc:
            #随机选择交叉点
            change_point=random.randint(1,4)
            temp1=[]
            temp2=[]
            temp1.extend(population[i][0:change_point])
            temp1.extend(population[i+1][change_point:])
            temp2.extend(population[i+1][0:change_point])
            temp2.extend(population[i][change_point:])
            population[i]=temp1
            population[i+1]=temp2

#进行变异
def mutation():
    for i in range(population_size):
        if random.random()<pm:
            mutation_point=random.randint(1,4)#随机变异点
            population[i][mutation_point]=random.uniform(-5,5)

#循环
init_population()
for step in range(generations):
    calculate_fitness()
    x_1,x_2,x_3,x_4=calculate_fitness()
    best_fitness=best_value()
    optimum_solution.append(best_fitness)
    selection()
    crossover()
    mutation()

#输出最优解和对应X数值
print('最优解：',best_fitness)
print('最优系数')
print('x1=',x_1)
print('x2=',x_2)
print('x3=',x_3)
print('x4=',x_4)
