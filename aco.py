# -*- coding: utf-8 -*-
import numpy as np
#import os
#from matplotlib import pyplot as plt
coordinates = np.array([[565.0, 575.0], [25.0, 185.0], [345.0, 750.0], [945.0, 685.0], [845.0, 655.0],
                        [880.0, 660.0], [25.0, 230.0], [525.0, 1000.0], [580.0, 1175.0], [650.0, 1130.0],
                        [1605.0, 620.0], [1220.0, 580.0], [1465.0, 200.0], [1530.0, 5.0], [845.0, 680.0],
                        [725.0, 370.0], [145.0, 665.0], [415.0, 635.0], [510.0, 875.0], [560.0, 365.0],
                        [300.0, 465.0], [520.0, 585.0], [480.0, 415.0], [835.0, 625.0], [975.0, 580.0],
                        [1215.0, 245.0], [1320.0, 315.0], [1250.0, 400.0], [660.0, 180.0], [410.0, 250.0],
                        [420.0, 555.0], [575.0, 665.0], [1150.0, 1160.0], [700.0, 580.0], [685.0, 595.0],
                        [685.0, 610.0], [770.0, 610.0], [795.0, 645.0], [720.0, 635.0], [760.0, 650.0],
                        [475.0, 960.0], [95.0, 260.0], [875.0, 920.0], [700.0, 500.0], [555.0, 815.0],
                        [830.0, 485.0], [1170.0, 65.0], [830.0, 610.0], [605.0, 625.0], [595.0, 360.0],
                        [1340.0, 725.0], [1740.0, 245.0]])

class acs:
    def __init__(self,coordinates, alpha, beta, Rho, pheromone0):
        self.coordinates=coordinates
        self.numant=10
        self.alpha=alpha#
        self.beta=beta#
        self.Rho=Rho#信息素全局蒸发
        self.numcity=coordinates.shape[0]
        self.pheromone0=pheromone0#
        self.pathtable = np.zeros((self.numant, self.numcity)).astype(int)#路径记录表，转化为int
        self.pheromone_mat = np.zeros((self.numcity , self.numcity)).astype(float)
        self.distmat=np.zeros((self.numcity,self.numcity)).astype(float)
        self.Tb=[]
        self.Xi =0.1#信息素局部蒸发


    def __get_dist_mat(self,coordinates,numcity,distmat):
        # 初始化生成52*52的矩阵distense matrix
        for i in range(numcity):
            for j in range(i, numcity):
                distmat[i][j] = distmat[j][i] = np.linalg.norm(coordinates[i] - coordinates[j])       

    def init_pheromone_mat(self,coordinates,numcity,pheromone_mat,pheromone0):
        for i in range(numcity):
            for j in range(i, numcity):
                pheromone_mat[i][j] = pheromone_mat[j][i] = pheromone0

    def choose_depart_city_randomly(self,coordinates,numcity,numant,pathtable):
        pathtable[:, 0] = np.random.permutation(range(numcity))[:numant]

    def _get_visiting(self, pathtable,i,numant):
        visiting = np.zeros(numant).astype(int)
        for k in range(numant):    
            visiting[k] = pathtable[k][i]
        return visiting

    def pheromone_globally_update(self,numcity,Tb,pheromone_mat,Rho,pheromone0):
        delta = pheromone0
        for i in numcity:
            for j in numcity:
                if ([i,j] in Tb):
                    pheromone_mat[i][j]=pheromone_mat[j][i]=(1-Rho)*pheromone_mat[i][j] + Rho*delta

    def get_path_length(self,k,numcity,pathtable,distmat):
        for i in range(1,numcity):
            length += distmat[pathtable[k][i-1]][pathtable[k][i]]
        return length

    def get_globally_best(self,Tb,numant,get_path_length,numcity,pathtable,distmat):#Tb源头
        length = []
        for k in range(numant):
            length.append(get_path_length(k,numcity,pathtable,distmat))
        lengthbest = length.index(min(length))
        Tb = pathtable[lengthbest]

    def status_move(self,numcity,q0,pathtable,i,numant,_get_visiting,distmat,beta,pheromone_mat,pheromone0,Xi):
        visiting = _get_visiting(pathtable,i,numant)
        for k in range(numant):
            Jk = [x for x in range(numcity) if x not in pathtable[k]]
            Posibility = np.zeros((numcity,numcity))
            q=np.random.random()
            if q<=q0:
                for j in Jk:#Develop
                    if j == np.argmax((np.power(1/distmat[visiting[k]][j],beta),pheromone_mat[visiting[k]][j])):
                        pathtable[k][i+1] = j                    
            else:           
                summ = []
                for u in Jk:
                    summ.append(np.power(1/distmat[visiting[k]][u],beta)*pheromone_mat[visiting[k]][u])
                for j in Jk:
                    if j == 0:
                        Posibility[visiting[k]][j] = summ[j]/sum(summ)
                    elif j>0:
                        Posibility[visiting[k]][j] = summ[j]/sum(summ) + Posibility[visiting[k]][j-1]
                #用模拟轮盘选择下个地点
                q=np.random.random()
                for j in range(1,numcity):
                    if (Posibility[visiting[k]][j-1]<q<=Posibility[visiting[k]][j]):
                        pathtable[k][i+1] = j
                        pheromone_mat[visiting[k]][j] = pheromone_mat[j][visiting[k]] = \
                            (1-Xi)*pheromone_mat[visiting[k]][j] + Xi*pheromone0#信息素局部更新 