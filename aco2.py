# -*- coding: utf-8 -*-
import numpy as np
#import os
from matplotlib import pyplot as plt
coordinates = np.array([[37.0, 52.0], [49.0, 49.0], [52.0, 64.0], [20.0, 26.0], [40.0, 30.0],
                        [21.0, 47.0], [17.0, 63.0], [31.0, 62.0], [52.0, 33.0], [51.0, 21.0],
                       [42.0, 41.0], [31.0, 32.0], [5.0, 25.0], [12.0, 42.0], [36.0, 16.0],
                        [52.0, 41.0], [27.0, 23.0], [17.0, 33.0], [13.0, 13.0], [57.0, 58.0],
                        [62.0, 42.0], [42.0, 57.0], [16.0, 57.0], [8.0, 52.0], [7.0, 38.0],
                        [27.0, 68.0], [30.0, 48.0], [43.0, 67.0], [58.0, 48.0],[58.0,27.0],
                        [37.0, 69.0], [38.0, 46.0], [46.0, 10.0], [61.0, 33.0], [62.0, 63.0],
                        [63.0, 69.0], [32.0, 22.0], [45.0, 35.0], [59.0, 15.0], [5.0, 6.0],
                        [10.0, 17.0], [21.0, 10.0], [5.0, 64.0], [30.0, 15.0], [39.0, 10.0],
                        [32.0, 39.0], [25.0, 32.0], [25.0, 55.0], [48.0, 28.0], [56.0, 37.0],
                        [30.0,40.0]])
    #def __init__(,coordinates, beta, Rho, pheromone0,q0,itermax):

numant=10
beta=4#
Rho=0.1 #信息素全局蒸发
numcity=coordinates.shape[0]
pheromone0=0.01
pathtable = [] #np.zeros((numant, numcity)).astype(int)#路径记录表，转化为int
pheromone_mat = np.zeros((numcity , numcity)).astype(float)
distmat=np.zeros((numcity,numcity)).astype(float)
TB=[]
Xi =0.01#信息素局部蒸
itermax = 5
q0 = 0.2
lengthBest = []
length = []

# 初始化生成52*52的矩阵distense matrix
for i in range(numcity):
    for j in range(i, numcity):
        distmat[i][j] = distmat[j][i] = np.linalg.norm(coordinates[i] - coordinates[j])       

def init_pheromone_mat():
    global pheromone_mat
    for i in range(numcity):
        for j in range(i, numcity):
            pheromone_mat[i][j] = pheromone_mat[j][i] = pheromone0

def choose_depart_city_randomly():
    global pathtable
    for i in range(numant):
        a=[np.random.randint(0,numcity)]
        pathtable.append(a) 

def _get_visiting(k):
    path = pathtable[k]
    visiting = path.pop()
    pathtable[k].append(visiting)
    return visiting

def get_globally_best():#Tb源头
    global TB,lengthBest,length
    for k in range(numant):
        g= _get_path_length(k)
        length.append(g)
        if g<=min(length):
            TB=k
    print("minlength",min(length))
    lengthBest.append(min(length))
    print ("minlengthBest",min(lengthBest))

def pheromone_globally_update():
    global pheromone_mat
    delta = pheromone0
    get_globally_best()
    for i in range(numcity):
        for j in range(numcity):
            if ([i,j] in TB):
                pheromone_mat[i][j]=pheromone_mat[j][i]=(1-Rho)*pheromone_mat[i][j] + Rho*delta

def _get_path_length(k):
    length = 0
    for i in range(len(pathtable[k])-1):
        length += distmat[i][i+1]
    return length

def status_move():#i 是各个蚂蚁当前到访的城市在pathtable中的编号
    global pathtable,pheromone_mat   
    for k in range(numant):
        visiting = _get_visiting(k) #shape(visiting)=10
        Jk = [x for x in range(numcity) if x not in pathtable[k]]
        summ = np.zeros(numcity)
        for u in Jk:
            summ[u]=(1/distmat[visiting][u] ** beta)*pheromone_mat[visiting][u]
        Eta = list(1/distmat[visiting] ** beta)
        Eta[np.argmax(Eta)]=0
        Posibility = np.zeros(numcity)
        q=np.random.random()
        if q<=q0:
            for j in Jk:#Develop
                if j == np.argmax(Eta) or (iter>=1 and j == np.argmax(pheromone_mat[visiting])):
                    pathtable[k].append(j)
                    pheromone_mat[visiting][j] = pheromone_mat[j][visiting] = (1-Xi)*pheromone_mat[visiting][j] + Xi*pheromone0#信息素局部更新
                    break
                    #print(j)                    
        else:           
            for j in Jk:
                if j == 0:
                    Posibility[j] = summ[j]/sum(summ)
                elif j>0:
                    Posibility[j] = summ[j]/sum(summ) + Posibility[Jk[Jk.index(j)-1]]
            #用模拟轮盘选择下个地点
            q1=np.random.random()
            for j in Jk:
                if (j==0 and q<Posibility[j]):
                    pathtable[k].append(j)
                    pheromone_mat[visiting][j] = pheromone_mat[j][visiting] = (1-Xi)*pheromone_mat[visiting][j] + Xi*pheromone0#信息素局部更新
                    break
                elif (Posibility[Jk[Jk.index(j)-1]]<q1<=Posibility[j]):
                    pathtable[k].append(j)
                    pheromone_mat[visiting][j] = pheromone_mat[j][visiting]=(1-Xi)*pheromone_mat[visiting][j] + Xi*pheromone0#信息素局部更新
                    break
        #print (pathtable)

init_pheromone_mat()
for iter in range(itermax):  
    choose_depart_city_randomly()    
    for i in range(numcity):
        status_move()
    pheromone_globally_update()
    print("iter",iter,"end")

print(TB)
print (sum(TB))
lenth = 
print(lenth)

plt.plot(range(itermax),lengthBest)
plt.title('title')
plt.ylabel('y')
plt.xlabel('x')
plt.xticks(list(range(itermax)))
plt.yticks(list(range(700,1600,20)))
plt.ylim([700,1500])
plt.xlim([0,itermax])

plt.grid(True)
plt.show()