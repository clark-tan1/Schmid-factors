# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 23:40:09 2017

@author: GM Zhu
"""
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import math
import numpy as np
import pandas as pd

data_file = pd.read_excel('MgCaE350-0-deletesomegrains(0-2)-20180328.xlsx')
data_file2 = pd.read_excel('Schmid_Factor_Slip_Systems.xlsx')
#print(data_file[0:3])
Euler = data_file.ix[:,'Numbers':'E3']
SlipSystem = data_file2.ix[:,'Numbers':'Direction4']
EulerMatrix = np.matrix(Euler)
SlipMatrix = np.matrix(SlipSystem)
#print(EulerMatrix[0,0])

a = 0
print('第二行为该滑移系编号,第三行开始每一行为不同滑移系的施密特因子')
print(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,52,53,54)
print(SlipMatrix)

while a < 430:
    Euler1 = float(EulerMatrix[a,1])
    Euler2 = float(EulerMatrix[a,2])
    Euler3 = float(EulerMatrix[a,3])
    a += 1
    
    #Euler1 = float(input("请输入Euler1\n"))
    #Euler2 = float(input("请输入Euler2\n"))
    #Euler3 = float(input("请输入Euler3\n"))
    Euler1rad = Euler1 * math.pi / 180
    Euler2rad = Euler2 * math.pi / 180
    Euler3rad = Euler3 * math.pi / 180
    Euler = [Euler1rad, Euler2rad, Euler3rad]
    gT11 = math.cos(Euler1rad) * math.cos(Euler3rad) - math.sin(Euler1rad) * math.sin(Euler3rad) * math.cos(Euler2rad)
    gT12 = -math.cos(Euler1rad) * math.sin(Euler3rad) - math.sin(Euler1rad) * math.cos(Euler3rad) * math.cos(Euler2rad)
    gT13 = math.sin(Euler1rad) * math.sin(Euler2rad)
    gT21 = math.sin(Euler1rad) * math.cos(Euler3rad) + math.cos(Euler1rad) * math.sin(Euler3rad) * math.cos(Euler2rad)
    gT22 = - math.sin(Euler1rad) * math.sin(Euler3rad) + math.cos(Euler1rad) * math.cos(Euler3rad) * math.cos(Euler2rad)
    gT23 = - math.cos(Euler1rad) * math.sin(Euler2rad)
    gT31 = math.sin(Euler3rad) * math.sin(Euler2rad)
    gT32 = math.cos(Euler3rad) * math.sin(Euler2rad)
    gT33 = math.cos(Euler2rad)
    gT = np.matrix([[gT11, gT12, gT13],  [gT21, gT22, gT23], [gT31, gT32, gT33]])
    g = gT.T
    sigmanormalized = np.matrix([[-1, 0, 0], [0, 0, 0], [0, 0, 0]])
    gsigma = g * sigmanormalized
    gsigmagT = gsigma * gT
    # Basal slip (0001)[2-1-10] as B1
    # (hkil) = hB1, kB1, iB1, lB1; [uvtw] = uB1, vB1, tB1, wB1
    #c/a=1.62
    #input (hkil)[uvtw]
    b = 0
    n1 = np.array([11.1,22.2,33.3,44.4,55.5,66.6,7.7,8.8,9.9,10.1,11.1,22.2,33.3,44.4,55.5,66.6,7.7,8.8,9.9,10.1,11.1,22.2,33.3,44.4,55.5,66.6,7.7,8.8,9.9,10.1,11.1,22.2,33.3,44.4,55.5,66.6,7.7,8.8,9.9,10.1,11.1,22.2,33.3,44.4,55.5,66.6,7.7,8.8,9.9,10.1,11.1,22.2,33.3,44.4])
    while b < 54:
        hB1 = float(SlipMatrix[b, 1])
        kB1 = float(SlipMatrix[b, 2])
        iB1 = float(SlipMatrix[b, 3])
        lB1 = float(SlipMatrix[b, 4])
        uB1 = float(SlipMatrix[b, 5])
        vB1 = float(SlipMatrix[b, 6])
        tB1 = float(SlipMatrix[b, 7])
        wB1 = float(SlipMatrix[b, 8])

        #calculate n and m
        n1B1 = hB1
        n2B1 = (2 * kB1 + hB1)/(np.sqrt(3))
        n3B1 = lB1/1.62
        n4B1 = np.sqrt(n1B1 * n1B1 + n2B1 * n2B1 + n3B1 * n3B1)
        m1B1 = 1.5 * uB1
        m2B1 = np.sqrt(3)/2 *(2 * vB1 + uB1)
        m3B1 = wB1 * 1.62
        m4B1 = np.sqrt(m1B1 * m1B1 + m2B1 * m2B1 + m3B1 * m3B1)
    
        #calculate normalized n and m
        n1norB1 = n1B1/n4B1
        n2norB1 = n2B1/n4B1
        n3norB1 = n3B1/n4B1
        m1norB1 = m1B1/m4B1
        m2norB1 = m2B1/m4B1
        m3norB1 = m3B1/m4B1
    
        #calculate minj and nimj
        minjB1 = np.matrix([[1.0,1.0,1.0],[1.0,1.0,1.0],[1.0,1.0,1.0]])
        minjB1[0,0] = n1norB1 * m1norB1
        minjB1[0,1] = n2norB1 * m1norB1
        minjB1[0,2] = n3norB1 * m1norB1
        minjB1[1,0] = n1norB1 * m2norB1
        minjB1[1,1] = n2norB1 * m2norB1
        minjB1[1,2] = n3norB1 * m2norB1
        minjB1[2,0] = n1norB1 * m3norB1
        minjB1[2,1] = n2norB1 * m3norB1
        minjB1[2,2] = n3norB1 * m3norB1
        nimjB1 = minjB1.T
    
        #calculate basal slip bs
        bsB1 = (minjB1 + nimjB1) * 0.5
    
        #calculate Schmid factor
        #Sfmatrix is the matrix fo Schmid factor
        sfmatrixB1 = np.multiply(bsB1, gsigmagT)
        SFB1 = np.sum(sfmatrixB1)
        n1[b] = SFB1
        b += 1
    print(n1[0], n1[1], n1[2], n1[3], n1[4], n1[5], n1[6], n1[7], n1[8], n1[9], n1[10], n1[11], n1[12], n1[13], n1[14], n1[15], n1[16], n1[17], n1[18], n1[19], n1[20], n1[21], n1[22], n1[23], n1[24], n1[25], n1[26], n1[27], n1[28], n1[29], n1[30], n1[31], n1[32], n1[33], n1[34], n1[35], n1[36], n1[37], n1[38], n1[39], n1[40], n1[41], n1[42], n1[43], n1[44], n1[45], n1[46], n1[47], n1[48], n1[49], n1[50], n1[51], n1[52], n1[53])
