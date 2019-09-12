# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 22:16:54 2018

@author: 朱高明
"""

import numpy as np
import pandas as pd
import math
import xlsxwriter

def in_Excel(file_path):
    data_file = pd.read_excel(file_path)
    dataF = data_file.ix[:, :]
    data_matrix = np.matrix(dataF)
    return(data_file, dataF, data_matrix)

def Schmid(E0, E1, E2, hh, kk, ii, ll, uu, vv, tt, ww):
    Euler1rad = E0 * math.pi / 180
    Euler2rad = E1 * math.pi / 180
    Euler3rad = E2 * math.pi / 180
#    Euler = [Euler1rad, Euler2rad, Euler3rad]
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
    sigmanormalized = np.matrix([[0, 0, 0], [0, 1, 0], [0, 0, 0]])
    gsigma = g * sigmanormalized
    gsigmagT = gsigma * gT

    hB1 = hh
    kB1 = kk
    iB1 = ii
    lB1 = ll
    uB1 = uu  
    vB1 = vv
    tB1 = tt
    wB1 = ww
    
    n1B1 = hB1
    n2B1 = (2 * kB1 + hB1)/(np.sqrt(3))
    n3B1 = lB1/1.62
    n4B1 = np.sqrt(n1B1 * n1B1 + n2B1 * n2B1 + n3B1 * n3B1)
    m1B1 = 1.5 * uB1
    m2B1 = np.sqrt(3)/2 *(2 * vB1 + uB1)
    m3B1 = wB1 * 1.62
    m4B1 = np.sqrt(m1B1 * m1B1 + m2B1 * m2B1 + m3B1 * m3B1)

    n1norB1 = n1B1/n4B1
    n2norB1 = n2B1/n4B1
    n3norB1 = n3B1/n4B1
    m1norB1 = m1B1/m4B1
    m2norB1 = m2B1/m4B1
    m3norB1 = m3B1/m4B1

    # burgers vectors
    Burg01 = gT11 * m1norB1 + gT12 * m2norB1 + gT13 * m3norB1
    Burg02 = gT21 * m1norB1 + gT22 * m2norB1 + gT23 * m3norB1
    Burg03 = gT31 * m1norB1 + gT32 * m2norB1 + gT33 * m3norB1
#-----------Schmid factors
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

    return(Burg01, Burg02, Burg03, SFB1)


data_file, data0F, data0 = in_Excel('gf_raw_45deg.xlsx')

out_18 = np.zeros((len(data0), 15))
Burgers = np.zeros((len(data0), 3))
basal_SF = np.zeros((len(data0), 3))
pri_SF = np.zeros((len(data0), 3))
pyr_a_SF = np.zeros((len(data0), 6))
pyr_ca_SF = np.zeros((len(data0), 12))
pyrII_ca_SF = np.zeros((len(data0), 6))
twin_SF = np.zeros((len(data0), 6))
#area_neighbor = np.zeros((len(data0), 1))
#cos_neighbor = np.zeros((len(data0), 1))

hkil = np.matrix([[0,0,0,1,2,-1,-1,0], #basal
                  [0,0,0,1,-1,2,-1,0],
                  [0,0,0,1,-1,-1,2,0], 
                  [0,1,-1,0,2,-1,-1,0], #pri
                  [1,0,-1,0,-1,2,-1,0],
                  [-1,1,0,0,-1,-1,2,0],
                  [0,1,-1,1,2,-1,-1,0], #pyr_a
                  [1,0,-1,1,-1,2,-1,0],
                  [-1,1,0,1,-1,-1,2,0],
                  [0,-1,1,1,2,-1,-1,0],
                  [-1,0,1,1,-1,2,-1,0],
                  [1,-1,0,1,-1,-1,2,0],
                  [1,0,-1,1,2,-1,-1,-3], #pyr_ca
                  [1,0,-1,1,1,1,-2,-3],
                  [0,1,-1,1,1,1,-2,-3],
                  [0,1,-1,1,-1,2,-1,-3],
                  [-1,1,0,1,-1,2,-1,-3],
                  [-1,1,0,1,-2,1,1,-3],
                  [-1,0,1,1,-2,1,1,-3],
                  [-1,0,1,1,-1,-1,2,-3],
                  [0,-1,1,1,-1,-1,2,-3],
                  [0,-1,1,1,1,-2,1,-3],
                  [1,-1,0,1,1,-2,1,-3],
                  [1,-1,0,1,2,-1,-1,-3],
                  [2,-1,-1,2,2,-1,-1,-3],#pyrII_ca
                  [1,1,-2,2,1,1,-2,-3],
                  [-1,2,-1,2,-1,2,-1,-3],
                  [-2,1,1,2,-2,1,1,-3],
                  [-1,-1,2,2,-1,-1,2,-3],
                  [1,-2,1,2,1,-2,1,-3],
                  [1,0,-1,2,-1,0,1,1], #twin
                  [0,1,-1,2,0,-1,1,1],
                  [-1,1,0,2,1,-1,0,1],
                  [-1,0,1,2,1,0,-1,1],
                  [0,-1,1,2,0,1,-1,1],
                  [1,-1,0,2,-1,1,0,1]])

for i in range(len(data0)):  #Schmid factors
    xxx1, xxx2, xxx3, xxx4 = Schmid(data0[i, 1], data0[i, 2], data0[i, 3], 0,0,0,0,0,0,0,1)
    Burgers[i, 0], Burgers[i, 1], Burgers[i, 2] = abs(xxx1), abs(xxx2), abs(xxx3)
    for j in range(3):
        basal_SF[i, j] = abs(Schmid(data0[i, 1], data0[i, 2], data0[i, 3], hkil[j, 0], hkil[j, 1], hkil[j, 2], hkil[j, 3], hkil[j, 4], hkil[j, 5], hkil[j, 6], hkil[j, 7])[3])
    for k in range(3):
        j = k + 3
        pri_SF[i, k] = abs(Schmid(data0[i, 1], data0[i, 2], data0[i, 3], hkil[j, 0], hkil[j, 1], hkil[j, 2], hkil[j, 3], hkil[j, 4], hkil[j, 5], hkil[j, 6], hkil[j, 7])[3])
    for k in range(6):
        j = k + 6
        pyr_a_SF[i, k] = abs(Schmid(data0[i, 1], data0[i, 2], data0[i, 3], hkil[j, 0], hkil[j, 1], hkil[j, 2], hkil[j, 3], hkil[j, 4], hkil[j, 5], hkil[j, 6], hkil[j, 7])[3])
    for k in range(12):
        j = k + 12
        pyr_ca_SF[i, k] = abs(Schmid(data0[i, 1], data0[i, 2], data0[i, 3], hkil[j, 0], hkil[j, 1], hkil[j, 2], hkil[j, 3], hkil[j, 4], hkil[j, 5], hkil[j, 6], hkil[j, 7])[3])
    for k in range(6):
        j = k + 24
        pyrII_ca_SF[i, k] = abs(Schmid(data0[i, 1], data0[i, 2], data0[i, 3], hkil[j, 0], hkil[j, 1], hkil[j, 2], hkil[j, 3], hkil[j, 4], hkil[j, 5], hkil[j, 6], hkil[j, 7])[3])
    for k in range(6):
        j = k + 30
        twin_SF[i, k] = Schmid(data0[i, 1], data0[i, 2], data0[i, 3], hkil[j, 0], hkil[j, 1], hkil[j, 2], hkil[j, 3], hkil[j, 4], hkil[j, 5], hkil[j, 6], hkil[j, 7])[3]
    basal_SF[i].sort()
    pri_SF[i].sort()
    pyr_a_SF[i].sort()
    pyr_ca_SF[i].sort()
    pyrII_ca_SF[i].sort()
    twin_SF[i].sort()


for i in range(len(data0)):
    for j in range(15):
         out_18[i, j] = data0[i, j]
   
c_z = np.zeros((len(data0), 6))
for i in range(len(data0)):
    c_z[i, 0] = math.acos(Burgers[i, 1])*180/math.pi

out_all = np.hstack((out_18, basal_SF, pri_SF, pyr_a_SF, pyr_ca_SF, pyrII_ca_SF, twin_SF, c_z))

titles = ('GrainNo', 'phi1', 'PHI', 'phi2', 'PositionX', 'PositionY', 'IQ', 'CI', 'Edge_grain', 'Area', 'Diameter', 'Aspect_ratio',
          'Orientation_spread', 'Misorientation', 'No_neighbor', 'Basal1', 'Basal2', 'Basal3', 'Prism1', 'Prism2', 'Prism3', 'PyrI<a>1', 'PyrI<a>2', 
          'PyrI<a>3', 'PyrI<a>4', 'PyrI<a>5', 'PyrI<a>6', 'PyrI<c+a>1', 'PyrI<c+a>2', 'PyrI<c+a>3', 'PyrI<c+a>4', 
          'PyrI<c+a>5', 'PyrI<c+a>6', 'PyrI<c+a>7', 'PyrI<c+a>8', 'PyrI<c+a>9', 'PyrI<c+a>10', 'PyrI<c+a>11', 'PyrI<c+a>12', 
          'PyrII<c+a>1', 'PyrII<c+a>2', 'PyrII<c+a>3', 'PyrII<c+a>4', 'PyrII<c+a>5', 'PyrII<c+a>6', 
          'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'c_x')
    
workbook = xlsxwriter.Workbook('SF_45deg_20190908.xlsx')
worksheet = workbook.add_worksheet()

i = 0
while i < len(data0):  #lines
    j = 0
    while j < len(titles):  
        worksheet.write(i + 1, j, out_all[i, j])
        if i == 0:
            worksheet.write(i, j, titles[j])
        j += 1
    i += 1

workbook.close()

    
    
    
    
    
    
    
    
    