# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 18:20:21 2019

@author: Shashank
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#filename='CHOH_training'
filename='H_training'

filename_all=[filename+'_vx',filename+'_vy',filename+'_tof']
### Simion x-> actual z, simion y-> actual x, simion z actual y
for j in range(0,3):
    data = np.loadtxt(filename_all[j]+'.csv',skiprows=9,delimiter=',')
    val_min=int(np.min(data[:,0]))
    val_max=int(np.max(data[:,0]))
    vx=[]
    vy=[]
    vz=[]
    KE=[]
    x=[]
    y=[]
    tof=[]
    for i in range(val_min,val_max+1):
        vx.append(data[2*(i-1),9])
        vy.append(data[2*(i-1),10])
        vz.append(data[2*(i-1),8])
        KE.append(data[2*(i-1),11])
        x.append(data[2*i-1,6])
        y.append(data[2*i-1,7])
        tof.append(data[2*i-1,2]*1000)   # converting tof from us to ns
        
    
#    vx=np.asarray(vx)
    x=np.asarray(x)
    y=np.asarray(y)
    tof=np.asarray(tof)
#    to
    if j==0:
        plt.figure()
        plt.plot(x,vx,'o')
        plt.xlabel('X(mm)')
        plt.ylabel('Vx (mm/us)')
        m,b = np.polyfit(x,vx,1)
        equation = r'$Vx = $' + str(round(m,4))+ r'X' ' + ' + str(round(b,4))
        plt.plot(x,m*x+b, '-',linewidth=2, color='r')
        plt.title(equation,fontsize=12)
    if j==1:
        plt.figure()
        plt.plot(y,vy,'o')
        plt.xlabel('Y(mm)')
        plt.ylabel('Vy (mm/us)')
        m,b = np.polyfit(y,vy,1)
        equation = r'$Vy = $' + str(round(m,4))+ r'Y' ' + ' + str(round(b,4))
        plt.plot(y,m*y+b, '-',linewidth=2, color='r')
        plt.title(equation,fontsize=12)
    if j==2:
        plt.figure()
        plt.plot(tof,vz,'o')
        plt.xlabel('TOF(ns)')
        plt.ylabel('Vz (mm/us)')
        m,b = np.polyfit(tof,vz,1)
        equation = r'$Vz = $' + str(round(m,4))+ r' TOF' ' + ' + str(round(b,4))
        plt.plot(tof,m*tof+b, '-',linewidth=2, color='r')
        plt.title(equation,fontsize=12)
