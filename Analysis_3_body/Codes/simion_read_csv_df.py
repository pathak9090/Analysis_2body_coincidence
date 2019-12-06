# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 18:20:21 2019

@author: Shashank
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#filename='CHOH_training'
#channel='CH3_C3H2_C3H3'
#filename='C3H2'
#Run_name='Toulene_400eV'

fragments=['CH2','C2H4','S']
Run_name='TMS_220eV'


channel=fragments[0]+'_'+fragments[1]+'_'+fragments[2]
#channel='CH2_C2H3_S'
basedir='C:/Users/patha/Documents/ALS_Nov_2019_Analysis/'
rep_dir=basedir+'simion_repository/'
#basedir1='C:/Users/patha/Documents/ALS_19th_March_start/'
basedir1=basedir+Run_name+'/'+channel+'/'
processed_dir=basedir1+'Processed/'

final_array=[]
for k in range(0,np.size(fragments)):
    filename=fragments[k]
    filename_all=['_training_vx','_training_vy','_training_tof']
    ### Simion x-> actual z, simion y-> actual x, simion z actual y
    for j in range(0,3):
        data = np.loadtxt(rep_dir+filename+filename_all[j]+'.csv',skiprows=8,delimiter=',')
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
            plt.plot(x,m*x+b, '-',linewidth=2, color='r',label=fragments[k])
            plt.title(equation,fontsize=12)
            plt.legend(frameon=False)
            plt.savefig(rep_dir+filename+'_Vx.png',bbox_inches='tight')
        if j==1:
            plt.figure()
            plt.plot(y,vy,'o')
            plt.xlabel('Y(mm)')
            plt.ylabel('Vy (mm/us)')
            m,b = np.polyfit(y,vy,1)
            equation = r'$Vy = $' + str(round(m,4))+ r'Y' ' + ' + str(round(b,4))
            plt.plot(y,m*y+b, '-',linewidth=2, color='r',label=fragments[k])
            plt.title(equation,fontsize=12)
            plt.legend(frameon=False)
            plt.savefig(rep_dir+filename+'_Vy.png',bbox_inches='tight')
        if j==2:
            plt.figure()
            plt.plot(tof,vz,'o')
            plt.xlabel('TOF(ns)')
            plt.ylabel('Vz (mm/us)')
            m,b = np.polyfit(tof,vz,1)
            equation = r'$Vz = $' + str(round(m,4))+ r' TOF' ' + ' + str(round(b,4))
            plt.plot(tof,m*tof+b, '-',linewidth=2, color='r',label=fragments[k])
            plt.title(equation,fontsize=12)
            plt.legend(frameon=False)
            plt.savefig(rep_dir+filename+'_vz.png',bbox_inches='tight')
        if (j==0) or (j==2):
            final_array.append(m)
            final_array.append(b)

np.savetxt(processed_dir+'vel_calc_fac_'+channel+'.txt',np.asarray(final_array))
