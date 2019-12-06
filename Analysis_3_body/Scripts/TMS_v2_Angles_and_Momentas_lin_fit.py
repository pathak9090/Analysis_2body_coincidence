# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:39:37 2019

@author: shashank
"""

import numpy as np
import matplotlib.pyplot as plt
from time import time
import h5py
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'C:/Users/patha/Documents/ALS_Nov_2019_Analysis/Codes/')
import useful_definitions
from useful_definitions import fhist1d,fhist2d,cosf1f2,NP,cosAB_keAB
plt.rcParams.update({'font.size': 16})


###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"Change these parameters for a new channel/run"
Run_name='TMS_220eV'
channel='CH2_C2H3_S'
label_species = ['CH$_2$$^+$','C$_2$H$_3$$^+$','S$^+$']
###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


basedir1='C:/Users/patha/Documents/ALS_Nov_2019_Analysis/'+Run_name+'/'+channel+'/'
basedir=basedir1+channel
processed_dir=basedir1+'Processed/'

hf1=h5py.File(processed_dir+Run_name+'_vel_KE_3body_v2.h5','r')
chunks=np.array(hf1['chunks'])
mass_species=np.array(hf1['masses'])
v1xall=[]
v1yall=[]
v1zall=[]
v2xall=[]
v2yall=[]
v2zall=[]
v3xall=[]
v3yall=[]
v3zall=[]
KE1all=[]
KE2all=[]
KE3all=[]
for i in range(1,int(chunks)+1):
    v1x=np.array(hf1['v1x/%d'%i])
    v1y=np.array(hf1['v1y/%d'%i])
    v1z=np.array(hf1['v1z/%d'%i])
    v2x=np.array(hf1['v2x/%d'%i])
    v2y=np.array(hf1['v2y/%d'%i])
    v2z=np.array(hf1['v2z/%d'%i])
    v3x=np.array(hf1['v3x/%d'%i])
    v3y=np.array(hf1['v3y/%d'%i])
    v3z=np.array(hf1['v3z/%d'%i])
#    KE1=np.array(hf1['KE1/%d'%i])
#    KE2=np.array(hf1['KE2/%d'%i])
#    KE3=np.array(hf1['KE3/%d'%i])
    
    # Append if more than 1 chunk present
    v1xall=np.append(v1xall,v1x)
    v1yall=np.append(v1yall,v1y)
    v1zall=np.append(v1zall,v1z)
    v2xall=np.append(v2xall,v2x)
    v2yall=np.append(v2yall,v2y)
    v2zall=np.append(v2zall,v2z)
    v3xall=np.append(v3xall,v3x)
    v3yall=np.append(v3yall,v3y)
    v3zall=np.append(v3zall,v3z)
#    KE1all=np.append(KE1all,KE1)
#    KE2all=np.append(KE2all,KE2)
#    KE3all=np.append(KE3all,KE3)
hf1.close()
print('Successfully read velocities from h5 file!')

#%% "REST FRAME ANALYSIS"


rest_frame_analysis=False


if rest_frame_analysis:
    "Rest frame of 12 -  subtract velocity of 3"
    
    v1xcom12=np.add(v1xall,v3xall)
    v1ycom12=np.add(v1yall,v3yall)
    v1zcom12=np.add(v1zall,v3zall)
    v2xcom12=np.add(v2xall,v3xall)
    v2ycom12=np.add(v2yall,v3yall)
    v2zcom12=np.add(v2zall,v3zall)
    
    
    P1xcom12 = 1e3 * 0.5E24 * mass_species[0] * v1xcom12
    P1ycom12 = 1e3 * 0.5E24 * mass_species[0] * v1ycom12
    P1zcom12 = 1e3 * 0.5E24 * mass_species[0] * v1zcom12
    
    
    P2xcom12 = 1e3 * 0.5E24 * mass_species[1] * v2xcom12
    P2ycom12 = 1e3 * 0.5E24 * mass_species[1] * v2ycom12
    P2zcom12 = 1e3 * 0.5E24 * mass_species[1] * v2zcom12
    
    "Rest frame of 13 -  subtract velocity of 2"
    
    #v1xcom13=np.add(v1xall,v3xall)
    #v1ycom13=np.add(v1yall,v3yall)
    #v1zcom13=np.add(v1zall,v3zall)
    #v2xcom13=np.add(v2xall,v3xall)
    #v2ycom13=np.add(v2yall,v3yall)
    #v2zcom13=np.add(v2zall,v3zall)
    v1xcom13=np.add(v1xall,v2xall)
    v1ycom13=np.add(v1yall,v2yall)
    v1zcom13=np.add(v1zall,v2zall)
    v2xcom13=np.add(v3xall,v2xall)
    v2ycom13=np.add(v3yall,v2yall)
    v2zcom13=np.add(v3zall,v2zall)
    
    P1xcom13 = 1e3 * 0.5E24 * mass_species[0] * v1xcom13
    P1ycom13 = 1e3 * 0.5E24 * mass_species[0] * v1ycom13
    P1zcom13 = 1e3 * 0.5E24 * mass_species[0] * v1zcom13
    
    
    P2xcom13 = 1e3 * 0.5E24 * mass_species[2] * v2xcom13
    P2ycom13 = 1e3 * 0.5E24 * mass_species[2] * v2ycom13
    P2zcom13 = 1e3 * 0.5E24 * mass_species[2] * v2zcom13
#%%
#Calculating Momentas

P1x = 1e3 * 0.5E24 * mass_species[0] * v1xall
P1y = 1e3 * 0.5E24 * mass_species[0] * v1yall
P1z = 1e3 * 0.5E24 * mass_species[0] * v1zall
P2x = 1e3 * 0.5E24 * mass_species[1] * v2xall
P2y = 1e3 * 0.5E24 * mass_species[1] * v2yall
P2z = 1e3 * 0.5E24 * mass_species[1] * v2zall
P3x = 1e3 * 0.5E24 * mass_species[2] * v3xall
P3y = 1e3 * 0.5E24 * mass_species[2] * v3yall
P3z = 1e3 * 0.5E24 * mass_species[2] * v3zall
Px=np.array((P1x,P2x,P3x))
Py=np.array((P1y,P2y,P3y))
Pz=np.array((P1z,P2z,P3z))
#KE=np.array((KE1,KE2,KE3))
Px_all=Px.sum(axis=0)
Py_all=Py.sum(axis=0)
Pz_all=Pz.sum(axis=0)

Psum_arr=np.array((Px_all,Py_all,Pz_all))
Psum=Psum_arr.sum(axis=0)
absPsum=np.abs(Psum)

P1sq=np.array((P1x**2,P1y**2,P1z**2))
P2sq=np.array((P2x**2,P2y**2,P2z**2))
P3sq=np.array((P3x**2,P3y**2,P3z**2))

KE1=((P1sq.sum(axis=0)/(0.25E48))/(2*mass_species[0]))*6.242E18            # Momentun square is converted to SI units and multipled by factor 6.242E18 to convert SI units -> eV
KE2=((P2sq.sum(axis=0)/(0.25E48))/(2*mass_species[1]))*6.242E18
KE3=((P3sq.sum(axis=0)/(0.25E48))/(2*mass_species[2]))*6.242E18

KER_arr=np.array((KE1,KE2,KE3))
KER=KER_arr.sum(axis=0)
#================================================================================================
#++++++++++++++++++++++++++++++++++++ PLOTS++++++++++++++++++++++++++++++++++++++++++++++++++++++
#================================================================================================
Plot_momentum_panel=True
Plot_angles=True
Plot_NP=True
Plot_KE=True
Psum_KER=True
Plot_KEcos=True
prange=[-400,400,2]

if Plot_momentum_panel:
    plt.rcParams.update({'font.size': 13})
    fig = plt.figure(figsize=(22,16))
    ax = fig.add_subplot(431)
    x,y=fhist1d(P1x,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Px (%s)'%label_species[0])
    ax.legend(frameon=False)
    ax = fig.add_subplot(432)
    x,y=fhist1d(P1y,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Py (%s)'%label_species[0])
    ax.legend(frameon=False)
    ax = fig.add_subplot(433)
    x,y=fhist1d(P1z,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Pz (%s)'%label_species[0])
    ax.legend(frameon=False)
    ax = fig.add_subplot(434)
    x,y=fhist1d(P2x,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Px (%s)'%label_species[1])
    ax.legend(frameon=False)
    ax = fig.add_subplot(435)
    x,y=fhist1d(P2y,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Py (%s)'%label_species[1])
    ax.legend(frameon=False)
    ax = fig.add_subplot(436)
    x,y=fhist1d(P2z,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Pz (%s)'%label_species[1])
    ax.legend(frameon=False)
    ax = fig.add_subplot(437)
    x,y=fhist1d(P3x,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Px (%s)'%label_species[2])
    ax.legend(frameon=False)
    ax = fig.add_subplot(438)
    x,y=fhist1d(P3y,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Py (%s)'%label_species[2])
    ax.legend(frameon=False)
    ax = fig.add_subplot(439)
    x,y=fhist1d(P3z,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Pz (%s)'%label_species[2])
    ax.legend(frameon=False)
    ax = fig.add_subplot(4,3,10)
    x,y=fhist1d(Px_all,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Px Sum')
    ax.legend(frameon=False)
    ax = fig.add_subplot(4,3,11)
    x,y=fhist1d(Py_all,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Py Sum')
    ax.legend(frameon=False)
    ax = fig.add_subplot(4,3,12)
    x,y=fhist1d(Pz_all,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Pz Sum')
    ax.legend(frameon=False)
    plt.savefig(basedir+'momentum_panel.png',bbox_inches='tight')


cosrange=[-1,1,0.05]
if Plot_angles:
    cos12=cosf1f2(P1x,P1y,P1z,P2x,P2y,P2z)
    cos23=cosf1f2(P2x,P2y,P2z,P3x,P3y,P3z)
    cos13=cosf1f2(P1x,P1y,P1z,P3x,P3y,P3z)
    fig = plt.figure(figsize=(16,4))
    ax = fig.add_subplot(131)
    x,y=fhist1d(cos12,cosrange[0],cosrange[1],cosrange[2])
    ax.plot(x,y)
    ax.set_title(r'Cos$\theta$(%s,%s)'%(label_species[0],label_species[1]))
    ax.set_xlabel(r'Cos$\theta$')
    ax = fig.add_subplot(132)
    x,y=fhist1d(cos23,cosrange[0],cosrange[1],cosrange[2])
    ax.plot(x,y)
    ax.set_title(r'Cos$\theta$(%s,%s)'%(label_species[1],label_species[2]))
    ax.set_xlabel(r'Cos$\theta$')
    ax = fig.add_subplot(133)
    x,y=fhist1d(cos13,cosrange[0],cosrange[1],cosrange[2])
    ax.plot(x,y)
    ax.set_title(r'Cos$\theta$(%s,%s)'%(label_species[0],label_species[2]))
    ax.set_xlabel(r'Cos$\theta$')
    plt.savefig(basedir+'angles.png',bbox_inches='tight')
    
    
if Plot_NP:
    NP(P1x,P1y,P1z,label_species[0],P2x,P2y,P2z,label_species[1],P3x,P3y,P3z,label_species[2],1,-10,10,0.08,basedir+'_NP1.png',0,0,Run_name)
    NP(P1x,P1y,P1z,label_species[0],P2x,P2y,P2z,label_species[1],P3x,P3y,P3z,label_species[2],2,-2,2,0.01,basedir+'_NP2.png',0,0,Run_name)
    NP(P1x,P1y,P1z,label_species[0],P2x,P2y,P2z,label_species[1],P3x,P3y,P3z,label_species[2],3,-2,2,0.01,basedir+'_NP3.png',0,0,Run_name)
    
    
KErange=[0,40,0.25]
if Plot_KE:
    x,y=fhist1d(KE1,KErange[0],KErange[1],KErange[2])
    fig, ax = plt.subplots()
    ax.plot(x,y,label='KE1')
    x,y=fhist1d(KE2,KErange[0],KErange[1],KErange[2])
    ax.plot(x,y,label='KE2')
    x,y=fhist1d(KE3,KErange[0],KErange[1],KErange[2])
    ax.plot(x,y,label='KE3')
    ax.legend(frameon=False)
    ax.set_xlabel('KE (eV)')
    ax.set_ylabel('Counts (arb units)')
    plt.savefig(basedir+'KE_all_frag.png',bbox_inches='tight')
    
PKER_range=[0,200,2,0,40,0.25]
if Psum_KER:
    cmap='jet'
    x,y,z=fhist2d(absPsum,KER,PKER_range[0],PKER_range[1],PKER_range[2],PKER_range[3],PKER_range[4],PKER_range[5])
    fig, ax = plt.subplots()
    f1=ax.pcolormesh(x,y,np.transpose(z),norm=LogNorm(),cmap=cmap)
    fig.colorbar(f1)
#    ax.pcolormesh(x,y,z,label='KE1')
    ax.set_xlabel('Abs P (a.u.)')
    ax.set_ylabel('KER (eV)')
    plt.savefig(basedir+'AbsP_KER.png',bbox_inches='tight')
    
cmap='jet'
keABcosABrange=[0,40,0.25,0,180,2]
mom_all=[P1x,P1y,P1z,P2x,P2y,P2z,P3x,P3y,P3z]
if Plot_KEcos:
    KE12,cosH12H3=cosAB_keAB(mom_all,mass_species,3,keABcosABrange,label_species,cmap,basedir,log=True)
    KE13,cosH13H2=cosAB_keAB(mom_all,mass_species,2,keABcosABrange,label_species,cmap,basedir,log=True)
    KE23,cosH23H1=cosAB_keAB(mom_all,mass_species,1,keABcosABrange,label_species,cmap,basedir,log=True)

#%% Tests
if rest_frame_analysis:
    "COM analysis -  Calculating momenta and energy of P1 and P3 in com frame of P13"
    
    #P1xcom=np.add(P1x,P2x)
    #P1ycom=np.add(P1y,P2y)
    #P1zcom=np.add(P1z,P2z)
    #P3xcom=np.add(P3x,P2x)
    #P3ycom=np.add(P3y,P2y)
    #P3zcom=np.add(P3z,P2z)
    
    P1comsq12=np.array((P1xcom12**2,P1ycom12**2,P1zcom12**2))
    P2comsq12=np.array((P2xcom12**2,P2ycom12**2,P2zcom12**2))
    
    P1comsq13=np.array((P1xcom13**2,P1ycom13**2,P1zcom13**2))
    P2comsq13=np.array((P2xcom13**2,P2ycom13**2,P2zcom13**2))
    
    KE1com12=((P1comsq12.sum(axis=0)/(0.25E48))/(2*mass_species[0]))*6.242E18            # Momentun square is converted to SI units and multipled by factor 6.242E18 to convert SI units -> eV
    KE2com12=((P2comsq12.sum(axis=0)/(0.25E48))/(2*mass_species[1]))*6.242E18
    KE1com13=((P1comsq13.sum(axis=0)/(0.25E48))/(2*mass_species[0]))*6.242E18
    KE2com13=((P2comsq13.sum(axis=0)/(0.25E48))/(2*mass_species[2]))*6.242E18
    
    KERcom_arr12=np.array((KE1com12,KE2com12))
    KERcom_arr13=np.array((KE1com13,KE2com13))
    KERcom12=KERcom_arr12.sum(axis=0)
    KERcom13=KERcom_arr13.sum(axis=0)
    Plot_Pcom=True
    Plot_KEcom=True
    Plot_coscom=True
    
    ptestx=np.add(P1xcom13,P2xcom13)
    ptesty=np.add(P1ycom13,P2ycom13)
    ptestz=np.add(P1zcom13,P2zcom13)
    #pt=np.add(ptestx,ptesty)
    #ptt=np.add(pt,ptestz)
    
    psum13=np.sum([P1xcom13,P2xcom13,P1ycom13,P2ycom13,P1zcom13,P2zcom13],axis=0)
    psum12=np.sum([P1xcom12,P2xcom12,P1ycom12,P2ycom12,P1zcom12,P2zcom12],axis=0)
    if Plot_Pcom:
        x,y=fhist1d(psum12,prange[0],prange[1],prange[2])
        fig, ax = plt.subplots()
        ax.plot(x,y,label='P (%s-%s) com'%(label_species[0],label_species[1]))
        x,y=fhist1d(psum13,prange[0],prange[1],prange[2])
        ax.plot(x,y,label='P (%s-%s) com'%(label_species[0],label_species[2]))
    #    x,y=fhist1d(KE2,KErange[0],KErange[1],KErange[2])
    #    ax.plot(x,y,label='KE2')
    #    x,y=fhist1d(KE3com,KErange[0],KErange[1],KErange[2])
    #    ax.plot(x,y,label='KE (HCHO+) com')
    ##    ax.plot(x,y,label='KE (C2H+) com')
        ax.legend(frameon=False)
        ax.set_xlabel('Pcom (a.u.)')
        ax.set_ylabel('Counts (arb units)')
        ax.set_title('Rest frame method')
        plt.savefig(basedir+'P_com',bbox_inches='tight')
        
    #    x,y=fhist1d(ptt,prange[0],prange[1],prange[2])
    #    fig, ax = plt.subplots()
    #    ax.plot(x,y,label='P')
    #    
    #    x,y=fhist1d(ptesty,prange[0],prange[1],prange[2])
    #    fig, ax = plt.subplots()
    #    ax.plot(x,y,label='P')
    
    
    KErange=[0,40,0.25]
    if Plot_KEcom:
        x,y=fhist1d(KERcom12,KErange[0],KErange[1],KErange[2])
        fig, ax = plt.subplots()
        ax.plot(x,y,label='KER (%s-%s) com'%(label_species[0],label_species[1]))
        x,y=fhist1d(KERcom13,KErange[0],KErange[1],KErange[2])
        ax.plot(x,y,label='KER (%s-%s) com'%(label_species[0],label_species[2]))
        ax.legend(frameon=False)
        ax.set_xlabel('KERcom (eV)')
        ax.set_ylabel('Counts (arb units)')
        ax.set_title('Rest frame KER method')
        plt.savefig(basedir+'KER_com',bbox_inches='tight')
        
    cos12com=cosf1f2(P1xcom12,P1ycom12,P1zcom12,P2xcom12,P2ycom12,P2zcom12)
    cos13com=cosf1f2(P1xcom13,P1ycom13,P1zcom13,P2xcom13,P2ycom13,P2zcom13)
    
    if Plot_coscom:
        x,y=fhist1d(cos12com,cosrange[0],cosrange[1],cosrange[2])
        fig, ax = plt.subplots()
        ax.plot(x,y,label=r'Cos $\theta$ (%s-%s) com'%(label_species[0],label_species[1]))
        x,y=fhist1d(cos13com,cosrange[0],cosrange[1],cosrange[2])
        ax.plot(x,y,label=r'Cos $\theta$ (%s-%s) com'%(label_species[0],label_species[2]))
        ax.legend(frameon=False)
        plt.title('Rest frame KER method')
        plt.savefig(basedir+'momenta_angle_second_step.png',bbox_inches='tight')


#    
#####################################################################################
##################################### GATING ########################################
#####################################################################################
Gating=False

"Put the conditions here"
condition=((absPsum<13))
if Gating:
    condition=condition
    P1xg=P1x[condition]
    P1yg=P1y[condition]
    P1zg=P1z[condition]
    P2xg=P2x[condition]
    P2yg=P2y[condition]
    P2zg=P2z[condition]
    P3xg=P3x[condition]
    P3yg=P3y[condition]
    P3zg=P3z[condition]
    KE1g=KE1[condition]
    KE2g=KE2[condition]
    KE3g=KE3[condition]
    Px_allg=Px_all[condition]
    Py_allg=Py_all[condition]
    Pz_allg=Pz_all[condition]
    absPsumg=absPsum[condition]
    KERg=KER[condition]
    
    
    mom_KE_arr=[P1xg,P1yg,P1zg,P2xg,P2yg,P2zg,P3xg,P3yg,P3zg,KE1g,KE2g,KE3g,Px_allg,Py_allg,Pz_allg,absPsumg,KERg]
    
    
    
    from basic_plots import basic_plots
    basic_plots(mom_KE_arr,prange,cosrange,KErange,PKER_range,keABcosABrange,mass_species,label_species,basedir,Run_name)