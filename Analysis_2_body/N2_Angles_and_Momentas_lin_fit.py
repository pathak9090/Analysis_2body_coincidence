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
from useful_definitions import fhist1d,fhist2d,cosf1f2,NP,cosAB_keAB
plt.rcParams.update({'font.size': 16})



#Run_name='Propargyl_alcohol_140ev_60C'
#Run_name='Propargyl_alcohol_140ev_160C'
#label_species = ['H$^+$','C$_2$H$^+$','HCHO$^+$']
#basedir1='C:/Users/patha/Documents/ALS_19th_March_start/'
#basedir=basedir1+Run_name+'/'

Run_name='N2_140eV'
basedir1='C:/Users/patha/Documents/ALS_19th_March_start/'
basedir=basedir1+Run_name+'/'

label_species = ['N$^+$','N$^+$']

hf1=h5py.File(basedir+Run_name+'_vel_2body_v2.h5','r')
chunks=np.array(hf1['chunks'])
mass_species=np.array(hf1['masses'])
v1xall=[]
v1yall=[]
v1zall=[]
v2xall=[]
v2yall=[]
v2zall=[]

KE1all=[]
KE2all=[]

for i in range(1,int(chunks)+1):
    v1x=np.array(hf1['v1x/%d'%i])
    v1y=np.array(hf1['v1y/%d'%i])
    v1z=np.array(hf1['v1z/%d'%i])
    v2x=np.array(hf1['v2x/%d'%i])
    v2y=np.array(hf1['v2y/%d'%i])
    v2z=np.array(hf1['v2z/%d'%i])

    
    # Append if more than 1 chunk present
    v1xall=np.append(v1xall,v1x)
    v1yall=np.append(v1yall,v1y)
    v1zall=np.append(v1zall,v1z)
    v2xall=np.append(v2xall,v2x)
    v2yall=np.append(v2yall,v2y)
    v2zall=np.append(v2zall,v2z)

#    KE1all=np.append(KE1all,KE1)
#    KE2all=np.append(KE2all,KE2)
#    KE3all=np.append(KE3all,KE3)
hf1.close()
print('Successfully read velocities from h5 file!')


#%%
#Calculating Momentas

P1x = 1e3 * 0.5E24 * mass_species[0] * v1xall
P1y = 1e3 * 0.5E24 * mass_species[0] * v1yall
P1z = 1e3 * 0.5E24 * mass_species[0] * v1zall
P2x = 1e3 * 0.5E24 * mass_species[1] * v2xall
P2y = 1e3 * 0.5E24 * mass_species[1] * v2yall
P2z = 1e3 * 0.5E24 * mass_species[1] * v2zall

Px=np.array((P1x,P2x))
Py=np.array((P1y,P2y))
Pz=np.array((P1z,P2z))

Px_all=Px.sum(axis=0)
Py_all=Py.sum(axis=0)
Pz_all=Pz.sum(axis=0)

Psum_arr=np.array((Px_all,Py_all,Pz_all))
Psum=Psum_arr.sum(axis=0)
absPsum=np.abs(Psum)

P1sq=np.array((P1x**2,P1y**2,P1z**2))
P2sq=np.array((P2x**2,P2y**2,P2z**2))


KE1=((P1sq.sum(axis=0)/(0.25E48))/(2*mass_species[0]))*6.242E18            # Momentun square is converted to SI units and multipled by factor 6.242E18 to convert SI units -> eV
KE2=((P2sq.sum(axis=0)/(0.25E48))/(2*mass_species[1]))*6.242E18


KER_arr=np.array((KE1,KE2))
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
    ax = fig.add_subplot(331)
    x,y=fhist1d(P1x,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Px (%s)'%label_species[0])
    ax.legend(frameon=False)
    ax = fig.add_subplot(332)
    x,y=fhist1d(P1y,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Py (%s)'%label_species[0])
    ax.legend(frameon=False)
    ax = fig.add_subplot(333)
    x,y=fhist1d(P1z,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Pz (%s)'%label_species[0])
    ax.legend(frameon=False)
    ax = fig.add_subplot(334)
    x,y=fhist1d(P2x,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Px (%s)'%label_species[1])
    ax.legend(frameon=False)
    ax = fig.add_subplot(335)
    x,y=fhist1d(P2y,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Py (%s)'%label_species[1])
    ax.legend(frameon=False)
    ax = fig.add_subplot(336)
    x,y=fhist1d(P2z,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Pz (%s)'%label_species[1])
    ax.legend(frameon=False)
    ax = fig.add_subplot(3,3,7)
    x,y=fhist1d(Px_all,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Px Sum')
    ax.legend(frameon=False)
    ax = fig.add_subplot(3,3,8)
    x,y=fhist1d(Py_all,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Py Sum')
    ax.legend(frameon=False)
    ax = fig.add_subplot(3,3,9)
    x,y=fhist1d(Pz_all,prange[0],prange[1],prange[2])
    ax.plot(x,y,label='Pz Sum')
    ax.legend(frameon=False)


cosrange=[-1,1,0.05]
if Plot_angles:
    cos12=cosf1f2(P1x,P1y,P1z,P2x,P2y,P2z)

    fig = plt.figure(figsize=(16,4))
    ax = fig.add_subplot(111)
    x,y=fhist1d(cos12,cosrange[0],cosrange[1],cosrange[2])
    ax.plot(x,y)
    ax.set_title(r'Cos$\theta$(%s,%s)'%(label_species[0],label_species[1]))
    ax.set_xlabel(r'Cos$\theta$')


KErange=[0,40,0.25]
if Plot_KE:
    x,y=fhist1d(KE1,KErange[0],KErange[1],KErange[2])
    fig, ax = plt.subplots()
    ax.plot(x,y,label='KE1')
    x,y=fhist1d(KE2,KErange[0],KErange[1],KErange[2])
    ax.plot(x,y,label='KE2')
    x,y=fhist1d(np.add(KE1,KE2),KErange[0],KErange[1],KErange[2])
    ax.plot(x,y,label='KER')
    ax.legend(frameon=False)
    ax.set_xlabel('KER (eV)')
    ax.set_ylabel('Counts (arb units)')
    
    
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