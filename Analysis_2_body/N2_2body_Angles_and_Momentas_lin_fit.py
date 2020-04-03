# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:39:37 2019

@author: Shashank Pathak
"""

import numpy as np
import matplotlib.pyplot as plt
from time import time
import h5py
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'C:/Users/patha/Documents/ALS_Nov_2019_Analysis/Codes/')
from useful_definitions import fhist1d,fhist2d,cosf1f2,NP,cosAB_keAB
plt.rcParams.update({'font.size': 16})



###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"Change these parameters for a new channel/run"
Run_name='N2_160eV'
channel='N_N'
label_species = ['N$^+$','N$^+$']
###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

basedir1='C:/Users/patha/Documents/ALS_Nov_2019_Analysis/'+Run_name+'/'+channel+'/'

basedir=basedir1+channel
processed_dir=basedir1+'Processed/'

hf1=h5py.File(processed_dir+Run_name+'_vel_2body_v2.h5','r')
chunks=np.array(hf1['chunks'])
mass_species=np.array(hf1['masses'])
pipico_range=np.array(hf1['pipico_range'])
pipicorot_range=np.array(hf1['pipicorot_range'])
tof_x_range=np.array(hf1['tof_x_range'])
tof_frag=np.array(hf1['tof_frag'])

Charge_positions_known=False
dist=3.229 # In Angstrom
if Charge_positions_known:
    KER_th=(1.6**2)*(8.987)*(0.624)/(dist)


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
x1all=[]
x2all=[]
x3all=[]
y1all=[]
y2all=[]
y3all=[]
t1all=[]
t2all=[]
t3all=[]
for i in range(1,int(chunks)+1):
    v1x=np.array(hf1['v1x/%d'%i])
    v1y=np.array(hf1['v1y/%d'%i])
    v1z=np.array(hf1['v1z/%d'%i])
    v2x=np.array(hf1['v2x/%d'%i])
    v2y=np.array(hf1['v2y/%d'%i])
    v2z=np.array(hf1['v2z/%d'%i])

    x1=np.array(hf1['x1/%d'%i])
    x2=np.array(hf1['x2/%d'%i])

    y1=np.array(hf1['y1/%d'%i])
    y2=np.array(hf1['y2/%d'%i])

    t1=np.array(hf1['t1/%d'%i])
    t2=np.array(hf1['t2/%d'%i])


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

    x1all=np.append(x1all,x1)
    x2all=np.append(x2all,x2)

    y1all=np.append(y1all,y1)
    y2all=np.append(y2all,y2)

    t1all=np.append(t1all,t1)
    t2all=np.append(t2all,t2)

    

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
# P_mid=np.add(P1sq,P2sq).sum(axis=0)
# P=np.sqrt(P_mid)

Pmid=np.array((Px_all**2,Py_all**2,Pz_all**2))
P=np.sqrt(np.sum(Pmid,axis=0))
cosz=np.divide(Pz_all,P)

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
KER_cosz=True
plot_det_images=True


par=np.loadtxt(processed_dir+Run_name+'_'+channel+'_final_parameters.txt')
# x1_c=par[0]
# x2_c=par[2]
# y1_c=par[1]
# y2_c=par[3]
x1_c=0
x2_c=0
y1_c=0
y2_c=0
xc = np.linspace(-50, 50, 100)
yc = np.linspace(-50, 50, 100)
X, Y = np.meshgrid(xc,yc)
R=800

# F2 = (X-0)**2 + (Y-0)**2 - 400
# plt.contour(X, Y, F1,[0], colors='r', linewidths=3)
cmap='jet'
if plot_det_images:
    fig = plt.figure(figsize=(18,4))

    x,y,z=fhist2d(x1all,y1all,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(121)
    f1=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f1)
    F1 = (X-x1_c)**2 + (Y-y1_c)**2 - R
    ax.contour(X, Y, F1,[0], colors='r', linewidths=3)
#    ax.set_title(label_species[0])
    ax.set_aspect('equal')

    x,y,z=fhist2d(x2all,y2all,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(122)
    f2=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f2)
    F1 = (X-x2_c)**2 + (Y-y2_c)**2 - R
    ax.contour(X, Y, F1,[0], colors='r', linewidths=3)
#    ax.set_title(label_species[1])
    ax.set_aspect('equal')








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


cosrange=[-1,1,0.01]
if Plot_angles:
    cos12=cosf1f2(P1x,P1y,P1z,P2x,P2y,P2z)

    fig = plt.figure(figsize=(16,4))
    ax = fig.add_subplot(111)
    x,y=fhist1d(cos12,cosrange[0],cosrange[1],cosrange[2])
    ax.plot(x,y)
    ax.set_title(r'Cos$\theta$(%s,%s)'%(label_species[0],label_species[1]))
    ax.set_xlabel(r'Cos$\theta$')
    ax.set_title('%s + %s'%(label_species[0],label_species[1]))



# KErange=[0,40,0.2] # last ppt
KErange=[0,40,0.4]
if Plot_KE:
    # x,y=fhist1d(KE1,KErange[0],KErange[1],KErange[2])
    fig, ax = plt.subplots()
    # ax.plot(x,y,label=label_species[0])
    x,y=fhist1d(KE2,KErange[0],KErange[1],KErange[2])
    # ax.plot(x,y,label=label_species[1])
    x,y=fhist1d(np.add(KE1,KE2),KErange[0],KErange[1],KErange[2])
    KERx=x
    KERy=y
    ax.step(x,y,label='KER')
    if Charge_positions_known:
        ax.axvline(x=KER_th,color='g',linestyle='--',label='KER (CEI model) ')
    ax.legend(frameon=False)
    ax.set_xlabel('KER (eV)')
    ax.set_ylabel('Counts (arb units)')
    ax.set_title('%s + %s'%(label_species[0],label_species[1]))
    plt.savefig(basedir1+'no_mom_gate_KER.png',bbox_inches='tight')

    x,y=fhist1d(KE1,KErange[0],KErange[1],KErange[2])
    fig, ax = plt.subplots()
    ax.plot(x,y,label=label_species[0])
    x,y=fhist1d(KE2,KErange[0],KErange[1],KErange[2])
    ax.plot(x,y,label=label_species[1])
    x,y=fhist1d(np.add(KE1,KE2),KErange[0],KErange[1],KErange[2])
    KERx=x
    KERy=y
    ax.step(x,y,label='KER')
    if Charge_positions_known:
        ax.axvline(x=KER_th,color='g',linestyle='--',label='KER (CEI model) ')
    ax.legend(frameon=False)
    ax.set_xlabel('KER (eV)')
    ax.set_ylabel('Counts (arb units)')
    ax.set_title('%s + %s'%(label_species[0],label_species[1]))

    
    
PKER_range=[0,200,0.5,0,40,0.1]
if Psum_KER:
    cmap='jet'
    x,y,z=fhist2d(absPsum,KER,PKER_range[0],PKER_range[1],PKER_range[2],PKER_range[3],PKER_range[4],PKER_range[5])
    fig, ax = plt.subplots()
    f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap)#,norm=LogNorm())
    fig.colorbar(f1)
#    ax.pcolormesh(x,y,z,label='KE1')
    ax.set_xlabel('Abs P (a.u.)')
    ax.set_ylabel('KER (eV)')
    ax.set_title('%s + %s'%(label_species[0],label_species[1]))

    
KER_cosz_range=[0,40,0.1,-1.1,1.1,0.05]
if KER_cosz:
    cmap='jet'
    x,y,z=fhist2d(KER,cosz,KER_cosz_range[0],KER_cosz_range[1],KER_cosz_range[2],KER_cosz_range[3],KER_cosz_range[4],KER_cosz_range[5])
    fig, ax = plt.subplots()
    f1=ax.pcolormesh(x,y,np.transpose(z),norm=LogNorm(),cmap=cmap)
    fig.colorbar(f1)
    #    ax.pcolormesh(x,y,z,label='KE1')
    ax.set_ylabel('cos theta')
    ax.set_xlabel('KER (eV)')
    ax.set_title('%s + %s'%(label_species[0],label_species[1]))
    

#####################################################################################
##################################### GATING ########################################
#####################################################################################
Gating=True

"Put the conditions here"
# condition=((absPsum<600))
# condition=((Py_all<5) & (Py_all>-5) & (Px_all<5) & (Px_all>-5) & (Pz_all<8) & (Pz_all>-8))#(Px_all<3) & (Px_all>-3) & 
# condition=((cos12<-0.99))
r1=np.sqrt(x1all**2 + y1all**2)
condition=((r1<5))#(Px_all<3) & (Px_all>-3) &
if Gating:
    condition=condition
    x1g=x1[condition]
    x2g=x2[condition]
    y1g=y1[condition]
    y2g=y2[condition]
    t1g=t1all[condition]
    t2g=t2all[condition]
    P1xg=P1x[condition]
    P1yg=P1y[condition]
    P1zg=P1z[condition]
    P2xg=P2x[condition]
    P2yg=P2y[condition]
    P2zg=P2z[condition]
    KE1g=KE1[condition]
    KE2g=KE2[condition]
    Px_allg=Px_all[condition]
    Py_allg=Py_all[condition]
    Pz_allg=Pz_all[condition]
    absPsumg=absPsum[condition]
    KERg=KER[condition]
    coszg=cosz[condition]
    
    
    mom_KE_arr=[P1xg,P1yg,P1zg,P2xg,P2yg,P2zg,KE1g,KE2g,Px_allg,Py_allg,Pz_allg,absPsumg,KERg]
    pos_tof_arr=[x1g,x2g,y1g,y2g,t1g,t2g]
    
    
    
    from basic_plots_2body import basic_plots
    KERgx,KERgy=basic_plots(mom_KE_arr,pos_tof_arr,coszg,pipico_range,pipicorot_range,tof_x_range,tof_frag,prange,cosrange,KErange,PKER_range,KER_cosz_range,mass_species,label_species,basedir,Run_name)


"Write final histogram arrays to h5 file for comparison later"

"-To add a new entry first create x and y arrays (like KERx and KERy) both in this file and in basic plots"

with h5py.File(processed_dir+Run_name+'_'+channel+'_final_hist_arrays_2body.h5','w') as f:
    if Plot_KE:
        f.create_dataset('label_species', data = np.string_(label_species), maxshape=(None), dtype='S20')  
        f.create_dataset('KER', data = np.transpose([KERx,KERy]), maxshape=(None), dtype='f4')
    # f.create_dataset('gating/condition', data = np.string_(condition), maxshape=(None), dtype='S10')
    f.create_dataset('gating/KER', data = np.transpose([KERgx,KERgy]), maxshape=(None), dtype='f4')
f.close()
