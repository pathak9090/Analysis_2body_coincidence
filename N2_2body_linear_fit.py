# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 22:19:34 2019

@author: Shashank
"""

import numpy as np
from useful_definitions_2body import fhist1d,fhist2d,twoConditions
import matplotlib.pyplot as plt
from time import time
#from sklearn.externals import joblib
import h5py
from matplotlib.colors import LogNorm
from optimise_parameters_2body import optimise_parameters
import tkinter as tk
plt.rcParams.update({'font.size': 16})


t1=time()
#Run_name='Propargyl_alcohol_140ev_60C'
Run_name='N2_140eV'
basedir1='C:/Users/patha/Documents/ALS_19th_March_start/'
basedir=basedir1+Run_name+'/'

path = basedir+'/Raw/'+Run_name+'_DAn.bin'
filename = path
print(filename)
f = open(filename)
#Read file all at once
print('Reading file...')
LineNum = 0
DataFile = np.fromfile(f, dtype = 'int32', count = -1)
f.close()
TotalNumLines = DataFile.size
print('File read successfully!')

mass_species = [14,14]
label_species = ['N$^+$','N$^+$']
mass_species = np.array(mass_species) * 1.66E-27
ionGate1 = [1690, 1800]
ionGate2 = [1760, 1840]
#ionGate2 = [2362,2412]
#
#ionGate3 = [2500,3000]



LineNum=0
ion1x=[]
ion1y=[]
ion2x=[]
ion2y=[]
#ion3x=[]
#ion3y=[]
ion1t=[]
ion2t=[]
#ion3t=[]

event_counter=0
all_evt_counter=0
while True:
#        if event_counter > 5000:
#            break
    #checks if end of file is reached
    if LineNum >= TotalNumLines:
        break

    #reading the hits in terms of x, y, tof
    numberOfHits = DataFile[LineNum]
    LineNum = LineNum +1
#        numberOfElectrons = ChunkData[LineNum]
#        LineNum = LineNum + 1

    ion = DataFile[LineNum:LineNum + numberOfHits * 3]
    LineNum = LineNum + numberOfHits * 3
    #dividing by factor since thats how the numbers were saved in terms of int32
    ion = ion.reshape((numberOfHits, 3 )).T/1000
    ionX = ion[0]
    ionY = ion[1]
    ionTOF = ion[2]

    checkCondition = twoConditions(ionTOF, ionGate1, ionGate2)
    all_evt_counter=all_evt_counter+1
    if checkCondition.sum() == 2:

        ionTOF = ionTOF[checkCondition]
        ionX = ionX[checkCondition]
        ionY = ionY[checkCondition]

#        ionX[0] = ionX[0] - pos_offset_1[0]
#        ionX[1] = ionX[1] - pos_offset_2[0]
#        ionX[2] = ionX[2] - pos_offset_3[0]
#
#        ionY[0] = ionY[0] - pos_offset_1[1]
#        ionY[1] = ionY[1] - pos_offset_2[1]
#        ionY[2] = ionY[2] - pos_offset_3[1]

        ion1x.append(ionX[0])
        ion1y.append(ionY[0])
        ion2x.append(ionX[1])
        ion2y.append(ionY[1])
        ion1t.append(ionTOF[0])
        ion2t.append(ionTOF[1])

        event_counter=event_counter+1
        if event_counter%10000 == 0:
            print(event_counter)

thbody_after_permutation=np.size(ion1x)
print('Total events = %d'%all_evt_counter)
print('Total 2 body events after permutation = %d'%thbody_after_permutation)
t2=time()

print("time for reading file and getting position and time arrays = %f sec"%(t2-t1))

ion1x=np.asarray(ion1x)
ion2x=np.asarray(ion2x)
ion1y=np.asarray(ion1y)
ion2y=np.asarray(ion2y)
ion1t=np.asarray(ion1t)
ion2t=np.asarray(ion2t)

#ion23t=np.add(ion2t,ion3t)
#trotplus23=np.add(ion23t,ion1t)
#trotminus23=np.add(ion23t,-ion1t)
#
#ion12t=np.add(ion1t,ion2t)
#trotplus12=np.add(ion12t,ion3t)
#trotminus12=np.add(ion12t,-ion3t)

#%%
cmap='jet'
## Histograms before selecting pipico channel
plot_det_images=True      # plots detector image for first, second hit
plot_all_pipico=True     # plots gated, non-gated pipico

# Gates for subplot 1 (tof1 vs tof2+tof3)
sp1_xmin=ionGate1[0]
sp1_xmax=ionGate1[1]
sp1_binsize=1
sp1_ymin=ionGate2[0]
sp1_ymax=ionGate2[1]


# Gates for subplot 2 (tof2+tof3-tof1 vs tof2+tof3+tof1)
#sp2_xmin=sp1_ymin-sp1_xmin
sp2_xmin=0
#sp2_xmax=sp1_ymax-sp1_xmax
sp2_xmax=90
sp2_binsize=1
sp2_ymin=sp1_ymin+sp1_xmin
sp2_ymax=sp1_ymax+sp1_xmax



## Time
if plot_all_pipico:
    ion1t=np.asarray(ion1t)
    ion2t=np.asarray(ion2t)

    fig = plt.figure(figsize=(18,4))
    ax = fig.add_subplot(121)
    x,y,z=fhist2d(ion1t,ion2t,sp1_xmin,sp1_xmax,sp1_binsize,sp1_ymin,sp1_ymax,sp1_binsize)
    f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    fig.colorbar(f1)
    ax.set_xlabel('TOF1 (ns)')
    ax.set_ylabel('TOF2 (ns)')
    ax.set_title('pipico')

    t12minus=np.add(ion2t,-ion1t)
    t12plus=np.add(ion2t,ion1t)
    
    ax = fig.add_subplot(122)
    x,y,z=fhist2d(t12minus,t12plus,sp2_xmin,sp2_xmax,sp2_binsize,sp2_ymin,sp2_ymax,sp2_binsize) 
    f2=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    fig.colorbar(f2)
    ax.set_xlabel('(TOF2-TOF1 (ns)')
    ax.set_ylabel('TOF2+TOF1 (ns)')
    ax.set_title('pipico rotated')
    plt.savefig(basedir+'pipico_rotated.png',bbox_inches='tight')


    ##Positions
if plot_det_images:
    fig = plt.figure(figsize=(18,4))

    x,y,z=fhist2d(ion1x,ion1y,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(121)
    f1=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f1)
#    ax.set_title(label_species[0])
    ax.set_aspect('equal')

    x,y,z=fhist2d(ion2x,ion2y,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(122)
    f2=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f2)
#    ax.set_title(label_species[1])
    ax.set_aspect('equal')


    plt.savefig(basedir+'det_image_non_gated_v1.png',bbox_inches='tight')


#%%
# Selection of tripico channel


t12rotGate = [30,62, 3555,3568]   # Select the channel of interest in rotated piipico!

condition=((t12minus>t12rotGate[0]) & (t12minus<t12rotGate[1]) & (t12plus>t12rotGate[2]) & (t12plus<t12rotGate[3]))

ion1tg=ion1t[condition]
ion2tg=ion2t[condition]

ion1xg=ion1x[condition]
ion2xg=ion2x[condition]

ion1yg=ion1y[condition]
ion2yg=ion2y[condition]

t12minusg=t12minus[condition]
t12plusg=t12plus[condition]
#%%
## Histograms after selecting pipico channel

plot_det_images_gated=True      # plots gated detector image 
plot_all_pipico_gated=True     # plots gated pipico
plot_tof_x_gated = True

## Time
if plot_all_pipico_gated:
    ion1tg=np.asarray(ion1tg)
    ion2tg=np.asarray(ion2tg)

    fig = plt.figure(figsize=(18,4))
    ax = fig.add_subplot(121)
    x,y,z=fhist2d(ion1tg,ion2tg,sp1_xmin,sp1_xmax,sp1_binsize,sp1_ymin,sp1_ymax,sp1_binsize)
    f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    fig.colorbar(f1)
    ax.set_xlabel('TOF1 (ns)')
    ax.set_ylabel('TOF2 (ns)')
    ax.set_title('pipico (gated)')
    
    ax = fig.add_subplot(122)
    x,y,z=fhist2d(t12minusg,t12plusg,sp2_xmin,sp2_xmax,sp2_binsize,sp2_ymin,sp2_ymax,sp2_binsize) 
    f2=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    fig.colorbar(f2)
    ax.set_xlabel('(TOF2-TOF1 (ns)')
    ax.set_ylabel('TOF2+TOF1 (ns)')
    ax.set_title('pipico rotated (gated)')
    plt.savefig(basedir+'pipico_rotated_gated.png',bbox_inches='tight')


    ##Positions
if plot_det_images_gated:
    fig = plt.figure(figsize=(18,4))

    x,y,z=fhist2d(ion1xg,ion1yg,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(121)
    f1=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f1)
    ax.set_title(label_species[0]+' (gated)')
    ax.set_aspect('equal')

    x,y,z=fhist2d(ion2xg,ion2yg,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(122)
    f2=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f2)
    ax.set_title(label_species[1]+' (gated)')
    ax.set_aspect('equal')
    plt.savefig(basedir+'det_image_gated.png',bbox_inches='tight')

if plot_tof_x_gated:
    all_x=np.concatenate((ion1yg,ion2yg),axis=0)
    all_tof=np.concatenate((ion1tg,ion2tg),axis=0)
    fig, ax = plt.subplots(figsize=(18,4))

    x,y,z=fhist2d(all_tof,all_x,0,4000,1,-50,50,0.5)
    f1=ax.pcolormesh(x,y,np.transpose(z),norm=LogNorm(),cmap=cmap)
    fig.colorbar(f1)
    plt.savefig(basedir+'tof_x_gated_v1.png',bbox_inches='tight')

#%% Loading functions for calculating velocities and KE

prange=[-500,500,1]

optimise_parameter=True  # Set time offset and position offset to make momentum sums peak at zero
save_vel_KE = False        # Save velocities and KE for all data in h5 file (Do it after optimizing parameters)

# Fit to get t0
tof = [459.1,3024.4,2131,2022,2527.7]
testmq = [1.008,40,20.1797,18.01,28]
_, t0 =  np.polyfit(np.sqrt(testmq), np.array(tof), 1)

fac=0.2
    #if tweaking in t0 is needed
# =============================================================================
t0 = t0 - fac*t0

# position offsets (should be adjusted based on momentum sums later!)
pos_offset_1 = [0.1,0.0]
pos_offset_2 = [2.1,0.3]


load_old_param=False
if load_old_param:
    par=np.loadtxt(basedir+Run_name+'final_parameters.txt')
    pos_offset_1[0],pos_offset_1[1],pos_offset_2[0],pos_offset_2[1],t0=par

param_list=[pos_offset_1[0],pos_offset_1[1],pos_offset_2[0],pos_offset_2[1],t0]


# Divide data in chunks if using Rbf else keep it 1

chunks=1
a=0
fac=np.int(np.size(ion1tg)/chunks)
if optimise_parameter:
        "Add GUI for optimization here"
        print("Now calculating velocities and KE for small chunck of data...")

        ion1xgc=np.asarray(ion1xg[a:fac])
        ion1ygc=np.asarray(ion1yg[a:fac])
        ion1tgc=np.asarray(ion1tg[a:fac])
        ion2xgc=np.asarray(ion2xg[a:fac])
        ion2ygc=np.asarray(ion2yg[a:fac])
        ion2tgc=np.asarray(ion2tg[a:fac])


        xyt_arr=[ion1xgc,ion1ygc,ion1tgc,ion2xgc,ion2ygc,ion2tgc]
        vel_calc_fac=[1.8399,-0.0707,-5.5784,2689.3848,0.368,-0.0141,-0.2231,537.8464]
        mx1,cx1,mz1,cz1,mx2,cx2,mz2,cz2=vel_calc_fac
        plt.rcParams.update({'font.size': 12})

        "Calling GUI..."

        top = tk.Tk()
        top.geometry("2000x2000")
        optimise_parameters(top,xyt_arr,vel_calc_fac,param_list,mass_species,basedir+Run_name)
        top.mainloop()
        plt.rcParams.update({'font.size': 16})


if save_vel_KE:
    print("Now calculating velocities and KE for all data...")
    par=np.loadtxt(basedir+Run_name+'final_parameters.txt')
    pos_offset_1[0],pos_offset_1[1],pos_offset_2[0],pos_offset_2[1],t0=par

    for i in range(1,chunks+1):
        b=i*fac
        if i==chunks:
            ion1xgc=np.asarray(ion1xg[a:])
            ion1ygc=np.asarray(ion1yg[a:])
            ion1tgc=np.asarray(ion1tg[a:])
            ion2xgc=np.asarray(ion2xg[a:])
            ion2ygc=np.asarray(ion2yg[a:])
            ion2tgc=np.asarray(ion2tg[a:])

        else:
            ion1xgc=np.asarray(ion1xg[a:b])
            ion1ygc=np.asarray(ion1yg[a:b])
            ion1tgc=np.asarray(ion1tg[a:b])
            ion2xgc=np.asarray(ion2xg[a:b])
            ion2ygc=np.asarray(ion2yg[a:b])
            ion2tgc=np.asarray(ion2tg[a:b])


        ion1xgc=np.add(ion1xgc,-pos_offset_1[0])
        ion2xgc=np.add(ion2xgc,-pos_offset_2[0])

        ion1ygc=np.add(ion1ygc,-pos_offset_1[1])
        ion2ygc=np.add(ion2ygc,-pos_offset_2[1])




        vx1 = (1.8399*ion1xgc)-0.0707
        vy1 = (1.8399*ion1ygc)-0.0707
        vz1 = (-5.5784*(ion1tgc-t0))+2689.3848

        vx2 = (0.368*ion2xgc)-0.0141
        vy2 = (0.368*ion2ygc)-0.0141
        vz2 = (-0.2231*(ion2tgc-t0))+537.8464





#
        a=b+1
        if i==1:
            with h5py.File(basedir+Run_name+'_vel_2body_v2.h5','w') as f:
                f.create_dataset('chunks', data = chunks, maxshape=(None), dtype='int16')
                f.create_dataset('masses', data = mass_species, maxshape=(None), dtype='f4')
                f.create_dataset('v1x/%d'%i, data = vx1, maxshape=(None), dtype='f4')
                f.create_dataset('v1y/%d'%i, data = vy1, maxshape=(None), dtype='f4')
                f.create_dataset('v1z/%d'%i, data = vz1, maxshape=(None), dtype='f4')
                f.create_dataset('v2x/%d'%i, data = vx2, maxshape=(None), dtype='f4')
                f.create_dataset('v2y/%d'%i, data = vy2, maxshape=(None), dtype='f4')
                f.create_dataset('v2z/%d'%i, data = vz2, maxshape=(None), dtype='f4')

                f.close()
        else:
            with h5py.File(basedir+Run_name+'_vel_2body_v2.h5','a') as f:
                f.create_dataset('v1x/%d'%i, data = vx1, maxshape=(None), dtype='f4')
                f.create_dataset('v1y/%d'%i, data = vy1, maxshape=(None), dtype='f4')
                f.create_dataset('v1z/%d'%i, data = vz1, maxshape=(None), dtype='f4')
                f.create_dataset('v2x/%d'%i, data = vx2, maxshape=(None), dtype='f4')
                f.create_dataset('v2y/%d'%i, data = vy2, maxshape=(None), dtype='f4')
                f.create_dataset('v2z/%d'%i, data = vz2, maxshape=(None), dtype='f4')



        print('chunks %d of %d'%(i,chunks))
    t3=time()
    print("time for calculating velocities = %f sec"%(t3-t2))
