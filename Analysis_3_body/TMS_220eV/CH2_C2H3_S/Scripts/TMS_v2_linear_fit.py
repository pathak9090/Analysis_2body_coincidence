# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 22:19:34 2019

@author: Shashank
"""

import numpy as np
import matplotlib.pyplot as plt
from time import time
#from sklearn.externals import joblib
import h5py
from matplotlib.colors import LogNorm
import tkinter as tk

import sys
sys.path.insert(0,'C:/Users/patha/OneDrive/Documents/GitHub/Double_VMI_coincidence_analysis/Analysis_3_body/Codes/')
import useful_definitions,optimise_parameters
from useful_definitions import fhist1d,fhist2d,threeConditions
from optimise_parameters import optimise_parameters
plt.rcParams.update({'font.size': 16})


t1=time()
#Run_name='Propargyl_alcohol_140ev_60C'

Run_name='TMS_220eV'
channel='CH2_C2H3_S'
basedir1='C:/Users/patha/OneDrive/Documents/GitHub/Double_VMI_coincidence_analysis/Analysis_3_body/'+Run_name+'/'+channel+'/'
basedir=basedir1+channel
processed_dir=basedir1+'Processed/'
path = basedir1+'Raw/'+Run_name+'_DAn.00000.bin'
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

mass_species = [14,27,32]
label_species = ['CH$_2$$^+$','C$_2$H$_3$$^+$','S$^+$']
ionGate1 = [1400,1600]
ionGate2 = [2000,2300]
ionGate3 = [2350,2550]


read_bin_file=True

# Fit to get t0
tof = [191,2733.66,1739.05,2239.94,4296.46]
testmq = [1.008,40,18.01,28,92.064]
slp, t0 =  np.polyfit(np.sqrt(testmq), np.array(tof), 1)

tof_frag=np.add(slp*(np.sqrt(mass_species)),t0)

mass_species = np.array(mass_species) * 1.66E-27


if read_bin_file:
    LineNum=0
    ion1x=[]
    ion1y=[]
    ion2x=[]
    ion2y=[]
    ion3x=[]
    ion3y=[]
    ion1t=[]
    ion2t=[]
    ion3t=[]
    
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
        if np.size(ion)<3: continue
        LineNum = LineNum + numberOfHits * 3
        #dividing by factor since thats how the numbers were saved in terms of int32
        ion = ion.reshape((numberOfHits, 3 )).T/1000
        ionX = ion[0]
        ionY = ion[1]
        ionTOF = ion[2]
    
        checkCondition = threeConditions(ionTOF, ionGate1, ionGate2, ionGate3)
        all_evt_counter=all_evt_counter+1
        if checkCondition.sum() == 3:
    
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
            ion3x.append(ionX[2])
            ion3y.append(ionY[2])
            ion1t.append(ionTOF[0])
            ion2t.append(ionTOF[1])
            ion3t.append(ionTOF[2])
    
            event_counter=event_counter+1
            if event_counter%10000 == 0:
                print(event_counter)
    
    thbody_after_permutation=np.size(ion1x)
    print('Total events = %d'%all_evt_counter)
    print('Total 3 body events after permutation = %d'%thbody_after_permutation)
    t2=time()
    
    print("time for reading file and getting position and time arrays = %f sec"%(t2-t1))
    
    ion1x=np.asarray(ion1x)
    ion2x=np.asarray(ion2x)
    ion3x=np.asarray(ion3x)
    ion1y=np.asarray(ion1y)
    ion2y=np.asarray(ion2y)
    ion3y=np.asarray(ion3y)
    ion1t=np.asarray(ion1t)
    ion2t=np.asarray(ion2t)
    ion3t=np.asarray(ion3t)
    
    
    
    with h5py.File(processed_dir+Run_name+'_raw_xyt.h5','w') as f:
        f.create_dataset('ion1x', data = ion1x, maxshape=(None), dtype='f4')
        f.create_dataset('ion2x', data = ion2x, maxshape=(None), dtype='f4')
        f.create_dataset('ion3x', data = ion3x, maxshape=(None), dtype='f4')
        f.create_dataset('ion1y', data = ion1y, maxshape=(None), dtype='f4')
        f.create_dataset('ion2y', data = ion2y, maxshape=(None), dtype='f4')
        f.create_dataset('ion3y', data = ion3y, maxshape=(None), dtype='f4')
        f.create_dataset('ion1t', data = ion1t, maxshape=(None), dtype='f4')
        f.create_dataset('ion2t', data = ion2t, maxshape=(None), dtype='f4')
        f.create_dataset('ion3t', data = ion3t, maxshape=(None), dtype='f4')
    f.close()

if not read_bin_file:
    hf1=h5py.File(processed_dir+Run_name+'_raw_xyt.h5','r')
    ion1x=np.array(hf1['ion1x'])
    ion2x=np.array(hf1['ion2x'])
    ion3x=np.array(hf1['ion3x'])
    ion1y=np.array(hf1['ion1y'])
    ion2y=np.array(hf1['ion2y'])
    ion3y=np.array(hf1['ion3y'])
    ion1t=np.array(hf1['ion1t'])
    ion2t=np.array(hf1['ion2t'])
    ion3t=np.array(hf1['ion3t'])
    hf1.close()
    
ion23t=np.add(ion2t,ion3t)
trotplus23=np.add(ion23t,ion1t)
trotminus23=np.add(ion23t,-ion1t)

ion12t=np.add(ion1t,ion2t)
trotplus12=np.add(ion12t,ion3t)
trotminus12=np.add(ion12t,-ion3t)

#%%
cmap='jet'
## Histograms before selecting tripico channel
plot_det_images=True      # plots detector image for first, second, third hit
plot_all_tripico=True     # plots gated, non-gated tripico
plot_tof=True


tof_x_range=[0,3000,1,-50,50,0.5]
# Gates for subplot 1 (tof1 vs tof2+tof3)
sp1_xmin=ionGate1[0]-100
sp1_xmax=ionGate1[1]+100
sp1_binsize=1
sp1_ymin=ionGate2[0]+ionGate3[0]
sp1_ymax=ionGate2[1]+ionGate3[1]

tripico23_range=[sp1_xmin,sp1_xmax,sp1_binsize,sp1_ymin,sp1_ymax,sp1_binsize]

# Gates for subplot 2 (tof2+tof3-tof1 vs tof2+tof3+tof1)
sp2_xmin=sp1_ymin-sp1_xmin-200
sp2_xmax=sp1_ymax-sp1_xmax+400
sp2_binsize=1
sp2_ymin=sp1_ymin+sp1_xmin-200
sp2_ymax=sp1_ymax+sp1_xmax+400

tripico23rot_range=[sp2_xmin,sp2_xmax,sp2_binsize,sp2_ymin,sp2_ymax,sp2_binsize]
# Gates for subplot 1 (tof3 vs tof1+tof2)
sp3_xmin=ionGate3[0]-200
sp3_xmax=ionGate3[1]+200
sp3_binsize=1
sp3_ymin=ionGate1[0]+ionGate2[0]
sp3_ymax=ionGate1[1]+ionGate2[1]

tripico12_range=[sp3_xmin,sp3_xmax,sp3_binsize,sp3_ymin,sp3_ymax,sp3_binsize]
# Gates for subplot 2 (tof1+tof2-tof3 vs tof1+tof2+tof3)
sp4_xmin=sp3_ymin-sp3_xmin-400
sp4_xmax=sp3_ymax-sp3_xmin+400
sp4_binsize=1
sp4_ymin=sp3_ymin+sp3_xmin-200
sp4_ymax=sp3_ymax+sp3_xmin+400
tripico12rot_range=[sp4_xmin,sp4_xmax,sp4_binsize,sp4_ymin,sp4_ymax,sp4_binsize]
## Time
if plot_all_tripico:
    ion1t=np.asarray(ion1t)
    ion2t=np.asarray(ion2t)
    ion3t=np.asarray(ion3t)
    ion23t=np.add(ion2t,ion3t)
    ion12t=np.add(ion1t,ion2t)
    tof_frag12=tof_frag[0]+tof_frag[1]
    tof_frag23=tof_frag[1]+tof_frag[2]
    
    fig = plt.figure(figsize=(18,4))
    ax = fig.add_subplot(121)
    x,y,z=fhist2d(ion1t,ion23t,sp1_xmin,sp1_xmax,sp1_binsize,sp1_ymin,sp1_ymax,sp1_binsize)   #" 2) Hardcoded -Needs to be changed"
    f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    ax.scatter(tof_frag[0], tof_frag23,marker='x',color='red',s=95)
    fig.colorbar(f1)
    ax.set_xlabel('TOF1 (ns)')
    ax.set_ylabel('TOF2+TOF3 (ns)')
    ax.set_title('Tripico23')

    tof_frag12plus=tof_frag12+tof_frag[2]
    tof_frag12minus=tof_frag12-tof_frag[2]
    tof_frag23plus=tof_frag23+tof_frag[0]
    tof_frag23minus=tof_frag23-tof_frag[0]
    ax = fig.add_subplot(122)
    x,y,z=fhist2d(trotminus23,trotplus23,sp2_xmin,sp2_xmax,sp2_binsize,sp2_ymin,sp2_ymax,sp2_binsize)   #" 3) Hardcoded -Needs to be changed"
    f2=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    ax.scatter(tof_frag23minus,tof_frag23plus,marker='x',color='red',s=95)
    fig.colorbar(f2)
    ax.set_xlabel('(TOF2+TOF3)-TOF1 (ns)')
    ax.set_ylabel('(TOF2+TOF3)+TOF1 (ns)')
    ax.set_title('Tripico23 rotated')
    plt.savefig(basedir+'tripico23_raw.png',bbox_inches='tight')

    fig = plt.figure(figsize=(18,4))
    ax = fig.add_subplot(121)
    x,y,z=fhist2d(ion3t,ion12t,sp3_xmin,sp3_xmax,sp3_binsize,sp3_ymin,sp3_ymax,sp3_binsize)   #" 2) Hardcoded -Needs to be changed"
    f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    ax.scatter(tof_frag[2], tof_frag12,marker='x',color='red',s=95)
    fig.colorbar(f1)
    ax.set_xlabel('TOF3 (ns)')
    ax.set_ylabel('TOF1+TOF2 (ns)')
    ax.set_title('Tripico12')

    ax = fig.add_subplot(122)
    x,y,z=fhist2d(trotminus12,trotplus12,sp4_xmin,sp4_xmax,sp4_binsize,sp4_ymin,sp4_ymax,sp4_binsize)
    f2=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    ax.scatter(tof_frag12minus,tof_frag12plus,marker='x',color='red',s=95)
    fig.colorbar(f2)
    ax.set_xlabel('(TOF1+TOF2)-TOF3 (ns)')
    ax.set_ylabel('(TOF1+TOF2)+TOF3 (ns)')
    ax.set_title('Tripico12 rotated')
    plt.savefig(basedir+'tripico12_raw.png',bbox_inches='tight')
    ##Positions
if plot_det_images:
    fig = plt.figure(figsize=(18,4))

    x,y,z=fhist2d(ion1x,ion1y,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(131)
    f1=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f1)
    ax.set_title('All 1$^{st}$ hits')
    ax.set_aspect('equal')

    x,y,z=fhist2d(ion2x,ion2y,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(132)
    f2=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f2)
    ax.set_title('All 2$^{nd}$ hits')
    ax.set_aspect('equal')

    x,y,z=fhist2d(ion3x,ion3y,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(133)
    f3=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f3)
    ax.set_title('All 3$^{rd}$ hits')
    ax.set_aspect('equal')
    plt.savefig(basedir+'det_image_non_gated_v1.png',bbox_inches='tight')
    
if plot_tof:
    fig,ax = plt.subplots(figsize=(18,4))
    x1,y1=fhist1d(ion1t,tof_x_range[0],tof_x_range[1],tof_x_range[2])
    x2,y2=fhist1d(ion2t,tof_x_range[0],tof_x_range[1],tof_x_range[2])
    x3,y3=fhist1d(ion3t,tof_x_range[0],tof_x_range[1],tof_x_range[2])
    ax.plot(x1,y1)
    ax.plot(x2,y2)
    ax.plot(x3,y3)
    ax.axvline(x=tof_frag[0],color='g',linestyle='--')
    ax.axvline(x=tof_frag[1],color='g',linestyle='--')
    ax.axvline(x=tof_frag[2],color='g',linestyle='--')
    val=np.max([y1,y2,y3])
    ax.text(tof_frag[0]+0.01*tof_frag[0],val-(0.05*val),label_species[0])
    ax.text(tof_frag[1]+0.01*tof_frag[1],val-(0.05*val),label_species[1])
    ax.text(tof_frag[2]++0.01*tof_frag[1],val-(0.05*val),label_species[2])
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('TOF (ns)')
    plt.savefig(basedir+'tof_all.png',bbox_inches='tight')
#%%
# Selection of tripico channel

#t123rotGate = [1306,1349,6165,6183]
t123rotGate = [1210,1348,6103,6120]  # Specify gating condition on tripico here


# Specify which tripico plot you are gating on
Gate_on_tripico23=False
Gate_on_tripico12=True
if Gate_on_tripico12:
    ## In addition to tripico gates you can also add ion position gates as shown below:
    condition=((trotminus12>t123rotGate[0]) & (trotminus12<t123rotGate[1]) & (trotplus12>t123rotGate[2]) & (trotplus12<t123rotGate[3]) & (ion2t>2159) & (ion2t<2233) & (ion1t>1470) & (ion1t<1530))

if Gate_on_tripico23:
    condition=((trotminus23>t123rotGate[0]) & (trotminus23<t123rotGate[1]) & (trotplus23>t123rotGate[2]) & (trotplus23<t123rotGate[3]))

ion1tg=ion1t[condition]
ion2tg=ion2t[condition]
ion3tg=ion3t[condition]
ion1xg=ion1x[condition]
ion2xg=ion2x[condition]
ion3xg=ion3x[condition]
ion1yg=ion1y[condition]
ion2yg=ion2y[condition]
ion3yg=ion3y[condition]

ion23tg=ion23t[condition]
ion12tg=ion12t[condition]
trotminus23g=trotminus23[condition]
trotplus23g=trotplus23[condition]
trotminus12g=trotminus12[condition]
trotplus12g=trotplus12[condition]
#%%
## Histograms after selecting tripico channel

plot_det_images_gated=True
plot_gated_channel=True
plot_tof_x_gated=True

if plot_gated_channel:

    fig = plt.figure(figsize=(18,4))
    ax = fig.add_subplot(121)
    x,y,z=fhist2d(ion1tg,ion23tg,sp1_xmin,sp1_xmax,sp1_binsize,sp1_ymin,sp1_ymax,sp1_binsize)   #" 2) Hardcoded -Needs to be changed"
    f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    fig.colorbar(f1)
    ax.set_xlabel('TOF1 (ns)')
    ax.set_ylabel('TOF2+TOF3 (ns)')
    ax.set_title('Tripico23 gated (%s+%s+%s)'%(label_species[0],label_species[1],label_species[2]))

    ax = fig.add_subplot(122)
    x,y,z=fhist2d(trotminus23g,trotplus23g,sp2_xmin,sp2_xmax,sp2_binsize,sp2_ymin,sp2_ymax,sp2_binsize)
    f2=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    fig.colorbar(f2)
    ax.set_xlabel('(TOF2+TOF3)-TOF1 (ns)')
    ax.set_ylabel('(TOF2+TOF3)+TOF1 (ns)')
    ax.set_title('Tripico23 gated + rotated (%s+%s+%s)'%(label_species[0],label_species[1],label_species[2]))
    plt.savefig(basedir+'tripico23_gated.png',bbox_inches='tight')

    fig = plt.figure(figsize=(18,4))
    ax = fig.add_subplot(121)
    x,y,z=fhist2d(ion3tg,ion12tg,sp3_xmin,sp3_xmax,sp3_binsize,sp3_ymin,sp3_ymax,sp3_binsize)   #" 2) Hardcoded -Needs to be changed"
    f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    fig.colorbar(f1)
    ax.set_xlabel('TOF3 (ns)')
    ax.set_ylabel('TOF1+TOF2 (ns)')
    ax.set_title('Tripico12 gated (%s+%s+%s)'%(label_species[0],label_species[1],label_species[2]))

    ax = fig.add_subplot(122)
    x,y,z=fhist2d(trotminus12g,trotplus12g,sp4_xmin,sp4_xmax,sp4_binsize,sp4_ymin,sp4_ymax,sp4_binsize)   #" 3) Hardcoded -Needs to be changed"
    f2=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
    fig.colorbar(f2)
    ax.set_xlabel('(TOF1+TOF2)-TOF3 (ns)')
    ax.set_ylabel('(TOF1+TOF2)+TOF3 (ns)')
    ax.set_title('Tripico12 gated + rotated (%s+%s+%s)'%(label_species[0],label_species[1],label_species[2]))
    plt.savefig(basedir+'tripico12_gated.png',bbox_inches='tight')

if plot_tof_x_gated:
    all_x=np.concatenate((ion1yg,ion2yg,ion3yg),axis=0)
    all_tof=np.concatenate((ion1tg,ion2tg,ion3tg),axis=0)
    fig, ax = plt.subplots(figsize=(18,4))

    x,y,z=fhist2d(all_tof,all_x,tof_x_range[0],tof_x_range[1],tof_x_range[2],tof_x_range[3],tof_x_range[4],tof_x_range[5])
    f1=ax.pcolormesh(x,y,np.transpose(z),norm=LogNorm(),cmap=cmap)
    fig.colorbar(f1)
    ax.axvline(x=tof_frag[0],color='g',linestyle='--')
    ax.axvline(x=tof_frag[1],color='g',linestyle='--')
    ax.axvline(x=tof_frag[2],color='g',linestyle='--')
    ax.text(tof_frag[0]+0.01*tof_frag[0],40,label_species[0])
    ax.text(tof_frag[1]+0.01*tof_frag[1],40,label_species[1])
    ax.text(tof_frag[2]++0.01*tof_frag[1],40,label_species[2])
    ax.set_title('X-TOF')
    ax.set_xlabel('TOF (ns)')
    ax.set_ylabel('X (mm)')
    plt.savefig(basedir+'tof_x_gated_v1.png',bbox_inches='tight')



if plot_det_images_gated:
    fig = plt.figure(figsize=(18,4))

    x,y,z=fhist2d(ion1xg,ion1yg,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(131)
    f1=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f1)
    ax.set_title(label_species[0])
    ax.set_aspect('equal')

    x,y,z=fhist2d(ion2xg,ion2yg,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(132)
    f2=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f2)
    ax.set_title(label_species[1])
    ax.set_aspect('equal')

    x,y,z=fhist2d(ion3xg,ion3yg,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(133)
    f3=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
    fig.colorbar(f3)
    ax.set_title(label_species[2])
    ax.set_aspect('equal')
    plt.savefig(basedir+'det_image_gated_v1.png',bbox_inches='tight')

if plot_tof:
    fig,ax = plt.subplots(figsize=(18,4))
    x1,y1=fhist1d(ion1tg,tof_x_range[0],tof_x_range[1],tof_x_range[2])
    x2,y2=fhist1d(ion2tg,tof_x_range[0],tof_x_range[1],tof_x_range[2])
    x3,y3=fhist1d(ion3tg,tof_x_range[0],tof_x_range[1],tof_x_range[2])
    ax.plot(x1,y1)
    ax.plot(x2,y2)
    ax.plot(x3,y3)
    ax.axvline(x=tof_frag[0],color='g',linestyle='--')
    ax.axvline(x=tof_frag[1],color='g',linestyle='--')
    ax.axvline(x=tof_frag[2],color='g',linestyle='--')
    val=np.max([y1,y2,y3])
    ax.text(tof_frag[0]+0.01*tof_frag[0],val-(0.05*val),label_species[0])
    ax.text(tof_frag[1]+0.01*tof_frag[1],val-(0.05*val),label_species[1])
    ax.text(tof_frag[2]++0.01*tof_frag[1],val-(0.05*val),label_species[2])
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('TOF (ns)')
    plt.savefig(basedir+'tof_gated.png',bbox_inches='tight')
#%% Loading functions for calculating velocities and KE


prange=[-500,500,1]

optimise_parameter=True  # Set time offset and position offset to make momentum sums peak at zero
save_vel_KE = False       # Save velocities and KE for all data in h5 file (Do it after optimizing parameters)



# position offsets (should be adjusted based on momentum sums later!)
pos_offset_1 = [0.1,0.0]
pos_offset_2 = [2.1,0.3]
pos_offset_3 = [2.6, 0.8]

load_old_param=True
if load_old_param:
    par=np.loadtxt(processed_dir+Run_name+'_'+channel+'_final_parameters.txt')
    pos_offset_1[0],pos_offset_1[1],pos_offset_2[0],pos_offset_2[1],pos_offset_3[0],pos_offset_3[1],t0=par

param_list=[pos_offset_1[0],pos_offset_1[1],pos_offset_2[0],pos_offset_2[1],pos_offset_3[0],pos_offset_3[1],t0]



#### Loading velocity calculation factors
vel_calc_fac=np.loadtxt(processed_dir+'vel_calc_fac_'+channel+'.txt')

#vel_calc_fac=[0.5154,-0.0221,-0.3883,694.5532,0.3711,-0.0159,-0.2013,500.1138,0.3409,-0.0146,-0.1699,459.418]
mx1,cx1,mz1,cz1,mx2,cx2,mz2,cz2,mx3,cx3,mz3,cz3=vel_calc_fac

# Divide data in chunks if using Rbf (for multidimentional fitting) else keep it 1
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
        ion3xgc=np.asarray(ion3xg[a:fac])
        ion3ygc=np.asarray(ion3yg[a:fac])
        ion3tgc=np.asarray(ion3tg[a:fac])

        xyt_arr=[ion1xgc,ion1ygc,ion1tgc,ion2xgc,ion2ygc,ion2tgc,ion3xgc,ion3ygc,ion3tgc]
#       vel_calc_fac=[1.8399,-0.0707,-5.5784,2689.3848,0.368,-0.0141,-0.2231,537.8464,0.3359,-0.0129,-0.1859,490.9972]
#        mx1,cx1,mz1,cz1,mx2,cx2,mz2,cz2,mx3,cx3,mz3,cz3=vel_calc_fac
        plt.rcParams.update({'font.size': 12})

        "Calling GUI..."

        top = tk.Tk()
        top.geometry("2000x2000")
        optimise_parameters(top,xyt_arr,vel_calc_fac,param_list,mass_species,processed_dir+Run_name,channel)
        top.mainloop()
        plt.rcParams.update({'font.size': 16})


if save_vel_KE:
    print("Now calculating velocities and KE for all data...")
    par=np.loadtxt(processed_dir+Run_name+'_'+channel+'_final_parameters.txt')
    pos_offset_1[0],pos_offset_1[1],pos_offset_2[0],pos_offset_2[1],pos_offset_3[0],pos_offset_3[1],t0=par

    for i in range(1,chunks+1):
        b=i*fac
        if i==chunks:
            ion1xgc=np.asarray(ion1xg[a:])
            ion1ygc=np.asarray(ion1yg[a:])
            ion1tgc=np.asarray(ion1tg[a:])
            ion2xgc=np.asarray(ion2xg[a:])
            ion2ygc=np.asarray(ion2yg[a:])
            ion2tgc=np.asarray(ion2tg[a:])
            ion3xgc=np.asarray(ion3xg[a:])
            ion3ygc=np.asarray(ion3yg[a:])
            ion3tgc=np.asarray(ion3tg[a:])
        else:
            ion1xgc=np.asarray(ion1xg[a:b])
            ion1ygc=np.asarray(ion1yg[a:b])
            ion1tgc=np.asarray(ion1tg[a:b])
            ion2xgc=np.asarray(ion2xg[a:b])
            ion2ygc=np.asarray(ion2yg[a:b])
            ion2tgc=np.asarray(ion2tg[a:b])
            ion3xgc=np.asarray(ion3xg[a:b])
            ion3ygc=np.asarray(ion3yg[a:b])
            ion3tgc=np.asarray(ion3tg[a:b])

        ion1xgc=np.add(ion1xgc,-pos_offset_1[0])
        ion2xgc=np.add(ion2xgc,-pos_offset_2[0])
        ion3xgc=np.add(ion3xgc,-pos_offset_3[0])
        ion1ygc=np.add(ion1ygc,-pos_offset_1[1])
        ion2ygc=np.add(ion2ygc,-pos_offset_2[1])
        ion3ygc=np.add(ion3ygc,-pos_offset_3[1])


#        vx1 = (1.8399*ion1xgc)-0.0707
#        vy1 = (1.8399*ion1ygc)-0.0707
#        vz1 = (-5.5784*(ion1tgc-t0))+2689.3848
#
#        vx2 = (0.368*ion2xgc)-0.0141
#        vy2 = (0.368*ion2ygc)-0.0141
#        vz2 = (-0.2231*(ion2tgc-t0))+537.8464
#
#        vx3 = (0.3359*ion3xgc)-0.0129
#        vy3 = (0.3359*ion3ygc)-0.0129
#        vz3 = (-0.1859*(ion3tgc-t0))+490.9972
        ion1tgc_s=ion1tgc-t0
        ion2tgc_s=ion2tgc-t0
        ion3tgc_s=ion3tgc-t0
        vx1 = (mx1*ion1xgc)+cx1
        vy1 = (mx1*ion1ygc)+cx1
        vz1 = (mz1*(ion1tgc_s))+cz1

        vx2 = (mx2*ion2xgc)+cx2
        vy2 = (mx2*ion2ygc)+cx2
        vz2 = (mz2*(ion2tgc_s))+cz2

        vx3 = (mx3*ion3xgc)+cx3
        vy3 = (mx3*ion3ygc)+cx3
        vz3 = (mz3*(ion3tgc_s))+cz3




#
        a=b+1
        if i==1:
            with h5py.File(processed_dir+Run_name+'_vel_KE_3body_v2.h5','w') as f:
                f.create_dataset('chunks', data = chunks, maxshape=(None), dtype='int16')
                f.create_dataset('masses', data = mass_species, maxshape=(None), dtype='f4')
                f.create_dataset('tripico23_range', data = tripico23_range, maxshape=(None), dtype='f4')
                f.create_dataset('tripico23rot_range', data = tripico23rot_range, maxshape=(None), dtype='f4')
                f.create_dataset('tripico12_range', data = tripico12_range, maxshape=(None), dtype='f4')
                f.create_dataset('tripico12rot_range', data = tripico12rot_range, maxshape=(None), dtype='f4')
                f.create_dataset('tof_x_range', data = tof_x_range, maxshape=(None), dtype='f4')
                f.create_dataset('tof_frag', data = tof_frag, maxshape=(None), dtype='f4')
                f.create_dataset('v1x/%d'%i, data = vx1, maxshape=(None), dtype='f4')
                f.create_dataset('v1y/%d'%i, data = vy1, maxshape=(None), dtype='f4')
                f.create_dataset('v1z/%d'%i, data = vz1, maxshape=(None), dtype='f4')
                f.create_dataset('v2x/%d'%i, data = vx2, maxshape=(None), dtype='f4')
                f.create_dataset('v2y/%d'%i, data = vy2, maxshape=(None), dtype='f4')
                f.create_dataset('v2z/%d'%i, data = vz2, maxshape=(None), dtype='f4')
                f.create_dataset('v3x/%d'%i, data = vx3, maxshape=(None), dtype='f4')
                f.create_dataset('v3y/%d'%i, data = vy3, maxshape=(None), dtype='f4')
                f.create_dataset('v3z/%d'%i, data = vz3, maxshape=(None), dtype='f4')
                f.create_dataset('x1/%d'%i, data = ion1xgc, maxshape=(None), dtype='f4')
                f.create_dataset('x2/%d'%i, data = ion2xgc, maxshape=(None), dtype='f4')
                f.create_dataset('x3/%d'%i, data = ion3xgc, maxshape=(None), dtype='f4')
                f.create_dataset('y1/%d'%i, data = ion1ygc, maxshape=(None), dtype='f4')
                f.create_dataset('y2/%d'%i, data = ion2ygc, maxshape=(None), dtype='f4')
                f.create_dataset('y3/%d'%i, data = ion3ygc, maxshape=(None), dtype='f4')
                f.create_dataset('t1/%d'%i, data = ion1tgc, maxshape=(None), dtype='f4')
                f.create_dataset('t2/%d'%i, data = ion2tgc, maxshape=(None), dtype='f4')
                f.create_dataset('t3/%d'%i, data = ion3tgc, maxshape=(None), dtype='f4')
#                f.create_dataset('KE1/%d'%i, data = rKE1, maxshape=(None), dtype='f4')
#                f.create_dataset('KE2/%d'%i, data = rKE2, maxshape=(None), dtype='f4')
#                f.create_dataset('KE3/%d'%i, data = rKE3, maxshape=(None), dtype='f4')
                f.close()
        else:
            with h5py.File(processed_dir+Run_name+'_vel_KE_3body_v2.h5','a') as f:
                f.create_dataset('v1x/%d'%i, data = vx1, maxshape=(None), dtype='f4')
                f.create_dataset('v1y/%d'%i, data = vy1, maxshape=(None), dtype='f4')
                f.create_dataset('v1z/%d'%i, data = vz1, maxshape=(None), dtype='f4')
                f.create_dataset('v2x/%d'%i, data = vx2, maxshape=(None), dtype='f4')
                f.create_dataset('v2y/%d'%i, data = vy2, maxshape=(None), dtype='f4')
                f.create_dataset('v2z/%d'%i, data = vz2, maxshape=(None), dtype='f4')
                f.create_dataset('v3x/%d'%i, data = vx3, maxshape=(None), dtype='f4')
                f.create_dataset('v3y/%d'%i, data = vy3, maxshape=(None), dtype='f4')
                f.create_dataset('v3z/%d'%i, data = vz3, maxshape=(None), dtype='f4')
                f.create_dataset('x1/%d'%i, data = ion1xgc, maxshape=(None), dtype='f4')
                f.create_dataset('x2/%d'%i, data = ion2xgc, maxshape=(None), dtype='f4')
                f.create_dataset('x3/%d'%i, data = ion3xgc, maxshape=(None), dtype='f4')
                f.create_dataset('y1/%d'%i, data = ion1ygc, maxshape=(None), dtype='f4')
                f.create_dataset('y2/%d'%i, data = ion2ygc, maxshape=(None), dtype='f4')
                f.create_dataset('y3/%d'%i, data = ion3ygc, maxshape=(None), dtype='f4')
                f.create_dataset('t1/%d'%i, data = ion1tgc_s, maxshape=(None), dtype='f4')
                f.create_dataset('t2/%d'%i, data = ion2tgc_s, maxshape=(None), dtype='f4')
                f.create_dataset('t3/%d'%i, data = ion3tgc_s, maxshape=(None), dtype='f4')
#                f.create_dataset('KE1/%d'%i, data = rKE1, maxshape=(None), dtype='f4')
#                f.create_dataset('KE2/%d'%i, data = rKE2, maxshape=(None), dtype='f4')
#                f.create_dataset('KE3/%d'%i, data = rKE3, maxshape=(None), dtype='f4')

        print('chunks %d of %d'%(i,chunks))
    if read_bin_file:
        t3=time()
        print("time for calculating velocities = %f sec"%(t3-t2))
