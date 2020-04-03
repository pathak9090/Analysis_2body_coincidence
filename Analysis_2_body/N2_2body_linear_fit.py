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
sys.path.insert(0,'C:/Users/patha/Documents/ALS_Nov_2019_Analysis/Codes/')
from useful_definitions_2body import fhist1d,fhist2d,twoConditions
from optimise_parameters_2body import optimise_parameters
plt.rcParams.update({'font.size': 16})


t1=time()

Run_name='N2_160eV'
channel='N_N'
basedir1='C:/Users/patha/Documents/ALS_Nov_2019_Analysis/'+Run_name+'/'+channel+'/'
basedir=basedir1+channel
processed_dir=basedir1+'Processed/'
# path = basedir1+'Raw/'+Run_name+'_DAn.00000.bin'
path = basedir1+'Raw/air_160eV_2016_1st_April_new_par_DAn.bin'
# path = basedir+'/Raw/'+Run_name+'_DAn.bin'
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

mass_species = [14.0067,14.0067]
label_species = ['N$^+$','N$^+$']

# ionGate1 = [1690, 1800]
# ionGate2 = [1760, 1840]
ionGate1 = [2460,2600]
ionGate2 = [2520,2700]
# ionGate1 = [3570,3680]
# ionGate2 = [1,8000]
read_bin_file=True
# Fit to get t0
testmq=[31.9988,28.014,14.0067,15.999,39.948,1.008]
tof=[3857.75,3608.5,2547,2723.25,4313.5,672.75]


# tof = [3858,3608,2546,2723,673]
# testmq = [31.9988,28.014,14,15.999,1]
slp, t0 =  np.polyfit(np.sqrt(testmq), np.array(tof), 1)

tof_frag=np.add(slp*(np.sqrt(mass_species)),t0)

mass_species = np.array(mass_species) * 1.66E-27

tof_x_range=[1000,5000,1,-50,50,0.25]

if read_bin_file:
    LineNum=0
    ion1x=[]
    ion1y=[]
    ion2x=[]
    ion2y=[]
    ion1t=[]
    ion2t=[]
    
    event_counter=0
    all_evt_counter=0
    
    hits=[]
    while True:
    #        if event_counter > 5000:
    #            break
        #checks if end of file is reached
        if LineNum >= TotalNumLines:
            break
    
        #reading the hits in terms of x, y, tof
        numberOfHits = DataFile[LineNum]
        hits.append(numberOfHits)
        
        
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
    
    with h5py.File(processed_dir+Run_name+'_raw_xyt.h5','w') as f:
        f.create_dataset('ion1x', data = ion1x, maxshape=(None), dtype='f4')
        f.create_dataset('ion2x', data = ion2x, maxshape=(None), dtype='f4')
        f.create_dataset('ion1y', data = ion1y, maxshape=(None), dtype='f4')
        f.create_dataset('ion2y', data = ion2y, maxshape=(None), dtype='f4')
        f.create_dataset('ion1t', data = ion1t, maxshape=(None), dtype='f4')
        f.create_dataset('ion2t', data = ion2t, maxshape=(None), dtype='f4')
    f.close()
    
    plt.figure()
    x,y=fhist1d(hits,-0.5,10.5,1)
    plt.bar(x,y)
    plt.xlim(left=0)
    plt.ylabel("No. of Events")
    plt.xlabel("No. of Hits")
    plt.savefig(basedir+'Hits_vs_Events.png',bbox_inches='tight')
    
    print("Total no of events = %d"%np.size(hits))
    print("Total no of lines in binary file = %d"%LineNum)
    print("Calculated lines from hits array = %d"%((np.sum(hits)*3)+(np.size(hits))))
    
t2=time()
if not read_bin_file:
    hf1=h5py.File(processed_dir+Run_name+'_raw_xyt.h5','r')
    ion1x=np.array(hf1['ion1x'])
    ion2x=np.array(hf1['ion2x'])
    ion1y=np.array(hf1['ion1y'])
    ion2y=np.array(hf1['ion2y'])
    ion1t=np.array(hf1['ion1t'])
    ion2t=np.array(hf1['ion2t'])
    hf1.close()
    
    
#### Dibromoethane swapping hit2 and hit3
swap_hit1_hit2=False
if swap_hit1_hit2:
    # ion1x_t=np.concatenate((ion1x,ion1x),axis=0)
    ion1x_t=np.concatenate((ion1x,ion2x),axis=0)
    ion2x_t=np.concatenate((ion2x,ion1x),axis=0)
    ion1y_t=np.concatenate((ion1y,ion2y),axis=0)
    ion2y_t=np.concatenate((ion2y,ion1y),axis=0)
    ion1t_t=np.concatenate((ion1t,ion2t),axis=0)
    ion2t_t=np.concatenate((ion2t,ion1t),axis=0)
    
    
    ion1x=ion1x_t
    ion2x=ion2x_t

    ion1y=ion1y_t
    ion2y=ion2y_t

    ion1t=ion1t_t
    ion2t=ion2t_t


#%%
cmap='viridis'
## Histograms before selecting pipico channel
plot_det_images=True      # plots detector image for first, second hit
plot_all_pipico=True     # plots gated, non-gated pipico
plot_tof=True

# Gates for subplot 1 (tof1 vs tof2+tof3)
sp1_xmin=ionGate1[0]
sp1_xmax=ionGate1[1]
sp1_binsize=1
sp1_ymin=ionGate2[0]
sp1_ymax=ionGate2[1]


# Gates for subplot 2 (tof2+tof3-tof1 vs tof2+tof3+tof1)
sp2_xmin=sp1_ymin-sp1_xmax
# sp2_xmin=0
sp2_xmax=sp1_ymax-sp1_xmin
# sp2_xmax=90
sp2_binsize=1
sp2_ymin=sp1_xmin+sp1_ymin
sp2_ymax=sp1_ymax+sp1_xmax

pipico_range=[sp1_xmin,sp1_xmax,sp1_binsize,sp1_ymin,sp1_ymax,sp1_binsize]
pipicorot_range=[sp2_xmin,sp2_xmax,sp2_binsize,sp2_ymin,sp2_ymax,sp2_binsize]

## Time
if plot_all_pipico:
    ion1t=np.asarray(ion1t)
    ion2t=np.asarray(ion2t)

    fig = plt.figure(figsize=(18,4))
    ax = fig.add_subplot(121)
    x,y,z=fhist2d(ion1t,ion2t,sp1_xmin,sp1_xmax,sp1_binsize,sp1_ymin,sp1_ymax,sp1_binsize)
    f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap)#,norm=LogNorm())
    fig.colorbar(f1)
    ax.set_xlabel('TOF1 (ns)')
    ax.set_ylabel('TOF2 (ns)')
    ax.set_title('pipico')

    t12minus=np.add(ion2t,-ion1t)
    t12plus=np.add(ion2t,ion1t)

    ax = fig.add_subplot(122)
    x,y,z=fhist2d(t12minus,t12plus,sp2_xmin,sp2_xmax,sp2_binsize,sp2_ymin,sp2_ymax,sp2_binsize)
    f2=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap)#,norm=LogNorm())
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
    f1=ax.pcolormesh(x,y,z,cmap=cmap)
    fig.colorbar(f1)
#    ax.set_title(label_species[0])
    ax.set_aspect('equal')

    x,y,z=fhist2d(ion2x,ion2y,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(122)
    f2=ax.pcolormesh(x,y,z,cmap=cmap)
    fig.colorbar(f2)
#    ax.set_title(label_species[1])
    ax.set_aspect('equal')


    plt.savefig(basedir+'det_image_non_gated_v1.png',bbox_inches='tight')



if plot_tof:
    
    fig,ax = plt.subplots(figsize=(18,4))
    x1,y1=fhist1d(ion1t,tof_x_range[0],tof_x_range[1],tof_x_range[2])
    x2,y2=fhist1d(ion2t,tof_x_range[0],tof_x_range[1],tof_x_range[2])
    ax.plot(x1,y1)
    ax.plot(x2,y2)
    ax.axvline(x=tof_frag[0],color='g',linestyle='--')
    ax.axvline(x=tof_frag[1],color='g',linestyle='--')
    val=np.max([y1,y2])
    ax.text(tof_frag[0]+0.01*tof_frag[0],val-(0.05*val),label_species[0])
    ax.text(tof_frag[1]+0.01*tof_frag[1],val-(0.05*val),label_species[1])
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('TOF (ns)')
    plt.savefig(basedir+'tof_all.png',bbox_inches='tight')
    
plot_tof_x=True
if plot_tof_x:
    all_x=np.concatenate((ion1x,ion2x),axis=0)
    all_tof=np.concatenate((ion1t,ion2t),axis=0)
    fig, ax = plt.subplots(figsize=(18,4))

    x,y,z=fhist2d(all_tof,all_x,tof_x_range[0],tof_x_range[1],tof_x_range[2],tof_x_range[3],tof_x_range[4],tof_x_range[5])
    f1=ax.pcolormesh(x,y,np.transpose(z),norm=LogNorm(),cmap=cmap)
    fig.colorbar(f1)
    ax.axvline(x=tof_frag[0],color='g',linestyle='--')
    ax.axvline(x=tof_frag[1],color='g',linestyle='--')
    ax.text(tof_frag[0]+0.01*tof_frag[0],40,label_species[0])
    ax.text(tof_frag[1]+0.01*tof_frag[1],40,label_species[1])
    ax.set_title('X-TOF (Channel gated)')
    ax.set_xlabel('TOF (ns)')
    ax.set_ylabel('X (mm)')
    plt.savefig(basedir+'tof_x_v1.png',bbox_inches='tight')

#%%
# Selection of tripico channel


# t12rotGate = [51,116,5091,5099]   # No rec

# t12rotGate = [0,116,5091,5102]   # 
t12rotGate = [0,116,5081,5106]   # 
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
    f1=ax.pcolormesh(x,y,z,cmap=cmap)#,norm=LogNorm()
    fig.colorbar(f1)
    ax.set_title(label_species[0]+' (gated)')
    ax.set_aspect('equal')

    x,y,z=fhist2d(ion2xg,ion2yg,-50,50,0.5,-50,50,0.5)
    ax = fig.add_subplot(122)
    f2=ax.pcolormesh(x,y,z,cmap=cmap)
    fig.colorbar(f2)
    ax.set_title(label_species[1]+' (gated)')
    ax.set_aspect('equal')
    plt.savefig(basedir+'det_image_gated.png',bbox_inches='tight')

if plot_tof_x_gated:
    all_x=np.concatenate((ion1yg,ion2yg),axis=0)
    all_tof=np.concatenate((ion1tg,ion2tg),axis=0)
    fig, ax = plt.subplots(figsize=(18,4))

    x,y,z=fhist2d(all_tof,all_x,tof_x_range[0],tof_x_range[1],tof_x_range[2],tof_x_range[3],tof_x_range[4],tof_x_range[5])
    f1=ax.pcolormesh(x,y,np.transpose(z),norm=LogNorm(),cmap=cmap)
    fig.colorbar(f1)
    ax.axvline(x=tof_frag[0],color='g',linestyle='--')
    ax.axvline(x=tof_frag[1],color='g',linestyle='--')
    ax.text(tof_frag[0]+0.01*tof_frag[0],40,label_species[0])
    ax.text(tof_frag[1]+0.01*tof_frag[1],40,label_species[1])
    ax.set_title('X-TOF (Channel gated)')
    ax.set_xlabel('TOF (ns)')
    ax.set_ylabel('X (mm)')
    plt.savefig(basedir+'tof_x_gated_v1.png',bbox_inches='tight')
    
if plot_tof:
    fig,ax = plt.subplots(figsize=(18,4))
    x1,y1=fhist1d(ion1tg,tof_x_range[0],tof_x_range[1],tof_x_range[2])
    x2,y2=fhist1d(ion2tg,tof_x_range[0],tof_x_range[1],tof_x_range[2])
    ax.plot(x1,y1)
    ax.plot(x2,y2)
    ax.axvline(x=tof_frag[0],color='g',linestyle='--')
    ax.axvline(x=tof_frag[1],color='g',linestyle='--')
    val=np.max([y1,y2])
    ax.text(tof_frag[0]+0.01*tof_frag[0],val-(0.05*val),label_species[0])
    ax.text(tof_frag[1]+0.01*tof_frag[1],val-(0.05*val),label_species[1])
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('TOF (ns)')
    plt.savefig(basedir+'tof_gated.png',bbox_inches='tight')
    
    



#%% Loading functions for calculating velocities and KE

prange=[-500,500,1]

optimise_parameter=True# Set time offset and position offset to make momentum sums peak at zero








if optimise_parameter:
    save_vel_KE = False
else:
    save_vel_KE = True
            # Save velocities and KE for all data in h5 file (Do it after optimizing parameters)

#fac=0.2
#    #if tweaking in t0 is needed
## =============================================================================
#t0 = t0 - fac*t0

# position offsets (should be adjusted based on momentum sums later!)
pos_offset_1 = [0.1,0.0]
pos_offset_2 = [2.1,0.3]


load_old_param=True
if load_old_param:
    par=np.loadtxt(processed_dir+Run_name+'_'+channel+'_final_parameters.txt')
    pos_offset_1[0],pos_offset_1[1],pos_offset_2[0],pos_offset_2[1],t0=par

param_list=[pos_offset_1[0],pos_offset_1[1],pos_offset_2[0],pos_offset_2[1],t0]

# "vel_calc_fac is the array of all factors we need to calculate velocities. These are obtained from simion. The order is mx1,cx1,mz1,cz1,mx2,cx2,mz2,cz2 where mx1 means slope of linear fit between x and vx for particle 1 and cx is the intercept. y is not included as x and y are identical due to cylindrical symmetry "
# vel_calc_fac=[1.8399,-0.0707,-5.5784,2689.3848,0.368,-0.0141,-0.2231,537.8464]
# mx1,cx1,mz1,cz1,mx2,cx2,mz2,cz2=vel_calc_fac
# # Divide data in chunks if using Rbf else keep it 1

#### Loading velocity calculation factors
vel_calc_fac=np.loadtxt(processed_dir+'vel_calc_fac_'+channel+'.txt')

#vel_calc_fac=[0.5154,-0.0221,-0.3883,694.5532,0.3711,-0.0159,-0.2013,500.1138,0.3409,-0.0146,-0.1699,459.418]
mx1,cx1,mz1,cz1,mx2,cx2,mz2,cz2=vel_calc_fac


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
    
    
    if plot_det_images_gated:
        x1_c=pos_offset_1[0]
        x2_c=pos_offset_2[0]
        y1_c=pos_offset_1[1]
        y2_c=pos_offset_2[1]
        xc = np.linspace(-50, 50, 100)
        yc = np.linspace(-50, 50, 100)
        X, Y = np.meshgrid(xc,yc)
        R=700
        cmap='viridis'
        
        fig = plt.figure(figsize=(18,4))
        x,y,z=fhist2d(ion1xg,ion1yg,-50,50,0.5,-50,50,0.5)
        ax = fig.add_subplot(121)
        f1=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
        fig.colorbar(f1)
        ax.set_title(label_species[0]+' (gated)')
        ax.set_aspect('equal')
        F1 = (X-x1_c)**2 + (Y-y1_c)**2 - R
        ax.contour(X, Y, F1,[0], colors='r', linewidths=3)
    
        x,y,z=fhist2d(ion2xg,ion2yg,-50,50,0.5,-50,50,0.5)
        ax = fig.add_subplot(122)
        f2=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
        fig.colorbar(f2)
        ax.set_title(label_species[1]+' (gated)')
        ax.set_aspect('equal')
        F1 = (X-x2_c)**2 + (Y-y2_c)**2 - R
        ax.contour(X, Y, F1,[0], colors='r', linewidths=3)
        plt.savefig(basedir+'det_image_gated_circled.png',bbox_inches='tight')
        
        plt.rcParams.update({'font.size': 12})
        x,y,z=fhist2d(ion1xg,ion1yg,-50,50,0.5,-50,50,0.5)
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(221)
        ax.plot(x,np.sum(z,axis=1))
        ax.vlines(x1_c,0,np.max(np.sum(z,axis=1)),linestyle='dashed')
        ax.set_xticklabels([])
        ax.set_xticks([])
        plt.subplots_adjust(wspace=0.05,hspace=0.05)
        ax.tick_params(labelbottom=True, labeltop= False, labelleft=True, labelright=False,bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()
        ax.tick_params(axis='x', which='minor', direction='inout',length=4,width=1,top=True)
        ax.tick_params(axis='y', which='minor', direction='inout',length=4,width=1,right=True)
        # ax.set_aspect(0.125)
        plt.title('Detector image - %s'%label_species[0])
        ax = fig.add_subplot(224)
       
        ax.plot(np.sum(z,axis=0),y)
        ax.hlines(y1_c,0,np.max(np.sum(z,axis=0)),linestyle='dashed')
        ax.tick_params(labelbottom=True, labeltop= False, labelleft=True, labelright=False,bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()
        ax.tick_params(axis='x', which='minor', direction='inout',length=4,width=1,top=True)
        ax.tick_params(axis='y', which='minor', direction='inout',length=4,width=1,right=True)
        # ax.set_aspect(1)
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax = fig.add_subplot(223)
        ax.set_aspect('equal')
        f1=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
        F1 = (X-x1_c)**2 + (Y-y1_c)**2 - R
        ax.contour(X, Y, F1,[0], colors='r', linewidths=2,alpha=0.7)
        ax.hlines(y1_c,-50,50,linestyle='dashed',color='red',linewidth=1)
        ax.vlines(x1_c,-50,50,linestyle='dashed',color='red',linewidth=1)
        ax.tick_params(direction='inout',axis='both', length=6, width=2, colors='black',grid_color='b', grid_alpha=0.5)
        ax.tick_params(labelbottom=True, labeltop= False, labelleft=True, labelright=False,bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()
        ax.tick_params(axis='x', which='minor', direction='inout',length=4,width=1,top=True)
        ax.tick_params(axis='y', which='minor', direction='inout',length=4,width=1,right=True)
        fig.tight_layout()
        plt.savefig(basedir+'det_image_center_ion1.png',bbox_inches='tight')

        
        x,y,z=fhist2d(ion2xg,ion2yg,-50,50,0.5,-50,50,0.5)
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(221)
        ax.plot(x,np.sum(z,axis=1))
        ax.vlines(x2_c,0,np.max(np.sum(z,axis=1)),linestyle='dashed')
        ax.set_xticklabels([])
        ax.set_xticks([])
        plt.subplots_adjust(wspace=0.05,hspace=0.05)
        ax.tick_params(labelbottom=True, labeltop= False, labelleft=True, labelright=False,bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()
        ax.tick_params(axis='x', which='minor', direction='inout',length=4,width=1,top=True)
        ax.tick_params(axis='y', which='minor', direction='inout',length=4,width=1,right=True)
        # ax.set_aspect(0.125)
        plt.title('Detector image - %s'%label_species[1])
        ax = fig.add_subplot(224)
       
        ax.plot(np.sum(z,axis=0),y)
        ax.hlines(y2_c,0,np.max(np.sum(z,axis=0)),linestyle='dashed')
        ax.tick_params(labelbottom=True, labeltop= False, labelleft=True, labelright=False,bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()
        ax.tick_params(axis='x', which='minor', direction='inout',length=4,width=1,top=True)
        ax.tick_params(axis='y', which='minor', direction='inout',length=4,width=1,right=True)
        # ax.set_aspect(1)
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax = fig.add_subplot(223)
        ax.set_aspect('equal')
        f1=ax.pcolormesh(x,y,z,norm=LogNorm(),cmap=cmap)
        F1 = (X-x2_c)**2 + (Y-y2_c)**2 - R
        ax.contour(X, Y, F1,[0], colors='r', linewidths=2,alpha=0.7)
        ax.hlines(y2_c,-50,50,linestyle='dashed',color='red',linewidth=1)
        ax.vlines(x2_c,-50,50,linestyle='dashed',color='red',linewidth=1)
        ax.tick_params(direction='inout',axis='both', length=6, width=2, colors='black',grid_color='b', grid_alpha=0.5)
        ax.tick_params(labelbottom=True, labeltop= False, labelleft=True, labelright=False,bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()
        ax.tick_params(axis='x', which='minor', direction='inout',length=4,width=1,top=True)
        ax.tick_params(axis='y', which='minor', direction='inout',length=4,width=1,right=True)
        plt.savefig(basedir+'det_image_center_ion2.png',bbox_inches='tight')
        plt.rcParams.update({'font.size': 16})


        ion1xgc=np.add(ion1xgc,-pos_offset_1[0])
        ion2xgc=np.add(ion2xgc,-pos_offset_2[0])

        ion1ygc=np.add(ion1ygc,-pos_offset_1[1])
        ion2ygc=np.add(ion2ygc,-pos_offset_2[1])
        
        
        #%% 
        "Adding correction to TOF based on polyfit of TOF vs X"
        correct_second_order=False
        if correct_second_order:
            r1g=np.sqrt(np.add(ion1xgc**2,ion1ygc**2))
            r2g=np.sqrt(np.add(ion2xgc**2,ion2ygc**2))
            def polynom(x,A,B,C,D):
                return A+(B*x)+(C*x**2)+(D*x**3)
            coefs=np.loadtxt(processed_dir+'N_parameter_X_tof_fit.txt')
            corr1=polynom(r1g,coefs[0],coefs[1],coefs[2],coefs[3])
            corr2=polynom(r2g,coefs[0],coefs[1],coefs[2],coefs[3])
            
            ion1tgc=np.add(ion1tgc,-corr1)
            ion2tgc=np.add(ion2tgc,-corr2)

        vx1 = (mx1*ion1xgc)+cx1
        vy1 = (mx1*ion1ygc)+cx1
        vz1 = (mz1*(ion1tgc-t0))+cz1

        vx2 = (mx2*ion2xgc)+cx2
        vy2 = (mx2*ion2ygc)+cx2
        vz2 = (mz2*(ion2tgc-t0))+cz2



#
        a=b+1
        if i==1:
            with h5py.File(processed_dir+Run_name+'_vel_2body_v2.h5','w') as f:
                f.create_dataset('chunks', data = chunks, maxshape=(None), dtype='int16')
                f.create_dataset('masses', data = mass_species, maxshape=(None), dtype='f4')
                f.create_dataset('pipico_range', data = pipico_range, maxshape=(None), dtype='f4')
                f.create_dataset('pipicorot_range', data = pipicorot_range, maxshape=(None), dtype='f4')

                f.create_dataset('tof_x_range', data = tof_x_range, maxshape=(None), dtype='f4')
                f.create_dataset('tof_frag', data = tof_frag, maxshape=(None), dtype='f4')
                f.create_dataset('v1x/%d'%i, data = vx1, maxshape=(None), dtype='f4')
                f.create_dataset('v1y/%d'%i, data = vy1, maxshape=(None), dtype='f4')
                f.create_dataset('v1z/%d'%i, data = vz1, maxshape=(None), dtype='f4')
                f.create_dataset('v2x/%d'%i, data = vx2, maxshape=(None), dtype='f4')
                f.create_dataset('v2y/%d'%i, data = vy2, maxshape=(None), dtype='f4')
                f.create_dataset('v2z/%d'%i, data = vz2, maxshape=(None), dtype='f4')

                f.create_dataset('x1/%d'%i, data = ion1xgc, maxshape=(None), dtype='f4')
                f.create_dataset('x2/%d'%i, data = ion2xgc, maxshape=(None), dtype='f4')

                f.create_dataset('y1/%d'%i, data = ion1ygc, maxshape=(None), dtype='f4')
                f.create_dataset('y2/%d'%i, data = ion2ygc, maxshape=(None), dtype='f4')

                f.create_dataset('t1/%d'%i, data = ion1tgc, maxshape=(None), dtype='f4')
                f.create_dataset('t2/%d'%i, data = ion2tgc, maxshape=(None), dtype='f4')


                f.close()
        else:
            with h5py.File(processed_dir+Run_name+'_vel_2body_v2.h5','a') as f:
                f.create_dataset('v1x/%d'%i, data = vx1, maxshape=(None), dtype='f4')
                f.create_dataset('v1y/%d'%i, data = vy1, maxshape=(None), dtype='f4')
                f.create_dataset('v1z/%d'%i, data = vz1, maxshape=(None), dtype='f4')
                f.create_dataset('v2x/%d'%i, data = vx2, maxshape=(None), dtype='f4')
                f.create_dataset('v2y/%d'%i, data = vy2, maxshape=(None), dtype='f4')
                f.create_dataset('v2z/%d'%i, data = vz2, maxshape=(None), dtype='f4')
                f.create_dataset('x1/%d'%i, data = ion1xgc, maxshape=(None), dtype='f4')
                f.create_dataset('x2/%d'%i, data = ion2xgc, maxshape=(None), dtype='f4')

                f.create_dataset('y1/%d'%i, data = ion1ygc, maxshape=(None), dtype='f4')
                f.create_dataset('y2/%d'%i, data = ion2ygc, maxshape=(None), dtype='f4')

                f.create_dataset('t1/%d'%i, data = ion1tgc, maxshape=(None), dtype='f4')
                f.create_dataset('t2/%d'%i, data = ion2tgc, maxshape=(None), dtype='f4')



        print('chunks %d of %d'%(i,chunks))
    t3=time()
    print("time for calculating velocities = %f sec"%(t3-t2))
