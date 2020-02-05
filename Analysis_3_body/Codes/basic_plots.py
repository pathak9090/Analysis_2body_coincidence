# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 13:09:30 2019

@author: patha
"""
import numpy as np
import matplotlib.pyplot as plt
from time import time
import h5py
from matplotlib.colors import LogNorm
from useful_definitions import fhist1d,fhist2d,cosf1f2,NP,cosAB_keAB,DP
plt.rcParams.update({'font.size': 16})

def basic_plots(mom_KE_arr,pos_tof_arr,tripico23_range,tripico23rot_range,tripico12_range,tripico12rot_range,tof_x_range,tof_frag,prange,cosrange,KErange,PKER_range,keABcosABrange,DPrange,mass_species,label_species,basedir,Run_name):
    P1x,P1y,P1z,P2x,P2y,P2z,P3x,P3y,P3z,KE1,KE2,KE3,Px_all,Py_all,Pz_all,absPsum,KER=mom_KE_arr
    x1,x2,x3,y1,y2,y3,t1,t2,t3=pos_tof_arr
    KE=[KE1,KE2,KE3]
    Plot_momentum_panel=True
    Plot_angles=True
    Plot_NP=True
    Plot_KE=True
    Psum_KER=True
    Plot_KEcos=True
    plot_tripico=True
    plot_tof_x_gated=True
    plot_tof=True
    Plot_DP=True
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
        plt.savefig(basedir+'momentum_panel_GATED.png',bbox_inches='tight')
    
    
    if Plot_angles:
#        cosrange=[-1,1,0.05]
        cos12=cosf1f2(P1x,P1y,P1z,P2x,P2y,P2z)
        cos23=cosf1f2(P2x,P2y,P2z,P3x,P3y,P3z)
        cos13=cosf1f2(P1x,P1y,P1z,P3x,P3y,P3z)
        fig = plt.figure(figsize=(16,4))
        ax = fig.add_subplot(131)
        x,y=fhist1d(cos12,cosrange[0],cosrange[1],cosrange[2])
        ax.plot(x,y)
        ax.set_title(r'Cos$\theta$(%s,%s)'%(label_species[0],label_species[1]))
        ax.set_xlabel(r'Cos$\theta$')
        ax.set_ylim(bottom=0)
        ax = fig.add_subplot(132)
        x,y=fhist1d(cos23,cosrange[0],cosrange[1],cosrange[2])
        ax.plot(x,y)
        ax.set_title(r'Cos$\theta$(%s,%s)'%(label_species[1],label_species[2]))
        ax.set_xlabel(r'Cos$\theta$')
        ax.set_ylim(bottom=0)
        ax = fig.add_subplot(133)
        x,y=fhist1d(cos13,cosrange[0],cosrange[1],cosrange[2])
        ax.plot(x,y)
        ax.set_title(r'Cos$\theta$(%s,%s)'%(label_species[0],label_species[2]))
        ax.set_xlabel(r'Cos$\theta$')
        ax.set_ylim(bottom=0)
        plt.savefig(basedir+'angles_GATED.png',bbox_inches='tight')
        
        
        
    if Plot_NP:
        NP(P1x,P1y,P1z,label_species[0],P2x,P2y,P2z,label_species[1],P3x,P3y,P3z,label_species[2],1,-4,4,0.05,basedir+'_NP1_GATED.png',0,0,Run_name)
        NP(P1x,P1y,P1z,label_species[0],P2x,P2y,P2z,label_species[1],P3x,P3y,P3z,label_species[2],2,-4,4,0.05,basedir+'_NP2_GATED.png',0,0,Run_name)
        NP(P1x,P1y,P1z,label_species[0],P2x,P2y,P2z,label_species[1],P3x,P3y,P3z,label_species[2],3,-4,4,0.05,basedir+'_NP3_GATED.png',0,0,Run_name)
        
    if Plot_KE:
#        KErange=[0,40,0.5]
        x,y=fhist1d(KE1,KErange[0],KErange[1],KErange[2])
        fig, ax = plt.subplots()
        ax.plot(x,y,label=label_species[0])
        x,y=fhist1d(KE2,KErange[0],KErange[1],KErange[2])
        ax.plot(x,y,label=label_species[1])
        x,y=fhist1d(KE3,KErange[0],KErange[1],KErange[2])
        ax.plot(x,y,label=label_species[2])
        ax.legend(frameon=False)
        ax.set_xlabel('KE (eV)')
        ax.set_ylabel('Counts (arb units)')
        plt.savefig(basedir+'KE_all_frag_GATED.png',bbox_inches='tight')
        
    if Psum_KER:
        cmap='jet'
#        PKER_range=[0,200,1,0,80,0.5]
        x,y,z=fhist2d(absPsum,KER,PKER_range[0],PKER_range[1],PKER_range[2],PKER_range[3],PKER_range[4],PKER_range[5])
        fig, ax = plt.subplots()
        f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap)
        fig.colorbar(f1)
    #    ax.pcolormesh(x,y,z,label='KE1')
        ax.set_xlabel('Abs P (a.u.)')
        ax.set_ylabel('KER (eV)')
        plt.savefig(basedir+'AbsP_KER_GATED.png',bbox_inches='tight')
        
    mom_all=[P1x,P1y,P1z,P2x,P2y,P2z,P3x,P3y,P3z]
    if Plot_KEcos:
        plt.rcParams.update({'font.size': 16})
        KE12,cosH12H3=cosAB_keAB(mom_all,mass_species,3,keABcosABrange,label_species,cmap,basedir,log=True,)
        KE13,cosH13H2=cosAB_keAB(mom_all,mass_species,2,keABcosABrange,label_species,cmap,basedir,log=True,)
        KE23,cosH23H1=cosAB_keAB(mom_all,mass_species,1,keABcosABrange,label_species,cmap,basedir,log=True,)
        
        fig, ax = plt.subplots()
        x,y=fhist1d(KE12,keABcosABrange[0],keABcosABrange[1],keABcosABrange[2])
        ax.plot(x,y,label='KER (%s-%s)'%(label_species[0],label_species[1]))
        x,y=fhist1d(KE13,keABcosABrange[0],keABcosABrange[1],keABcosABrange[2])
        ax.plot(x,y,label='KER (%s-%s)'%(label_species[0],label_species[2]))
        ax.set_title('Native frames method')
        ax.legend(frameon=False)
        ax.set_xlabel('KER (eV)')
        ax.set_ylabel('Counts (arb units)')
        plt.savefig(basedir+'KER_NF_GATED.png',bbox_inches='tight')
    
    ion1t=np.asarray(t1)
    ion2t=np.asarray(t2)
    ion3t=np.asarray(t3)
    if plot_tripico:

        ion23t=np.add(ion2t,ion3t)
        ion12t=np.add(ion1t,ion2t)
        
        trotplus23=np.add(ion23t,ion1t)
        trotminus23=np.add(ion23t,-ion1t)

        trotplus12=np.add(ion12t,ion3t)
        trotminus12=np.add(ion12t,-ion3t)
    
        fig = plt.figure(figsize=(18,4))
        ax = fig.add_subplot(121)
        x,y,z=fhist2d(ion1t,ion23t,tripico23_range[0],tripico23_range[1],tripico23_range[2],tripico23_range[3],tripico23_range[4],tripico23_range[5])   #" 2) Hardcoded -Needs to be changed"
        f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
        fig.colorbar(f1)
        ax.set_xlabel('TOF1 (ns)')
        ax.set_ylabel('TOF2+TOF3 (ns)')
        ax.set_title('Tripico23 (momentum gated)')
    
    
        ax = fig.add_subplot(122)
        x,y,z=fhist2d(trotminus23,trotplus23,tripico23rot_range[0],tripico23rot_range[1],tripico23rot_range[2],tripico23rot_range[3],tripico23rot_range[4],tripico23rot_range[5])   #" 3) Hardcoded -Needs to be changed"
        f2=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
        fig.colorbar(f2)
        ax.set_xlabel('(TOF2+TOF3)-TOF1 (ns)')
        ax.set_ylabel('(TOF2+TOF3)+TOF1 (ns)')
        ax.set_title('Tripico23 rotated (momentum gated)')
        plt.savefig(basedir+'tripico23_gated_momenta.png',bbox_inches='tight')
    
        fig = plt.figure(figsize=(18,4))
        ax = fig.add_subplot(121)
        x,y,z=fhist2d(ion3t,ion12t,tripico12_range[0],tripico12_range[1],tripico12_range[2],tripico12_range[3],tripico12_range[4],tripico12_range[5])   #" 2) Hardcoded -Needs to be changed"
        f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
        fig.colorbar(f1)
        ax.set_xlabel('TOF3 (ns)')
        ax.set_ylabel('TOF1+TOF2 (ns)')
        ax.set_title('Tripico12 (momentum gated)')
    
        ax = fig.add_subplot(122)
        x,y,z=fhist2d(trotminus12,trotplus12,tripico12rot_range[0],tripico12rot_range[1],tripico12rot_range[2],tripico12rot_range[3],tripico12rot_range[4],tripico12rot_range[5])
        f2=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
        fig.colorbar(f2)
        ax.set_xlabel('(TOF1+TOF2)-TOF3 (ns)')
        ax.set_ylabel('(TOF1+TOF2)+TOF3 (ns)')
        ax.set_title('Tripico12 rotated (momentum gated)')
        plt.savefig(basedir+'tripico12_gated_momenta.png',bbox_inches='tight')
        
        
    if plot_tof_x_gated:
        all_x=np.concatenate((y1,y2,y3),axis=0)
        all_tof=np.concatenate((ion1t,ion2t,ion3t),axis=0)
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
        ax.set_title('X-TOF (momentum gated)')
        ax.set_xlabel('TOF (ns)')
        ax.set_ylabel('X (mm)')
        plt.savefig(basedir+'tof_x_momentum_gated.png',bbox_inches='tight')
        
    if plot_tof:
        fig,ax = plt.subplots(figsize=(18,4))
        x1,y1=fhist1d(t1,tof_x_range[0],tof_x_range[1],tof_x_range[2])
        x2,y2=fhist1d(t2,tof_x_range[0],tof_x_range[1],tof_x_range[2])
        x3,y3=fhist1d(t3,tof_x_range[0],tof_x_range[1],tof_x_range[2])
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
        ax.set_title('TOF (channel + momentum gated)')
        ax.set_xlabel('TOF (ns)')
        plt.savefig(basedir+'tof_momentum_gated.png',bbox_inches='tight')
        
        

    if Plot_DP:
        cmap='jet'
        DP(KE,DPrange,label_species,cmap,basedir,log=True)