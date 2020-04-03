# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 13:09:30 2019

@author: Shashank Pathak
"""
import numpy as np
import matplotlib.pyplot as plt
from time import time
import h5py
from matplotlib.colors import LogNorm
from useful_definitions import fhist1d,fhist2d,cosf1f2,NP,cosAB_keAB,DP
plt.rcParams.update({'font.size': 16})

def basic_plots(mom_KE_arr,pos_tof_arr,cosz,pipico_range,pipicorot_range,tof_x_range,tof_frag,prange,cosrange,KErange,PKER_range,KER_cosz_range,mass_species,label_species,basedir,Run_name):
    P1x,P1y,P1z,P2x,P2y,P2z,KE1,KE2,Px_all,Py_all,Pz_all,absPsum,KER=mom_KE_arr
    x1,x2,y1,y2,t1,t2=pos_tof_arr
    KE=[KE1,KE2]
    Plot_momentum_panel=True
    Plot_angles=True
    Plot_NP=True
    Plot_KE=True
    Psum_KER=True
    plot_pipico=True
    plot_tof_x_gated=True
    plot_tof=True
    Plot_DP=True
    KER_cosz=True
    plot_det_img=True
    if Plot_momentum_panel:
        plt.rcParams.update({'font.size': 13})
        # fig = plt.figure(figsize=(22,16))
        # ax = fig.add_subplot(331)
        # x,y=fhist1d(P1x,prange[0],prange[1],prange[2])
        # ax.plot(x,y,label='Px (%s)'%label_species[0])
        # ax.legend(frameon=False)
        # ax = fig.add_subplot(332)
        # x,y=fhist1d(P1y,prange[0],prange[1],prange[2])
        # ax.plot(x,y,label='Py (%s)'%label_species[0])
        # ax.legend(frameon=False)
        # ax = fig.add_subplot(333)
        # x,y=fhist1d(P1z,prange[0],prange[1],prange[2])
        # ax.plot(x,y,label='Pz (%s)'%label_species[0])
        # ax.legend(frameon=False)
        # ax = fig.add_subplot(334)
        # x,y=fhist1d(P2x,prange[0],prange[1],prange[2])
        # ax.plot(x,y,label='Px (%s)'%label_species[1])
        # ax.legend(frameon=False)
        # ax = fig.add_subplot(335)
        # x,y=fhist1d(P2y,prange[0],prange[1],prange[2])
        # ax.plot(x,y,label='Py (%s)'%label_species[1])
        # ax.legend(frameon=False)
        # ax = fig.add_subplot(336)
        # x,y=fhist1d(P2z,prange[0],prange[1],prange[2])
        # ax.plot(x,y,label='Pz (%s)'%label_species[1])
        # ax.legend(frameon=False)
        # ax = fig.add_subplot(3,3,7)
        # x,y=fhist1d(Px_all,prange[0],prange[1],prange[2])
        # ax.plot(x,y,label='Px Sum')
        # ax.legend(frameon=False)
        # ax = fig.add_subplot(3,3,8)
        # x,y=fhist1d(Py_all,prange[0],prange[1],prange[2])
        # ax.plot(x,y,label='Py Sum')
        # ax.legend(frameon=False)
        # ax = fig.add_subplot(3,3,9)
        # x,y=fhist1d(Pz_all,prange[0],prange[1],prange[2])
        # ax.plot(x,y,label='Pz Sum')
        # ax.legend(frameon=False)
        # plt.savefig(basedir+'momentum_panel_GATED.png',bbox_inches='tight')
        fig = plt.figure(figsize=(22,16))
        ax = fig.add_subplot(331)
        x,y=fhist1d(P1x,prange[0],prange[1],prange[2])
        ax.step(x,y,label='Px (%s)'%label_species[0])
        ax.legend(frameon=False)
        ax = fig.add_subplot(332)
        x,y=fhist1d(P1y,prange[0],prange[1],prange[2])
        ax.step(x,y,label='Py (%s)'%label_species[0])
        ax.legend(frameon=False)
        ax = fig.add_subplot(333)
        x,y=fhist1d(P1z,prange[0],prange[1],prange[2])
        ax.step(x,y,label='Pz (%s)'%label_species[0])
        ax.legend(frameon=False)
        ax = fig.add_subplot(334)
        x,y=fhist1d(P2x,prange[0],prange[1],prange[2])
        ax.step(x,y,label='Px (%s)'%label_species[1])
        ax.legend(frameon=False)
        ax = fig.add_subplot(335)
        x,y=fhist1d(P2y,prange[0],prange[1],prange[2])
        ax.step(x,y,label='Py (%s)'%label_species[1])
        ax.legend(frameon=False)
        ax = fig.add_subplot(336)
        x,y=fhist1d(P2z,prange[0],prange[1],prange[2])
        ax.step(x,y,label='Pz (%s)'%label_species[1])
        ax.legend(frameon=False)
        ax = fig.add_subplot(3,3,7)
        x,y=fhist1d(Px_all,-50,50,0.1)
        ax.step(x,y,label='Px Sum')
        ax.legend(frameon=False)
        ax = fig.add_subplot(3,3,8)
        x,y=fhist1d(Py_all,-50,50,0.1)
        ax.step(x,y,label='Py Sum')
        ax.legend(frameon=False)
        ax = fig.add_subplot(3,3,9)
        x,y=fhist1d(Pz_all,-50,50,0.1)
        ax.step(x,y,label='Pz Sum')
        ax.legend(frameon=False)
        plt.savefig(basedir+'momentum_panel_GATED.png',bbox_inches='tight')
    
    if Plot_angles:
#        cosrange=[-1,1,0.05]
        cos12=cosf1f2(P1x,P1y,P1z,P2x,P2y,P2z)
    
        fig = plt.figure(figsize=(16,4))
        ax = fig.add_subplot(111)
        x,y=fhist1d(cos12,cosrange[0],cosrange[1],cosrange[2])
        ax.plot(x,y)
        ax.set_title(r'Cos$\theta$(%s,%s)'%(label_species[0],label_species[1]))
        ax.set_xlabel(r'Cos$\theta$')
        ax.set_title('%s + %s'%(label_species[0],label_species[1]))
        plt.savefig(basedir+'angles_GATED.png',bbox_inches='tight')
        

        
    if Plot_KE:
#        KErange=[0,40,0.5]
        x,y=fhist1d(KE1,KErange[0],KErange[1],KErange[2])
        fig, ax = plt.subplots()
        # ax.plot(x,y,label=label_species[0])
        x,y=fhist1d(KE2,KErange[0],KErange[1],KErange[2])
        # ax.plot(x,y,label=label_species[1])
        x,y=fhist1d(np.add(KE1,KE2),KErange[0],KErange[1],KErange[2])
        ax.step(x,y,label='KER')
        KERgx=x
        KERgy=y
        ax.legend(frameon=False)
        ax.set_xlabel('KER (eV)')
        ax.set_ylabel('Counts (arb units)')
        ax.set_title('%s + %s'%(label_species[0],label_species[1]))
        plt.savefig(basedir+'KE_all_frag_GATED.png',bbox_inches='tight')
        
    if Psum_KER:
        cmap='jet'
        x,y,z=fhist2d(absPsum,KER,PKER_range[0],PKER_range[1],PKER_range[2],PKER_range[3],PKER_range[4],PKER_range[5])
        fig, ax = plt.subplots()
        f1=ax.pcolormesh(x,y,np.transpose(z),norm=LogNorm(),cmap=cmap)
        fig.colorbar(f1)
    #    ax.pcolormesh(x,y,z,label='KE1')
        ax.set_xlabel('Abs P (a.u.)')
        ax.set_ylabel('KER (eV)')
        ax.set_title('%s + %s'%(label_species[0],label_species[1]))
        plt.savefig(basedir+'AbsP_KER_GATED.png',bbox_inches='tight')
        
    mom_all=[P1x,P1y,P1z,P2x,P2y,P2z]

    
    ion1t=np.asarray(t1)
    ion2t=np.asarray(t2)

    if plot_pipico:

        t1=np.asarray(t1)
        t2=np.asarray(t2)
        t12minusg=np.add(t2,-t1)
        t12plusg=np.add(t2,t1)
    
        fig = plt.figure(figsize=(18,4))
        ax = fig.add_subplot(121)
        x,y,z=fhist2d(t1,t2,pipico_range[0],pipico_range[1],pipico_range[2],pipico_range[3],pipico_range[4],pipico_range[5])
        f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
        fig.colorbar(f1)
        ax.set_xlabel('TOF1 (ns)')
        ax.set_ylabel('TOF2 (ns)')
        ax.set_title('pipico (momentum gated)')
    
        ax = fig.add_subplot(122)
        x,y,z=fhist2d(t12minusg,t12plusg,pipicorot_range[0],pipicorot_range[1],pipicorot_range[2],pipicorot_range[3],pipicorot_range[4],pipicorot_range[5])
        f2=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap,norm=LogNorm())
        fig.colorbar(f2)
        ax.set_xlabel('(TOF2-TOF1 (ns)')
        ax.set_ylabel('TOF2+TOF1 (ns)')
        ax.set_title('pipico rotated (momentum gated)')
        plt.savefig(basedir+'pipico_rotated_momentum_gated.png',bbox_inches='tight')
        
        
    if plot_tof_x_gated:
        all_x=np.concatenate((y1,y2),axis=0)
        all_tof=np.concatenate((t1,t2),axis=0)
        fig, ax = plt.subplots(figsize=(18,4))
    
        x,y,z=fhist2d(all_tof,all_x,tof_x_range[0],tof_x_range[1],tof_x_range[2],tof_x_range[3],tof_x_range[4],tof_x_range[5])
        f1=ax.pcolormesh(x,y,np.transpose(z),norm=LogNorm(),cmap=cmap)
        fig.colorbar(f1)
        ax.axvline(x=tof_frag[0],color='g',linestyle='--')
        ax.axvline(x=tof_frag[1],color='g',linestyle='--')
        ax.text(tof_frag[0]+0.01*tof_frag[0],40,label_species[0])
        ax.text(tof_frag[1]+0.01*tof_frag[1],40,label_species[1])
        ax.set_title('X-TOF (Momentum gated)')
        ax.set_xlabel('TOF (ns)')
        ax.set_ylabel('X (mm)')
        plt.savefig(basedir+'tof_x_momentum_gated.png',bbox_inches='tight')

    if plot_det_img:
        cmap='jet'
        fig = plt.figure(figsize=(18,4))
    
        x,y,z=fhist2d(x1,y1,-50,50,0.5,-50,50,0.5)
        ax = fig.add_subplot(121)
        f1=ax.pcolormesh(x,y,z,cmap=cmap)
        fig.colorbar(f1)
        # F1 = (X-x1_c)**2 + (Y-y1_c)**2 - R
        # ax.contour(X, Y, F1,[0], colors='r', linewidths=3)
        # ax.set_title(label_species[0])
        # ax.set_aspect('equal')
    
        x,y,z=fhist2d(x2,y2,-50,50,0.5,-50,50,0.5)
        ax = fig.add_subplot(122)
        f2=ax.pcolormesh(x,y,z,cmap=cmap)
        fig.colorbar(f2)
        # F1 = (X-x2_c)**2 + (Y-y2_c)**2 - R
        # ax.contour(X, Y, F1,[0], colors='r', linewidths=3)
        #ax.set_title(label_species[1])
        # ax.set_aspect('equal')
        
    if plot_tof:
        fig,ax = plt.subplots(figsize=(18,4))
        x11,y11=fhist1d(t1,tof_x_range[0],tof_x_range[1],tof_x_range[2])
        x22,y22=fhist1d(t2,tof_x_range[0],tof_x_range[1],tof_x_range[2])
        ax.plot(x11,y11)
        ax.plot(x22,y22)
        ax.axvline(x=tof_frag[0],color='g',linestyle='--')
        ax.axvline(x=tof_frag[1],color='g',linestyle='--')
        val=np.max([y11,y22])
        ax.text(tof_frag[0]+0.01*tof_frag[0],val-(0.05*val),label_species[0])
        ax.text(tof_frag[1]+0.01*tof_frag[1],val-(0.05*val),label_species[1])
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.set_xlabel('TOF (ns)')
        plt.savefig(basedir+'tof_momentum_gated.png',bbox_inches='tight')


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
        plt.savefig(basedir+'KER_cosz_momentum_gated.png',bbox_inches='tight')
    return KERgx,KERgy

    

