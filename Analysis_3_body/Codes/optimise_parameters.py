# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 02:44:49 2019

@author: Shashank (credit-Nathan)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import math as ma
from scipy.interpolate import Rbf
from sklearn.externals import joblib
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import tkinter as tk
import numpy as np
from scipy.optimize import curve_fit
import tkinter
from useful_definitions import fhist1d,fhist2d
#
#x=np.random.random([20])
#y=np.random.random([20])

#top = tkinter.Tk()
#top.geometry("2000x2000")
class optimise_parameters:
    def __init__(self,top,xyt_array,vel_calc_fact,param_list,mass,Run_name,channel):
#        # Code to add widgets will go here...
#        self.x=x
#        self.y=y
        pos_offset_1x,pos_offset_1y,pos_offset_2x,pos_offset_2y,pos_offset_3x,pos_offset_3y,t0=param_list
        self.vel_calc_fact=vel_calc_fact
        self.mass=mass
        self.xyt_array=xyt_array
        self.Run_name=Run_name
        self.channel=channel
        self.param_list=param_list
        top.title('Parameter Optimisation GUI')
        top.bind('<Return>', self.refresh)

        self.label1=tk.Label(top,text='t0')
        self.label1.grid(row=0,column=0)
        self.entrywin1=tk.Entry(top,width=15)
        self.entrywin1.grid(row=1,column=0)
        self.entrywin1.insert(0, t0)

        self.label2=tk.Label(top,text='x1')
        self.label2.grid(row=0,column=1)
        self.entrywin2=tk.Entry(top,width=15)
        self.entrywin2.grid(row=1,column=1)
        self.entrywin2.insert(0, pos_offset_1x)

        self.label3=tk.Label(top,text='y1')
        self.label3.grid(row=0,column=2)
        self.entrywin3=tk.Entry(top,width=15)
        self.entrywin3.grid(row=1,column=2)
        self.entrywin3.insert(0, pos_offset_1y)

        self.label4=tk.Label(top,text='x2')
        self.label4.grid(row=0,column=3)
        self.entrywin4=tk.Entry(top,width=15)
        self.entrywin4.grid(row=1,column=3)
        self.entrywin4.insert(0, pos_offset_2x)

        self.label5=tk.Label(top,text='y2')
        self.label5.grid(row=0,column=4)
        self.entrywin5=tk.Entry(top,width=15)
        self.entrywin5.grid(row=1,column=4)
        self.entrywin5.insert(0, pos_offset_2y)

        self.label6=tk.Label(top,text='x3')
        self.label6.grid(row=0,column=5)
        self.entrywin6=tk.Entry(top,width=15)
        self.entrywin6.grid(row=1,column=5)
        self.entrywin6.insert(0, pos_offset_3x)

        self.label7=tk.Label(top,text='y3')
        self.label7.grid(row=0,column=6)
        self.entrywin7=tk.Entry(top,width=15)
        self.entrywin7.grid(row=1,column=6)
        self.entrywin7.insert(0, pos_offset_3y)


        self.button1=tk.Button(top,command=self.write_parameters,text='Write Parameters')
        self.button1.grid(row=1,column=7)

        self.fig=Figure(figsize=(18,9),dpi=100)
        self.ax1 = plt.subplot2grid((4,3),(0,0), fig=self.fig)
        self.ax2 = plt.subplot2grid((4,3),(0,1), fig=self.fig)
        self.ax3 = plt.subplot2grid((4,3),(0,2), fig=self.fig)
        self.ax4 = plt.subplot2grid((4,3),(1,0), fig=self.fig)
        self.ax5 = plt.subplot2grid((4,3),(1,1), fig=self.fig)
        self.ax6 = plt.subplot2grid((4,3),(1,2), fig=self.fig)
        self.ax7 = plt.subplot2grid((4,3),(2,0), fig=self.fig)
        self.ax8 = plt.subplot2grid((4,3),(2,1), fig=self.fig)
        self.ax9 = plt.subplot2grid((4,3),(2,2), fig=self.fig)
        self.ax10 = plt.subplot2grid((4,3),(3,0), fig=self.fig)
        self.ax11 = plt.subplot2grid((4,3),(3,1), fig=self.fig)
        self.ax12 = plt.subplot2grid((4,3),(3,2), fig=self.fig)
        self.canvas=FigureCanvasTkAgg(self.fig,master=top)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=6,columnspan=8)
        self.toolbar_frame = tk.Frame(top)
        self.toolbar_frame.grid(row=10, column=1,columnspan=6)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbar_frame)
        self.toolbar.update()

        self.plots()
#        self.button1=tk.Button(top,command=self.write_file,text='submit')
#        self.button1.grid(row=5,column=0)

    def write_file(self):
        print(self.entrywin1.get())

    def plots(self):
#        x=np.random.random([20])
#        y=np.random.random([20])
#        x=self.x
#        y=self.y
        t0=float(self.entrywin1.get())
        pos_offset_1x=float(self.entrywin2.get())
        pos_offset_1y=float(self.entrywin3.get())
        pos_offset_2x=float(self.entrywin4.get())
        pos_offset_2y=float(self.entrywin5.get())
        pos_offset_3x=float(self.entrywin6.get())
        pos_offset_3y=float(self.entrywin7.get())

        self.param_list=[pos_offset_1x,pos_offset_1y,pos_offset_2x,pos_offset_2y,pos_offset_3x,pos_offset_3y,t0]

        mx1,cx1,mz1,cz1,mx2,cx2,mz2,cz2,mx3,cx3,mz3,cz3=self.vel_calc_fact
        ion1xgc,ion1ygc,ion1tgc,ion2xgc,ion2ygc,ion2tgc,ion3xgc,ion3ygc,ion3tgc=self.xyt_array
        m1,m2,m3=self.mass
#        pos_offset_1x,pos_offset_1y,pos_offset_2x,pos_offset_2y,pos_offset_3x,pos_offset_3y,t0=self.param_list
        ion1xgc=np.add(ion1xgc,-pos_offset_1x)
        ion2xgc=np.add(ion2xgc,-pos_offset_2x)
        ion3xgc=np.add(ion3xgc,-pos_offset_3x)
        ion1ygc=np.add(ion1ygc,-pos_offset_1y)
        ion2ygc=np.add(ion2ygc,-pos_offset_2y)
        ion3ygc=np.add(ion3ygc,-pos_offset_3y)

        vx1 = (mx1*ion1xgc)+cx1
        vy1 = (mx1*ion1ygc)+cx1
        vz1 = (mz1*(ion1tgc-t0))+cz1

        vx2 = (mx2*ion2xgc)+cx2
        vy2 = (mx2*ion2ygc)+cx2
        vz2 = (mz2*(ion2tgc-t0))+cz2

        vx3 = (mx3*ion3xgc)+cx3
        vy3 = (mx3*ion3ygc)+cx3
        vz3 = (mz3*(ion3tgc-t0))+cz3

        P1x = 1e3 * 0.5E24 * m1 * vx1
        P1y = 1e3 * 0.5E24 * m1 * vy1
        P1z = 1e3 * 0.5E24 * m1 * vz1
        P2x = 1e3 * 0.5E24 * m2 * vx2
        P2y = 1e3 * 0.5E24 * m2 * vy2
        P2z = 1e3 * 0.5E24 * m2 * vz2
        P3x = 1e3 * 0.5E24 * m3 * vx3
        P3y = 1e3 * 0.5E24 * m3 * vy3
        P3z = 1e3 * 0.5E24 * m3 * vz3

        Px=np.array((P1x,P2x,P3x))
        Py=np.array((P1y,P2y,P3y))
        Pz=np.array((P1z,P2z,P3z))
#        KE=np.array((KE1,KE2,KE3))
        prange=[-500,500,1]
        Px_all=Px.sum(axis=0)
        Py_all=Py.sum(axis=0)
        Pz_all=Pz.sum(axis=0)
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.ax4.clear()
        self.ax5.clear()
        self.ax6.clear()
        self.ax7.clear()
        self.ax8.clear()
        self.ax9.clear()
        self.ax10.clear()
        self.ax11.clear()
        self.ax12.clear()
#        self.ax1.plot(x,y)
        x,y=fhist1d(P1x,prange[0],prange[1],prange[2])
        self.ax1.plot(x,y,label='P1x')
        self.ax1.plot(-x,y,label='-P1x')
        self.ax1.legend(frameon=False)
        x,y=fhist1d(P1y,prange[0],prange[1],prange[2])
        self.ax2.plot(x,y,label='P1y')
        self.ax2.plot(-x,y,label='-P1y')
        self.ax2.legend(frameon=False)
        x,y=fhist1d(P1z,prange[0],prange[1],prange[2])
        self.ax3.plot(x,y,label='P1z')
        self.ax3.plot(-x,y,label='-P1z')
        self.ax3.legend(frameon=False)
        x,y=fhist1d(P2x,prange[0],prange[1],prange[2])
        self.ax4.plot(x,y,label='P2x')
        self.ax4.plot(-x,y,label='-P2x')
        self.ax4.legend(frameon=False)
        x,y=fhist1d(P2y,prange[0],prange[1],prange[2])
        self.ax5.plot(x,y,label='P2y')
        self.ax5.plot(-x,y,label='-P2y')
        self.ax5.legend(frameon=False)
        x,y=fhist1d(P2z,prange[0],prange[1],prange[2])
        self.ax6.plot(x,y,label='P2z')
        self.ax6.plot(-x,y,label='-P2z')
        self.ax6.legend(frameon=False)
        x,y=fhist1d(P3x,prange[0],prange[1],prange[2])
        self.ax7.plot(x,y,label='P3x')
        self.ax7.plot(-x,y,label='-P3x')
        self.ax7.legend(frameon=False)
        x,y=fhist1d(P3y,prange[0],prange[1],prange[2])
        self.ax8.plot(x,y,label='P3y')
        self.ax8.plot(-x,y,label='-P3y')
        self.ax8.legend(frameon=False)
        x,y=fhist1d(P3z,prange[0],prange[1],prange[2])
        self.ax9.plot(x,y,label='P3z')
        self.ax9.plot(-x,y,label='-P3z')
        self.ax9.legend(frameon=False)
        x,y=fhist1d(Px_all,prange[0],prange[1],prange[2])
        self.ax10.plot(x,y,label='Px sum')
        self.ax10.plot(-x,y,label='-Px sum')
        self.ax10.legend(frameon=False)
        x,y=fhist1d(Py_all,prange[0],prange[1],prange[2])
        self.ax11.plot(x,y,label='Py sum')
        self.ax11.plot(-x,y,label='-Py sum')
        self.ax11.legend(frameon=False)
        x,y=fhist1d(Pz_all,prange[0],prange[1],prange[2])
        self.ax12.plot(x,y,label='Pz sum')
        self.ax12.plot(-x,y,label='-Pz sum')
        self.ax12.legend(frameon=False)
        self.canvas.draw()
    def refresh(self,event):
        self.plots()

    def write_parameters(self):
        with open(self.Run_name+'_'+self.channel+'_final_parameters.txt', mode='w') as file:
            for param in self.param_list:
                file.write(str(param))
                file.write('\n')

#optimise_parameters(top,x,y)
#top.mainloop()
