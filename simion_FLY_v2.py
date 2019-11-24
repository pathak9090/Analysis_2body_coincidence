# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 17:36:53 2019

@author: Shashank Pathak (credit: Razib Obaid)
"""

import numpy as np
### Simion x-> actual z, simion y-> actual x, simion z actual y
def flywrite(filename):
    filename_all=[filename+'_vx',filename+'_vy',filename+'_tof']
    for j in range(0,3):
        f = open(filename_all[j]+'.fly2', 'w')
        if j==0:
            vec='(0,1,0)'
        elif j==1:
            vec='(0,0,1)'
        elif j==2:
            vec='(1,0,0)'
        #Writing the Simion FLY file from here
        foreText = "particles {coordinates = 0,"
        f.write(foreText)
    
        for i in np.arange(0.1, 10, 0.1):
            midText = "standard_beam {n = 1, tob = 0, mass = 39, charge = 1, \
            ke ="+ str(i) + ", cwf = 1, color = 0, direction = vector"+vec+", position = vector(0, 0, 0)},"
            
            f.write(midText)
    
        f.write( "}")
        f.close()
    
if __name__=='__main__':
    filename = 'C3H3v2'
    flywrite(filename)
