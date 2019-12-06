# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:47:14 2019

@author: shashank
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from fast_histogram import histogram1d,histogram2d
plt.rcParams.update({'font.size': 16})

def fhist1d(x1,xmin,xmax,xbinsize):
    xbins=int((xmax-xmin)/xbinsize)
    x=np.linspace(xmin,xmax,xbins+1)
    centersx=(x[:-1]+x[1:])/2
    test_z=histogram1d(x1,bins=xbins,range=[xmin,xmax])
    return centersx,test_z


def fhist2d(x1,y1,xmin,xmax,xbinsize,ymin,ymax,ybinsize):
    xbins=int((xmax-xmin)/xbinsize)
    ybins=int((ymax-ymin)/ybinsize)
    x=np.linspace(xmin,xmax,xbins+1)
    y=np.linspace(ymin,ymax,ybins+1)
    centersx=(x[:-1]+x[1:])/2
    centersy=(y[:-1]+y[1:])/2
    test_z=histogram2d(x1,y1,bins=[xbins,ybins],range=[[xmin,xmax],[ymin,ymax]])
    return centersx,centersy,test_z

def PA(ma,mb,mc,M,mab,muab,pxABC,pyABC,pzABC,pxABcC,pyABcC,pzABcC,pxAB,pyAB,pzAB):
    t1=mab*mb*pxAB
    t2=mb*mc*pxAB
    t3=-mab*pxABC*muab
    t4=M*pxABcC*muab
    arr_int=np.array((t1,t2,t3,t4))
    PAx=(arr_int.sum(axis=0))*(-ma/((ma+mb)*(mab+mc)*muab))

    t1=mab*mb*pyAB
    t2=mb*mc*pyAB
    t3=-mab*pyABC*muab
    t4=M*pyABcC*muab
    arr_int=np.array((t1,t2,t3,t4))
    PAy=(arr_int.sum(axis=0))*(-ma/((ma+mb)*(mab+mc)*muab))

    t1=mab*mb*pzAB
    t2=mb*mc*pzAB
    t3=-mab*pzABC*muab
    t4=M*pzABcC*muab
    arr_int=np.array((t1,t2,t3,t4))
    PAz=(arr_int.sum(axis=0))*(-ma/((ma+mb)*(mab+mc)*muab))
    return PAx,PAy,PAz



def PB(ma,mb,mc,M,mab,muab,pxABC,pyABC,pzABC,pxABcC,pyABcC,pzABcC,pxAB,pyAB,pzAB):
    t1=ma*(mab+mc)*pxAB
    t2=mab*pxABC*muab
    t3=-muab*M*pxABcC
    arr_int=np.array((t1,t2,t3))
    PBx=(arr_int.sum(axis=0))*(mb/((ma+mb)*(mab+mc)*muab))

    t1=ma*(mab+mc)*pyAB
    t2=mab*pyABC*muab
    t3=-muab*M*pyABcC
    arr_int=np.array((t1,t2,t3))
    PBy=(arr_int.sum(axis=0))*(mb/((ma+mb)*(mab+mc)*muab))

    t1=ma*(mab+mc)*pzAB
    t2=mab*pzABC*muab
    t3=-muab*M*pzABcC
    arr_int=np.array((t1,t2,t3))
    PBz=(arr_int.sum(axis=0))*(mb/((ma+mb)*(mab+mc)*muab))
    return PBx,PBy,PBz

def PC(mc,M,mab,pxABC,pyABC,pzABC,pxABcC,pyABcC,pzABcC):
    PCx=np.divide(np.add((mc*pxABC),(M*pxABcC)),(mab+mc))    ### Remember AB -> 13 for this case
    PCy=np.divide(np.add((mc*pyABC),(M*pyABcC)),(mab+mc))
    PCz=np.divide(np.add((mc*pzABC),(M*pzABcC)),(mab+mc))
    return PCx,PCy,PCz


#def cos12KE12(p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z,m1,m2,m3,m12,M,mu12):
#    p123x=np.add(((m12/M)*p3x),-((m3/M)*np.add(p1x,p2x)))
#    p123y=np.add(((m12/M)*p3y),-((m3/M)*np.add(p1y,p2y)))
#    p123z=np.add(((m12/M)*p3z),-((m3/M)*np.add(p1z,p2z)))
#
#    p12x=mu12*(np.add(np.divide(p2x,m2),-np.divide(p1x,m1)))
#    p12y=mu12*(np.add(np.divide(p2y,m2),-np.divide(p1y,m1)))
#    p12z=mu12*(np.add(np.divide(p2z,m2),-np.divide(p1z,m1)))
#
#    #Dot
#    dx=np.multiply(p123x,p12x)
#    dy=np.multiply(p123y,p12y)
#    dz=np.multiply(p123z,p12z)
#    dot=np.add(np.add(dx,dy),dz)
#    ptt=np.array((p12x**2,p12y**2,p12z**2))
#    P12 = np.sqrt(ptt.sum(axis=0))
#    pttt=np.array((p123x**2,p123y**2,p123z**2))
#    P123 = np.sqrt(pttt.sum(axis=0))
#
#    cosH12H3=np.divide(dot,(P123*P12))
#    test=np.arccos(cosH12H3)
#    cosH12H3=np.multiply(test,np.divide(180,np.pi))
#    KE12=np.divide((P12**2),(2*mu12*66.4))
#
#    return KE12,cosH12H3


def cosAB_keAB(pall,mass,single_frag,keABcosABrange,label_species,cmap,basedir,log=None):
    "if o/p cos12 and KE12, single_frag=3 // if o/p cos13 and KE13, single_frag=2"
    if single_frag==1:
        p_all=np.roll(pall,-1*3,axis=0)
        m_all=np.roll(mass,-1,axis=0)
        label=np.roll(label_species,-1,axis=0)
    elif single_frag==2:
        p_all=np.roll(pall,-2*3,axis=0)
        m_all=np.roll(mass,-2,axis=0)
        label=np.roll(label_species,-2,axis=0)
    elif single_frag==3:
        p_all=pall
        m_all=mass
        label=label_species
    p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z=p_all
    m_all=m_all/(1.6E-27)
    m12=m_all[0]+m_all[1]
    M=m12+m_all[2]
    mu12=(m_all[0]*m_all[1])/(m12)
    p123x=np.add(((m12/M)*p3x),-((m_all[2]/M)*np.add(p1x,p2x)))
    p123y=np.add(((m12/M)*p3y),-((m_all[2]/M)*np.add(p1y,p2y)))
    p123z=np.add(((m12/M)*p3z),-((m_all[2]/M)*np.add(p1z,p2z)))

    p12x=mu12*(np.add(np.divide(p2x,m_all[1]),-np.divide(p1x,m_all[0])))
    p12y=mu12*(np.add(np.divide(p2y,m_all[1]),-np.divide(p1y,m_all[0])))
    p12z=mu12*(np.add(np.divide(p2z,m_all[1]),-np.divide(p1z,m_all[0])))

    #Dot
    dx=np.multiply(p123x,p12x)
    dy=np.multiply(p123y,p12y)
    dz=np.multiply(p123z,p12z)
    dot=np.add(np.add(dx,dy),dz)
    ptt=np.array((p12x**2,p12y**2,p12z**2))
    P12 = np.sqrt(ptt.sum(axis=0))
    pttt=np.array((p123x**2,p123y**2,p123z**2))
    P123 = np.sqrt(pttt.sum(axis=0))

    cosAB=np.divide(dot,(P123*P12))
    test=np.arccos(cosAB)
    cosAB=np.multiply(test,np.divide(180,np.pi))
#    KE12=np.divide((P12**2),(2*mu12*66.4))
    keAB=((P12**2/(0.25E48))/(2*mu12*(1.6E-27)))*6.242E18

    x,y,z=fhist2d(keAB,cosAB,keABcosABrange[0],keABcosABrange[1],keABcosABrange[2],keABcosABrange[3],keABcosABrange[4],keABcosABrange[5])
    fig, ax = plt.subplots()
    if log:
        f1=ax.pcolormesh(x,y,np.transpose(z),norm=LogNorm(),cmap=cmap)
    else:
        f1=ax.pcolormesh(x,y,np.transpose(z),cmap=cmap)
    fig.colorbar(f1)
    ax.set_xlabel('KE(%s-%s) (eV)'%(label[0],label[1]))
    ax.set_ylabel(r'$\theta$(%s-%s,%s) (deg)'%(label[0],label[1],label[2]))
    plt.savefig(basedir+'KER_costh_%s_gated.png'%str(single_frag),bbox_inches='tight')
#    fig, ax = plt.subplots()
#    x,y=fhist1d(keAB,keABcosABrange[0],keABcosABrange[1],keABcosABrange[2])
#    ax.plot(x,y)
#    ax.set_xlabel('KE(%s-%s) (eV)'%(label[0],label[1]))
    return keAB,cosAB

" Newton Plot with 3rd fragment fixed with momentum 1 along x axis "

def NP3(p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z):

    pt = np.array((p1x**2,p1y**2,p1z**2))
    P1 = np.sqrt(pt.sum(axis=0))
    ptt = np.array((p2x**2,p2y**2,p2z**2))
    P2 = np.sqrt(ptt.sum(axis=0))
    pttt = np.array((p3x**2,p3y**2,p3z**2))
    P3 = np.sqrt(pttt.sum(axis=0))

    p1xyz=np.array((p1x,p1y,p1z))
    p2xyz=np.array((p2x,p2y,p2z))
    p3xyz=np.array((p3x,p3y,p3z))

    p1mod=np.sqrt((p1xyz**2).sum(axis=0))
    p2mod=np.sqrt((p2xyz**2).sum(axis=0))
    p3mod=np.sqrt((p3xyz**2).sum(axis=0))

    dot23=(p2xyz*p3xyz).sum(axis=0)
    mod23=np.multiply(p2mod,p3mod)
    cosH2H3=np.divide(dot23,mod23)

    dot13=(p1xyz*p3xyz).sum(axis=0)
    mod13=np.multiply(p1mod,p3mod)
    cosH1H3=np.divide(dot13,mod13)

    x1=(np.divide(-P1,-P3))*cosH1H3
    y1=-(np.divide(P1,-P3))*np.sqrt(np.add(1,-cosH1H3**2))
    x2=(np.divide(-P2,-P3))*(cosH2H3)
    y2=(np.divide(P2,-P3))*np.sqrt(np.add(1,-cosH2H3**2))

    return x1,y1,x2,y2


" Plot Normalised Newton Plot - more general "

def NP(p1x,p1y,p1z,label1,p2x,p2y,p2z,label2,p3x,p3y,p3z,label3,fragment_fixed,rangemin,rangemax,binsize,filename,clim_min,clim_max,run_name):
    "input are momenta of all three fragments, order is important! 1->2->3"
    "fragment_fixed=>which fragment you want to fix with momentum 1 along x axis "
    "rangemax,rangemin,binsize are for plotting 2d histogram"
    "Plot is normalised by no of coincident events i.e. summing the histograms 2D matrix, s then each bin value is divided by s/2(counts have values from 2 hits) "
    "clim -> specify lower and upper limit in z for histogram"
    mom1x=p1x
    mom1y=p1y
    mom1z=p1z
    mom2x=p2x
    mom2y=p2y
    mom2z=p2z
    mom3x=p3x
    mom3y=p3y
    mom3z=p3z
    l1=label1
    l2=label2
    l3=label3

    if fragment_fixed==1:   #2 in top quadrant, # 3 in bottom quadrant
        p1x=mom2x
        p1y=mom2y
        p1z=mom2z
        p2x=mom3x
        p2y=mom3y
        p2z=mom3z
        p3x=mom1x
        p3y=mom1y
        p3z=mom1z
        label1=l2
        label2=l3
        label3=l1
    elif fragment_fixed==2:   #1 in top quadrant, # 3 in bottom quadrant
        p1x=mom1x
        p1y=mom1y
        p1z=mom1z
        p2x=mom3x
        p2y=mom3y
        p2z=mom3z
        p3x=mom2x
        p3y=mom2y
        p3z=mom2z
        label1=l1
        label2=l3
        label3=l2
    elif fragment_fixed==3:  #1 in top quadrant, # 2 in bottem quadrant
        pass
    else:
        raise NameError('fragment_fixed takes value from 1-3')

    pt = np.array((p1x**2,p1y**2,p1z**2))
    P1 = np.sqrt(pt.sum(axis=0))
    ptt = np.array((p2x**2,p2y**2,p2z**2))
    P2 = np.sqrt(ptt.sum(axis=0))
    pttt = np.array((p3x**2,p3y**2,p3z**2))
    P3 = np.sqrt(pttt.sum(axis=0))

    p1xyz=np.array((p1x,p1y,p1z))
    p2xyz=np.array((p2x,p2y,p2z))
    p3xyz=np.array((p3x,p3y,p3z))

    p1mod=np.sqrt((p1xyz**2).sum(axis=0))
    p2mod=np.sqrt((p2xyz**2).sum(axis=0))
    p3mod=np.sqrt((p3xyz**2).sum(axis=0))

    dot23=(p2xyz*p3xyz).sum(axis=0)
    mod23=np.multiply(p2mod,p3mod)
    cosH2H3=np.divide(dot23,mod23)

    dot13=(p1xyz*p3xyz).sum(axis=0)
    mod13=np.multiply(p1mod,p3mod)
    cosH1H3=np.divide(dot13,mod13)

    x1=(np.divide(-P1,-P3))*cosH1H3
    y1=-(np.divide(P1,-P3))*np.sqrt(np.add(1,-cosH1H3**2))
    x2=(np.divide(-P2,-P3))*(cosH2H3)
    y2=(np.divide(P2,-P3))*np.sqrt(np.add(1,-cosH2H3**2))


    #Plot graph
#    NPRange = []
#    NPh = hist2d(*init_hist2d(*NPRange))
#    NPh.fill(x1,y1)
#    NPh.fill(x2,y2)
    a1,b1,c1=fhist2d(x1,y1,rangemin,rangemax,binsize,rangemin,rangemax,binsize)
    a2,b2,c2=fhist2d(x2,y2,rangemin,rangemax,binsize,rangemin,rangemax,binsize)
#    NPhx, NPhy, NPhz = NPh.centerData
    NPhz=c1+c2
    counts=sum2D(NPhz)
#    NPhz=np.divide(NPhz,(counts/2))    # Normalizing
    NPhz=np.divide(NPhz,np.max(NPhz))
    fac=1+binsize
    plt.figure()
    plt.pcolormesh(a1,b1,np.transpose(NPhz), cmap='jet')#,norm=LogNorm())
#    plt.arrow( 0, 0, 1, 0, fc='white', ec='white',lw=1,head_width=0.05, head_length=0.1)
    plt.arrow( 0, 0, 1, 0, fc='black', ec='black',lw=1,head_width=0.05, head_length=0.1)
    plt.text(rangemin+fac,rangemax-fac,label1,fontweight='bold',color='red')
    plt.text(rangemin+fac,rangemin+fac-0.35,label2,fontweight='bold',color='red')
    plt.text(rangemax-1.25*fac,0.08,label3,fontweight='bold',color='red')
#    if (run_name[-4:]=='105C') or (run_name[-4:]=='160C'):
#        plt.title('Newton plot (T = '+run_name[-4:]+')')
#    else:
    plt.title('Newton plot')
    plt.colorbar()
    if (clim_min==0) & (clim_max==0):
        pass
    else:
        plt.clim(clim_min,clim_max)
    plt.xlabel('P(arb units)')
    plt.ylabel('P(arb units)')
    plt.savefig(filename, bbox_inches="tight")
#    plt.grid(True)


    return x1,y1,x2,y2,counts

"Calculating Cos theta between any fragments"

def cosf1f2(p1x,p1y,p1z,p2x,p2y,p2z):
    p1xyz=np.array((p1x,p1y,p1z))
    p2xyz=np.array((p2x,p2y,p2z))
#    p3xyz=np.array((p3x,p3y,p3z))
    p1mod=np.sqrt((p1xyz**2).sum(axis=0))
    p2mod=np.sqrt((p2xyz**2).sum(axis=0))
#    p3mod=np.sqrt((p3xyz**2).sum(axis=0))

    dot12=(p1xyz*p2xyz).sum(axis=0)
    mod12=np.multiply(p1mod,p2mod)
    cos12=np.divide(dot12,mod12)
    return cos12

#def cos13(p1x,p1y,p1z,p3x,p3y,p3z):
#    p1xyz=np.array((p1x,p1y,p1z))
##    p2xyz=np.array((p2x,p2y,p2z))
#    p3xyz=np.array((p3x,p3y,p3z))
#    p1mod=np.sqrt((p1xyz**2).sum(axis=0))
##    p2mod=np.sqrt((p2xyz**2).sum(axis=0))
#    p3mod=np.sqrt((p3xyz**2).sum(axis=0))
#
#    dot13=(p1xyz*p3xyz).sum(axis=0)
#    mod13=np.multiply(p1mod,p3mod)
#    cosH1H3=np.divide(dot13,mod13)
#    return cosH1H3
#
#def cos23(p2x,p2y,p2z,p3x,p3y,p3z):
##    p1xyz=np.array((p1x,p1y,p1z))
#    p2xyz=np.array((p2x,p2y,p2z))
#    p3xyz=np.array((p3x,p3y,p3z))
##    p1mod=np.sqrt((p1xyz**2).sum(axis=0))
#    p2mod=np.sqrt((p2xyz**2).sum(axis=0))
#    p3mod=np.sqrt((p3xyz**2).sum(axis=0))
#
#    dot23=(p2xyz*p3xyz).sum(axis=0)
#    mod23=np.multiply(p2mod,p3mod)
#    cosH2H3=np.divide(dot23,mod23)
#    return cosH2H3

"This function sums total counts of a 2D histogram (sums 2D array)"

def sum2D(z2d):
    z1d=np.sum(z2d,axis=0)
    z=np.sum(z1d)
    return z

def threeConditions(ionTOF, cond1, cond2, cond3):
    "credit-Razib Obaid(UConn)"
    indices = ionTOF.size
    checkList = np.zeros(indices, dtype=bool)
    for i in range(0, indices):
        if cond1[0] < ionTOF[i] < cond1[1]:
            checkList[i] = True
            for j in range(i+1, indices):
                if cond2[0] < ionTOF[j] < cond2[1]:
                    checkList[j] = True
                    for k in range(j+1, indices):
                        if cond3[0] < ionTOF[k] < cond3[1]:
                            checkList[k] = True
                            break
                    break
            break

    return checkList
