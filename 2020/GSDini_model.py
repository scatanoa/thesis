#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def readList2array(fn):
#    Xo = open(fn)
#    Xr = Xo.readlines()
#    Xo.close()
##    Xr.split(',')
#    print(f'Xr={Xr}')
#    Xr.replace('\n','')    
#    X = np.array(Xr)
    X = np.loadtxt(fn, delimiter=',')
    return X

def readFlumeFlood(which_flume, flume_else_basin):
    fp = 'doct_fluvial/hg_model_flume/calib_flood_pitlick13/pp13_gsd_'
    D = readList2array(fp + 'D.txt')
    P = readList2array(fp + 'P.txt')
    print(f'9apr20, P=\n{P}')
    print(f'9apr20, D=\n{D}')
    dP = np.diff(P)
    fb = [P[0]]
    fb.extend(list(dP))
    Dmm = 1e3*D
    Dc_mm = getDcent(Dmm)
    plot(which_flume, flume_else_basin, P, fb, Dmm, Dc_mm)    
    return fb, D, Dc_mm    

def calc(nD, flume_else_basin, fl_gsd, which_flume):
#see table3 in nonlinear bedload paper of recking13
#mar24, 2020.
#    if flume_else_basin==0: #ref sklar17 for basin?
#        Fs = .2 #higher than bed GSD curves in recking13, as here it is intended to produce supply GSD.
#        D84 = 1e-1
#        D50 = D84/2
#    else:
    Fs = fl_gsd[0]  #.05 #higher than real bed GSD curves in recking13, as here it is intended to produce supply GSD.
    D50 = fl_gsd[1]
    D84 = fl_gsd[2]
    
    Dmn = 2e-3 #portion of sizes below Dmn is Fs = sand fraction.
    cn = [-99, 16, 3.3, 1.9, 1.3, -99, 5.9, 2.3, -99, 1.3, 2.5, 5.1]
    cn = np.array(cn)
    P = np.zeros((12))
    P[0] = 100*Fs
    P[1:5] = 100*Fs + (50-100*Fs)/cn[1:5]
    P[5:] = [50, 60, 70, 84, 90, 98, 100]
    D = np.zeros((12))
    D[0] = Dmn
    D[1:5] = (D50 - Dmn)/cn[1:5] + Dmn
    D[5] = D50
    D[6:8] = (D84 - D50)/cn[6:8] + D50
    D[8] = D84
    D[9:] = cn[9:]*D84
    P = P/100
    if flume_else_basin: D[10:] = [.029, .032]      #done to fit with real ubcTrench gsd.
    P, D = get_GSDnD(nD, P, D)
#    print(f'26mar, after resampling for nD, P={P}')
#    print(f'D={D}')    
    dP = np.diff(P)
    fb = [P[0]]
    fb.extend(list(dP))
    Dmm = 1e3*D
    Dc_mm = getDcent(Dmm)
    plot(which_flume, flume_else_basin, P, fb, Dmm, Dc_mm)
    return fb, D, Dc_mm

def get_GSDnD(nD, P, D):
    if nD==5:
        ip=[0,5,8,9,11]
    elif nD==6:
        ip=[0,3,5,8,9,11]
    elif nD==7:
        ip=[5,6,7,8,9,10,11]
    elif nD==8:
        ip=[0,1,3,5,8,9,10,11]
    elif nD==9:
        ip=[0,1,3,5,7,8,9,10,11]
    elif nD==10:
        ip=[0,1,2,3,5,7,8,9,10,11]        
    elif nD==11:
        ip=[0,1,2,3,5,6,7,8,9,10,11]
    elif nD==12:
        ip=[0,1,2,3,4,5,6,7,8,9,10,11]
    P2 = []; D2 = []
#    print(f'26mar, ip={ip}')     
    for i in ip:
        P2.append(P[i]); D2.append(D[i])
    return np.array(P2), np.array(D2)

def plot(which_flume, flume_else_basin, P,fb,Dmm,Dc_mm):
    fig, ax = plt.subplots()
    l1=ax.plot(Dmm, P, 'k', marker=11, label = 'Fraction finer')
    ax.set_xscale('log')
    ax.set_ylim(bottom= 0, top=1)
    if flume_else_basin==1: ax.set_xlim(left= 1e-1, right=1e2)
    ax.set_ylabel('Fraction finer')
    ax.set_xlabel('Grain size (mm)')
    ax1 = ax.twinx()
    l2=ax1.plot(Dc_mm, fb, 'r', ls='--', marker=6, label = 'Frequency') 
    ax1.set_ylabel(f'Frequency, whose sum is {sum(fb)}')
    ax.grid(which='both')
    
    ls=l1+l2
    labs = [l.get_label() for l  in ls]
    leg = ax.legend(ls, labs, loc='best', ncol=1, fancybox=True, fontsize=10) #ncol=2, 
    leg.get_frame().set_alpha(0.7)
    flumName = ['Trench', 'Flood']
    fn = 'basin' if flume_else_basin==0 else 'flume'
    spec = '' if flume_else_basin==0 else flumName[which_flume-1]
    plt.tight_layout()
    plt.savefig(f'doct_fluvial/GSDini_{fn}{spec}.png')
    plt.close('all')

def getDcent(D):
    Dc = [np.sqrt(D[0]*.062)] #mean sand size (mm)
    nD = len(D)
    for k in range(2,nD+1):
        Dc.append(np.sqrt(D[k-1]*D[k-2]))
    Dc = np.array(Dc)
    return Dc
