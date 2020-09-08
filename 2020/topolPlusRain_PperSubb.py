    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#dec20, 2019

#pre-existing modules
import numpy as np
import matplotlib.pyplot as plt
import time

import mapPlot as mp
#import bin_io as io

def printCumRumTime(t1):
    t2 = time.time()
    print(f'run_t[s] = {t2-t1}')


def P_perSubbasin(): #each day, get matrix: P(t= days, x= reaches):  
    t1 = time.time()
    X = np.loadtxt('doct_topol/Topol_X.txt')
    t_ev = np.loadtxt('doct_rain/t_ev.txt')
    Px0 = np.loadtxt('doct_rain/Px_day' + str(int(t_ev[0])) + '.txt')
#    for y in range(1, sizeNorth+1):
#        Px0[y-1, :], res = io.read_float_arr(fileunit, path_list[pos_ini_var[0]-1], y, sizeEast)
            #take path_list and pos_ini_var from txt file.
            #sizeNorth and sizeEast must have been saved in another txt file when rain fields were created.            
    Px0_sh = Px0.shape #1st and 2nd components say rows/y/northing and columns/x/easting lenght respectively.
    sizeNorth = Px0_sh[0];  sizeEast = Px0_sh[1]
    print('sides: ', sizeNorth, sizeEast)   
    n_reaches = len(X[:, 0])
    n_days_ev = len(t_ev)    
    Pday_subb = np.zeros((sizeNorth, sizeEast, n_days_ev))    
    P_perDayperReach = np.zeros((n_days_ev, n_reaches))    
    Masks_bin = np.zeros((sizeNorth, sizeEast, n_reaches))    
    for j in range(1, (n_reaches-1)+1): #because last X row is not upstream point of a reach.
        reachID = j
        Mask_bin_currentSubb = np.loadtxt('doct_supply/geom/Mask_bin_currentSubb_reachID' + str(reachID)  + '.txt')                
        Masks_bin[:, :, j-1] = Mask_bin_currentSubb
    Pmean_casc=[]        
    for i in range(1, n_days_ev+1):
        print(f'masking fields, i={i}')
        Px = np.zeros((sizeNorth, sizeEast))
        if i>1: 
            Px = np.loadtxt('doct_rain/Px_day' + str(int(t_ev[i-1])) + '.txt') #Px of the day.            
#            for y in range(1, sizeNorth+1):
#                Px[y-1, :], res = io.read_float_arr(fileunit, path_list[pos_ini_var[i-1]-1], y, sizeEast)
        else:
            Px = Px0    #If i=0 Px had been already load to get size of spatial matrix.
#        print()
#        print()        
#        print('--------------NEW DAY-------------')
        Pmean_casc.append(int(np.mean(Px)))
#        print()
        for j in range(1, (n_reaches-1)+1):
            Pday_localSubb = Px*Masks_bin[:, :, j-1]
            sumP_subb = np.sum(Pday_localSubb)
            cells_subb = np.sum(Masks_bin[:, :, j-1])
#            print('i= ', i, ', j= ', j)
#            print('sumP_subb= ', sumP_subb, ', cells_subb= ', cells_subb)
            P_perDayperReach[i-1, j-1] = sumP_subb/cells_subb #np.mean(Pday_localSubb).
            Pday_subb[:, :, i-1] = Pday_subb[:, :, i-1] + Pday_localSubb
#            print('P_perDayperReach= ', P_perDayperReach[i-1, j-1])
#            print()
        print('day= ', t_ev[i-1], ' , Pmean_subbs[mm]= ', int(np.mean(P_perDayperReach[i-1, :])))
        print('day= ', t_ev[i-1], ' , Pmean_casc[mm]= ', Pmean_casc[i-1])
        printCumRumTime(t1)
    mp.mapPlotter(P_perDayperReach,  'doct_supply/PperSubb', 'P_tx', 
        'P(mm)', '# subbasin', 'time (days)', mn=1e0, log_cbar=True)
    np.savetxt('doct_supply/PperSubb/P_tx.txt', P_perDayperReach, fmt='%d')
#    last_days = 0
    lastDays2PlotMaskedField = 5
    for i in range(n_days_ev - lastDays2PlotMaskedField, n_days_ev+1):
        print(f'masking fields, i={i}')
        y = Pday_subb[:, :, i-1]
        isSubb = np.sum(Masks_bin, axis = 2)
        y[isSubb==0] = -10
#        if n_days_ev - i < 5: last_days = 1   #22dec2019: plot fields4only (effic) last 5 events, to review spatial peaks.
#        print(f'{n_days_ev}n_days_evi, i: {i}, last_days: {last_days}')
#        if last_days ==1: mp.mapPlotter(y, 'doct_supply/PperSubb', 'P_subb_day' + str(int(t_ev[i-1])))
        mp.mapPlotter(y, 'doct_supply/PperSubb', 'P_subb_day' + str(int(t_ev[i-1])),
            'P(mm)', 'easting', 'northing')
        printCumRumTime(t1)
