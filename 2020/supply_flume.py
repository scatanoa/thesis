#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#5mar2020 
 
import numpy as np
import sys
import join_pp13exp_database as pp13

def fn(sufs, which_flume):
    X = np.loadtxt('doct_topol/Topol_X_flume.txt')
    n_nodes = len(X[:,0])
    if which_flume==1:
        Q2reach = [4, 6, 8, 10, 14, 18, 24, 28, 32, 32, 14, 8, 7, 6, 4] #l/s
        Q2reach = np.array(Q2reach)
        #5mar20: check len(Q2reach) = dT_days that is in param_dict.txt
    else:
        fd = pp13.get_df()
        dat = 'Q_half_lps'
        Q2reach = 2*fd[dat].to_numpy() #l/s for twice halfchannel
#        print(f'for floodFlume, Q2reach=\n{Q2reach}')

    fo = 'doct_supply/supply_wat_sed/'
    dims = (len(Q2reach), n_nodes) #(n_coarse_dt,n_nodes)
    Qhill = np.zeros(dims)
    slid_sh=np.zeros(dims)

    Q2reach = Q2reach*1e-3 #m3/s
    for j in range(1,n_nodes+1):
    #    if j==1:    #only water input in upstream boundary of flume #5mar20: same Q in all nodes, to avoid hydrology transient calcs. 
        Qhill[:,j-1]=Q2reach
    #    else:
    #        Qhill[:,j-1]=Q2reach/1e6 #to initialize water transit in flume, <1%err if 10 spatial nodes if /1e3
    if which_flume==1:
        c_feed = 0    #27ago20; null bro, I lost time calibrating unreal (no expUBC) feed!. 1e-3: deprec
        print(f'in supply2Fluvial.py; c_feed: {c_feed} ')
#        sys.exit()
            #5mar2020: 'concentration-like' units [m3sed/m3water]    
        slid_sh[:, 0] = c_feed*Q2reach
    else:        
        dat = 'Qsi_half_gps'
        slid_sh[:, 0] = 2*fd[dat].to_numpy()/2.65e6 #g/s to m3/s for twice halfchannel
        print(f'for floodFlume, slid_sh[0:5, 0]=\n{slid_sh[0:5, 0]}')
#        sys.exit('29ago20, slid_sh ok?')
    wash_sh=np.zeros(dims)    
    r = [Qhill, wash_sh, slid_sh]; rn = ['Qhill', 'wash', 'slid']
    for i, ri in enumerate(r): np.savetxt(fo + rn[i] + sufs[which_flume-1], ri)
    #    return 'holita' #Qhill, wash_sh, slid_sh
