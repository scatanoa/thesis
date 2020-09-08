#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pandas as pd
import numpy as np
import sys

def get_df():
    #read names of variables:
    fp = 'doct_fluvial/hg_model_flume/calib_flood_pitlick13/'
    fn = 'pp13vars'
    f = open(fp + fn)
    vn1 = f.readlines()
    f.close()
    vn = [[],[]]
    vt = []
    for i, vi in enumerate(vn1):
        vl = vi.split(' ')
        vl[-1] = vl[-1].replace('\n','')
        vt.append(vl[0])
        vn[i] = vl[1:]
    #    print(f'santi, at i={i}, vl[1:]=\n{vl[1:]}')
    #    print(f'santi, vn={vn}')
    #    print(f'vt={vt}')

    #read values of variables:
    vt = [vti.replace(':','') for vti in vt]
    vals = [[],[]]
    for i,vti in enumerate(vt):
        x = np.loadtxt(fp + 'pp13_' + vti)
        vals[i] = x.T #.tolist()
    vals = np.array(vals)
    #print(f'vals={vals}')

    #join coarse and fine variables in pandas dataframe (like database).
    cdi = {}; fdi = {}
    dil = [fdi, cdi]
    for i,vni in enumerate(vn):
        dil[i] = {vn_ij : vals[i][j].tolist()  for j, vn_ij in enumerate(vni)}
    #    print(f'dil[{i}]= \n {dil[i]}')
    fdf = pd.DataFrame(dil[0]) #.from_dict
    fdf = fdf.replace(-99, np.nan)
    #print(f'fdf:\n{fdf}')
    #print(f'cdf:\n{cdf}')
    cdf = pd.DataFrame(dil[1]) #.from_dict
    cdf = cdf.replace(-99, np.nan)
    #print(f'cdf:\n{cdf}')
    vnc = vn[1] [1:]    #do not include 1st column (time) which is redundant in coarse and fine info.
    for i,vni in enumerate(vnc): fdf[vni] = np.nan
    #print(f'with new NaN columns, fdf:\n{fdf}')
    t_fdf = fdf['t_min']
    t_cdf = cdf['t_min']
    #print(f'vnc={vnc}')
    for tc in t_cdf:
    #    print(f'\n\ntc{tc}')
        for i, vnci in enumerate(vnc):
    #        print('i, vnci: ', i, vnci)
            newval = cdf.loc[t_cdf == tc, [str(vnci)]]
    #        print(f'df newval = \n {newval}')
            newval = newval.to_numpy()
    #        print(f'np newval = \n {newval} \n')        
            fdf.loc[t_fdf == tc, [vnci]] = newval
    #as I know sediment supply was constant.
    fdf['Qsi_half_gps'] = 3.0
#    print(f'final fdf:\n{fdf}')
#    sys.exit('29ago: paro pa ver pp13 data')

    return fdf
