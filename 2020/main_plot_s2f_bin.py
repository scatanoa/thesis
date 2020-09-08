#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np

import readParamFromDict as rpfd
import list_str2float as lsf
import plot_s2f_bin as psfb
import GSDini_model as gsd
import read_model_results as rmr
from sloper import *
import sys

def sci(var_list, nd):
    st = '{:0.' + str(nd) + 'e}'
    z = []
    for i, v in enumerate(var_list):
        vsc = st.format(v)
        z.append(vsc)
    return z

def recalcVars(v, h, vf, w, c, h_s, h_sf, w_c, D50c, D84c, D50c_ss, D84c_ss, w_v, sa, ta, flagSo_smooth, base_level, SoMin):
    #all coming variables were computed by f90, but were not written as binary output for PCspeed. For details, see the f90 code.
    
#    print(f'28apr20; just IN recalcVars, check w_v; wv_mean={np.mean(w_v)}, v_median={np.median(w_v)}')
#    print(f'28apr20; just IN recalcVars, check w_c; wc_mean={np.mean(w_c)}, wc_median={np.median(w_c)}')
#    print(f'28apr20; just IN recalcVars;  w_v \n: {w_v}')
    hbf = h_sf - h_s
    w_b = w_c - 2*hbf/ta
    ww = .5*(w_b+w_c)
    jup, jdown = sloper.get_neighbors(Rd, nR)
#    print(f'27apr20; in recalcVars, jup \n {jup}')
#    print(f'27apr20; in recalcVars, jdown \n {jdown}')
        #27apr20: New f90 module sloper was needed because, if original subroutine were to be used, 
        #then no py2f array lagger was available. 
        #Lagger is needed to traslate array indexes from C/py, that start at 0, to f90 form, that start at 1.
        #Note array delag from f to py is not needed. This seems to be done by something in f2py/numpy already.
    Ac = np.zeros((nt, nR)); Pc = np.zeros((nt, nR)); Rp = np.zeros((nt, nR))
    diagUnitTwic = 2*np.sqrt(1 + 1/ta**2)
#    print(f'28apr20; diagUnitTwic= {diagUnitTwic}')
#    print(f'28apr20; median&mean(hbf): {np.median(hbf)}, {np.mean(hbf)};  median&mean(h): {np.median(h)}, {np.mean(h)}')
    for i in range(1, nt): #28apr20: nt+1 -1, with -1 as f90 calcs.
        if i%1000==0: print(f'27apr20; method recalcVars, i{i}')
        for j in range(1, nR): #nR+1-1; -1 as f90 calc Fr(x=L)=1, where sediment is bypassed via rigid crosssection.
#            print('------------------------------------------')
#            print(f'27apr20; recalcVars; i{i}, j{j}, just before calling sloper.update_geom()')
#            print(f'{nR}, {j}, {jup[j-1, 0]}, {jdown[j-1]}, {base_level}, {flagSo_smooth}, {SoMin}, \
#                {Sor[j-1]}, {dx[j-1]}')
#            print(f'{h_s[i-1, :]}')
            h2 = h[i-1, j-1]; hbf2 = hbf[i-1, j-1]; vc2 = v[i-1, j-1] ; wc2 = w_c[i-1, j-1] ; ww2 = ww[i-1, j-1]
            wv2 = w_v[j-1]; w2 = w[i-1, j-1]; wb2 = w_b[i-1, j-1]
#            print('just in new i,j {i},{j}------------------------')
            So[i-1, j-1] = sloper.update_geom(j, jup[j-1, 0], jdown[j-1], base_level, flagSo_smooth,
                SoMin, Sor[j-1], dx[j-1], h_s[i-1, :], nR) #27apr20: j is passed without lag, as array lagger will occur.
                #subscript 2 means nothing; it is used to make the name of 0Dvariable different from the array variable. 
            hfld2 = h2-hbf2 if h2>hbf2 else 0
            if hfld2==0:
                Pc[i-1, j-1] = wb2 + h2*diagUnitTwic
                Ac2 = .5*(ww2+wb2)*h2
                Ac[i-1, j-1] = Ac2
                Q[i-1, j-1] = vc2*Ac2
                hh = Ac2/w2
                Fr[i-1, j-1] = vc2/np.sqrt(10*hh)
                aspect[i-1, j-1] = w2/hh
                asp_b_inv = h2/wb2
                Rp[i-1, j-1] = h2*(1 + asp_b_inv/ta) / (1 + 2*tbt/sa*asp_b_inv)
            else:
                Pc[i-1, j-1] = wb2 + hbf2*diagUnitTwic #note interface with floodplain flow in not accounted.
                Ac2 = (wc2 - hbf2/ta)*hbf2
                Ac[i-1, j-1] = Ac2
                Af = (wv2 - wc2)*hfld2
                Q2 = vc2*Ac2 + vf[i-1, j-1]*Af
                Q[i-1, j-1] = Q2
                A = Ac2+Af
                hh = A/wv2
#                print(f'28apr20; in flood; wv2={round(wv2, 3)}, wc2={round(wc2, 3)}, hfld2={round(hfld2, 3)}')
#                print(f'28apr20; in flood; Ac2={round(Ac2, 3)}, Af={round(Af, 3)}, A={round(A, 3)}')        
                Fr[i-1, j-1] = Q2/(A*np.sqrt(10*hh)) #27apr20; here Fr is for composite channel, while f90 calc is just for ch.
                aspect[i-1, j-1] = wv2/hh
                w_sf = wv2-wc2
                asp_b_inv = hbf2/wb2
                Rp[i-1, j-1] = ( hbf2*(1 + asp_b_inv/ta) + wv2/wb2*hfld2 ) / ( 
                    1 + 2*tbt*( hfld2/h2*(hfld2+w_sf)/wb2 + asp_b_inv/sa) )
        subm[:-1, :-1] = Ac[:-1, :-1]/Pc[:-1, :-1]/D84c[:-1, :-1]
        ts_c[:-1, :-1] = Rp[:-1, :-1]*So[:-1, :-1]/(1.65*D50c[:-1, :-1])
            #stress in channel bed. For details, see f90 code, where Rp was also computed.
    tc_c = calib*(1.32*So + .037) #eq 23 recking09 (D50)
    return Q, So, w_b, ts_c, tc_c, Fr, aspect, subm, hbf

pd = rpfd.readParam('param_dict.txt')
plot_series = int(pd['plot_series'])
which_flume = int(pd['which_flume'])    
flags = lsf.f(pd['flags'])
flume_else_basin = flags[0]
flagSo_smooth = flags[3]
flagSupplyWash = flags[1]; flagSupplySlide = flags[2]; 
flb = lsf.f(pd['fluv_basin'])
base_level = flb[1]
fl_fg = lsf.f(pd['fluv_geom']) 
SoMin = fl_fg[4-1]
dry_repos_ang = fl_fg[1-1]; drp_rad = 3.14/180*dry_repos_ang
ta = np.tan(drp_rad); sa = np.sin(drp_rad)
Lf = float(pd['L_flume'])
calib = float(pd['calib'])
dt_ind = lsf.f(pd['temp_resol'])[3]  if flume_else_basin==0 else  float(pd['dt_flume'])   #prl: pts_rising_limb
caseNum = int(pd['caseNum'])
fl_fwqsy = lsf.f(pd['fluv_wqsy'])
tbt = fl_fwqsy[1]
prt_ddt = lsf.f(pd['prt_rg_ddt'])

fn_top = 'Topol_X' if flume_else_basin ==0 else 'Topol_X_flume'
pthtop= f'doct_topol/{fn_top}.txt'      
X = np.loadtxt(pthtop); 
nR = len(X[:, 0])
w_v = X[:, 9-1]
Rd = X[:, 4-1]
dx = X[:, 7-1]
Sor = X[:, 6-1]
A = X[:, 5-1]

    #path builder to get data saved when producing binaries.        
caseNum = str(int(caseNum))
st_dt_ind = str(int(dt_ind))
if flume_else_basin==1:    
    pwe = 'hg_model_flume/' #path for which environment
    st_dx_flum = str(int(round(10*dx[0], 0)))
    if which_flume==1: #ubc trench
        pwf = 'trench/' #'path for which flume'
        case = '_10wv' + str(int(round(10*w_v[0],0))) + '_Lf' + str(int(Lf)) + '_10dx' + \
            st_dx_flum + '_dt' + st_dt_ind 
    else:
        pwf = 'flood/'
        case = '_10dx' + st_dx_flum + '_dt' + st_dt_ind
else:
    pwe = 'hg_model/'    
    pwf = ''
    case = '_case' + caseNum + '_dtPrl' + st_dt_ind #4jun20; always dx is reach lenght associated to Ahill.
        #cases intend to differ in supply regime, High&Pulsed(onlySlid), moderate&contin(onlyWash). p240u.
    #build here the name of sensitivity case with parameters that are varied from param_dict. #see phdnotes p172dr.
fo = 'doct_fluvial/' + pwe
ofpb = fo + pwf + 'bins' + case
#if not os.path.exists(ofpb): os.makedirs(ofpb)

fn= f'{ofpb}/dt_serie.txt'
dt_serie = np.loadtxt(fn)
fn= f'{ofpb}/cdt.txt'
cdt = np.loadtxt(fn)
fn= f'{ofpb}/nt_numChar.txt'
n = np.loadtxt(fn)
nt = int(n[0]); 
num_char = n[1]
fn= f'{ofpb}/pos_ini_var.txt'
pos_ini_var = np.loadtxt(fn, dtype='int')
fn= f'{ofpb}/ooVarNames_list.txt'
x = 'U' + str(int(num_char))
path_list_fnames = np.loadtxt(fn, dtype=(x))
path_list = []
num_paths = len(path_list_fnames)
for q in range (1, num_paths + 1):    
    fn = path_list_fnames[q-1]
    fn = str(fn) #from numpy.bytes_ to str
#    print(f'26apr20: path{q}, fn: \n {fn}')
    path_list.append(ofpb + '/' + fn)
#print(f'26apr20: path_list= \n {path_list}')
    
fl_fw = lsf.f(pd['fluv_wlogist'])

if ( flume_else_basin==0 or (flume_else_basin==1 and which_flume==1) ):
    nD = fl_fg[5]
    fl_gsd = lsf.f(pd['gsd'])
    fb, Dfbr, Dc = gsd.calc(nD, flume_else_basin, fl_gsd, which_flume)
else: #which_flume=2
    fb, Dfbr, Dc  = gsd.readFlumeFlood(which_flume, flume_else_basin)
#else?? basin..        

#initialize variables for plot method 'plt_hc' of class instance 'aa' created.

prt_wdw = lsf.f(pd['print_range'])

nD = len(Dc)
Qh = np.zeros((nt,nR)); Q = np.zeros((nt,nR)); wash_sh = np.zeros((nt,nR)); slid_sh = np.zeros((nt,nR))
aspect = np.zeros((nt,nR)); subm = np.zeros((nt,nR)); Mgtm = np.zeros((nt,nR)); gam = np.zeros((nt,nR))
v = np.zeros((nt,nR)); h = np.zeros((nt,nR)); vf = np.zeros((nt,nR)); w = np.zeros((nt,nR)); 
wproy = np.zeros((nt,nR)); c = np.zeros((nt,nR)); hbf = np.zeros((nt,nR))
So = np.zeros((nt,nR)); h_s = np.zeros((nt,nR)); h_sf = np.zeros((nt,nR)); w_c = np.zeros((nt,nR)); 
w_b = np.zeros((nt,nR)); ts_c = np.zeros((nt,nR)); tc_c = np.zeros((nt,nR)); tsm_c = np.zeros((nt,nR))
Fr = np.zeros((nt,nR)); D50c = np.zeros((nt,nR)); D84c = np.zeros((nt,nR)); D50c_ss = np.zeros((nt,nR))
D84c_ss = np.zeros((nt,nR)); pi = np.zeros((nt,nR, nD)); 
Qbf = np.zeros((nt,nR)); Qsbf = np.zeros((nt,nR));  Q_f90 = np.zeros((nt,nR))
volIn = np.zeros((nt,nR)); dvol = np.zeros((nt,nR))
colm = np.zeros((nt,nR))
fileunit=10
for i in range(1, nt + 1 -1): #28apr20: nt + 1 -1, with -1 as from f90 last time did not writte binary output.
    prt_cnd = 1 if (i%1000==0 or i>nt-2) else 0
    if prt_cnd==1: print('----------------')
    if prt_cnd==1: print(f'28apr20; reading binaries; i{i} of {nt+1}')
#    print(f'26apr20; i{i}')
    Qh[i-1,:],res = rmr.read_float_arr(fileunit, path_list[pos_ini_var[1-1] - 1], i, nR)
#    Q[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[2-1] - 1], i, nR)  
    if flagSupplyWash==1: #3jun20; otherwise array remains null when to be plotted
        wash_sh[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[2-1] - 1], i, nR)
    if flagSupplySlide==1:        
        slid_sh[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[3-1] - 1], i, nR) 
#    Mgtm[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[7-1] - 1], i, nR) 
#    gam[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[8-1] - 1], i, nR) 
    v[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[4-1] - 1], i, nR) 
    h[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[5-1] - 1], i, nR) 
    vf[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[6-1] - 1], i, nR) 
    w[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[7-1] - 1], i, nR)
    wproy[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[8-1] - 1], i, nR)     
    c[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[9-1] - 1], i, nR) 
    h_s[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[10-1] - 1], i, nR) 
    h_sf[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[11-1] - 1], i, nR) 
    w_c[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[12-1] - 1], i, nR) 
    D50c[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[13-1] - 1], i, nR) 
    D84c[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[14-1] - 1], i, nR) 
    D50c_ss[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[15-1] - 1], i, nR) 
    D84c_ss[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[16-1] - 1], i, nR)
    Qbf[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[17-1] - 1], i, nR)
    Qsbf[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[18-1] - 1], i, nR)
    Q_f90[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[19-1] - 1], i, nR)            
    So[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[20-1] - 1], i, nR)    
    w_b[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[21-1] - 1], i, nR) 
    ts_c[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[22-1] - 1], i, nR) 
    tc_c[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[23-1] - 1], i, nR) 
    Fr[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[24-1] - 1], i, nR)    
    aspect[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[25-1] - 1], i, nR) 
    subm[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[26-1] - 1], i, nR) 
    hbf[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[27-1] - 1], i, nR)    
    volIn[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[28-1] - 1], i, nR)
    dvol[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[29-1] - 1], i, nR)    
    colm[i-1,:],res= rmr.read_float_arr(fileunit, path_list[pos_ini_var[30-1] - 1], i, nR)
    if prt_cnd==1: print('   all binaries were read for this time i')
#    for k in range(1,nD+1):
#        pi[i-1,:, k-1],res= rmr.read_float_arr(fileunit, path_list[ pos_ini_var[27-1] - 1 + k - 1 ], i, nR)

print(f'26apr20; leyó todo binario a npArrays pa plotearlos.')
#Q, So, w_b, ts_c, tc_c, Fr, aspect, subm, hbf = recalcVars(v, h, vf, w, c, h_s, h_sf, w_c, D50c, D84c, D50c_ss, D84c_ss, w_v, 
#    sa, ta, flagSo_smooth, base_level, SoMin)      
#       #30may20, 11:20: variables were ready to be binsaved from f90; 
#           #recalc here was deprec due to Q errors mostly when incised.
         
    #26apr20: faster than binary R/W. See phdSCnotes p174r.
print('27apr20; just went out from recalcVars. It will enter plotter plt_hc')
print(f'27apr20; before, check varType1 v; v_mean={np.mean(v)}, v_median={np.median(v)}')

print(f'27apr20; before, check varType1 v; c_mean={np.mean(c)}, c_median={np.median(c)}')
print(f'15may20; before, check varType1 v; dt_serie_mean={np.mean(dt_serie)}, dt_serie_median={np.median(dt_serie)}')

print(f'27apr20; before, check varType2; Q_mean={np.mean(Q_f90)}, Q_median={np.median(Q_f90)}')

#print( '29may20; Qsbf_mean, Qsbf5, Qsbf50, Qsbf95',
#    np.mean(Qsbf), np.percentile(Qsbf,5), np.median(Qsbf),np.percentile(Qsbf,95) )
#print( '29may20; Qbf_mean, Qbf5-10-50-90-95',
#    np.mean(Qbf), np.percentile(Qbf,5),np.percentile(Qbf,10), np.median(Qbf),np.percentile(Qbf,90),
#    np.percentile(Qbf,95) )
#cbf = Qsbf/Qbf
#print( '29may20; cbf_mean, cbf5, cbf50, cbf95',
#    np.mean(cbf), np.percentile(cbf,5), np.median(cbf),np.percentile(cbf,95) )

print(f'27apr20; before, check out f2pySloper; Smean={np.mean(So)}, Smedian={np.median(So)}')
print(f'27apr20; before, check varType2; Fr_mean={np.mean(Fr)}, Fr_median={np.median(Fr)}')
print(f'27apr20; before, check varType2; wb_mean={np.mean(w_b)}, wb_median={np.median(w_b)}')
print(f'27apr20; before, check varType2; tsc_mean={np.mean(ts_c)}, tsc_median={np.median(ts_c)}')
print(f'27apr20; before, check varType2; tcc_mean={np.mean(tc_c)}, tcc_median={np.median(tc_c)}')
print(f'27apr20; before, check varType2; aspect_mean={np.mean(aspect)}, aspect_median={np.median(aspect)}')
print(f'27apr20; before, check varType2; subm_mean={np.mean(subm)}, subm_median={np.median(subm)}')
print(f'27apr20; before, check varType2; wproy_mean={np.mean(wproy)}, wproy_median={np.median(wproy)}')

#parker07 dimensionless variables
k1 = 1/np.sqrt(10*D50c**5)
k2 = 10**.2/Qbf**.4
Qbf_dl = Qbf*k1  #dl: dimensionless
Qsbf_dl = Qsbf*k1
wc_dl = w_c*k2
hbf_dl = hbf*k2


aa = psfb.c(plot_series, fo, fl_fw, Dc, fb, Lf, dx[0])


#30APR20: DDT PLOT ONLY FOR BASIN
#15may20: recovering ddt plot:
##plot downscaled supply:
plot_tseries = 0 #16may20
pad = 0 #29may20: 'plot also dailyTimeseries' active only if len(var) = 3 and var contains only supply variables.
    #if pad=1 to see effect of ddt downscaling, add vdkj arg to aa.plt_ddt() instantiation, in code block below.
if flume_else_basin==0:
    n_days = int(pd['dT_days'])
    ntc = n_days
    Vtx = np.zeros((ntc, nR, 3))
    fi = 'doct_supply/supply_wat_sed'
    if pad==0:
#        var = ['Qhill', 'wash', 'slid', 'Qbf', 'Q_f90', 'Q_py']
        var = ['Qhill', 'wash', 'slid', 'Qbf', 'Q_f90', 'c', 'w', 'v', 'D84c',
            'Qbf_dl', 'wbf_dl', 'hbf_dl', 'Qsbf_dl']
        dim = ['(m3/s)', '(m3/s)', '(m3/s)', '(m3/s)', '(m3/s)', '', '(m)', '(m/s)', '(m)',
            '','','','']
        vmn2 = [1e-3, 1e-3, 1e-3, 1e-3, 1e-2, 1e-7, 1e-1, 1e-2, 1e-3,  1e-2,1e-2,1e-2,1e-2] 
        plog = [1,1,1,1,1,1,0,1,1,  1,1,1,1]
#        needPlt = [1,1,1,1,1,1,1,1,1,  1,1,1,1]
        typ = [1,1,1,1,1,1,1,1,1,  0,1,1,1] #2jun20; type 2 was deprec, as hztal axis (scaling var) bust be static.
            #ref values for min tick in y axis, per variable        
#        vmx2 = [1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-7, 1e-1]
    else:
        var = ['Qhill', 'wash', 'slid']                    
    print('16may20; prev to load supply')
    for k2 in range(1, 3+1): #3 kinds of supply to loda from file.
#            print(f'6apr20; before supply change, k2={k2}, sar.shape[1]={sar.shape[1]}, nR_now={nR_now}')
#            print(f'6may20; will read this file: {fi}/{f[k2-1]}{suf}')
        sar = np.loadtxt(f'{fi}/{var[k2-1]}_tx.txt')
#            print(f'6apr20; k2={k2}, after supply change, sar.shape[1]={sar.shape[1]}, nR_now={nR_now}')
        Vtx[:, :, k2-1] = sar
    d2pi = int(prt_ddt[0]);  d2pf = int(prt_ddt[1]) #days to plot
    x2pi = int(prt_ddt[2]);  x2pf = int(prt_ddt[3])
    t_fine = np.cumsum(dt_serie) + 1
    pct = [3e-3, 5e-2, .25, .5, .75, .95] #17jun20; 1st value changed from 3e-5 to 3e-3.
    v_pct = np.zeros( (nR, len(var), len(pct)) )
    for j in range(x2pi, x2pf + 1):
        print('----------------------------------')                    
        suptt = f'Reach{j}';  print(suptt)
        Q_coarse = Vtx[:, :, 0];  wash_coarse = Vtx[:, :, 1];  slid_coarse = Vtx[:, :, 2]
        Qj = Qh[:, j-1];  wshj = wash_sh[:, j-1];  sj = slid_sh[:, j-1]      #downscaled at intradaily timescale.
        if pad==0: 
            Qbfj = Qbf[:, j-1]; Q_f90_j = Q_f90[:, j-1]; cj = c[:, j-1]; wj = w[:, j-1]  #; Q_py_j = Q[:, j-1]
            vj = v[:, j-1]; D84cj = D84c[:, j-1];
            Qbf_dl_j = Qbf_dl[:, j-1]; Qsbf_dl_j = Qsbf_dl[:, j-1]; 
            wc_dl_j = wc_dl[:, j-1]; hbf_dl_j = hbf_dl[:, j-1];
            
            vj = np.stack( (Qj, wshj, sj, Qbfj, Q_f90_j, cj, wj, vj, D84cj, 
                Qbf_dl_j, wc_dl_j, hbf_dl_j, Qsbf_dl_j), axis = 1 )     #, Q_py_j
        else:
            vj = np.stack((Qj, wj, sj), axis = 1)        
        #Qm = Qdj.mean() #;  ymn = 1e-3*Qm                
        print(f'mostrame Qhill en py tras f90: {Qj}')
    #                Qj = Qj/Qj.mean()*Qm;  wj = wj/wj.mean()*wdj.mean()
            # mean-preserving rescaling; don't needed by slides as they don't suffer triangular pulse, see ddt.f90.
        if pad==1:
            Qdj = Q_coarse[:, j-1];  wdj = wash_coarse[:, j-1];  sdj = slid_coarse[:, j-1]        
            vdj = np.stack((Qdj, wdj, sdj), axis = 1)
            #stack process comes by columns, which is axis 1. Then, each vector runs by rows.                
        for k, vr in enumerate(var): #plot each of the 3 supply variables: Q, wash and slide.
            vkj = vj[:, k];   
            l = [f'{var[k]}', f'{var[k]}_day', 'dt']
            print('----------------------------------')
            print(f'6mar20: estamos en j{j}, var: {vr}') #j{x2pi + j -1}
            print('6mar20: len(dt_serie)= ', len(dt_serie))
            vkjm = np.average(vkj, weights = dt_serie)
            vl = [min(vkj), vkjm, max(vkj)]
            z = sci(vl, 1)     #to print values list in scientific notation.
            if pad==1:
                vdkj = vdj[:, k]
                vld = [min(vdkj), np.mean(vdkj), max(vdkj)]                
                zd = sci(vld, 1)
                print(f'{vr}_MinMeanMax {z[0]}/{z[1]}/{z[2]}. \
                    DAILY, {vr}_MinMeanMax {zd[0]}/{zd[1]}/{zd[2]}')
                vkj_pct = aa.plt_ddt(flume_else_basin, plot_tseries, cdt, t_fine, vkj, dt_serie, l, d2pi, d2pf, 
                    suptt, j, pct, plog[k], dim[k], pad, vdkj)
            else:                
                vkj_pct = aa.plt_ddt(flume_else_basin, plot_tseries, cdt, t_fine, vkj, dt_serie, l, d2pi, d2pf, 
                    suptt, j, pct, plog[k], dim[k], pad)
            v_pct[j-1, k, :] = vkj_pct
#    sys.exit('31may20; all v_pct read?')
    nv = len(var)
    xx = np.zeros( ( nR, nv ) )
    Qbf_dl_p50 = v_pct[:, 9, 3]
    for i in range(1, nv+1):
        if typ[i-1]==1:
            xx[:, i-1] = A
        elif typ[i-1]==2:
            xx[:, i-1] = Qbf_dl_p50
else: #27ago2020; flume. just to fill arg required in plt_hc().
    xx=1; v_pct=1; var=1; vmn2=1; dim=1; typ=1; A=1; 

    #3sep2020; save time series for ch4phd plot with special code (to compare num simuls of both flume exps). See p280dl.
    if which_flume==1:
        flumeName = 'ubc'
        xl = [2, 16] #1m from upstream and downstream boundary condition.
    else:
        flumeName = 'fld'
        xl = [2, 23] #1m from upstream and downstream boundary condition.
    of = 'doct_fluvial/hg_model_flume/ch4phd_bothFlumes' 
    for j in xl:
        jj = j-1
        np.savetxt(f'{of}/{flumeName}_D84c_j{j}.txt', D84c[:,jj])
        np.savetxt(f'{of}/{flumeName}_h_s_j{j}.txt', h_s[:,jj])
        np.savetxt(f'{of}/{flumeName}_w_j{j}.txt', w[:,jj])
        np.savetxt(f'{of}/{flumeName}_Q_j{j}.txt', Q_f90[:,jj])
        np.savetxt(f'{of}/{flumeName}_c_j{j}.txt', c[:,jj])
        np.savetxt(f'{of}/{flumeName}_ts_c_j{j}.txt', ts_c[:,jj])
        np.savetxt(f'{of}/{flumeName}_tc_c_j{j}.txt', tc_c[:,jj])
    np.savetxt(f'{of}/{flumeName}_nt.txt', len(w[:, 0]), fmt='%i')
    sys.exit('3sep20; check if ok vars txt saved.')
aa.plt_hc(flume_else_basin, which_flume, case, flags[7], flags[4], cdt, prt_wdw, dt_serie, Q_f90, v, h, vf, w, wproy,
    c, D50c,  D84c, D50c_ss, D84c_ss, So, h_s, h_sf, w_c, w_b, w_v, ts_c, tc_c, Fr, aspect, subm, colm, Qbf_dl, 
    Qsbf_dl, hbf_dl, wc_dl, xx, v_pct, var, vmn2, dim, typ, A, volIn, dvol, wash_sh, slid_sh)
        #, Mgtm, gam, pi..., tsm_c
        #tc_c, tsm_c are maybe deprecated, as flood calibration of pitlick13exp was done with D50 instead D84.     
