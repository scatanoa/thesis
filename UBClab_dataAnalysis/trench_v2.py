    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#oct11, 2018

#pre-existing modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from osgeo import gdal
from osgeo import gdal_array
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import linregress
from scipy import stats
import math
#from pandas.tools.plotting import autocorrelation_plot
#from hurst import compute_Hc, random_walk

rcParams.update({'figure.autolayout': True})    #auto adjust plot size, to make axis labels visible

#functions: output is variable or plot

def fun(x, t, y):   #self
    return  y - x[0] + x[1]*(t+x[2])**2   

def get_dt(f,f_ini):
    t=np.zeros(2)
    ff=[f_ini,f]
    for i in range(1,2+1):
        f3=ff[i-1]
        t_aux=int(f3[5:7+1])
        t[i-1]=1*(t_aux if t_aux<100 else t_aux/100*60) 
    dt = 1*(t[1] - t[0] if f[1:2+1] == f_ini[1:2+1] else 15)  #same discharge?, [minutes]
    dt = 60*dt
    return dt   

def c_profile_Exner(Qs_outFlume, delta_t, delta_x, dz_h, dz_c, Wx, W_flume, steps_per_m):        
    poros_complem = 0.7
    Qsi_ant = 0
    Qs_x = []
    dst_c_x = []
    dst_h_x = []
#    print('reach averaged, dz_c= ', np.mean(dz_c))
    for jj in range(1,steps_per_m*9+1):
        Qso = 1*(Qsi_ant if jj>1 else Qs_outFlume)      #/W_sx_long[0]
        dst_ch = dz_c[jj-1]*Wx[jj-1]*delta_x
        dst_hs = dz_h[jj-1]*(W_flume-Wx[jj-1])*delta_x
        Qsi= Qso + poros_complem*(dst_hs + dst_ch)/delta_t #as dzh<0 always
        Qs_x.append(Qsi)
        dst_c_x.append(dst_ch)
        dst_h_x.append(dst_hs)
        Qsi_ant = Qsi
    dst_c_x.append(dst_c_x[-1])   
    Qs_x.append(Qs_x[-1])
    Qs_x = np.array(Qs_x)
    return Qs_x, dst_c_x, dst_h_x
    
#def c_profile_D_Exner(Qs_outFlume, delta_t, delta_x, dz_h, dz_c, Wx, steps_per_m, pi, nD):      
#NO SE PUDO PUES NO SE CONOCEN GRANULOMETRÍAS EN EL TIEMPO PARA EL ALUVIÓN (textura superficial vía photo de lecho sería engañosa, pues la idea es asumir que toodo el aluvión es capa activa)  
#    poros_complem = 0.7
#    Qsi_ant = 0
#    Qs_x = []
#    Qs_outFlume_frac = Qs_outFlume*pi
#    for d in range(nD):
#        for jj in range(1,steps_per_m*9+1):
##		                    F(i+1,j,k,kk) = F(i,j,k,kk) &
##		                    + dt(i)/za*dLa_dt(kk)*(F(i,j,k,kk) - fI_ok(k)) &
##                           + dt(i)/za/(1 - poros)*((Is_kk_k/dt(i) - Qs(i,j,k,kk))/(dx(j)*Wflow(j,kk)) &
##                             -  (1 - poros)*fI_ok(k)*dzdt(kk))        

#            Qso = 1*(Qsi_ant if jj>1 else Qs_outFlume_frac[p])      #/W_sx_long[0]
#            Qsi= Qso  +  poros_complem*delta_x/delta_t*(dz_c[jj-1]*Wx[jj-1] + dz_h[jj-1]*(1 - Wx[jj-1]))
#            Qs_x.append(Qsi)
#            Qsi_ant = Qsi
#    Qs_x.append(Qs_x[-1])
#    Qs_x = np.array(Qs_x)
#    return Qs_x_D

def lump_stor_and_Qs_per_stretch(Qsed, vol_c_x, Wx, dt, dx, n_str, steps_by_m):   
    #define number and lenght of stretches.
    W_flume= 1
    xf_list = [167, 284, 451-steps_by_m]
    xi_list = [steps_by_m+1, xf_list[0]+1, xf_list[1]+1]
    n_stretches=len(xf_list)
    if n_str!=n_stretches: 
        exit()
        print('stopped because of inconsistent number of stretches')
    #Lump Stores (LS) per stretches.
#    dst_lump = []    
    vol_c_lump = []
    Qs_lump = []
    for i in range(1, n_stretches + 1):
        xi = xi_list[i-1]
        xf = xf_list[i-1]
#        dzh_mean = np.mean(dz_h[xi-1 : xf+1])
#        L_stretch = dx*(xf - xi)
#        L_hill_mean = W_flume - np.mean(Wx[xi-1 : xf+1])
#        Vol_sed_from_hill = -dzh_mean*L_stretch*L_hill_mean
#        dst = Vol_sed_from_hill + dt*(Qsed[xf-1] - Qsed[xi-1])    #mass balance: Qs enters at xf and exits at xi. '+'=depos.
#        dst = np.sum(dst_x[xi-1 : xf-1]) 
#        dst_lump.append(dst)
        vol_c_lump.append(np.sum(vol_c_x[xi-1 : xf-1]))
        Qs_lump.append(Qsed[xi-1])
#    dst_lump = np.array(dst_lump)
    vol_c_lump = np.array(vol_c_lump)
    Qs_lump = np.array(Qs_lump)
    Qs_lump[Qs_lump < 0] = np.nan
    return Qs_lump, vol_c_lump
    
def plot_cx(output_sed_kg, dz_h, dz_c, c_sx, f, xv_long):
    path='/media/santiago/Datos/trench_code/'    
    tit = f + '_' + str(round(output_sed_kg, 3))  + 'kg out'    
    fig, ax1= plt.subplots()
    ax1.plot(xv_long, c_sx, '.k')
    ax2=ax1.twinx()
    ax2.plot(xv_long, dz_h, '.r', xv_long, dz_c, '.b')
    plt.title(tit)
    plt.xlabel('distance from downstream [m]')
    ax1.set_ylabel('c=Qs/Q [m3/m3]')
    ax2.set_ylabel('dz [m]')
    img_file = path + 'cx_' + f + '.png'
    plt.savefig(img_file)	#pcU: /media/santiago/Datos
    plt.clf() 
    plt.plot(-dz_h, c_sx, '.k') #note that erosion here is plotted as positive value, allowing log scale in axis.
    plt.title(tit)
    plt.xlabel('dzh [m]')
    plt.ylabel('c=Qs/Q [m3/m3]')
    plt.xscale('log')
    plt.yscale('log')    
    img_file = path + 'cx_vs_dzh_' + f + '.png'
    plt.savefig(img_file)	#pcU: /media/santiago/Datos
    plt.clf()
    #save data for Eric (UBC figures editor)
    data = np.vstack((xv_long, dz_h, dz_c, c_sx))
    filename = 'cx_' + f + '.txt'
    np.savetxt(filename, data)        
    
    
def plot_c_xt_per_run(c, f): #runs: A, B, C
    path='/media/santiago/Datos/trench_code/'    
    plt.imshow(c, cmap = 'gray', aspect = 'auto')
    plt.colorbar()
    tit = 'run' + f
    plt.title(tit)
    img_file=path + 'c_xt_' + tit + '.png'
    plt.savefig(img_file)	#pcU: /media/santiago/Datos
    plt.clf()

def plot_flume_stores(Q, cum_output, vol_flume, run):
    path='/media/santiago/Datos/trench_code/'    
    fig, ax = plt.subplots()
    n=len(vol_flume)   
    vol = np.array(vol_flume)
    poros_complem = 1.-.3
    vol_from_warming = vol*poros_complem*2650
    ax.plot(Q[0:n], -vol_from_warming,'-k', label= 'alluv.erosion')
    ax.set_xlabel('Q(m3/s)')
    ax.set_ylabel('kg')
    ax.plot(Q[0:n], cum_output, '+k', label= 'output')
    ax.plot(Q[0:n], -vol_from_warming - cum_output, label= 'errBalance')    #alluv change explains output
    leg = ax.legend(loc='best',fancybox=True, fontsize=10)
    leg.get_frame().set_alpha(0.3)    
    tit = 'run' + run
    plt.title(tit)
    img_file=path + 'singleStore_vs_Q_' + tit + '.png'
    plt.savefig(img_file)	#pcU: /media/santiago/Datos
    plt.clf()
    
def plot_flume_stores_h_c(Q, cum_output, vol_c_flume, vol_h_flume, run):
    path='/media/santiago/Datos/trench_code/'    
    fig, ax = plt.subplots()
    n=len(vol_c_flume)   
    vc = np.array(vol_c_flume)
    vh = np.array(vol_h_flume)
    poros_complem = 1.-.3
    vol_c_from_warming = vc*poros_complem*2650
    vol_h_from_warming = vh*poros_complem*2650    #- vol_h_flume[0]
    ax.plot(Q[0:n], vol_c_from_warming,'-k', label= 'channel store')
    ax.plot(Q[0:n], vol_h_from_warming,'--k', label= 'bank supply')
    ax.set_xlabel('Q(m3/s)')
    ax.set_ylabel('kg')
    ax.plot(Q[0:n], cum_output, '+k', label= 'output')
    ax.plot(Q[0:n], -vol_h_from_warming - vol_c_from_warming - cum_output, label= 'errBalance')    #alluv change explains output
    leg = ax.legend(loc='best',fancybox=True, fontsize=10)
    leg.get_frame().set_alpha(0.3)
    tit = 'run' + run
    plt.title(tit)
    img_file=path + 'store_vs_Q_' + tit + '.png'
    plt.savefig(img_file)	#pcU: /media/santiago/Datos
    plt.clf()  
    #save data for Eric (UBC figures editor)
    flume_stores_h_c_dat = np.vstack((vol_c_from_warming, vol_h_from_warming, cum_output))
    filename = 'store_vs_Q_run' + run + '.txt'
    np.savetxt(filename, flume_stores_h_c_dat)
    
    
    #plot slope 2 check
    
def plot_r2_c_vs_x_per_run(Q, name_exp_list, r2_t, r2_t_c_Fr, r2_t_aspect_Fr, f):
    path='/media/santiago/Datos/trench_code/'
    fig, ax = plt.subplots(nrows=1, ncols=3)
#    t = list(range(len(r2_t)))
    minutes = [15, 30, 60]
    col_per_run = {'A':'b', 'B':'r', 'C':'g'}
    ls = [':', '--', '-']
    name_exp_list = name_exp_list.T
    name_exp_list = name_exp_list.tolist()
    name_exp_list = name_exp_list[0]
    data = []    
    for i, mn in enumerate(minutes):        
        y = []
        y_c_Fr = []
        y_aspect_Fr = []
        Qlist = []
        for j, qq in enumerate(Q):
            for k, fn in enumerate(name_exp_list):
                t_str2num = int(fn[5:])
                indQ_str2num = int(fn[1:2+1])
                Q_current = Q[indQ_str2num-1]
                minu_f = 1*(t_str2num if t_str2num<46 else 60)
                if minu_f==mn and Q_current==qq:
                    val2add = r2_t[k]  
                    y.append(val2add[0])    #so scalar, instead 1-element, is appended. 
                    val2add = r2_t_c_Fr[k]  
                    y_c_Fr.append(val2add[0])
                    val2add = r2_t_aspect_Fr[k]  
                    y_aspect_Fr.append(val2add[0])
                    Qlist.append(qq)
        t_label = str(mn) + ' min' 
        Qlist = np.array(Qlist)                            
        ax[0].plot(1e3*Qlist, y, 'o',  markersize=3, linestyle = ls[i], color = col_per_run[f], label = t_label)
        ax[1].plot(1e3*Qlist, y_c_Fr, 'o',  markersize=3, linestyle = ls[i], color = col_per_run[f])
        ax[2].plot(1e3*Qlist, y_aspect_Fr, 'o',  markersize=3, linestyle = ls[i], color = col_per_run[f])
        #falta d dato en Qini pa exp B:
        bug = 0
        if f=='B' and len(y)<9:
            bug = 1
            nanValue = [np.nan]
            y = nanValue + y
            y_c_Fr = nanValue + y_c_Fr
            y_aspect_Fr = nanValue + y_aspect_Fr

        print('run, min ', f, mn)        
#        print('r2= ', r2_t, r2_t_c_Fr, r2_t_aspect_Fr)
        print('len y= ', len(y), len(y_c_Fr), len(y_aspect_Fr))
            
        print('y, prearray= ', y)
        print('y_c_Fr, prearray= ', y_c_Fr)
        print('y_aspect_Fr, prearray= ', y_aspect_Fr) 
        y = np.array(y)
        y_c_Fr = np.array(y_c_Fr)
        y_aspect_Fr = np.array(y_aspect_Fr)
        if bug==1: 
            print('toy dentro de bug')
            y = y.T
        print('y, POSarray= ', y)
        print('dims y:', y.shape)
        if len(data)==0:
            data = np.vstack((y.T, y_c_Fr.T, y_aspect_Fr.T))
            data = np.array(data) #to apply vstack the next time
        else:            
            data = np.vstack((data, y.T, y_c_Fr.T, y_aspect_Fr.T))
        print('data= ', data)
        print('')
    ax[0].set_xlabel('Q [m3/s]')
    ax[0].set_ylabel('r2_cx')
    leg = ax[0].legend()   #(loc=1,fancybox=True, fontsize=8)
#    leg.get_frame().set_alpha(0.3)    
    tit = 'run' + f
    plt.title(tit)
    img_file = path + 'r2_cx_' + tit + '.png'
    plt.savefig(img_file)	#pcU: /media/santiago/Datos
    plt.clf()
    #data4Eric UBC fig editor
    txtfilename = 'r2_x_' + tit + '.txt'
    np.savetxt(txtfilename, data)        
    
    
#    Qso_t = Qso_t*2.65e6
#    fig, ax = plt.subplots()    
#    ax.semilogx(Qso_t, r2_t,'.')
#    ax.set_xlabel('Qs [g/s]')
#    ax.set_xlim([1e-2, 1e2])
#    ax.set_ylabel('r2_cx')
#    tit = 'run' + f
#    plt.title(tit)
#    img_file = path + 'r2cx_vs_Qso_' + tit + '.png'
#    plt.savefig(img_file)	#pcU: /media/santiago/Datos
#    plt.clf()    

def plot_lumped_stores_per_run(lowQ_vars, f, singleBalanc):
    path='/media/santiago/Datos/trench_code/'
    store_lowQ = lowQ_vars[0]
#    store_highQ = highQ_vars[0]
    delta_store_lowQ = lowQ_vars[1]
#    delta_store_highQ = highQ_vars[1]
    Qs_lowQ = lowQ_vars[2]*2.65e6
#    print('Qs_lowQ= ', Qs_lowQ)
#    Qs_highQ = highQ_vars[2]*2.65e6
    np.where(Qs_lowQ<0, np.nan, Qs_lowQ)
#    np.where(Qs_highQ<0, np.nan, Qs_highQ)
    n_stretches = 3    
    col_per_run = {'A':'b', 'B':'r', 'C':'g'}
    ls = ['-', '--', ':']   #x started downstream when reach-averaging, then this list also follows such direction.
    lbl = ['downstream', 'middle', 'upstream']
    if singleBalanc==0:
        fig_name = ['Qs_vs_ChStore', 'deltaSt_vs_ChStore']
    else:         
        fig_name = ['Qs_vs_AlluvStore', 'deltaSt_vs_AlluvStore']
    for k in range(2): #plot store related to both Qs and delta_store
        fig, ax = plt.subplots()
        if k==0:    #k=0 plots Qs=f(stores)
            for i in range(1, n_stretches + 1):
                ax.semilogy(store_lowQ[:, i-1], Qs_lowQ[:, i-1], 'o',  markersize=3, linestyle = ls[i-1], \
                color = col_per_run[f], label = lbl[i-1])   #compare with PriorLisle.
            ax.set_ylabel('Qs [g/s]')
            ax.set_ylim(bottom = 1e-1, top = 1e2)            
        else: #k=1 plots delta_store=f(stores)
            for i in range(1, n_stretches + 1):
                ax.plot(store_lowQ[:, i-1], delta_store_lowQ[:, i-1], 'o',  markersize=3, linestyle = ls[i-1], \
                color = col_per_run[f], label = lbl[i-1])   #compare with ReidHassan19.
            ax.set_ylabel('d_st [m3]')                
            ax.set_ylim(bottom = 1*(-14e-3 if singleBalanc==0 else -18e-3), top = 3e-3)
        ax.set_xlabel('st [m3]')
        leg = ax.legend(loc='best',fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.3)   
#        ax.set_xlim([-.09, .01])
        tit = 'run' + f
        fig.suptitle(tit)
        img_file=path + fig_name[k] + tit + '.png'
        plt.savefig(img_file)	#pcU: /media/santiago/Datos
        txtfilename = fig_name[k] + '_' + tit + '.txt'
        data = np.vstack((store_lowQ.T, Qs_lowQ.T, delta_store_lowQ.T))
        np.savetxt(txtfilename, data)
        plt.clf()


class output:

    def __init__(self):
        self.steps_by_m = 50       

    def c_x_t_run_singleBalance(self):    #21june2019: deprecated because sepparating hill (dry alluvium) and channel (wet) is no worthy and lacks cumulBalance.
    #takes processed DEM diffs and sediment output in trap to calculate mass balance along channel.            
        singleBalanc = 1
        manualData_and_dem = np.genfromtxt('has_manualData_and_DEM.txt',dtype='str')
        kg_out_sed = np.genfromtxt('kg_out.txt')
        mdd_vs_systematicRuns = np.genfromtxt('mdd_vs_systematicRuns.txt')       
        pi = np.genfromtxt('md_pi_t_D.txt')        
        dzh = np.genfromtxt('pd_dzh_t_xFine.txt')
        dzc = np.genfromtxt('pd_dzc_t_xFine.txt')
        W = np.genfromtxt('pd_width_t_xFine.txt')
        d = np.genfromtxt('pd_depth_t_xFine.txt')
        v = np.genfromtxt('pd_veloc_t_xFine.txt')
        steps_by_m = self.steps_by_m
        W_flume = 1
        n_stretches = 3
        #c(x,t):
        dx = 1/steps_by_m  
        xv_long = np.linspace(0, 9, (9-0)*steps_by_m + 1)        
        Q = [2, 4, 6, 8, 10, 14, 18, 24, 28, 32]
        Q_transic = 40
        Q = np.array(Q)
        Q = 1e-3*Q
        indQ = 1           
        dimW = W.shape
        c_xt=[]
        output_given_Q_list = []
        Qrun = []       
        ind_end_exp = []
        vol_flume = []
        zc_ini = 0
        zh_ini = 0        
        n_x_fine = len(xv_long)
        zc = zc_ini*np.ones((n_x_fine))
        zh = zh_ini*np.ones((n_x_fine))
        st_lump = 0
        st_lump_prev = np.zeros((n_stretches))
        t_reset_new_exp = 0
        for i, fname in enumerate(manualData_and_dem): #81 instead 84
            i_aux = i - 1*(1 if fname[0]=='A' else 2 if fname[0]=='B' else 3)
            read = 0
            if i>0 and fname[0]!=fname_ini[0]:   # and fname[6:7+1]=='15':  #pa no plotear c2.100 sino desde c3.15
                zc = zc_ini*np.ones((n_x_fine))
                zh = zh_ini*np.ones((n_x_fine))
                changeExp= 1
                st_lump = 0
                c_xt = c_xt.T
                a = c_xt.shape
                ind_end_exp.append(a[1])
                plot_c_xt_per_run(c_xt, fname_ini[0])   
                lowQ_vars = [st_t_stretch_lowQ, dst_t_stretch_lowQ, Qs_t_stretch_lowQ]
#                highQ_vars = [st_t_stretch_highQ, dst_t_stretch_highQ, Qs_t_stretch_highQ]                 
                plot_lumped_stores_per_run(lowQ_vars, fname_ini[0], singleBalanc)
                plot_r2_c_vs_x_per_run(Q, name_exp_list, r2_t, r2_t_c_Fr, r2_t_aspect_Fr, fname_ini[0])
#                print('will call plot_flume_stores')
                output_given_Q_list = np.array(output_given_Q_list)
                plot_flume_stores(Qrun, np.cumsum(output_given_Q_list), vol_flume, fname_ini[0])
                vol_flume = []
                output_given_Q_list = []
                Qrun = []             
            if i>0 and fname[0]==fname_ini[0] and mdd_vs_systematicRuns[i]==1:   #2nd cond: same experiment either A, B or C                 
                read = 1
                dt = get_dt(fname, fname_ini)
                Qso = kg_out_sed[i]/2650/dt
                Q_changed = 1*(1 if (fname[1:2+1]!=fname_ini[1:2+1] or fname_ini[0:2+1]=='B02') else 0) #simplif: output no cumul in 1st Q at expB
                output_given_Q = kg_out_sed[i] + 1*(output_given_Q if Q_changed==0 else 0)
#                print('Q_changed= ', Q_changed)
#                print('output in this last measure (kg)', kg_out_sed[i])
#                print('output_given_Q (kg)= ', output_given_Q)
                Qs_x, dst_c_x, dst_h_x = c_profile_Exner(Qso, dt, dx, dzh[i_aux,:], dzc[i_aux,:], W[i_aux,:], W_flume, steps_by_m)   #Qs(x) per exp.
                indQ = int(fname[2])
                Q_current = Q[indQ-1]
                cx = Qs_x/Q_current
                slope, intercept, r_value, p_value, std_err = stats.linregress(xv_long[steps_by_m : 8*steps_by_m], cx[steps_by_m : 8*steps_by_m])
                r2 = r_value**2
                if fname=='C03_0015' or fname=='C08_0015':  #choosen examples to show longitudinal profile of concentration for paper2 (trenchExp,sed).
                    plot_cx(kg_out_sed[i], dzh[i_aux,:], dzc[i_aux,:], cx, fname, xv_long)
                d_now = d[i_aux, steps_by_m : 8*steps_by_m]
                v_now = v[i_aux, steps_by_m : 8*steps_by_m]
                W_now = W[i_aux, steps_by_m : 8*steps_by_m]
                cx_data = cx[steps_by_m : 8*steps_by_m]
                Fr = v_now / np.sqrt(10*d_now)
                aspect = W_now / d_now                                                    
                if fname=='B05_0030' or fname=='B08_0100': #23ago2019
                    x_data = xv_long[steps_by_m : 8*steps_by_m]                    
                    data_hc_sed = np.vstack((x_data, Fr, cx_data, aspect))
                    np.savetxt('hc_vs_c_' + fname, data_hc_sed)    #saved longitudinal arrays of x, Fr, c and W/d.
#                print(cx.shape)
                slope, intercept, r_value, p_value, std_err = stats.linregress(Fr, cx_data)
                r2_c_Fr = r_value**2
                slope, intercept, r_value, p_value, std_err = stats.linregress(Fr, aspect)
                r2_aspect_Fr = r_value**2
#                print('all r2 for this exp: ', r2, r2_c_Fr, r2_aspect_Fr)   #ok, aqui no hay redund, cada r2 da diferente.            
                Wh = W_flume - W[i_aux,:]
                zh = zh + dzh[i_aux,:]
                zc = zc + dzc[i_aux,:]
                zh[-50:] = 0    # do not include upstream ramp zone
                zc[-50:] = 0
#                if fname=='A05_0015':
#                    print('fname= ', fname)
#                    print('zc= ', zc)
                vol_x = dx*(np.multiply(zc, W[i_aux,:]) + np.multiply(zh, Wh))  # *hypoth: starts at 0
#                vol_c_x[-steps_by_m:] = 0   #null values in ramp (upstream) zone
                t_meas_isEnd= 1*(0 if fname[0:2+1]=='C03' else 1)
                if ((fname[5:7+1]=='100' and t_meas_isEnd==1) or (fname[6:7+1]=='30' and t_meas_isEnd==0)): 
                #as c3100 flowMask has high offset -> W ~6 instead .3m of manualMeasure. 
                    vol_flume.append(np.sum(vol_x))  #cumulate if last measurement of given discharge.
                    output_given_Q_list.append(output_given_Q)
                    Qrun.append(Q_current)
                Qs_lump, vol_c_lump = lump_stor_and_Qs_per_stretch(Qs_x, vol_x, W[i_aux,:], dt, dx, n_stretches, steps_by_m) # dst_lump, 
                st_lump_prev = st_lump
                st_lump = vol_c_lump    #B05_0045 adds to this cumulative balance, but its store data is not saved or 'stacked' for plot.
                dst_lump = st_lump - st_lump_prev            
                if (len(c_xt)<1 or changeExp==1):   #initial measurement of given experiment (A, B or C)                
                    c_xt = cx                    
                    st_t_stretch_lowQ = st_lump
                    dst_t_stretch_lowQ = dst_lump
                    Qs_t_stretch_lowQ = Qs_lump
                    r2_t = r2
                    r2_t_c_Fr = r2_c_Fr
                    r2_t_aspect_Fr = r2_aspect_Fr
                    name_exp_list = fname 
                    Qso_t = Qso 
                elif  fname!='B05_0045':                
                    c_xt = np.vstack((c_xt, cx))                
                    r2_t = np.vstack((r2_t, r2))
                    r2_t_c_Fr = np.vstack((r2_t_c_Fr, r2_c_Fr))
                    r2_t_aspect_Fr = np.vstack((r2_t_aspect_Fr, r2_aspect_Fr))
                    name_exp_list = np.vstack((name_exp_list, fname))
                    Qso_t = np.vstack((Qso_t, Qso))
                    if Q_current < 1e-3*Q_transic: #10: #'subtransitional' flow
                        st_t_stretch_lowQ = np.vstack((st_t_stretch_lowQ, st_lump))
                        dst_t_stretch_lowQ = np.vstack((dst_t_stretch_lowQ, dst_lump))
                        Qs_t_stretch_lowQ = np.vstack((Qs_t_stretch_lowQ, Qs_lump))                       
                    else:
                        transi_flag = 1*(1 if (Q_current==1e-3*Q_transic and fname[6:7+1]=='15') else 0)
                        if transi_flag==1: print('transi_flag= ', transi_flag)
                        st_t_stretch_highQ = 1*(st_lump if transi_flag==1 else np.vstack((st_t_stretch_highQ, st_lump)))
                        dst_t_stretch_highQ = 1*(dst_lump if transi_flag==1 else np.vstack((dst_t_stretch_highQ, dst_lump)))
                        Qs_t_stretch_highQ = 1*(Qs_lump if transi_flag==1 else np.vstack((Qs_t_stretch_highQ, Qs_lump)))                
                changeExp = 0
            fname_ini = fname
        c_xt = c_xt.T   #this way distance is plotted is vertical axis.
        a = c_xt.shape
        ind_end_exp.append(a[1])        
        plot_c_xt_per_run(c_xt, fname_ini[0])
        lowQ_vars = [st_t_stretch_lowQ, dst_t_stretch_lowQ, Qs_t_stretch_lowQ]
#        highQ_vars = [st_t_stretch_highQ, dst_t_stretch_highQ, Qs_t_stretch_highQ]
        plot_lumped_stores_per_run(lowQ_vars, fname_ini[0], singleBalanc)
        plot_r2_c_vs_x_per_run(Q, name_exp_list, r2_t, r2_t_c_Fr, r2_t_aspect_Fr, fname_ini[0])
        print('will call plot_flume_stores')        
        plot_flume_stores(Qrun, np.cumsum(output_given_Q_list), vol_flume, fname_ini[0])
    

    def c_x_t_run(self):    #21june2019: deprecated because sepparating hill (dry alluvium) and channel (wet) is no worthy and lacks cumulBalance.
    #takes processed DEM diffs and sediment output in trap to calculate mass balance along channel.            
        singleBalanc = 0
        manualData_and_dem = np.genfromtxt('has_manualData_and_DEM.txt',dtype='str')
        kg_out_sed = np.genfromtxt('kg_out.txt')
        mdd_vs_systematicRuns = np.genfromtxt('mdd_vs_systematicRuns.txt')       
        pi = np.genfromtxt('md_pi_t_D.txt')        
        dzh = np.genfromtxt('pd_dzh_t_xFine.txt')
        dzc = np.genfromtxt('pd_dzc_t_xFine.txt')
        W = np.genfromtxt('pd_width_t_xFine.txt')
        steps_by_m = self.steps_by_m
        W_flume = 1
        n_stretches = 3
        #c(x,t):
        dx = 1/steps_by_m  
        xv_long = np.linspace(0, 9, (9-0)*steps_by_m + 1)        
        Q = [2, 4, 6, 8, 10, 14, 18, 24, 28, 32]
        Q_transic = 40
        Q = np.array(Q)
        Q = 1e-3*Q
        zc_ini = 0
        zh_ini = 0
        indQ = 1           
        dimW = W.shape
        c_xt=[]
        vol_c = 0
        vol_h = 0
        vol_c_flume = []
        vol_h_flume = []
        vol_h_flume_all_t = []
        vol_c_flume_all_t = []        
        output_given_Q_list = []
        Qrun = []       
        ind_end_exp = []
        st_lump = 0
        st_lump_prev = np.zeros((n_stretches))
        n_x_fine = len(xv_long)
        zc = zc_ini*np.ones((n_x_fine))
        zh = zh_ini*np.ones((n_x_fine))
        t_reset_new_exp = 0
        for i, fname in enumerate(manualData_and_dem): #81 instead 84
            i_aux = i - 1*(1 if fname[0]=='A' else 2 if fname[0]=='B' else 3)
            print('')
            read = 0            
            if i>0 and fname[0]!=fname_ini[0]:   # and fname[6:7+1]=='15':  #pa no plotear c2.100 sino desde c3.15
                zc = zc_ini*np.ones((n_x_fine))
                zh = zh_ini*np.ones((n_x_fine))
                changeExp= 1
                st_lump = 0
                c_xt = c_xt.T
                a = c_xt.shape
                ind_end_exp.append(a[1])
                plot_c_xt_per_run(c_xt, fname_ini[0])   
                lowQ_vars = [st_t_stretch_lowQ, dst_t_stretch_lowQ, Qs_t_stretch_lowQ]
#                highQ_vars = [st_t_stretch_highQ, dst_t_stretch_highQ, Qs_t_stretch_highQ]                 
                plot_lumped_stores_per_run(lowQ_vars, fname_ini[0], singleBalanc)
#                plot_r2_c_vs_x_per_run(Q, name_exp_list, r2_t, fname_ini[0])
                print('will call plot_flume_stores_h_c')
                output_given_Q_list = np.array(output_given_Q_list)
                plot_flume_stores_h_c(Qrun, np.cumsum(output_given_Q_list), vol_c_flume, vol_h_flume, fname_ini[0])
                vol_c_flume = []
                vol_h_flume = []
                output_given_Q_list = []
                Qrun = []
                #save figs data for Eric UBC to edit             
                np.savetxt('stores_h_tseries_exp' + fname_ini[0], np.array(vol_h_flume_all_t)*(1-.3)*2650)
                np.savetxt('stores_c_tseries_exp' + fname_ini[0], np.array(vol_c_flume_all_t)*(1-.3)*2650)
                vol_h_flume_all_t = []
                vol_c_flume_all_t = []                
            if i>0 and fname[0]==fname_ini[0] and mdd_vs_systematicRuns[i]==1:   #2nd cond: same experiment either A, B or C                 
                read = 1
                dt = get_dt(fname, fname_ini)
                Qso = kg_out_sed[i]/2650/dt
                Q_changed = 1*(1 if (fname[1:2+1]!=fname_ini[1:2+1] or fname_ini[0:2+1]=='B02') else 0) #simplif: output no cumul in 1st Q at expB
                output_given_Q = kg_out_sed[i] + 1*(output_given_Q if Q_changed==0 else 0)
#                print('fname= ', fname)
#                print('Q_changed= ', Q_changed)
#                print('output in this last measure (kg)', kg_out_sed[i])
#                print('output_given_Q (kg)= ', output_given_Q)
                Qs_x, dst_c_x, dst_h_x = c_profile_Exner(Qso, dt, dx, dzh[i_aux,:], dzc[i_aux,:], W[i_aux,:], W_flume, steps_by_m)   #Qs(x) per exp.
                indQ = int(fname[2])
                Q_current = Q[indQ-1]
                cx = Qs_x/Q_current
                slope, intercept, r_value, p_value, std_err = stats.linregress(xv_long[steps_by_m : 8*steps_by_m], cx[steps_by_m : 8*steps_by_m])
                r2 = r_value**2
#                plot_cx(kg_out_sed[i], dzh[i_aux,:], dzc[i_aux,:], cx, fname, xv_long)
#                print(cx.shape)
                Wh = W_flume - W[i_aux,:]
                zc = zc + dzc[i_aux,:]
                zh = zh + dzh[i_aux,:]
                vol_c_x = dx*np.multiply(zc, W[i_aux,:])  # *hypoth: starts at 0
                vol_h_x = dx*np.multiply(zh, Wh)       #*
                vol_c_x[-steps_by_m:] = 0   #null values in ramp (upstream) zone
                vol_h_x[-steps_by_m:] = 0
                                
#                d_vol_c = 9*np.mean(np.multiply(dzc[i_aux,:], W[i_aux,:]))  #flume lenght is 9m
#                d_vol_h = 9*np.mean(np.multiply(dzh[i_aux,:], Wh))
                wh= np.mean(Wh[:-50])
                wc= np.mean(W[i_aux,:-50])
                dzcm = np.mean(dzc[i_aux,:])
                dzhm = np.mean(dzh[i_aux,:])
                d_vol_c = 9*wc*dzcm  #flume lenght is 9m
                d_vol_h = 9*wh*dzhm
                m3out = kg_out_sed[i]/2650
                err_relat_balanc = ((1-.3)*(-d_vol_h - d_vol_c) - m3out)/m3out                
                print(fname)
                print('kg_out= ', kg_out_sed[i])
                print('wh, wc= ', wh, wc)
                print('zhm, zcm ', np.mean(zh), np.mean(zc))
                print('dzhm, dzcm ', dzhm, dzcm)
                print('mean hill modif layer (m)', d_vol_h/(9*wh))   
                print('mean ch modif_layer (m)', d_vol_c/(9*wc))
                print('err_relat_balanc = ', err_relat_balanc)                                
                if fname=='A05_0015': 
                    print('zh.W= ', vol_h_x/dx)
                    print('zc.W= ', vol_c_x/dx)
                    fig, ax1 = plt.subplots()
                    ax2  = ax1.twinx()
                    ax2.plot(W[i_aux,:],'.k')
                    ax1.plot(dzh[i_aux,:],'.r')
                    ax1.plot(dzc[i_aux,:],'.b')
                    plt.show()                              
                
                t_meas_isEnd= 1*(0 if fname[0:2+1]=='C03' else 1)
                vol_h_aggr = np.sum(vol_h_x)
                vol_c_aggr = np.sum(vol_c_x)
                if ((fname[5:7+1]=='100' and t_meas_isEnd==1) or (fname[6:7+1]=='30' and t_meas_isEnd==0)): 
                #as c3100 flowMask has high offset -> W ~6 instead .3m of manualMeasure. 
#                    print('W_m= ', np.mean(W[i_aux,:]))                    
#                    print('zc_m= ', np.mean(zc))
                    vol_h_flume.append(np.sum(vol_h_x))  #cumulate if last measurement of given discharge.
                    vol_c_flume.append(np.sum(vol_c_x))  #cumulate if last measurement of given discharge.
#                    vol_h_flume.append(d_vol_h + 1*(0 if len(vol_h_flume)==0 else vol_h_flume[-1]))  #cum if last measure of given Q.
#                    vol_c_flume.append(d_vol_c + 1*(0 if len(vol_c_flume)==0 else vol_c_flume[-1]))
                    output_given_Q_list.append(output_given_Q)
                    Qrun.append(Q_current)
#                    print('vol_c_flume= ', np.mean(vol_c_flume))                    
                vol_h_flume_all_t.append(np.sum(vol_h_x))  #cumulate if last measurement of given discharge.
                vol_c_flume_all_t.append(np.sum(vol_c_x))  #cumulate if last measurement of given discharge.
                Qs_lump, vol_c_lump = lump_stor_and_Qs_per_stretch(Qs_x, vol_c_x, W[i_aux,:], dt, dx, n_stretches, steps_by_m) # dst_lump, 
                st_lump_prev = st_lump
                st_lump = vol_c_lump    #B05_0045 adds to this cumulative balance, but its store data is not saved or 'stacked' for plot.
                dst_lump = st_lump - st_lump_prev            
                if (len(c_xt)<1 or changeExp==1):   #initial measurement of given experiment (A, B or C)                
                    c_xt = cx                    
                    st_t_stretch_lowQ = st_lump
                    dst_t_stretch_lowQ = dst_lump
                    Qs_t_stretch_lowQ = Qs_lump
                    r2_t = r2
                    name_exp_list = fname 
                    Qso_t = Qso 
                elif  fname!='B05_0045':                
                    c_xt = np.vstack((c_xt, cx))                
                    r2_t = np.vstack((r2_t, r2))
                    name_exp_list = np.vstack((name_exp_list, fname))
                    Qso_t = np.vstack((Qso_t, Qso))
                    if Q_current < 1e-3*Q_transic: #10: #'subtransitional' flow
                        st_t_stretch_lowQ = np.vstack((st_t_stretch_lowQ, st_lump))
                        dst_t_stretch_lowQ = np.vstack((dst_t_stretch_lowQ, dst_lump))
                        Qs_t_stretch_lowQ = np.vstack((Qs_t_stretch_lowQ, Qs_lump))                       
                    else:
                        transi_flag = 1*(1 if (Q_current==1e-3*Q_transic and fname[6:7+1]=='15') else 0)
                        if transi_flag==1: print('transi_flag= ', transi_flag)
                        st_t_stretch_highQ = 1*(st_lump if transi_flag==1 else np.vstack((st_t_stretch_highQ, st_lump)))
                        dst_t_stretch_highQ = 1*(dst_lump if transi_flag==1 else np.vstack((dst_t_stretch_highQ, dst_lump)))
                        Qs_t_stretch_highQ = 1*(Qs_lump if transi_flag==1 else np.vstack((Qs_t_stretch_highQ, Qs_lump)))                
                changeExp = 0
            fname_ini = fname
        c_xt = c_xt.T   #this way distance is plotted is vertical axis.
        a = c_xt.shape
        ind_end_exp.append(a[1])        
        plot_c_xt_per_run(c_xt, fname_ini[0])
        lowQ_vars = [st_t_stretch_lowQ, dst_t_stretch_lowQ, Qs_t_stretch_lowQ]
#        highQ_vars = [st_t_stretch_highQ, dst_t_stretch_highQ, Qs_t_stretch_highQ]
        plot_lumped_stores_per_run(lowQ_vars, fname_ini[0], singleBalanc)
#        plot_r2_c_vs_x_per_run(Q, name_exp_list, r2_t, fname_ini[0])
        plot_flume_stores_h_c(Qrun, np.cumsum(output_given_Q_list), vol_c_flume, vol_h_flume, fname_ini[0])
        #save figs data for Eric UBC to edit             
        np.savetxt('stores_h_tseries_exp' + fname_ini[0] + '.txt', np.array(vol_h_flume_all_t)*(1-.3)*2650)
        np.savetxt('stores_c_tseries_exp' + fname_ini[0] + '.txt', np.array(vol_c_flume_all_t)*(1-.3)*2650)

class process:
    def __init__(self):
        dem_files=np.genfromtxt('filenames.txt',dtype='str')
        sheets=np.genfromtxt('sheets_manualData.txt',dtype='str')    #this is why following lines are deprecated       
        has_DEM=np.genfromtxt('has_manualData_and_DEM.txt',dtype='str')
        has_DEM_bin=np.genfromtxt('has_DEM_bin.txt',dtype='int')
        y_center=np.genfromtxt('ycenter_dem.txt')
        width=np.genfromtxt('w.txt')
        height=np.genfromtxt('h.txt')
        kg_out_sed=np.genfromtxt('kg_out.txt')
        manualData_and_dem=np.genfromtxt('has_manualData_and_DEM.txt',dtype='str')    
#        y_center.astype(int)   #does not work               
        self.dem_files=dem_files
        self.sheets=sheets
        self.has_DEM=has_DEM
        self.has_DEM_bin=has_DEM_bin
        self.y_center=y_center
        self.width=width
        self.height=height
        self.kg_out_sed=kg_out_sed
        self.manualData_and_dem=manualData_and_dem
        self.steps_by_m=50

    def which_manualData_and_dem(self):        
        #get list of sheets with manual data of experiments
        dem_files=self.dem_files
        xl=pd.ExcelFile('trench_flume_actAgo2018.xlsx')
        sheets=xl.sheet_names        
        len_sheets=len(sheets)
        sheets=sheets[9:len_sheets-10]  #manual data
        sheets=np.array(sheets)
        np.savetxt('sheets_manualData.txt',sheets,fmt="%s")  
        print('sheets= ', sheets)
        len_sheets=len(sheets)  #refreshed after filtering manual data
        print('len_sheets= ', len_sheets)        
        
        #evaluate which experiments has corresponding sfmDEM         
        has_DEM=np.zeros(len_sheets)
        has_DEM_complete=[]
        for i,sh in enumerate(sheets):
            print('')
            print(sh)
            tomatch=sh[4:]+'_dem.asc'
            print('tomatch= ', tomatch)
            ind,=np.where(dem_files==tomatch)            
            if len(ind)>0: 
                print('ind[0]= ', ind[0])
                has_DEM[i]=1
                has_DEM_complete.append(sh[4:])
        has_DEM_complete=np.array(has_DEM_complete)   
        np.savetxt('has_manualData_and_DEM.txt',has_DEM_complete,fmt='%s')   #%i
        np.savetxt('has_DEM_bin.txt',has_DEM,fmt='%i')   #%i              
    
    def get_ycenter_dem(self):
        dem_files=self.dem_files
        cum_y0=[]
        for i,fname in enumerate(dem_files):
            dem=np.loadtxt(fname,skiprows=6)
            dem[dem<0]=0    #because NaN from Agisoft is -32767!!            
            dem_long_mean=dem.mean(axis=1)    
            y0,=np.where(dem_long_mean==max(dem_long_mean))
            y0=min(y0)  #only to convert list to scalar
            print(fname)
            print(y0)
            cum_y0.append(y0)                    
        print('cum_y0=', cum_y0)
        np.savetxt('ycenter_dem.txt',cum_y0)            
    
    def jupyter(self):
        print('ey, me leíste desde web, uawww!!')
        plt.plot([1, 2, 3, 4])
        plt.ylabel('some few numbers')
        plt.show()
        
    def read_fractional_transport(self): 
        #marach13/2019: this method saves 84-6= 78 files. 6 are experiments not considered in sequence of flows #2 to #10, ...
        #...in which rows in xls file were set
        manualData_and_dem=self.manualData_and_dem    #abbreviation: mdd
        f=pd.ExcelFile('trench_flume_actAgo2018.xlsx')
        sh='fractionalTransport_okL2mm'
        run_names=['A','B','C']
        minutes=['015','030','045','100']
        cols=['G:O','Q:Y','AA:AI']
        i=0     
        mdd_vs_systematicRuns=np.zeros(len(manualData_and_dem)) #used by method c_xt to multiply Qs 1D array by pi 2D array.
        for run in range(1,3+1):    
            number_Q=2  #ordinal numbers for discharge, e.g. 1 for 2l/s, 2 for 4 l/s and so on.    
            measures_given_Q=0
            for row in range(1,36+1):     #36 is total amount of rows with fractional data in excel.   
                read=0
                i+=1
                measures_given_Q+=1
                zero_char='0'
                if number_Q>9: zero_char=''
                exp_name=run_names[run-1] + zero_char + str(number_Q) + '_0' + minutes[measures_given_Q-1]                
                ind,=np.where(manualData_and_dem==exp_name)
                if len(ind)>0:
                    mdd_vs_systematicRuns[ind[0]]=1                
                    read=1
                    row_read_df=pd.read_excel(f,sheet_name=sh,skiprows=list(range(85+row-1)),usecols=cols[run-1],nrows=1)
                    row_read_array=row_read_df.values
                    fracts_qsb=1*(row_read_array if i==1 else np.vstack((fracts_qsb,row_read_array)))                    
                if measures_given_Q==4: 
                    number_Q+=1   #after 4 quarter hours flow is changed.                    
                    measures_given_Q=0                    
                print('')                    
                print('i= ', i)
                print('exp_name= ', exp_name)
                print('read= ', read)
                print('check_sum1_pi_rowRead: sum= ', np.sum(row_read_array, axis=1)) #array has still floating numbers
#                print('row_read_df= ', row_read_df)                                                    
#                print('row_read_array= ', row_read_array)
                print('N_experiment_read= ', fracts_qsb.shape)
                print('fracts_qsb= ', np.around(fracts_qsb))            
        np.savetxt('md_pi_t_D.txt', fracts_qsb)  #remember that header 'md' and 'pd' refers to measured and processed data respectively. 
        np.savetxt('mdd_vs_systematicRuns.txt', mdd_vs_systematicRuns, fmt='%i')
        
    def read_xls_trench(self):
        """29oct2018: reads widths and Qs from xlsx and save each variable in txt file. This serves for coarse_balance"""
        dem_files=self.dem_files            
        sheets=self.sheets
        has_DEM_bin=self.has_DEM_bin
        f='trench_flume_actAgo2018.xlsx'                                        
        cum_w=[]
        cum_kg_out=[]
        cum_h=[]
        for i,sh in enumerate(sheets):
#            if i>4: break
            if has_DEM_bin[i]==1:
                print(sh)
#                xls=pd.read_excel(f,sheet_name=sh,skiprows=list(range(11)),usecols='A:M',nrows=15) #whole set of variables(x)
                w=pd.read_excel(f,sheet_name=sh,skiprows=list(range(11)),usecols='F',nrows=15) #widths                 
                w=w['width (cm)'].tolist()
                cum_w.append(w)
                h=pd.read_excel(f,sheet_name=sh,skiprows=list(range(11)),usecols='G',nrows=15) #water depth in middle of flume
                h=h['h (cm)'].tolist()
                cum_h.append(h)                
                kg_out=pd.read_excel(f,sheet_name=sh,skiprows=list(range(29)),usecols='A:B',nrows=2)
                kg_out=kg_out['Unnamed: 1'].tolist()
                cum_kg_out.append(kg_out[1])                          
        cum_w=np.array(cum_w) 
        cum_h=np.array(cum_h)       
        cum_kg_out=np.array(cum_kg_out)
        np.savetxt('w.txt',cum_w)        #,delimiter=','
        np.savetxt('h.txt',cum_h)
        np.savetxt('kg_out.txt',cum_kg_out)      
        print('cum_w= ',cum_w)   
        print('cum_h= ',cum_h)        
        print('cum_kg_out= ',cum_kg_out)
        
    def count_measures_per_experiment(self, mData_and_DEM):
#        manualData_and_dem=self.manualData_and_dem
        meas_per_exp=np.zeros((3))  # #A, #B, #C
#        print(manualData_and_dem[61][0:8])
        for i,fname in enumerate(mData_and_DEM):
            if fname[0]=='A': meas_per_exp[0]+=1
            if fname[0]=='B': meas_per_exp[1]+=1
            if fname[0]=='C': meas_per_exp[2]+=1
#            indA,=np.where(manualData_and_dem[:][0]=='A')
#            indB,=np.where(manualData_and_dem[:][0]=='B')
#            indC,=np.where(manualData_and_dem[:][0]=='C')                    
        return(meas_per_exp)
        
                                               
    def balance_tif_dem_difs(self):
        #30oct18: given I needed to calculate DEM differences manually in Qgis.
        dem_files=self.dem_files
        y_center=self.y_center
        width=self.width
        height=self.height
        kg_out_sed=self.kg_out_sed
        manualData_and_dem=self.manualData_and_dem
        meas_per_exp= self.count_measures_per_experiment(manualData_and_dem)
        print('meas_per_exp= ', meas_per_exp)
#        print('sum(kg_out_sed)= ', sum(kg_out_sed))

#        print(manualData_and_dem)
        path='/media/santiago/Datos/trench_code/'
        poros_complem=.7
        dx_ave=1.3e-3
        half_flume=.5   #[m]
        dy=round(1.05*half_flume/dx_ave)   
        xmin=500
        xmax=9000
        dz_max=.1
        param_fit=np.array([0,1/500,-5])
        Q=[2,4,6,8,10,14,18,24,28,32]
        Q=np.array(Q)
        Q=Q*1e-3
        denud_out=1e-10     #any value greater than 0        
        y_lsq_ch_long_calib=0
        y_lsq_hs_long=0
        y_lsq_ch_long_calib=np.array(y_lsq_ch_long_calib)
        y_lsq_hs_long=np.array(y_lsq_hs_long)
        lump_flume_storChange=[]
        d_ave_manual_per_exp=[]
        nCols_SOC_evid=45
            #cols: 0= Q, 1-3: r2_linRegress_c_vs_x, 4-6: Fr_ave, 7-9: Fr_pct95, 10-12: Fr_pct05, 13-15: cumMass_hill[kg], 16-18: cumMass_chan[kg],
            #19-21: dw_R2, 22-24: dw_slop, 25-27: vw_R2, 28-30: cFr_R2, 31-33: wdFr_R2, 34-36: wdFr_slop, 37-39: wdFr_intpto.
            #40-42: d_mx/d_av in ch and hill, 43-45: ind_t_relax.
        SOC_evid_A=np.zeros((8-1+1, nCols_SOC_evid))
        SOC_evid_B=np.zeros((10-2+1, nCols_SOC_evid))
        SOC_evid_C=np.zeros((9-3+1, nCols_SOC_evid))
        n_exp_givenQ=0
#        n_exp_to_pos=[0,1,0,2]  #indexes are possible values for n_exp_givenQ
        SOC_evid_row=0
        indQ=1
        changeExp=0
        n_exp_givenQ_list=[]
        steps_by_m=self.steps_by_m        
#        dz=[]   #to validate correction in DEM differences, distributed in hillslope and channel
#        x_w=np.linspace(1.5,8.5,15)
        for i,fname in enumerate(manualData_and_dem):
#            if i>9: break
            if i>0 and fname[0]!=fname_ini[0]: changeExp=1
            if i>0 and fname[0]==fname_ini[0]:   #last cond: same experiment either A, B or C
                name_dif=fname + '_minus_' + fname_ini + '.tif'
                indQant=indQ
                indQ=int(fname[2])  #10 Q values --> max_index=9<10. So, only 1 int is needed. When filename has '>9' are special (e.g. Q hystheresis) and unnecesary in this analysis.
                
                print(indQ)
                #load array of DEMchange
                f=path+name_dif
                dif=gdal_array.LoadFile(f)
                dif[abs(dif)>dz_max]=0
                yc=int(y_center[i])
                ymin=yc-dy
                ymax=yc+dy                
                dif=dif[ymin:ymax,xmin:xmax] 
                dif_sh=dif.shape          
                #clean extreme values in rigid walls of flume. Do not consider bricks (x>6000 pixels)
                for j in range(round(.95*dif_sh[0]),dif_sh[0]+1):
                    cell_ext,=np.where(dif[j-1,0:6000]>0)   #deposition 5cm near wall is not coherent.
                    if len(cell_ext)>0: dif[j-1,:]=0                                       
#                #clean extreme values in rigid walls of flume. Do not consider bricks (x>6000 pixels)
#                for j in range(1,dif_sh[0]+1):
#                    cell_ext,=np.where(abs(dif[j-1,0:6000])==dz_max) 
#                    if len(cell_ext)>0: dif[j-1,:]=0    

                t_aux=int(fname[5:7+1])
                t=1*(t_aux if t_aux<100 else t_aux/100*60) 
                t_aux=int(fname_ini[5:7+1]) 
                t_ini=1*(t_aux if t_aux<100 else t_aux/100*60)
                dt=1*(t-t_ini if fname[1:2+1]==fname_ini[1:2+1] else 15)  #same discharge?, [minutes]
                print('')  
                print('dt= ', dt)
                print(name_dif)
                #mask flow in channel:                
                dem_fin=np.loadtxt(fname+'_dem.asc', skiprows=6)                
                dem_fin_sh=dem_fin.shape
                mask_flow=np.zeros((dem_fin_sh[0],dem_fin_sh[1]))
                dem_fin[dem_fin<0]=0    #because NaN from Agisoft is -32767!!            
                dem_long_mean=dem_fin.mean(axis=1)    
#                y0,=np.where(dem_long_mean==max(dem_long_mean))
#                y0=min(y0)  #only to convert list to scalar 
#                print('y0= ', y0)
                dem_fin=10-dem_fin  #deeper near center of flume's Xsection
#                dem_fin[dem_fin==10]=0
                x0,=np.where(dem_fin.mean(axis=0)!=10)   
                x_ini=min(x0)   #[pixels]
                x_fin=int(x_ini+(9-1)/dx_ave)  #prismatic channel from abscissa 9m. Sfm from abscissa 1m  [pixels]
                print('x_ini= ', x_ini)
                print('x_fin= ', x_fin)                
                d=np.zeros((dem_fin_sh[0],dem_fin_sh[1]))
                for j3 in range(1,15+1): #each place where width was measured
                    ind_x_mn=x_ini+int((j3-1)/dx_ave*.5)  #measurements each .5m
                    ind_x_mx=x_ini+int(j3/dx_ave*.5)
                    zw=dem_fin[yc,ind_x_mx]+height[i,j3-1]/100  #[height]=cm
                    mask_flow[dem_fin<zw]=1
                    d=zw-dem_fin
                mask_flow=mask_flow[ymin:ymax,xmin:xmax]
                d=d[ymin:ymax,xmin:xmax]
                d[d<0]=0
                dif_ch=dif*mask_flow
                dif_hs=dif*(1-mask_flow)
                print('w= ', width[i,:])
                print('h= ', height[i,:])                                
#                print('t= ', t)                
#                print('t_ini= ', t_ini)
#                print('dt= ', dt)                
                #check mass balance per subsystem: hill or channel                
#                x0,=np.where(dif.mean(axis=0)!=0)   
#                x_ini=min(x0)   #[pixels]
#                x_fin=int(x_ini+(9-1)/dx_ave)  #prismatic channel from abscissa 9m. Sfm from abscissa 1m  [pixels]
                eros_layer_by_m_ch=[0]
                eros_layer_by_m_hs=[0]
                A_sx=[]
                W_sx=[]
                d_sx=[]
                d_ave_manual=[]
                x_manual=.5*np.linspace(1,13,(13-1)+1)  #x_manual needs ofset given DEM starts at x=1
                jj_manual=0
                v_sx=[]
                Fr_sx=[]
                Q_current=Q[indQ-1]
                for jj in range(1,steps_by_m*8+1):  #save erosion per each abscissa [m]: 1,2,..9m. each 1/5 m
                    ind_x_mn=x_ini+int((jj-1)/dx_ave/steps_by_m)
                    ind_x_mx=x_ini+int(jj/dx_ave/steps_by_m)
#                    if np.mean(dif_ch[:,ind_x_mn:ind_x_mx])==0:
#                        eros_layer_by_m_ch.append(0)
#                    else:
##                        eros_layer_by_m_ch.append(-np.mean(dif_ch[:,ind_x_mn:ind_x_mx]))
                    
#                    eros_layer_by_m_ch.append(np.average(dif_ch[:,ind_x_mn:ind_x_mx],weights=(dif_ch[:,ind_x_mn:ind_x_mx]>0))) #error
                    
                    dif_hs_wdw=dif_hs[:,ind_x_mn:ind_x_mx]
                    dif_ch_wdw=dif_ch[:,ind_x_mn:ind_x_mx]
                    mask_wdw=mask_flow[:,ind_x_mn:ind_x_mx]
                    eros_layer_by_m_ch.append(np.mean(dif_ch_wdw[mask_wdw>0]))          #dif_ch_wdw>0
                    eros_layer_by_m_hs.append(np.mean(dif_hs_wdw[mask_wdw==0]))    #dif_hs_wdw>0
                    
#                    if np.mean(dif_hs[:,ind_x_mn:ind_x_mx])==0:
#                        eros_layer_by_m_hs.append(0)
#                    else:
##                        eros_layer_by_m_hs.append(-np.mean(dif_hs[:,ind_x_mn:ind_x_mx]))

#                    eros_layer_by_m_hs.append(np.average(dif_hs[:,ind_x_mn:ind_x_mx],weights=(dif_hs[:,ind_x_mn:ind_x_mx]>0)))
                        
#                    eros_layer_by_m_hs.append(-np.mean(dif_hs[:,ind_x_mn:ind_x_mx]))  




                    d_jj=np.average(d[:,ind_x_mn:ind_x_mx],weights=(d[:,ind_x_mn:ind_x_mx]>0))
                    d_sx.append(d_jj)
                    if jj_manual < len(x_manual) and jj == x_manual[jj_manual]*steps_by_m:
                        print('d_ave_manual= ', d_ave_manual)
                        d_ave_manual.append(d_jj)   
                        jj_manual+=1                 
                        
                        
                    W_sx.append(sum(mask_flow[:,ind_x_mx])*dx_ave)  #one random width in strip                  
                    A_sx.append(d_sx[-1]*W_sx[-1])
                    v_sx.append(Q_current/A_sx[-1])
                    Fr_sx.append(v_sx[-1]/np.sqrt(10*d_sx[-1]))                    
                if i==1: print ('W_sx= ', W_sx)                    
#                print ('d_sx= ', d_sx)
#                print ('v_sx= ', v_sx)
#                print ('Fr_sx= ', Fr_sx)
                eros_layer_by_m_ch = np.nan_to_num(eros_layer_by_m_ch)
                eros_layer_by_m_hs = np.nan_to_num(eros_layer_by_m_hs)                    
                mean_eros_layer= np.mean(dif[:, x_ini:x_fin+1])
#                print('eros_layer_by_m= \n', eros_layer_by_m)
#                print('eros_layer_by_m_ch= ', eros_layer_by_m_ch)
#                print('eros_layer_by_m_hs= ', eros_layer_by_m_hs)
                print('mean_eros_layer= ', mean_eros_layer)
                m3_removed_flume= poros_complem * mean_eros_layer * 1 * 8 * 9/8 #9/8 extrapolating with mean to 1st meter with no photo
                print('m3_removed_flume= ', m3_removed_flume)
                m3_out_sed= kg_out_sed[i]/2650
                print('m3_out_sed= ', m3_out_sed)   #*60*dt
                rel_err_demDif= 1*(m3_removed_flume/m3_out_sed if m3_removed_flume>0 else -99)
                print('outputDEM/outputQs= ', rel_err_demDif)
                #plot mask flow of channel
#                plt.imshow(mask_flow)
#                plt.colorbar()
#                plt.title(name_dif[:-4])
#                img_file=path+'flowMask_'+name_dif[:-4]+'.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()
                #plot depth field
#                plt.imshow(d)
#                plt.colorbar()
#                tit='flowDepth[m]_' + name_dif[:-4]
#                plt.title(tit)
#                img_file=path+'depth_'+name_dif[:-4]+'.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()                  
                #plot 1D hydraulics
                extrap_ones=np.ones((steps_by_m))   #fill 1m downstream due to lack of DEM
                steps2extrap=3*steps_by_m                                
#                fig,ax=plt.subplots(nrows=4,ncols=1)     
                xv=np.linspace(0,9,(9-0)*steps_by_m+1)      
                W_sx.append(W_sx[0])
                W_sx=np.array(W_sx)
                d_sx.append(d_sx[0])
                d_sx=np.array(d_sx)                
                v_sx.append(v_sx[0])
                v_sx=np.array(v_sx)
                Fr_sx.append(Fr_sx[0])
                Fr_sx=np.array(Fr_sx)                                
                W_sx_extrapol=W_sx[:steps2extrap].mean()
                d_sx_extrapol=d_sx[:steps2extrap].mean()                
                v_sx_extrapol=v_sx[:steps2extrap].mean()
                Fr_sx_extrapol=Fr_sx[:steps2extrap].mean()                
                W_sx_extrapol=W_sx_extrapol*extrap_ones
                d_sx_extrapol=d_sx_extrapol*extrap_ones
                v_sx_extrapol=v_sx_extrapol*extrap_ones
                Fr_sx_extrapol=Fr_sx_extrapol*extrap_ones                
                W_sx_long=np.concatenate((W_sx_extrapol,W_sx))
                d_sx_long=np.concatenate((d_sx_extrapol,d_sx))
                v_sx_long=np.concatenate((v_sx_extrapol,v_sx))
                Fr_sx_long=np.concatenate((Fr_sx_extrapol,Fr_sx))
                indx_8half_m=int(len(W_sx_long)-.5*steps_by_m)
                W_mean=W_sx_long[0:indx_8half_m].mean()
                A_ch=W_mean*9.
                A_hs=9.*1-A_ch
                print('A_ch= ', A_ch)                
#                ax[0].plot(xv,W_sx_long,'.b',xv,W_mean*np.ones((len(xv))),'r')        #to solve bug and couple sizes of arrays to plot (err<1/400)
#                ax[0].set_ylim(top=np.percentile(W_sx,98))
#                ax[0].set_ylabel('W [m]')  
#                ax[1].plot(xv,d_sx_long,'.b',x_manual+1,d_ave_manual,'or')        
#                ax[1].set_ylim(top=np.percentile(d_sx,98))
#                ax[1].set_ylabel('d [m]')                 
#                ax[2].plot(xv,v_sx_long,'.b')        #to solve bug and couple sizes of arrays to plot (err<1/400)
#                ax[2].set_ylim(top=np.percentile(v_sx,98))
#                ax[2].set_ylabel('v [m/s]')                 
#                ax[3].plot(xv,Fr_sx_long,'.b')        #to solve bug and couple sizes of arrays to plot (err<1/400)
#                ax[3].set_ylim(top=np.percentile(Fr_sx,98))
#                ax[3].set_ylabel('Fr')
#                ax[3].set_xlabel('distance from downstream [m]')                                                
#                img_file=path+'hydrau1d_'+name_dif[:-4]+'.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()
                
                W_defi=1*(W_sx_long if i==1 else np.vstack((W_defi,W_sx_long)))
                print('W_defi_shape= ', W_defi.shape)        
                d_defi=1*(d_sx_long if i==1 else np.vstack((d_defi,d_sx_long)))
                print('d_defi_shape= ', d_defi.shape)                
                v_defi=1*(v_sx_long if i==1 else np.vstack((v_defi,v_sx_long)))
                print('v_defi_shape= ', v_defi.shape)                
                
                #plot spatial correlations of hydraulics to explore riffle-pool mechanics (as HassanShawnFerrer18). Nice secondary axis
                W_sx_pureData=W_sx_long[steps_by_m:8*steps_by_m]
                d_sx_pureData=d_sx_long[steps_by_m:8*steps_by_m]  
                v_sx_pureData=v_sx_long[steps_by_m:8*steps_by_m]                                
#                fig,ax1=plt.subplots()
                slope_wd, intercept, r_value, p_value, std_err = stats.linregress(W_sx_pureData,d_sx_pureData) 
                R2_dw=r_value**2                                
                y_linFit=intercept+slope_wd*W_sx_pureData
#                ax1.plot(W_sx_pureData, d_sx_pureData, '.b',label='d')
#                ax1.plot(W_sx_pureData,y_linFit,'-b',label='d_fit')
                slope, intercept, r_value, p_value, std_err = stats.linregress(W_sx_pureData,v_sx_pureData) 
                R2_vw=r_value**2
                y_linFit=intercept+slope*W_sx_pureData                
                tit='v and d, vs W. R2_dw='+str(round(R2_dw,2))+'. R2_vw='+str(round(R2_vw,2))                
#                ax1.set_title(tit)                
#                ax1.set_ylabel('d[m]')
#                ax1.set_xlabel('W[m]')                
#                ax2=ax1.twinx()
#                ax2.plot(W_sx_pureData, v_sx_pureData, '.r',label='v')
#                ax2.plot(W_sx_pureData,y_linFit,'-r',label='v_fit')                
#                ax2.set_ylabel('v[m/s]')
#                img_file=path+'hydrau1d_'+name_dif[:-4]+'correls.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()                    
                #plot eros_layer_by_m to define extrapolation:     
                xv=np.linspace(1,9,(9-1)*steps_by_m+1)                 
                lsq_hs=least_squares(fun, param_fit, loss='soft_l1', f_scale=1e-4, args=(xv,eros_layer_by_m_hs))   #loss & fscale<< ~omit outliers
                lsq_ch=least_squares(fun, param_fit, loss='soft_l1', f_scale=1e-4, args=(xv,eros_layer_by_m_ch))
                print('hs_optimal_parameters= ', lsq_hs.x)
                print('ch_optimal_parameters= ', lsq_ch.x)
                y_lsq_hs=fun(lsq_hs.x,xv,eros_layer_by_m_hs)
                y_lsq_ch=fun(lsq_ch.x,xv,eros_layer_by_m_ch)
#            fig,ax=plt.subplots(nrows=2,ncols=1)
                erosHill_by_m_trendRemoved=-1*(y_lsq_hs-eros_layer_by_m_hs)
#            ax[0].plot(xv,eros_layer_by_m_ch,'.b',xv,eros_layer_by_m_hs,'.r',xv,erosHill_by_m_trendRemoved,'-r',\
#                xv,-1*(y_lsq_ch-eros_layer_by_m_ch),'-b')
#                plt.plot(xv,eros_layer_by_m_ch,'.b',xv,eros_layer_by_m_hs,'.r')  ,xv,y_lsq,'-r')
                y_lsq_hs[y_lsq_hs>0]=0  #only survive negative (erosion) anomalies in hillslope
#                y_lsq_hs[abs(y_lsq_hs)>np.percentile(abs(y_lsq_hs),50)]=0
#                y_lsq_ch[abs(y_lsq_ch)>np.percentile(abs(y_lsq_ch),50)]=0
                length_unchanged=1      #[m]
                y_lsq_hs[-length_unchanged*steps_by_m:]=0    #0<x1 upstream was ~inactive in experiment. I'm removing false big DEM changes there.
                y_lsq_ch[-length_unchanged*steps_by_m:]=0
                y_lsq_hs_extrapol=y_lsq_hs[:steps2extrap].mean()
                y_lsq_hs_extrapol=y_lsq_hs_extrapol*extrap_ones
                y_lsq_ch_extrapol=y_lsq_ch[:steps2extrap].mean()
                y_lsq_ch_extrapol=y_lsq_ch_extrapol*extrap_ones
                y_lsq_hs_long=np.concatenate((y_lsq_hs_extrapol,y_lsq_hs))
                y_lsq_ch_long=np.concatenate((y_lsq_ch_extrapol,y_lsq_ch))
                y_lsq_ch_long_anom=y_lsq_ch_long-y_lsq_ch_long.mean()
                mean_eros_layer_hs=y_lsq_hs_long.mean()
                denud_out=m3_out_sed/(1*(9-length_unchanged))
                print('denud_out= ', denud_out)
#                mean_eros_layer_ch=denud_out/poros_complem-mean_eros_layer_hs      #are you crazy? balance is done by volume, not by layer!
                mean_eros_layer_ch=(-m3_out_sed/poros_complem-mean_eros_layer_hs*A_hs)/A_ch
                y_lsq_ch_long_calib=mean_eros_layer_ch+y_lsq_ch_long_anom
                print('mean_eros_layer_ch= ', mean_eros_layer_ch)
                print('mean_minus_dzChInX9= ', mean_eros_layer_ch-y_lsq_ch_long_calib[-1])
                print('upstr_dzCh= ', y_lsq_ch_long_calib[-length_unchanged*steps_by_m:])
                len_noZero=len(y_lsq_ch_long_calib)-length_unchanged*steps_by_m
                dilus_factor=(length_unchanged*steps_by_m)/len_noZero
                print('dilus_factor= ', dilus_factor)                
                
                y_lsq_ch_long_calib[:len_noZero]= y_lsq_ch_long_calib[:len_noZero] + \
                y_lsq_ch_long_calib[-length_unchanged*steps_by_m:].mean()*dilus_factor
                
                y_lsq_ch_long_calib[-length_unchanged*steps_by_m:]=0                
#                mean_amplif=len_long/(len_long-length_unchanged*steps_by_m)
#                print('mean_amplif= ', mean_amplif)
#                y_lsq_ch_long_calib=y_lsq_ch_long_calib*mean_amplif
#                err_relat_balanc=(poros_complem*(mean_eros_layer_ch+mean_eros_layer_hs)-denud_out)/denud_out
                y_lsq_ch_long_calib_mean=y_lsq_ch_long_calib.mean()
                volChang_hs=mean_eros_layer_hs*A_hs
                volChang_ch=y_lsq_ch_long_calib_mean*A_ch
                err_relat_balanc=(-poros_complem*(volChang_ch + volChang_hs) - m3_out_sed) / m3_out_sed
                lump_flume_storChange.append([volChang_hs, volChang_ch])    #hill,channel [m3]. If erosion then negative
                print('lump_flume_storChange= ', lump_flume_storChange)
                print('fname= ', fname)
                print('err_relat_balanc= ', err_relat_balanc)
                print('Wc= ', A_ch/9)
                print('Wh= ', A_hs/9)
                print('kg_out= ', kg_out_sed[i])
                print('dzhm, dzcm ', np.mean(y_lsq_hs_long), np.mean(y_lsq_ch_long_calib))
                
                
#                So_calib.append(So+)
                
                xv_long=np.linspace(0,9,(9-0)*steps_by_m+1)
#            ax[1].plot(xv_long,y_lsq_hs_long,'.r',label='hillslope')
#            ax[1].plot(xv_long,y_lsq_ch_long_calib,'.b',label='channel')  #stats.mode(y_lsq)[0][0]
#            ax[0].set_title(name_dif[:-4] + ', ' + str(round(kg_out_sed[i],3)) + 'kg out' + ', relErr=' + str(round(err_relat_balanc,2)))
#            img_file=path+'eros_by_m_'+name_dif[:-4]+'.png'                                                               
#            ax[1].set_xlabel('distance from downstream [m]')
#            ax[1].set_ylabel('depos_thickness [m]')
#            leg = ax[1].legend(loc=1,fancybox=True, fontsize=8)
#            leg.get_frame().set_alpha(0.3)                                  
#            plt.savefig(img_file)	#pcU: /media/santiago/Datos
#            plt.clf()  
                
                dzh=1*(y_lsq_hs_long if i==1 else np.vstack((dzh,y_lsq_hs_long)))
                print('dzh_shape= ', dzh.shape)
                dzc=1*(y_lsq_ch_long_calib if i==1 else np.vstack((dzc,y_lsq_ch_long_calib)))
                print('dzc_shape= ', dzc.shape)                
                                              
                #plot hydraulics vs dz
#                y_ch_calib_short=y_lsq_ch_long_calib[steps_by_m:]
#                plt.plot(Fr_sx_long[steps_by_m:8*steps_by_m],y_lsq_ch_long_calib[steps_by_m:8*steps_by_m],'.')
#                plt.xlabel('Fr')
#                plt.ylabel('channel_depos_thickness[m]')
#                img_file=path+'hc_vs_dz_'+name_dif[:-4]+'.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()        
                #plot 1d Exner
                dx_exner=1/steps_by_m
                c_sx=[]
                Qsi_ant=0
                dt=60*dt
                for jj in range(1,steps_by_m*9+1):
                    Qso=1*(Qsi_ant if jj>1 else m3_out_sed/dt)      #/W_sx_long[0]
                    print('jj= ', jj)
                    print('Qso= ', Qso)
#                    Qshill=poros_complem*((1-)*dx_exner/dt)  #'y'<0 for erosion. However exner balance implies qsi_hill>0 for slide.
#                    print('qshill= ', qshill)
#                    print('term_dz= ', term_dz)
                    Qsi= Qso + poros_complem*dx_exner/dt*(y_lsq_ch_long_calib[jj-1]*W_sx_long[jj-1] + y_lsq_hs_long[jj-1]*(1-W_sx_long[jj-1]))
                        #y_hs<0 if slide so inCh
                    Qsi_ant=Qsi
                    c=Qsi/Q_current
                    c_sx.append(c)
                c_sx.append(c_sx[-1])
                print('c_sx= ', c_sx)
                c_sx=np.array(c_sx)
                slope, intercept, r_value, p_value, std_err = stats.linregress(xv_long[steps_by_m:8*steps_by_m],c_sx[steps_by_m:8*steps_by_m]) 
                R2=r_value**2
                y_linFit=intercept+slope*xv_long
#                plt.plot(xv_long,c_sx,'.k',xv_long[steps_by_m:8*steps_by_m],y_linFit[steps_by_m:8*steps_by_m],'-r')
##                plt.yscale('log')   #avoid: neg values

##                c_sx_H=c_sx+abs(min(c_sx))+1e-10     #vertical translation of series such that >0, cos only variability matters for Hurst
##                H, c, data = compute_Hc(c_sx_H, kind='price', simplified=True)
##                tit='H=' + str(round(H,2)) + '_' + name_dif[:-4]
#                tit='R2=' + str(round(R2,2)) + '_' + name_dif[:-4]
#                plt.title(tit)
#                plt.xlabel('distance from downstream [m]')
#                plt.ylabel('c=Qs/Q [m3/m3]')
#                img_file=path+'exner1d_'+name_dif[:-4]+'.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()
                #plot autocorrelation_plot for longitudinal transport (nice plot, but not useful to show Qcr~ 10-14 l/s via memory in c(x))
                #so try Hurst few lines above, in same plot of c(x)
#                c_sx=np.array(c_sx)
#                c_sx=pd.Series(c_sx)
#                plt.figure()                                                                    
#                autocorrelation_plot(c_sx)
#                plt.title(name_dif[:-4])
#                img_file=path+'exner_cx_autocorrel_'+name_dif[:-4]+'.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()                    
                #plot improved DEMchange                 
#                plt.imshow(dif)
#                plt.colorbar()
#                plt.title(name_dif[:-4])
#                img_file=path+name_dif[:-4]+'.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()
                #plot improved DEMchange in channel
#                plt.imshow(dif_ch)
#                plt.colorbar()
#                plt.title(name_dif[:-4])
#                img_file=path+'channel_'+name_dif[:-4]+'.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()
                #plot improved DEMchange in hillslopes
#                plt.imshow(dif_hs)
#                plt.colorbar()
#                plt.title(name_dif[:-4])
#                img_file=path+'hillslope_'+name_dif[:-4]+'.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()
                #plot hydraulics vs transportExner
                Fr_sx_pureData=Fr_sx_long[steps_by_m:8*steps_by_m]
                c_sx_pureData=c_sx[steps_by_m:8*steps_by_m]
                w_d_ratio=W_sx_pureData/d_sx_pureData
#                fig,ax1=plt.subplots()
                slope, intercept, r_value, p_value, std_err = stats.linregress(Fr_sx_pureData,c_sx_pureData) 
                R2_c_Fr=r_value**2                                
                y_linFit=intercept+slope*Fr_sx_pureData                
#                ax1.plot(Fr_sx_pureData,c_sx_pureData,'.k')
#                ax1.plot(Fr_sx_pureData,y_linFit,'-k')
#                ax1.set_xlabel('Fr')
#                ax1.set_ylabel('c=Qs/Q [m3/m3]')                                                
#                ax2=ax1.twinx()
                slope_wdFr, intercept_wdFr, r_value, p_value, std_err = stats.linregress(Fr_sx_pureData,w_d_ratio) 
                R2_wd_Fr=r_value**2                                
                y_linFit=intercept_wdFr + slope_wdFr*Fr_sx_pureData                
#                ax2.plot(Fr_sx_pureData,w_d_ratio,'.g')
#                ax2.plot(Fr_sx_pureData,y_linFit,'-g')
#                ax2.set_ylabel('W/d')
#                tit='R2_c_Fr=' + str(round(R2_c_Fr,2)) + '. R2_wd_Fr=' + str(round(R2_wd_Fr,2))
#                ax1.set_title(tit)                                                
#                img_file=path+'hc_vs_c_'+name_dif[:-4]+'.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()                
                #plot calibrated DEMchange
                dif_sh=dif_ch.shape
                dif_ch_calib=np.zeros((dif_sh[0],dif_sh[1]))
                dif_hs_calib=np.zeros((dif_sh[0],dif_sh[1]))
#                dif_hs_anom=np.zeros((dif_sh[0],dif_sh[1]))
#                dif_ch_anom=np.zeros((dif_sh[0],dif_sh[1]))                
                for jj in range(1,steps_by_m*8+1):  #save erosion per each abscissa [m]: 1,2,..9m. each 1/5 m
                    ind_x_mn=x_ini+int((jj-1)/dx_ave/steps_by_m)
                    ind_x_mx=x_ini+int(jj/dx_ave/steps_by_m)
                    dif_hs_anom=dif_hs[:,ind_x_mn:ind_x_mx] - dif_hs[:,ind_x_mn:ind_x_mx].mean()  # -eros_layer_by_m_hs[jj-1]
                    dif_ch_anom=dif_ch[:,ind_x_mn:ind_x_mx] - dif_ch[:,ind_x_mn:ind_x_mx].mean()     #-eros_layer_by_m_ch[jj-1]
                    
#                    zz np.average(d[:,ind_x_mn:ind_x_mx],weights=(d[:,ind_x_mn:ind_x_mx]>0)))
                    
                    print('anom_mean_sum= ', dif_ch_anom.mean()+dif_hs_anom.mean())                    
                    dif_hs_calib[:,ind_x_mn:ind_x_mx]= dif_hs_anom - (-1*(y_lsq_hs[jj-1]-eros_layer_by_m_hs[jj-1]))     #substraction given y_lsq is positive for erosion
                    print('y_lsq_ch_long_calib.shape=, ', y_lsq_ch_long_calib.shape)
                    dif_ch_calib[:,ind_x_mn:ind_x_mx]= dif_ch_anom - (-1*(y_lsq_ch_long_calib[steps_by_m+jj-1]-eros_layer_by_m_ch[jj-1]))
#                fig,(ax0,ax1)=plt.subplots(nrows=2)
#                dif_hs_calib=dif_hs_calib*(1-mask_flow)
#                hs_plot=ax0.imshow(dif_hs_calib)
#                fig.colorbar(hs_plot,ax=ax0)
#                tit='deposThickness[m]_' + name_dif[:-4]    
#                fig.suptitle(tit)
#                dif_ch_calib=dif_ch_calib*mask_flow
#                ch_plot=ax1.imshow(dif_ch_calib)
#                fig.colorbar(ch_plot,ax=ax1)    
#                img_file=path+'calib_hs_ch_'+name_dif[:-4]+'.png'
#                plt.savefig(img_file)	#pcU: /media/santiago/Datos
#                plt.clf()
                
#                resname='pd_dzh2D_yF_xF_'+name_dif[:-4]+'.txt'
#                np.savetxt(resname,dif_hs_calib)
#                resname='pd_dzc2D_yF_xF_'+name_dif[:-4]+'.txt'
#                np.savetxt(resname,dif_ch_calib)            
                
                #save lumped results to make evidence of Qcrit (SOC?)
                SOC_evid_row+=1 if (indQ!=indQant) else 0
                min45=0
                if indQ!=indQant:
                    n_exp_givenQ=0 
                elif int(fname[5:7+1])!=45: #2nd cnd given expC starts with c3,15-c2,100
                    'mrk, hasta comenzando c me desplaza 1 la posic!'
                    n_exp_givenQ+=1                   
                    print('entiendo q no soy 45')
                elif int(fname[5:7+1])==45: min45=1
                if changeExp==1:                        
                    SOC_evid_row=0  #reset row counter if new exp
                    if int(fname[5:7+1])>15:   #because expC starts DEMdiffs with C3,15; which implies n_exp_givenQ must be 0 not 1
                        n_exp_givenQ=1  #because this 'for' loops initializes with 1st Q and 2nd time of measure, because it's about DEMdiffs
                changeExp=0    
#                n_exp_givenQ=n_exp_givenQ+1 if fname[0:3]==fname_ini[0:3] else 0    #n_exp_givenQ can be 0,1,2,3 for 15,30,45 or 60min.
                    #next line only if data for 15,30 and 60min is already stored.
                print(fname, '. i= ', i, '. n_exp_givenQ=', n_exp_givenQ)
                print('changeExp= ', changeExp)
                if (n_exp_givenQ<3 and min45==0):  # and n_exp_givenQ!=2  #Save neither extraExp nor measure of 45min because it is only available for ExpA.
                    Fr_wdw=Fr_sx_long[steps_by_m:8*steps_by_m]
                    Fr_ave=Fr_wdw.mean()
                    Fr95=np.percentile(Fr_wdw,95)
                    Fr05=np.percentile(Fr_wdw,5)
                    print('Fr_ave,95,05= ',Fr_ave,' ', Fr95,' ', Fr05)
                    if fname[0]=='A': 
                        SOC_evid_A[SOC_evid_row, 1+n_exp_givenQ]=R2
                        SOC_evid_A[SOC_evid_row, 0]=Q[indQ-1]
                        SOC_evid_A[SOC_evid_row, 4+n_exp_givenQ]=Fr_ave
                        SOC_evid_A[SOC_evid_row, 7+n_exp_givenQ]=Fr95
                        SOC_evid_A[SOC_evid_row, 10+n_exp_givenQ]=Fr05
                        
                        SOC_evid_A[SOC_evid_row, 19+n_exp_givenQ]=R2_dw
                        SOC_evid_A[SOC_evid_row, 22+n_exp_givenQ]=slope_wd
                        SOC_evid_A[SOC_evid_row, 25+n_exp_givenQ]=R2_vw
                        SOC_evid_A[SOC_evid_row, 28+n_exp_givenQ]=R2_c_Fr
                        SOC_evid_A[SOC_evid_row, 31+n_exp_givenQ]=R2_wd_Fr
                        SOC_evid_A[SOC_evid_row, 34+n_exp_givenQ]=slope_wdFr
                        SOC_evid_A[SOC_evid_row, 37+n_exp_givenQ]=intercept_wdFr
                    if fname[0]=='B': 
                        SOC_evid_B[SOC_evid_row, 1+n_exp_givenQ]=R2
                        SOC_evid_B[SOC_evid_row, 0]=Q[indQ-1]
                        SOC_evid_B[SOC_evid_row, 4+n_exp_givenQ]=Fr_ave
                        SOC_evid_B[SOC_evid_row, 7+n_exp_givenQ]=Fr95
                        SOC_evid_B[SOC_evid_row, 10+n_exp_givenQ]=Fr05
                        
                        SOC_evid_B[SOC_evid_row, 19+n_exp_givenQ]=R2_dw
                        SOC_evid_B[SOC_evid_row, 22+n_exp_givenQ]=slope_wd
                        SOC_evid_B[SOC_evid_row, 25+n_exp_givenQ]=R2_vw
                        SOC_evid_B[SOC_evid_row, 28+n_exp_givenQ]=R2_c_Fr
                        SOC_evid_B[SOC_evid_row, 31+n_exp_givenQ]=R2_wd_Fr
                        SOC_evid_B[SOC_evid_row, 34+n_exp_givenQ]=slope_wdFr
                        SOC_evid_B[SOC_evid_row, 37+n_exp_givenQ]=intercept_wdFr                                                
                    if fname[0]=='C': 
                        SOC_evid_C[SOC_evid_row, 1+n_exp_givenQ]=R2
                        SOC_evid_C[SOC_evid_row, 0]=Q[indQ-1]
                        SOC_evid_C[SOC_evid_row, 4+n_exp_givenQ]=Fr_ave
                        SOC_evid_C[SOC_evid_row, 7+n_exp_givenQ]=Fr95
                        SOC_evid_C[SOC_evid_row, 10+n_exp_givenQ]=Fr05
                        
                        SOC_evid_C[SOC_evid_row, 19+n_exp_givenQ]=R2_dw
                        SOC_evid_C[SOC_evid_row, 22+n_exp_givenQ]=slope_wd
                        SOC_evid_C[SOC_evid_row, 25+n_exp_givenQ]=R2_vw
                        SOC_evid_C[SOC_evid_row, 28+n_exp_givenQ]=R2_c_Fr
                        SOC_evid_C[SOC_evid_row, 31+n_exp_givenQ]=R2_wd_Fr
                        SOC_evid_C[SOC_evid_row, 34+n_exp_givenQ]=slope_wdFr
                        SOC_evid_C[SOC_evid_row, 37+n_exp_givenQ]=intercept_wdFr                         
                print('SOC_evid_row= ', SOC_evid_row)
                print('SOC_evid_A= ', SOC_evid_A)
                print('SOC_evid_B= ', SOC_evid_B)
                print('SOC_evid_C= ', SOC_evid_C)
                n_exp_givenQ_list.append(n_exp_givenQ)
                d_ave_manual_per_exp.append(d_ave_manual)
                print('d_ave_manual_per_exp= ', d_ave_manual_per_exp)                                                    
            fname_ini=fname
        np.savetxt('pd_dzh_t_xFine.txt',dzh)                    
        np.savetxt('pd_dzc_t_xFine.txt',dzc)
        np.savetxt('pd_depth_t_xFine.txt',d_defi)
        np.savetxt('pd_veloc_t_xFine.txt',v_defi)
        np.savetxt('pd_width_t_xFine.txt',W_defi)
        lump_flume_storChange=np.array(lump_flume_storChange)
        d_ave_manual_per_exp=np.array(d_ave_manual_per_exp)
        np.savetxt('lump_flume_storChange.txt',lump_flume_storChange)
        np.savetxt('d_ave_manual_per_exp.txt', d_ave_manual_per_exp)
        np.savetxt('SOC_evid_A.txt', SOC_evid_A)
        np.savetxt('SOC_evid_B.txt', SOC_evid_B)
        np.savetxt('SOC_evid_C.txt', SOC_evid_C)
        n_exp_givenQ_list=np.array(n_exp_givenQ_list)
        np.savetxt('n_exp_givenQ_list.txt', n_exp_givenQ_list)                
            
#            mean_eros_qs_ok=m3_out_sed
#            print('mean_eros_layer_qs_ok= ', mean_eros_qs_ok)
#            mean_eros_ch_ok=
#            print('mean_eros_layer_ch_ok= ', mean_eros_ch_ok)
#            mean_eros_hs_ok=y_lsq_hs_long.mean()
#            print('mean_eros_layer_hs_ok= ', mean_eros_hs_ok)
#            err_relat_balanc_ok = (poros_complem*(mean_eros_layer_hs_ok*A_hs + mean_eros_layer_ch_ok) - mean_eros_qs_ok) / mean_eros_qs_ok
#            print('RelErr_qs_toPrint= ', err_relat_balanc_ok)

#            if fname[0]==fname_ini[0]:  #same experiment

    def soc_evid(self):
        manualData_and_dem=self.manualData_and_dem
        n_exp_givenQ_list=np.genfromtxt('n_exp_givenQ_list.txt')
        SOC_evid_A=np.genfromtxt('SOC_evid_A.txt')
        SOC_evid_B=np.genfromtxt('SOC_evid_B.txt')
        SOC_evid_C=np.genfromtxt('SOC_evid_C.txt')
        cum_volChange_hs, cum_volChange_ch= self.balance_analyzer()
        ii=0
        SOC_evid_row=0
        indQ=1
        changeExp=0
        #organize sediment balance:
        for i,fname in enumerate(manualData_and_dem):
            if i>0 and fname[0]!=fname_ini[0]: changeExp=1
            if i>0 and fname[0]==fname_ini[0]:
                ii+=1
                indQant=indQ
                indQ=int(fname[2])
                SOC_evid_row+=1 if (indQ!=indQant) else 0
                if changeExp==1: SOC_evid_row=0
                changeExp=0
#                print('ii= ', ii)
#                print('n_exp_givenQ_list[ii-1]= ', n_exp_givenQ_list[ii-1])
#                print('int(fname[5:7+1])= ', int(fname[5:7+1]))
#                print('exp',fname[0])                
                if n_exp_givenQ_list[ii-1]<3 and int(fname[5:7+1])!=45:
                    ind=int(n_exp_givenQ_list[ii-1])
                    if fname[0]=='A': 
                        SOC_evid_A[SOC_evid_row, 13+ind]=cum_volChange_hs[ii-1]
                        SOC_evid_A[SOC_evid_row, 16+ind]=cum_volChange_ch[ii-1]
                    if fname[0]=='B': 
                        SOC_evid_B[SOC_evid_row, 13+ind]=cum_volChange_hs[ii-1]
                        SOC_evid_B[SOC_evid_row, 16+ind]=cum_volChange_ch[ii-1]
                    if fname[0]=='C': 
                        SOC_evid_C[SOC_evid_row, 13+ind]=cum_volChange_hs[ii-1]
                        SOC_evid_C[SOC_evid_row, 16+ind]=cum_volChange_ch[ii-1]                                                
            fname_ini=fname
        print('SOC_evid_A= ', SOC_evid_A)
        print('SOC_evid_B= ', SOC_evid_B)
        print('SOC_evid_C= ', SOC_evid_C)                      
        path='/media/santiago/Datos/trench_code/'        
        Qa=1e3*SOC_evid_A[:,0]
        Qb=1e3*SOC_evid_B[:,0]
        Qc=1e3*SOC_evid_C[:,0]
        evid_name=['R2_c_vs_x', 'Fr', 'cum_depos [kg]','d_vs_w','v_vs_w','c_vs_Fr','W_d_ratio_vs_Fr']
        stat_name=[':R2',':slope',':intercept']
        ylim_mx=[1,5]
        ylim_mn=[0,0]
        n_evid=len(evid_name)
        evid_pos=[1,4,13,19,25,28,31]
        evid_n_x_stats=[1,3,2,2,1,1,3]  #number of indicators of a given the evid_name (e.g. Fr has 3: Fr_ave, pct95 and pct05)
        lin_type_stats=['-','--',':']
        for i in range(1, n_evid+1):
            fig,ax=plt.subplots()
            print('i=',i)
            for j in range(1, evid_n_x_stats[i-1]+1):            
                evid_pos_ok=evid_pos[i-1] if evid_n_x_stats[i-1]==1 else evid_pos[i-1]+3*(j-1)
                ya=SOC_evid_A[:,evid_pos_ok:evid_pos_ok+2+1]
                yb=SOC_evid_B[:,evid_pos_ok:evid_pos_ok+2+1]    
                yc=SOC_evid_C[:,evid_pos_ok:evid_pos_ok+2+1]
                print('i,j:',i,j)
                print('evid_pos_ok= ', evid_pos_ok)                
                print('ya=',ya)
                if i!=3:
                    yam=np.mean(ya,axis=1)
                    ybm=np.mean(yb,axis=1)
                    ycm=np.mean(yc,axis=1)
                else:       
                    yam=ya[:,-1]  #get last value of given Q, cos kg_out is cumulative info   
                    ybm=yb[:,-1]
                    ycm=yc[:,-1]
                    print('yam=',yam)                    
                    print('ybm=',ybm)
                    print('ycm=',ycm)                    
                    yam=yam-yam[0]
                    ybm=ybm-ybm[0]
                    ycm=ycm-ycm[0]
                    print('yam_ok0=',yam)                    
                    print('ybm_ok0=',ybm)
                    print('ycm_ok0=',ycm)
                    if j==1:    #hillslope
                        dvol_hs_lump_a=yam
                        dvol_hs_lump_b=ybm
                        dvol_hs_lump_c=ycm
                    else:   #channel
                        dvol_ch_lump_a=yam
                        dvol_ch_lump_b=ybm
                        dvol_ch_lump_c=ycm                        
#                for k in range(1,nc+1):
                lab_a='expA'
                lab_b='expB'
                lab_c='expC'
                if i==2: 
                    lab_a= lab_a + ('pct5/95' if j>1 else 'average')
                    lab_b= lab_b + ('pct5/95' if j>1 else 'average')
                    lab_c= lab_c + ('pct5/95' if j>1 else 'average')
                if i==3: 
                    lab_a= lab_a + ('channel' if j>1 else 'hillslope')
                    lab_b= lab_b + ('channel' if j>1 else 'hillslope')
                    lab_c= lab_c + ('channel' if j>1 else 'hillslope')
                if i>3:
                    lab_a= lab_a + stat_name[j-1]
                    lab_b= lab_b + stat_name[j-1]
                    lab_c= lab_c + stat_name[j-1]
                if (i==7 or i==4) and j==1: #only 2dary axis for WdFr_R2 and dw_R2
                    ax2=ax.twinx()
                    ax2.plot(Qa[1:],yam[1:],lin_type_stats[j-1]+'b',label=lab_a)
                    ax2.plot(Qb[1:],ybm[1:],lin_type_stats[j-1]+'r',label=lab_b)
                    ax2.plot(Qc[1:],ycm[1:],lin_type_stats[j-1]+'g',label=lab_c)
                    tit2=evid_name[i-1]+'_R2'
                    ax2.set_ylabel(tit2)
                else:
                    ax.plot(Qa[1:],yam[1:],lin_type_stats[j-1]+'b',label=lab_a)
                    ax.plot(Qb[1:],ybm[1:],lin_type_stats[j-1]+'r',label=lab_b)
                    ax.plot(Qc[1:],ycm[1:],lin_type_stats[j-1]+'g',label=lab_c)
                    ax.set_ylabel(evid_name[i-1])
#                if i<n_evid: ax.set_ylim(bottom=ylim_mn[i-1], top=ylim_mx[i-1])
                ax.grid(which='major',axis='y')
                ax.set_xlabel('Q [l/s]')
#                if i<n_evid: ax.get_xaxis().set_visible(False)
                if i==2: ax.set_ylim(bottom=0,top=3.5) 
                leg = ax.legend(loc='best',fancybox=True, fontsize=7,frameon=True)
                leg.get_frame().set_alpha(0.3)            
            img_file=path+'SOC_evid_'+evid_name[i-1]+'.png'
            plt.savefig(img_file)	#pcU: /media/santiago/Datos
            plt.clf()
        #better way to show phase transition in sediment store in hill compared to channel:
        fig,ax=plt.subplots()
        ddvol_a=dvol_ch_lump_a-dvol_hs_lump_a
        ddvol_b=dvol_ch_lump_b-dvol_hs_lump_b
        ddvol_c=dvol_ch_lump_c-dvol_hs_lump_c
        ax.plot(Qa[1:],ddvol_a[1:],lin_type_stats[0]+'b',label='expA')
        ax.plot(Qb[1:],ddvol_b[1:],lin_type_stats[0]+'r',label='expB')
        ax.plot(Qc[1:],ddvol_c[1:],lin_type_stats[0]+'g',label='expC')        
        ax.set_xlabel('Q [l/s]')
        ax.set_ylabel('cum_depos_ch - cum_depos_hs [kg]')
        leg = ax.legend(loc='best',fancybox=True, fontsize=7,frameon=True)
        leg.get_frame().set_alpha(0.3)        
        img_file=path+'SOC_evid_'+evid_name[3-1]+'_difs.png'
        plt.savefig(img_file)	#pcU: /media/santiago/Datos
        plt.clf()        
        
                    

    def balance_analyzer(self):
        lump_flume_storChange=np.genfromtxt('lump_flume_storChange.txt')
        poros_complem=.7
        lump_flume_storChange=lump_flume_storChange*poros_complem*2650 #[kg_sed]
        dvol_hs=lump_flume_storChange[:,0]
        dvol_ch=lump_flume_storChange[:,1]
        plt.figure()
        plt.plot(dvol_hs,dvol_ch,'.')
        plt.xlabel('hillslope_VolChange[kg_sed]')
        plt.ylabel('channel_VolChange[kg_sed]')
        cum_volChange_hs=[]
        cum_volChange_ch=[]
        for i,dv_hs in enumerate(dvol_hs):
            if i==0: 
                cum_volChange_hs.append(dv_hs)
            else:                       
                cum_volChange_hs.append(cum_volChange_hs[-1]+dv_hs)
        cum_volChange_hs=np.array(cum_volChange_hs)
        for i,dv_ch in enumerate(dvol_ch):
            if i==0: 
                cum_volChange_ch.append(dv_ch)
            else:                       
                cum_volChange_ch.append(cum_volChange_ch[-1]+dv_ch)
        cum_volChange_ch=np.array(cum_volChange_ch)
        t=np.arange(len(cum_volChange_ch))        
        plt.figure()
        plt.plot(t,cum_volChange_hs,'.r',t,cum_volChange_ch,'.b')        
        return cum_volChange_hs, cum_volChange_ch
         
    def coarse_balance(self):
        dem_files=self.dem_files        
        print('how many DEMs?: ', len(dem_files))
        cum_geog_hdr=[]
        cum_shapes=[]
        cum_x0=[]   #initial abscissa of lightable
        cum_y0=[]
#        offset_y=0        
#        offset_x=0
        cum_x0_pix=[]
        cum_y0_pix=[]
        ibreak=4
        for i,fname in enumerate(dem_files):
            if i>ibreak: break

#            dx=.00125644
    #        fname='A01_0015_dem.asc'

    #        with open(filename,'r') as f:
    #           dem=f.readlines()    #all lines
#            print(fname)
            dem_db=gdal.Open(fname)
            print(fname)
            geog_hdr=dem_db.GetGeoTransform()
#            print(geog_hdr)

#            hdr=np.genfromtxt(fname,max_rows=6)
#            print('header=',hdr)
            
            dem=np.loadtxt(fname,skiprows=6)
            dem[dem<0]=0    #because NaN from Agisoft is -32767!!  
            dem_sh=dem.shape        
            ymin=400   #round(dem_sh[0]*.25)   #+45
            ymax=1700   #round(dem_sh[0]*.75)   #-45
            xmin=500
            xmax=8000
            dem=dem[ymin:ymax,xmin:xmax]
            dem_sh=dem.shape                                   

            cum_shapes.append([dem_sh[0],dem_sh[1]])
            x0=geog_hdr[0]
            y0=geog_hdr[3]
            dx=geog_hdr[1]            
            cum_geog_hdr.append([x0,y0,dx])   #upleft x0,y0,dx=dy
            
#            x0,=np.where(dem.mean(axis=0)>0)
#            x0=min(x0)
#            dem_long_mean=dem.mean(axis=1)    
#            y0,=np.where(dem_long_mean==max(dem_long_mean))
#            y0=min(y0)  #only to convert list to scalar

#            if i==0: 
#                x0_ref=x0
#                y0_ref=y0        
#            else:
#                offset_x=round((x0-x0_ref)/dx)
#                offset_y=round((y0-y0_ref)/dx)
#                print('x0,y0= ',x0,',',y0)                    
#                print('offset_x=',offset_x)                
#                print('offset_y=',offset_y)                
        
            x0_pix=round(x0/dx)  #[pixels]
            y0_pix=round(y0/dx)  #[pixels]
            cum_x0.append(x0)
            cum_x0_pix.append(x0_pix)
            cum_y0_pix.append(y0_pix)            
        x0_ref=max(cum_x0_pix)  #upleft corner: left
        y0_ref=max(cum_y0_pix)  #upleft corner: up
        print('x0_ref= ', x0_ref)
        print('y0_ref= ', y0_ref)        
        
        
            
        for i,fname in enumerate(dem_files):                                
            if i>ibreak: break      
            offset_x=x0_ref-cum_x0_pix[i]
            offset_y=y0_ref-cum_y0_pix[i]
            print(fname)            
            print('offset_x= ', offset_x)
            print('offset_y= ', offset_y)            

            dem_trasl=np.zeros((dem_sh[0],dem_sh[1]))
            for x in range(offset_x,dem_sh[1]-offset_x):
                print('x=', x)
                dem_trasl[:,x+offset_x]=dem[:,x]
            for y in range(offset_y,dem_sh[0]-offset_y):
                print('y=', y)
                dem_trasl[y+offset_y,:]=dem[y,:]
            dem=dem_trasl
           
            plt.imshow(dem)
            plt.colorbar()
            img_file='/media/santiago/Datos/trench_code/'+fname[:-4]+'.png'
            plt.savefig(img_file)	#pcU: /media/santiago/Datos
            plt.clf()            
            
            if i==0: 
                dem_ini=dem                                      
            if i>0: 
                dem_dif=dem-dem_ini
                dem_dif[abs(dem_dif)>.2]=.2
                plt.imshow(dem_dif)
                plt.colorbar()                
                img_file='/media/santiago/Datos/trench_code/' + fname[:-8] + '_minus_' + fname_ini[:-8] + '.png'
                plt.savefig(img_file)
                plt.clf()
                dem_ini=dem 
            fname_ini=fname                      
            
        cum_geog_hdr=np.array(cum_geog_hdr)
        cum_x0=np.array(x0)
#        cum_offset_x=np.array(cum_offset_x)
#        cum_offset_y=np.array(cum_offset_y)        
        print('cum_geog_hdr',cum_geog_hdr) 
        print('cum_shapes',cum_shapes)        
        print('cum_x0',cum_x0)
#        print('cum_offset_x',cum_offset_x)
#        print('cum_offset_y',cum_offset_y)                    
            
      
      
      
      
      
      
      
      
      

##            #read DEM
##            
##            #clip DEM to leave only alluvium cells
##        
##            #read

