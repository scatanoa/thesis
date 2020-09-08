    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#jan2, 2020

import numpy as np
import sys

def refine_dt_flmFlood(nt_coarse, dt_ev, t_perQ): #9abr20: Qsteps NO equispaced in time.
    t_perQ = np.append(t_perQ, 1080*60)    #as per fig5 in pitlick13 I assume run 6.2 lasts 1080min.    
    dt_dQ_crs_arr = np.diff(t_perQ)
#    print(f'9abr20, dt_ev={dt_ev}, t_perQ[s]=\n{t_perQ}')
#    print(f'9abr20, len(dt_dQ_crs_arr) = {len(dt_dQ_crs_arr)}, dt_dQ_crs_arr = \n {dt_dQ_crs_arr}')
    n_copies_ev_flume = dt_dQ_crs_arr/dt_ev #float()
#    print(f'9abr20, nt_coarse={nt_coarse}, len(n_copies_ev_flume)={len(n_copies_ev_flume)}; \
#        n_copies_ev_flume=\n{n_copies_ev_flume}')
    coarse_dt=np.zeros((nt_coarse,2))  #to have relation [day,time_step*] *(at end of day) for plotting series with variable dt
    for i in range(1,nt_coarse+1):    #here i means time chunks with Qconstant in flume.
        coarse_dt[i-1,0] = i
        coarse_dt[i-1,1] = 1*(0 if i==1 else coarse_dt[i-2,1]) + n_copies_ev_flume[i-1]
#    print(f'9abr20, coarse_dt=\n{coarse_dt}')
#    nt = int(nt_coarse*n_copies_ev_flume)
    nt = int(coarse_dt[-1,1])
#    print(f'9abr20, nt=\n{nt}')
    dt_fine = np.array([dt_ev]*nt)
    n_copies_ev=0; n_copies_evt=0; n_copies_roe=0
    return dt_dQ_crs_arr, nt, dt_fine, n_copies_ev, n_copies_evt, n_copies_roe, coarse_dt
  
def refine_dt_flmTrench(nt_coarse, dt_ev, dt_dQ_coarse): #Qsteps equispaced in time.
    coarse_dt=np.zeros((nt_coarse,2))  #to have relation [day,time_step*] *(at end of day) for plotting series with variable dt
    n_copies_ev_flume = float(dt_dQ_coarse)/dt_ev
    for i in range(1,nt_coarse+1):    #here i means time chunks with Qconstant in flume.
        coarse_dt[i-1,0] = i
        coarse_dt[i-1,1] = 1*(0 if i==1 else coarse_dt[i-2,1]) + n_copies_ev_flume    
    nt = int(nt_coarse*n_copies_ev_flume)
    dt_fine = np.array([dt_ev]*nt)
    n_copies_ev=0; n_copies_evt=0; n_copies_roe=0
    return nt, dt_fine, n_copies_ev, n_copies_evt, n_copies_roe, coarse_dt

def refine_dt(nt_coarse, nR, flume_else_basin, Qhill, slid_sh_sum, dt_ev, dt_evt, dt_roe_prelim, dt_dQ_coarse, \
    tb_hs, t_mainstream, refine_pctile):
    #fine timesteps according to Courant criteria when high event (exced.prob~2%).
    if (flume_else_basin==1 and nt_coarse>1):
        Qref=.9*np.amin(Qhill)  #trick: set Qref such that code decide to produce finest dt always 
    elif (flume_else_basin==0):   
#        refine_pctile=1        #jan11, 2020: deprecated. Now it comes from list temp_resol of parameters dictionary.
        Qref=np.percentile(Qhill,100-refine_pctile)
    t_ev_spatial=[];      #t_slid=[]; xt_slid[]
    np.set_printoptions(precision=3, suppress=False)    
    print("Qref= ", Qref)
    slid_sh_mx = np.amax(slid_sh_sum, axis=0)
    print("slid_sh_mx= ", slid_sh_mx)
    for j in range(1, nR + 1): 
        t_ev,=np.where(Qhill[:,j-1] >= Qref)
        t_ev2,=np.where(slid_sh_sum[:,j-1] > 1e-9)   #jan10, 2020: I know min slid is 5e-10 (nD=5) when no one occurs.
#        t_slid.extend(t_ev2)
            #deprecated in jan10, 2020:   > slid_sh_sum_LocalMean[j-1])
        t_ev_spatial.extend(t_ev)
        t_ev_spatial.extend(t_ev2)                   
    t_ev_spatial = np.array(t_ev_spatial)
    t_ev_spatial = np.unique(t_ev_spatial)
    print("t_ev_spatial", t_ev_spatial)
    dt_fine=[]        

    tb_hs_mx = max(tb_hs)   #refine dt for subbasin with the most spread event hydrograph.
    n_copies_ev = np.around(tb_hs_mx/dt_ev);    n_copies_ev = n_copies_ev.astype(int)
    dt_ev_list= [dt_ev]*n_copies_ev    
    n_copies_evt = np.around((t_mainstream - tb_hs_mx) / dt_evt);  n_copies_evt = n_copies_evt.astype(int)
    dt_evt_list = [dt_evt]*n_copies_evt
        #evt: max event travel time through basin, along mainstream.
    n_copies_roe_float = (86400 - dt_ev*n_copies_ev - dt_evt*n_copies_evt)/dt_roe_prelim
    n_copies_roe = round(n_copies_roe_float, 0).astype(int)     #dt_roe is to be repaired, due to this rounding.

    dt_roe = (86400 - dt_ev*n_copies_ev - dt_evt*n_copies_evt)/n_copies_roe
    print(f'22mar, dt_roe{dt_roe}, n_copies_roe{n_copies_roe}')
    dt_roe_list = [dt_roe]*n_copies_roe
        #roe: restOfEventCoarseDt (e.g. rest of day, after event pulse traveled through basin).
    ss = f'tb_hs_mx/t_mainstream/dt_dQ_coarse must increase: {tb_hs_mx}/{t_mainstream}/{dt_dQ_coarse}'
    if (flume_else_basin == 0) and ( t_mainstream < tb_hs_mx or dt_dQ_coarse < t_mainstream): 
        sys.exit(ss)
    print(ss)
    print(f'The 3 dt resolutions in event day are: dte{dt_ev}, dtt{dt_evt}, dtr{dt_roe}')
    print(f'n copies in day event: n_e{n_copies_ev}, n_et{n_copies_evt}, n_roe{n_copies_roe}')
    print(f'must sum 86400, i.e. seconds/day: {dt_ev*n_copies_ev + dt_evt*n_copies_evt + dt_roe*n_copies_roe}')
    
    coarse_dt=np.zeros((nt_coarse,2))  #to have relation [day,time_step*] *(at end of day) for plotting series with variable dt
    n_copies_ev_flume = float(dt_dQ_coarse)/dt_ev #6may20; is this used, given flume has different refine_dt methods in this file?
    cont_fine = 0;  cont_coarse = 0     #cont_interm = 0; 
    for i in range(1,nt_coarse+1):    #here i means coarse time chunks (days except for trench experiment -hours)
        coarse_dt[i-1,0]=i
        any_reach_ev,=np.where(t_ev_spatial==i-1)
        if (i==1 or len(any_reach_ev>0)):   #1st day was refined anyway for easy debug of dt downscaling (ddt.f90)
            cont_fine += 1
#            i_last_ev=i
            dt_fine.extend(dt_ev_list);  dt_fine.extend(dt_evt_list);  dt_fine.extend(dt_roe_list)            
            coarse_dt[i-1,1] = 1*(0 if i==1 else coarse_dt[i-2,1]) + n_copies_ev + n_copies_evt + n_copies_roe
#        elif (i>1 and i==i_last_ev+1):
#            cont_interm += 1
#            dt_fine.extend(dtinterm_list) 
#            coarse_dt[i-1,1]=1*(0 if i==1 else coarse_dt[i-2,1])+len_dtinterm_list
        else: #current Qhill > Q
            cont_coarse += 1
            dt_fine.append(dt_dQ_coarse)
            coarse_dt[i-1,1]=1*(0 if i==1 else coarse_dt[i-2,1])+1
    pf = round(cont_fine/nt_coarse, 2); pc = round(cont_coarse/nt_coarse, 2)    #  ; pi = round(cont_interm/n, 2);
    print(f'p_days_event{pf}, p_days_noEv{pc}') 
    print('coarse_dt[day, t_step_endDay]= ', coarse_dt)
    dt_fine=np.array(dt_fine)
#        dt=dt.T    !no considered to test if fortran read this new array as argument (12Apr2018)
    print('6mar20, dtFromQ.py, dt_ev=', dt_ev, 'dt_dQ_coarse=', dt_dQ_coarse, 'nt_coarse=', nt_coarse)
    nt = len(dt_fine)
    print(f'in basin, must be 86400: {sum(dt_fine)/nt_coarse}')
    return nt, dt_fine, n_copies_ev, n_copies_evt, n_copies_roe, coarse_dt
