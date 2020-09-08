#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import os
import os.path as op

import FolderCleaner as fc #25apr20: forbidden bro; files don't  go even to trash bin, so this tool is dangerous.
    #to recover files:
    #https://itsfoss.com/recover-deleted-files-linux/
import readParamFromDict as rpfd
import list_str2float as lsf

import var_paths as vp
import dtFromQ as dtfq
import plot_s2f as ps2f
import supply2fluvial_model as s2fm
import GSDini_model as gsd
import main_topol as mt
import supply_flume as sf
import join_pp13exp_database as pp13

class a:
    def __init__(self):
        pd = rpfd.readParam('param_dict.txt')
        flags = lsf.f(pd['flags']);   prt_wdw = lsf.f(pd['print_range']); prt_ddt = lsf.f(pd['prt_rg_ddt'])
        self.fl_fw = lsf.f(pd['fluv_wlogist'])
        self.fl_fwqsy = lsf.f(pd['fluv_wqsy'])
        self.HGelseFract = flags[6]
        self.flume_else_basin = flags[0]
        
        if self.flume_else_basin ==0:
            tr = lsf.f(pd['temp_resol'])
            self.dt_dQ_coarse = int(tr[0])
        else:
            self.dt_dQ_coarse = pd['dt_dQ_coarse_flume']
        ti = time.time()
        self.runFluvialModel(prt_wdw, prt_ddt, flags, pd)
        tf = time.time();  t = self.sci([tf-ti], 1)[0]
        print('----------------------------')
        print(f'in supply2Fluvial.py, method runFluvialModel lasts {t}s')

    def sci(self, var_list, nd):
        st = '{:0.' + str(nd) + 'e}'
        z = []
        for i, v in enumerate(var_list):
            vsc = st.format(v)
            z.append(vsc)
        return z                 
      
        
    def topol2lenghtMx(self, X):
        #jan7, 2020. 
        #repaired by 6may20: z column in X in not the last one anymore. The col is set according to readme_topol.
        z = X[:, 7]
        max_z = max(z)
        rchID_mxZ,= np.where(z == max_z)  #indexes from where function are pythonic, i.e. start at 0.
        IDdown = int(X[rchID_mxZ[0], 3])     
        L_mainstr_topol = X[rchID_mxZ[0], 6]
        while IDdown!=0:    #while not arriving to basin outlet. 
            L_mainstr_topol += X[IDdown-1, 6]  #indexes at X array are fortranic, i.e. start at 1.
            IDdown = int(X[IDdown-1, 3])
        return L_mainstr_topol
            
            
    def pickBasicParams(self, flags, prt_wdw, output_folder, pd):
        #user defined indexes:
        dt_dQ_coarse = self.dt_dQ_coarse
        flume_else_basin = self.flume_else_basin
        which_flume = 1 #just to fill argument to create topology, even if basin does not need it.
        if flume_else_basin==1: 
            which_flume = int(pd['which_flume'])
            self.which_flume = which_flume
        fi = 'doct_supply/supply_wat_sed'
        n = ['fluv_subbasin', 'geotech2', 'A_hs', 'temp_resol', 'dT_days', 'fluv_basin', 'fluv_geom', \
            'fluv_hc_powlaws', 'fluv_qs_powlaw', 'fluv_wlogist', 'fluv_gtm', 'gsd'];
        fls = lsf.f(pd[n[0]]);  gt = lsf.f(pd[n[1]]);  tr = lsf.f(pd[n[3]]);  n_days = int(pd[n[4]]);  
        A_hs = float(pd[n[2]])
        flb = lsf.f(pd[n[5]]); fl_fg = lsf.f(pd[n[6]]); fl_fh = lsf.f(pd[n[7]]); fl_fqs = lsf.f(pd[n[8]]);
        fl_gtm = lsf.f(pd[n[10]])
        fl_gsd = lsf.f(pd[n[11]])
        fl_fw = self.fl_fw  #lsf.f(pd[n[9]])
        fl_fwqsy = self.fl_fwqsy
        msg2 = 'Min/max number of grain size classes is 5/12, to account for sand, D50, D84, D90 and Dmx; \
            avoiding exceed number of grain sizes estimated by Recking13 GSDini_model.py.'
        nD = fl_fg[5]
        if (nD<5 or nD>12): msg2
        if fl_gsd[2]<fl_gsd[1]: sys.exit(f'error: D84={fl_gsd[2]}<D50={fl_gsd[1]}')
            #see examples of reading both single and list parameters in main_Ptx2supply.py.
        nDh = gt[5]
        vh = fls[0];  v_mn_mainstr = fls[1];  v_mx_mainstr = fls[2]
        if v_mn_mainstr > v_mx_mainstr: sys.exit('v_mn_mainstr > v_mx_mainstr !!')        
        dt_roe = tr[1]; refine_pctile = tr[2]; prl = tr[3]
        h_s_ini = flb[1];  w_b_ini = flb[2]
#        nW = flb[0]        

        #drainage network topology and hillslope base time tb_hs of SCS unit hydrograph (see my msc thesis, Cataño2015):
        fn_top = 'Topol_X' if flume_else_basin ==0 else 'Topol_X_flume'
#        if not os.path.exists(pthtop): mt.fn()
        mt.fn(flume_else_basin, which_flume)
        pthtop= f'doct_topol/{fn_top}.txt'      
        Xprev = np.loadtxt(pthtop); n_reaches = len(Xprev[:, 0])
            #Xprev means before cheking is domain size (L, Wv) was not updated in param_dict.txt.
        if flume_else_basin==1:
            self.Lf = float(pd['L_flume'])
            self.dx = float(pd['dx_flume'])
            nR_now = int(self.Lf/self.dx)
            print(f'nR_now=', nR_now)
            if n_reaches!=nR_now: mt.fn(flume_else_basin, which_flume)      #abs(n_reaches - nR_now) > 1e-2: mt.fn()
#            execfile('main_topol.py')
##            sys.exit('update topol by running main_topol.py')
            if abs(Xprev[0, 8]- float(pd['Wflume'])) >1e-2: mt.fn(flume_else_basin, which_flume)
        X = np.loadtxt(f'doct_topol/{fn_top}.txt')
        if flume_else_basin==0:
            nR_now = len(X[:,0])
#        print(f'in py, 11apr20, last col that is wc: \n {X[:,-1]}')
        L_reach = X[:, 6]   #2ene2019: 1D array; elements refer to each subbasin.
        L_hs = np.round(A_hs/(2*L_reach))   #might be greater due to topology of drainage network in subbasin.
        tc_hs = L_hs/vh
        dur_rain = tc_hs;  lag = .6*tc_hs;  tp_hs = dur_rain/2 + lag;  tb_hs = 2.7*tp_hs
#        prl = 10 
        dt_ev = int(round(min(tp_hs)/prl))  if flume_else_basin == 0 else float(pd['dt_flume'])
            #min is used as tc = f(Lreach).        
        print(f'dt_ev{dt_ev}, dt_dQ_coarse{dt_dQ_coarse}')
#        dt_ev: check if it is smaller than base time of event hydrograph, as this must include several timesteps. 
#            For flume, 0.5 minutes works for fractional model calibration of UBCtrechExp.
#            For Tokunaga basin, if A_hs is around 1km2, tc~12min, to tb_ev (base time of event hydrograph) is 2.7tc:    
#            30-35min. Thus, 5-10min seems ok. As dt_ev holds for all reaches and only one day per event is downscaled, then    
#            flood transit time < 1 day.
        Lr_min = min(L_reach)
        dt_evt = int(round(Lr_min/v_mx_mainstr))  #Courant condition for flow velocity around 3m/s while event travels 
            #through mainstream toward the basin outlet.      
        L_mainstr_topol = self.topol2lenghtMx(X)
        L_mainstr_Hack = 1e3*2.02*X[-1,4]**.49 #Hacks law [km, km2], after Montgomery&Dietrich92.
        t_mainstream = L_mainstr_topol/v_mn_mainstr
        if flume_else_basin==0:
            print(f'tc [seg]{tc_hs}')        
            print(f'Mainstream in whole basin A{X[-1,4]}km2. \
                Lhack/Ltopol: {round(L_mainstr_Hack,0)}/{round(L_mainstr_topol,0)}, t_mainstream<={round(t_mainstream, 1)}')        
            print(f'dx/dt_ev = {round(Lr_min/dt_ev, 1)} inside triangular subbasin event')
            print(f'dx/dt_evt = {Lr_min/dt_evt} while event travel in basin')
            print(f'dx/dt_roe = {Lr_min/dt_roe} while event day finishes, after flood traveled in basin.')
        else:
            print(f'Lr_min{Lr_min}')
            print(f'dx/dt_ev = dx/dt = {round(Lr_min/dt_ev, 4)} for all time steps in flume')
            
        #supply flows:
        #Reads in Vtx each supply set (water, sedWash and sedSlide) of timeseries for all reaches x.
        f = ['Qhill', 'wash', 'slid']
        sufs = ['_tx.txt', '_tx_flmTrench.txt', '_tx_flmFlood.txt']
        n_supplyItems = len(f)
        if flume_else_basin==0:
            suf = sufs[0]
        elif which_flume==1:
            suf = sufs[1]
        else:
            suf = sufs[2]
            df = pp13.get_df()
            print(f'11apr20: Pitlick13 flood database: \n {df}')
            nQsteps = df.shape[0]
            t_perQ_min = df['t_min'].to_numpy()
#            print(f'9apr20, nQsteps={nQsteps}')
#        suf = sufs[0] if flume_else_basin==0 else sufs[1]        
        fp = f'{fi}/{f[0]}{suf}'
        cnd = op.isfile(fp)
#        print(f'9apr20: cnd={cnd}')
        if cnd:
            sar = np.loadtxt(f'{fi}/{f[0]}{suf}')
#        else
#            sys.exit(f'No file whose name ends in {suf} exists.')
        cnd2 = flume_else_basin==0 or (flume_else_basin==1 and which_flume==1)        
        ntc = n_days if cnd2==True else nQsteps #22ago20; nQsteps only for pp13 flume?
        Vtx = np.zeros((ntc, nR_now, n_supplyItems))
#        k2i = 2     
           
#        if (flume_else_basin==1 and (cnd==False or sar.shape[1] != nR_now)): #27ago2020, p276d; quito if: update always.
        sf.fn(sufs[1:], which_flume) #update size of input file due to change of problem dimensions from param_dict.txt.
            
#            k2i+=-1
#        else:
#            Vtx[:, :, 1] = sar
        for k2 in range(1, n_supplyItems+1): #k2i, 
#            print(f'6apr20; before supply change, k2={k2}, sar.shape[1]={sar.shape[1]}, nR_now={nR_now}')
#            print(f'6may20; will read this file: {fi}/{f[k2-1]}{suf}')
            fname = f'{fi}/{f[k2-1]}{suf}'
            print('25ago20; fname: ', fname)
            sar = np.loadtxt(fname)
#            print(f'6apr20; k2={k2}, after supply change, sar.shape[1]={sar.shape[1]}, nR_now={nR_now}')
            Vtx[:, :, k2-1] = sar
            
        #produce timeseries of timesteps dt(t)
        if flume_else_basin==0:
            nt, dt_serie, n_copies_ev, n_copies_evt, n_copies_roe, cdt = dtfq.refine_dt(n_days, nR_now, flume_else_basin, \
                Vtx[:, :, 0], Vtx[:, :, 2], dt_ev, dt_evt, dt_roe, dt_dQ_coarse, tb_hs, t_mainstream, refine_pctile)
            dt_dQ_crs_arr = [dt_dQ_coarse]*ntc
            dt_dQ_crs_arr = np.array(dt_dQ_crs_arr)
        elif which_flume==1:
            nt, dt_serie, n_copies_ev, n_copies_evt, n_copies_roe, cdt = dtfq.refine_dt_flmTrench(n_days, dt_ev, dt_dQ_coarse)
            dt_dQ_crs_arr = [dt_dQ_coarse]*ntc
            dt_dQ_crs_arr = np.array(dt_dQ_crs_arr)            
        else: #if flood flume, note that coarse dt (dt_dQ_crs_arr) is variable as per experimental database.
            dt_dQ_crs_arr, nt, dt_serie, n_copies_ev, n_copies_evt, n_copies_roe, cdt = dtfq.refine_dt_flmFlood(
            nQsteps, dt_ev, 60*t_perQ_min)
        ddt_args = [dt_ev, dt_evt, dt_roe, n_copies_ev, n_copies_evt, n_copies_roe] 
        #, dt_dQ_coarse: deprecated, as used only for flume ddt, and flume can have variable dt_dQ_coarse e.g. pp13 flood.
            
        #produce initial alluvium GSD to get mean gravel and sand size, and initial sand fraction.
        if (flume_else_basin==0 or (flume_else_basin==1 and which_flume==1) ):
            fb, Dfbr, self.Dc = gsd.calc(nD, flume_else_basin, fl_gsd, which_flume)
        else:
            fb, Dfbr, self.Dc  = gsd.readFlumeFlood(which_flume, flume_else_basin)
        if nD!=len(Dfbr): sys.exit('read gsd and nD differ')
        self.fb = fb
        m1= 'Check if current initial GSD is the one showed in figure at /doct_fluvial. \n'
        m2= '  Quit comment in call of plot method at GSDini_model.py to update the figure.'
        msg = m1+m2
        print(msg)
#        sys.exit('6may20, paro pa debug de gsd based on its prints')
        
        #output paths:
        w_v = X[:, 9-1]
        dt_ind = lsf.f(pd['temp_resol'])[3]  if flume_else_basin==0 else  float(pd['dt_flume'])   #prl: pts_rising_limb        
        caseNum = int(pd['caseNum'])
        st_dt_ind = str(int(dt_ind))
        caseNum = str(int(caseNum))
        if flume_else_basin==1:
            dx = self.dx
            if which_flume==1: #ubc trench
                Lf = self.Lf
                pwf = 'trench/' #'path for which flume'
                case = '_10wv' + str(int(round(10*w_v[0],0))) + '_Lf' + str(int(Lf)) + '_10dx' + \
                    str(int(round(10*dx, 0))) + '_dt' + st_dt_ind
            else:
                pwf = 'flood/'
                case = '_10dx' + str(int(round(10*dx, 0))) + '_dt' + st_dt_ind
        else:
            pwf = ''
            case = '_case' + caseNum + '_dtPrl' + st_dt_ind #6may20: prl means points in rising limb of hydrograph.
            #build here the name of sensitivity case with parameters that are varied from param_dict. #see phdnotes p172dr.
        ofpb = output_folder + pwf + 'bins' + case  #output_folder is already affected by 'environment': either flume or basin.
        if not os.path.exists(ofpb): os.makedirs(ofpb)
            #definitive output_folder for paths. Note all folder name strings end with '/', 
                #except for the smaller hierarchy subfolder, as filename is added with its '/' in f90 code that writes bin output.
        b = vp.a(self.HGelseFract)                
        ooVarNames_list, path_list2f90, pos_ini_var, nvars, num_char, num_paths = b.path_list_builder(1, int(nD), ofpb)
        print(f'25apr20: output_folder for paths of binaries: \n {ofpb}')
        print(f'25apr20: pathList names length: {num_char}, num_paths: {num_paths}')
        np.savetxt(f'{ofpb}/pos_ini_var.txt', pos_ini_var, fmt='%i')
        np.savetxt(f'{ofpb}/ooVarNames_list.txt', ooVarNames_list, delimiter='', fmt='%s')
        np.savetxt(f'{ofpb}/dt_serie.txt', dt_serie)
        np.savetxt(f'{ofpb}/cdt.txt', cdt, fmt = ' '.join(['%i']*2))
        np.savetxt(f'{ofpb}/nt_numChar.txt', np.array([nt, num_char])) 
#        sys.exit('6may20: paro aquí para debug de path pa bin outputs')

#        print(f'26apr20; OEELO, tuvo q haber guardado aux vars en esta ofpb: \n {ofpb}')
#        sys.exit(f'26apr20; ooVarNames_list: \n {ooVarNames_list}')

        #pick parameters (arguments) for fluvial model:
        print(f'nt{nt}, n_days{n_days}')
#        print(f'Are nt/n_days int?: {isinstance(nt, int)}/{isinstance(nt, int)}')
#        print(f'pathList: {path_list2f90}')

        calib = float(pd['calib']) 
        calibv = float(pd['calibv'])        
        print(f'31mar; calib{calib}')
        basic_args = (prt_wdw[0], prt_wdw[1], prt_wdw[2], prt_wdw[3], \
        	ofpb, pos_ini_var, \
        	which_flume, \
            flags[0], flags[1], flags[2], flags[3], flags[4], flags[5], flags[8], \
        	ddt_args, dt_dQ_crs_arr, dt_serie, tb_hs, \
            path_list2f90, \
        	h_s_ini, w_b_ini, fl_fg, fl_fh, fl_fqs, fl_fw, fl_fwqsy, fl_gtm, calib,calibv, \
        	X, fb, Dfbr, Vtx[:, :, 2], Vtx[:, :, 1], Vtx[:, :, 0])
        nD = int(nD), #see 'nd' instead 'nD' in keywords below.                    	
#        print(f'26mar, before basic_kwargs, nD={nD}, nt={nt}')
        basic_kwargs = {'nd':nD, 'nt':nt, 'nr':nR_now, 'n_topol_vars':len(X[0,:]), \
            'num_paths':num_paths, 'nvars': nvars, 'num_char':num_char, 'n_coarse_dt':ntc}
        return cdt, dt_serie, Vtx[:, :, 2], Vtx[:, :, 1], Vtx[:, :, 0], basic_args, basic_kwargs
        

#    def pickFractionalParams(self, flags, pd):
#        #output paths:
##        b = vp.a(HGelseFract)        
##        path_list2f90, pos_ini_var, nvars, num_char, num_paths = b.path_list_builder(nW, nDh)        
##        return fr_args


    def runFluvialModel(self, prt_wdw, prt_ddt, flags, pd):
        fl_fw = self.fl_fw
        dt_dQ_coarse = self.dt_dQ_coarse
        fo1 = 'doct_fluvial/'
        flume_else_basin = self.flume_else_basin
        if flume_else_basin == 1:
            sf = 'hg_model_flume/'
        else:
            sf = 'hg_model/'
        fo = fo1 + sf
        cdt, dt_serie, slid_day, wash_day, Q_day, basic_args, basic_kwargs = self.pickBasicParams(flags, prt_wdw, fo, pd)
        flag_day_else_tStep = flags[4]
#        d2pi = 1;  d2pf = 10
        d2pi = int(prt_ddt[0]);  d2pf = int(prt_ddt[1])
        x2pi = int(prt_ddt[2]);  x2pf = int(prt_ddt[3])
        if self.HGelseFract == 1:
#            fc.cleanFolder(fo) #25apr20: dont use cleanFolder any more; erased input files accidentally!. Overwrite files instead.
            t_fine = np.cumsum(dt_serie) + 1 #/float(dt_dQ_coarse) #*(1 if flume_else_basin ==0 else 0)
                #so 1<t<2 for fine steps of day 1.
#            t_fine_days = np.append(1, t_fine_days[:-1])    #left-offset for timeseries, so plot day value at morning if no event.

#            print ('vea doc d f90 fluvial module pues: ')
#            print(hg_model_class.model_sc.__doc__)
#            print ('vea q TERMINÓ D MOSTRAR doc d f90 fluvial module pues')

            ti = time.time()            
            s2fma= s2fm.a()
            Qhill, Q, wash, slid, v, h, vf, w, c, D50, D84, D50c_ss, D84c_ss, S, h_s, h_sf, w_c, w_b, w_v, ts, tc84, \
                tsm, Fr, aspect, subm, Mgtm, gam, pi  \
                = s2fma.calc(basic_args, basic_kwargs)
            tf = time.time();  t = tf-ti #self.sci([tf-ti], 1)[0]
            print('----------------------------')
            print(f'in supply2Fluvial.py, main f90 calc (s2fma) lasts {t}s')               
            plot_series = int(pd['plot_series'])
            aa = ps2f.c(plot_series, fo, fl_fw, self.Dc, self.fb, self.Lf, self.dx)
            aa.plt_hc(self.which_flume, flags[7], prt_wdw, dt_serie, Q, v, h, vf, w, c, D50, D84, D50c_ss, D84c_ss, S, 
                h_s, h_sf, w_c, w_b, w_v, ts, tc84, tsm, Fr, aspect, subm, Mgtm, gam, pi)
            print('in py, from f90, lets print downscaled supply!!')
            var = ['Qhill', 'wash', 'slid']
#            if flag_day_else_tStep == 1:    #as this plot mainly compare fine and daily resolutions.
            for j in range(x2pi, x2pf + 1):
                print('----------------------------------')                    
                suptt = f'Reach{j}';  print(suptt)
                Qdj = Q_day[:, j-1];  wdj = wash_day[:, j-1];  sdj = slid_day[:, j-1]
                Qj = Qhill[:, j-1];  wj = wash[:, j-1];  sj = slid[:, j-1]      #downscaled at intradaily timescale.
                Qm = Qdj.mean()
                print(f'downscaled Qh_m{Qm}')
                print(f'mostrame Qhill en py tras f90: {Qj}')
                sys.exit(f'16may20; stop here to rev Qh')

#                Qj = Qj/Qj.mean()*Qm;  wj = wj/wj.mean()*wdj.mean()
                    # mean-preserving rescaling; don't needed by slides as they don't suffer triangular pulse, see ddt.f90.
                vj = np.stack((Qj, wj, sj), axis = 1);  vdj = np.stack((Qdj, wdj, sdj), axis = 1)
                    #stack process comes by columns, which is axis 1. Then, each vector runs by rows.                
                for k, v in enumerate(var): #plot each of the 3 supply variables: Q, wash and slide.
                    vkj = vj[:, k];   vdkj = vdj[:, k]
                    l = [f'{var[k]}', f'{var[k]}_day', 'dt']
                    print('----------------------------------')
                    print(f'6mar20: estamos en j{x2pi + j -1}, var: {v}')
                    print('6mar20: len(dt_serie)= ', len(dt_serie))
                    vkjm = np.average(vkj, weights = dt_serie)
                    vl = [min(vkj), vkjm, max(vkj), min(vdkj), np.mean(vdkj), max(vdkj)]
                    z = self.sci(vl, 1)     #to print values list in scientific notation.
                    print(f'{v}_MinMeanMax {z[0]}/{z[1]}/{z[2]}. \
                        DAILY, {v}_MinMeanMax {z[3]}/{z[4]}/{z[5]}')                                            
                    aa.plt_ddt(flume_else_basin, cdt, t_fine, vkj, vdkj, dt_serie, l, d2pi, d2pf, suptt, j)
        else:
            fo = fo1 + 'fract_model'; 
#            fc.cleanFolder(fo)
#           fr_args = join list of basic_args and FractParams output
#            fract_model_class.model_sc(*basic_args)
        print("review if effectively results were saved as pngs")
