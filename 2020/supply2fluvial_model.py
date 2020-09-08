    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#jan14, 2020

import numpy as np

#import supply2fluvial_model_varxt as s2fmv


#from fract_model_class import *
#from ddt_class import *
#import ddt_class
from hg_model_fldbin_ubc_qt import *       #4ene2020:note that only this syntax works, instead of simply 'import <module>'.

class a:
    def __init__(self): #, args, kargs
        xx=1
#        par_w = calc_params_wat()       #par is tuple, len(par) = 2; one per cross section zone: channel or floodplain.    

#        self.calc(args, kargs) #, par_w

    def calc(self, args, kargs):  #, par_w
        Qhill, Q, wash_sh, slid_sh, v, h, vf, w, c, D50, D84, D50c_ss, D84c_ss, \
            S, h_s, h_sf, w_c, w_b, w_v, ts, tc, tsm, Fr, aspect, subm, Mgtm, gam, pi = \
            hg_model_fldbin_ubc_qt.model_sc(*args, **kargs)
            #mar5, 2020: *basic_args, **basic_kwargs
            #jan14, 2020: only applied to get ddt.
            #jan14, 2020: too many args were passed when I though physics does not need numerical solvers for eqns, but it DOES.
                #Given scipy.optimize.root call minpack from .f file, efficiency must be decent.
                
#        explode params                
#        z[1-1, :] = z_ini
#        Sor= X[:,6-1]; dx = X[:,7-1]; Wv= X[:,8-1]  #jan14, 2020: column 8 still does not exist. Update Xtopol calc to produce it.
#        jup, jdown = self.get_neighbors(X[:,4-1], nR)
#        for i in range(1, nt + 1):
#            i_plot =  i if flag_day_else_tStep == 0 else day_of_calc[i-1]
#            for j in range(1, nR + 1):
#                St+1, Wt+1, vol_bank_fail = s2fmv.update_geom(z_next[0], Q-Q2, D50[0])
#                Q, Q2, v1, v2 = s2fmv.wat_balance(Qhill i, store_prev, dxj, dti, Wvj, Wij, Sij)
#                d, D50, D84 = s2fmv.get_flow_traits(Q, Q2, v1, v2, Wvj, Wij, Sij)
#                qsb1, qsb2 = s2fmv.get_sed_tpt_capac(d_bf, dij, Sij)
#                z_next = s2fmv.sed_balance(vol_bank_fail, W, qb1, qb2)

        return Qhill, Q, wash_sh, slid_sh, v, h, vf, w, c, D50, D84, D50c_ss, D84c_ss, \
            S, h_s, h_sf, w_c, w_b, w_v, ts, tc, tsm, Fr, aspect, subm, Mgtm, gam, pi

        
#    def get_neighbors(Rd, nR):
#        jup = np.zeros(nR);  jdown = np.zeros(nR)
#        Rdi = Rd.astype(int)
#		for j, rd in enumerate(Rdi):
#		    if (Rdi[j] > 0):  #if it is not basin outlet
#    		    jup[Rdi[j]] = j
#    		else:
#    		    jup[j] = nR-1
#    		jdown[j] = Rdi[j]
#    return jup, jdown 
        
        
    def ini_geom(zr, A_basin_km2):
        D50ini = 1e-2
        Wvj_de_X    
        return z_next_ini, Qbf_reg, D50_ini

#    def get_hydraulics():
#        solve reck 11 given S&q
#        return d, D50, D84
        
        
#    def get_sed_tpt_capac(d_bf, d, S):
#        solve reck 11 given S&q
#        return qsb1, qsb2        

        
#    def get_flow_traits(Q, Q2, v1, v2, Wvj, Wij, Sij):
#        q = Q/W
#        d, D50, D84 = self.get_hydraulics(q):  #D50 & D84 have shape 2; one per cross section zone: channel or floodplain.
#        qsb1, qsb2 = self.get_sed_tpt_capac(q, d_bf, d, S)        

        
#    def calc_params_w():
#        #channel:
#        pdict = {}
#        calc coefs y expons pa NoFlood
#        
#        #floodplain:
#        pdict = {}        
#        calc coefs y expons pa flood
#        
#        
#    def wat_balance(self, par_w):
#        solve eqn
#        if d > d_bf:    #flood
#            solve eqn
#        return Q, Q2, v1, v2

#            
#    def sed_balance(vol_bank_fail, W, qb1, qb2):
