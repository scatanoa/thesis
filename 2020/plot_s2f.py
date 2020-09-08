    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#jan11, 2020.

import numpy as np
import matplotlib.pyplot as plt
import FolderCleaner as fc
import matplotlib.gridspec as grsp
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker
import os


class c:

    def __init__(self, plot_series, of, par_wlog, Dc, fb, Lf, dx):
        self.dx = dx
        self.plot_series = plot_series
        self.of = of
        self.par_wlog = par_wlog
        self.Dc = Dc
        self.fb = fb
        self.Lf = Lf

    def seriesSed(self, fon, dt, t, cumOut, Q, braidInd, S, c, ts, tc, tsm, pi, Mgtm, gam, j):
#        fig,ax=plt.subplots(nrows = 2, ncols = 1)
        print(f'In plots2f.py/seriesSed; printing j={j}')
        col = 'tab:red'        
        Dc = self.Dc
        fb = self.fb
        ms= 2
        gs = grsp.GridSpec(4, 1, height_ratios = [3, 2, 3, 3])
        fig = plt.figure()
        axd = fig.add_subplot(gs[0,:])
        axb = fig.add_subplot(gs[1,:])
        axa = fig.add_subplot(gs[2,:])                
        axc = fig.add_subplot(gs[3,:])
        axb.yaxis.set_major_locator(ticker.MultipleLocator(100))
        l1 = axb.plot(t, cumOut[:-1], '--k', label = 'cum', markersize=ms)
        axb.set_ylabel('cumSedTpt\n(kg)', fontsize=9)
        ax1 = axb.twinx()
        ax1.set_ylim(bottom=1e-6, top=1e-3)
        ax1.set_yscale('log')        
        majors = [1e-6, 1e-5, 1e-4, 1e-3]
        ax1.yaxis.set_major_locator(ticker.FixedLocator(majors)) #1apr20: LogLocator did not work for 2nd axis.        
        l2 = ax1.plot(t, c, '-r', label = 'c', markersize=ms)
        ax1.set_ylabel('$c (m^3/m^3)$', color=col)
        ax1.tick_params(axis = 'y', labelcolor=col)  
#        ax1.yaxis.grid(which='both')        
#        ax2 = axb.twinx()
#        l3 = ax2.plot(t, D84, '.g', label = '$D_{84}$', markersize=ms)
#        ax2.set_ylabel('$D_{84}$(m) or S(m/m)')
#        ax2.set_yscale('log')

#        ax2.set_ylim(bottom=1e-3, top=1e-1)

#        axb.set_ylim(bottom=0, top=1e-3)
        ls=l1+l2
        labs = [l.get_label() for l  in ls]        
        leg = axb.legend(ls, labs, loc='best', ncol=1, fontsize=9) #ncol=2, , fancybox=True

        l1= axd.plot(t, ts, '-b', markersize=ms, label = r'$\tau^*$')
        l2= axd.plot(t, tc, ':k', markersize=ms, label = r'$\tau^*_{c84}$')
        axd.set_ylabel(r'$\tau^*$ ; $\tau^*_{c84}$')
        axd.set_ylim(bottom=1e-3, top=1e-1)
        axd.yaxis.grid(which='both', color='0.8')
        axd.set_yscale('log')
        ax8=axd.twinx()
        lbtc = r'$\tau^*/\tau^*_{c84}$'
        l3 = ax8.plot(t, ts/tc, '-.r', markersize=ms, label = lbtc)
        lbtsm = r'$\tau^*/\tau^*_{m}$'
        l4 = ax8.plot(t, ts/tsm, '--r', markersize=ms, label = lbtsm)
        ax8.tick_params(axis = 'y', labelcolor=col)
        ls = l1+l2+l3+l4
        labs = [l.get_label() for l  in ls]        
        leg = axd.legend(ls, labs, loc='best', ncol=4, fancybox=True, fontsize=8) #ncol=2,        
        ax8.set_yscale('log')
        ax8.set_ylim(bottom=1e-1, top=1e1)
        lbtc2 = lbtc + r' ; $\tau^*/\tau^*_{m}$'
        ax8.set_ylabel(lbtc2, color=col)
        
        nD = len(pi[0,:])
        nD2 = nD+1
        kl = np.arange(nD)
        axa.set_yscale('log')
        axa.set_ylim(bottom=1e-2, top=1e1)        
        majors = [1e-2, 1e-1, 1, 10]
        axa.yaxis.set_major_locator(ticker.FixedLocator(majors)) #1apr20: LogLocator did not work for 2nd axis.        
        for k in kl:
            cl = (kl[k]/nD2, kl[k]/nD2, kl[k]/nD2)
            Dck = Dc[k]
            lbn = np.around(Dck, decimals=2) if Dck<1 else int(Dck)
            lb = str(lbn) + 'mm'
            pik = pi[:,k]
            rat = pik/fb[k]
            rat[rat<1e-2] = np.nan
            axa.plot(t, rat, '-', c=cl, markersize=ms, label = lb)
        axa.legend(title='D(mm)',fontsize=6,ncol=6, fancybox=True, title_fontsize=7)
#        axa.set_ylim(bottom=0, top=pi_mx)
        axa.set_ylabel('$p_i / f_{i,b}$')
        ax1 = axa.twinx()
        Mgtm[-1]=Mgtm[-2] #to avoid 0
        gam[-1]=gam[-2] #to avoid 0
        majors = [1e-1, 1, 1e1, 1e2]
        ax1.set_yscale('log')
        ax1.set_ylim(bottom=1e-1, top=1e2)
        ax1.yaxis.set_major_locator(ticker.FixedLocator(majors))                                            
        ax1.plot(t, Mgtm, '--r')
        ax1.plot(t, gam, ':r')
        ax1.set_ylabel(r'(--)$M_{gtm}$; (..)$\gamma_{gtm}$', color=col)
        ax1.tick_params(axis = 'y', labelcolor=col)
        
        l1= axc.plot(t, Q, ':k', label = 'Q')
        axc.set_ylabel('$Q(m^3/s)$')
#        axc.set_xlabel('time (hours)') #dont work: , ha = 'left', va = 'top'
        axc.xaxis.set_label_coords(.95,-.1)
        axc.set_xlabel('time (hours)', ha='left')
        ax1 = axc.twinx()
        ax1.set_ylim(0, 1)
#        ax1.set_yscale('log')        
#        ax1.set_ylim(1e-3, 1e-1)
#        majors = [1e-3, 1e-2, 1e-1]
#        ax1.yaxis.set_major_locator(ticker.FixedLocator(majors))        
        lbi = 'i'
        ons = np.ones(len(t))   #reference braidInd for AB (anabranching) and B (braiding)
        mc=.4   #after Eaton10, threshold braidInd for multichannel mc is .28 or .4 for eqns 4 or 7 respect/. 
            #In this eq 7, threshold grows with bank resistance miu.
        ana_arr = mc*ons; bra_arr = 1.82*ana_arr
        l2 = ax1.plot(t, braidInd[:-1], '-r', label = lbi)
        l3 = ax1.plot(t, ana_arr, '-.r', label = 'i_ana')
        l4 = ax1.plot(t, bra_arr, ':r', label = 'i_bra')
        lbi2 = lbi + r'$=SQ^{*0.44}$'
        ax1.set_ylabel(lbi2, color=col)
        ax1.tick_params(axis = 'y', labelcolor=col)        
        ls=l1+l2+l3+l4
        labs = [l.get_label() for l  in ls]        
        leg = axc.legend(ls, labs, loc=8, ncol=2, fancybox=True, fontsize=8) #ncol=2,

        axc.get_shared_x_axes().join(axb, axc, axd, axa)
        axa.set_xticklabels([])
        axb.set_xticklabels([])
        axd.set_xticklabels([])                
        
        fn = '/sed_x' + str(j) + '.png'
        gs.tight_layout(fig)
        gs.update(hspace=.13)
        fnfull = fon + fn
        plt.savefig(fnfull)
        plt.close(fig)

        
    def series(self, fon, t, c, D50, D84, jam, D50ss, D84ss, Q, h, v, vf, w, S, h_s, h_sf, w_c, w_b, w_v,  
        Fr, asp, subm, j):        
        col = 'tab:red'
#        fig,ax=plt.subplots(nrows = 2, ncols = 1)
        print(f'In plots2f.py/series; printing j={j}')
        par_wlog = self.par_wlog
        ms= 2
        gs = grsp.GridSpec(4, 1, height_ratios = [1, 1, 1, 1])
            #https://stackoverflow.com/questions/52560527/how-to-fully-customize-subplot-size-in-matplotlib
        fig = plt.figure()       
        axa = fig.add_subplot(gs[1,:])
        axaa = fig.add_subplot(gs[0,:])
        axb = fig.add_subplot(gs[2,:])
        axc = fig.add_subplot(gs[3,:])
        
        l1 = axa.plot(t, h_s + h, '-b', label = '$h+h_s$', markersize=ms)
        l4 = axa.plot(t, h_sf, '-.', c='.6',  label = '$h_{sf}$', markersize=ms)
        l5 = axa.plot(t, h_s, '--k', label = '$h_s$', markersize=ms)
        axa.set_ylabel('h (m)')        
        ax12 = axa.twinx()        
        l11 = ax12.plot(t, subm, ':r', markersize=ms, label= '$R/D_{84}$')
        ax12.set_ylabel('subm$=R/D_{84}$', color = col)
        ax12.tick_params(axis = 'y', labelcolor=col)
        ax12.set_ylim(bottom=0, top=10)
        l2 = axaa.plot(t, w, '-b', label = 'w', markersize = ms)
        l6 = axaa.plot(t, w_c, ls= '-.', c='.6', label = '$w_c$', markersize = ms) #, linestyle= 'None', marker='.'
        l7 = axaa.plot(t, w_b, '--k', label = '$w_b$')
        l8 = axaa.plot(t, par_wlog[1]*Q**par_wlog[4], ':g', label = '$w_s$')
        w_v_t = w_v*np.ones(len(t))
#        l9 = axaa.plot(t, w_v_t, c='0.2', label = 'w_v', linewidth=.3)
        ax11 = axaa.twinx()
        ax11.set_ylim(bottom=1e0, top=1e2)
        ax11.set_yscale('log')        
        majors = [1e0, 1e1, 1e2]
        ax11.yaxis.set_major_locator(ticker.FixedLocator(majors))
#        l10 = ax11.plot(t, par_wlog[0]*Q**par_wlog[3], c='0.6', label = 't_sat')
        lba = 'a:$w/d$'
        l10 = ax11.plot(t, asp, ':r', markersize=ms, label = 'a')  #this python d is f90 hh: hydraulic depth (A/w)
        lbj = 'j:$w/D_{84}$'
        l12 = ax11.plot(t, jam, '-.r', markersize=ms, label = 'j')
        axaa.set_ylabel('w (m)')
        y2lb = lba + '. ' + lbj
        ax11.set_ylabel(y2lb, color=col)
        ax11.tick_params(axis = 'y', labelcolor=col)
        ls=l1+l5+l4+l11
        labs = [l.get_label() for l  in ls]
        leg = axa.legend(ls, labs, loc='best', ncol=4, fancybox=True, fontsize=9) #ncol=2, 
        leg.get_frame().set_alpha(0.7)
        ls = l2 + l6 + l7 + l8 + l10 + l12 #+ l9
        labs = [l.get_label() for l  in ls]
        leg = axaa.legend(ls, labs, loc='best', ncol=6, fancybox=True, fontsize=9) #ncol=2, 
        leg.get_frame().set_alpha(0.7)
        
        l2 = axb.plot(t, D84, '-k', label = '$D_{84}$', markersize=ms)
        l3 = axb.plot(t, D50, '-g', label = '$D_{50}$', markersize=ms)
        l5 = axb.plot(t, D84ss, ':k', label = '$D_{ss,84}$', markersize=ms)
        l6 = axb.plot(t, D50ss, ':g', label = '$D_{ss,50}$', markersize=ms)        
        axb.set_ylabel('$D$(m)')
        ax2 = axb.twinx()
        l4 = ax2.plot(t, S, '--r', label = 'S', markersize=ms)
        ax2.set_ylabel('S(m/m)', color=col)
        ax2.tick_params(axis = 'y', labelcolor=col)        
#        axb.set_yscale('log')        
        ax2.set_yscale('log')
        axb.set_ylim(bottom=.005, top=.035) #bottom=1e-3, top=1e-1
        ax2.set_ylim(bottom=1e-3, top=1e-1) #bottom=0, top=.08
        ls=l2+l3+l4+l5+l6
        labs = [l.get_label() for l  in ls]        
        leg = axb.legend(ls, labs, loc='best', ncol=5, fancybox=True, fontsize=9) #ncol=2,         


        l1= axc.plot(t, Q, ':k', label = 'Q')
        axc.set_ylabel('$Q(m^3/s)$')
#        axc.set_xlabel('time (hours)')
        axc.xaxis.set_label_coords(.95,-.1)
        axc.set_xlabel('time (hours)', ha='left')        
        ax3 = axc.twinx()
        l2 = ax3.plot(t, Fr, '--r', label = 'Fr', markersize = ms)
        l3 = ax3.plot(t, v, '-b', label = 'v', markersize = ms)
        l4 = ax3.plot(t, vf, '-.b', label = '$v_f$', markersize = ms)
        ax3.set_ylim(bottom=0, top=1.2)
        ax3.set_ylabel('v(m/s); Fr')
        ls=l1+l3+l2+l4
        labs = [l.get_label() for l  in ls]        
        leg = axc.legend(ls, labs, loc=1, ncol=2, fancybox=True, fontsize=8) #ncol=2,

        axa.get_shared_x_axes().join(axaa, axa, axb, axc)
        axa.set_xticklabels([])
        axaa.set_xticklabels([])        
        axb.set_xticklabels([])
        
        fn = '/hc_x' + str(j) + '.png'
        gs.tight_layout(fig)
        gs.update(hspace=.13)        
        fnfull = fon + fn        
        plt.savefig(fnfull)
        plt.close(fig)
                      
            
    def get_t_list(self, dt_serie):
        t = [] 
        t.append(dt_serie[0])
        dts = list(dt_serie)
        dts.pop(0)  #pop(m) erase element in position m.
#        print(f'in get_t_list, dts={dts}')
        for dti in dts:
            t.append(t[-1] + dti)
        return t
        
    def mapPlotter(self, matrix, relative_folder_path, filename, mn, mx, xf, tf, log_cbar=False, tit=None):
        fig=plt.figure()
        plt.title(tit)
        cm = plt.get_cmap('magma')
        cm.set_bad(color='r')
        if log_cbar:
            plt.imshow(matrix, interpolation='nearest', cmap=cm, aspect = 'auto', \
                norm = LogNorm(vmin=mn, vmax=mx), extent=[0,xf,tf,0])
        else:
            plt.imshow(matrix, interpolation='nearest', cmap=cm, aspect = 'auto', \
                vmin=mn, vmax=mx, extent=[0,xf,tf,0])        
        plt.xlabel('downstream distance (m)')
        plt.ylabel('time (hours)', rotation=270)
        plt.colorbar()
        fig.tight_layout()        
        plt.savefig(relative_folder_path + '/' + filename + '.png')
        plt.clf
        
    def mapSaver(self, fon1, case_n, var, var_n, zmn, zmx, xf, tf, log=False, tit=None):
        fi_n = var_n + case_n 
        fon = fon1 + 'tx'
        if not os.path.exists(fon): os.makedirs(fon)
        self.mapPlotter(var, fon, fi_n, zmn, zmx, xf, tf, log, tit)
        
    def Xanomalizer(self, var, nR, nt):
        var_mi = np.mean(var, axis=1)
        varXanom = np.zeros((nt, nR))        
        for j in range(1, nR+1): #each abscissa
            varXanom[:, j-1] = var[:, j-1] / var_mi - 1
        return varXanom
        
    def plotXcorr(self, xn, Qgoal, d, case_n, fon1): #, firstElseLast
        mk = ['xr', 'xm', '+r', '+m', '.r', '.m', '.b', '.k']
        Qgoal = np.array(Qgoal)
        Qlps = Qgoal*1e3
        fon = fon1 + 'tx'
        if not os.path.exists(fon): os.makedirs(fon)        

        fig,ax = plt.subplots(nrows = 2, ncols = 1)
        anom_mn=-.5; anom_mx=.5
        for i, Q in enumerate(Qlps):
            ax[0].plot(d['wXanom'][i], d['hsXanom'][i], mk[i])
            ax[0].set_xlabel(r'$w_{anom}$'); ax[0].set_ylabel(r'$h_{s, anom}$')
            ax[0].set_xlim(left=anom_mn, right=anom_mx); ax[0].set_ylim(bottom=anom_mn, top=anom_mx)
            lb1= f'Q={int(Q)}lps_'
            lb2= 'ini' if i%2==0 else 'fin'
            ax[1].plot(d['wXanom'][i], d['vXanom'][i], mk[i], label=lb1+lb2)
            ax[1].set_xlabel(r'$w_{anom}$'); ax[1].set_ylabel(r'$v_{anom}$')
            ax[1].legend(loc='best', ncol=1, fancybox=True, fontsize=8)
            ax[1].set_xlim(left=anom_mn, right=anom_mx); ax[1].set_ylim(bottom=anom_mn, top=anom_mx)
        fig.tight_layout()        
        fn = 'ej1Xcorrel' + case_n        
        plt.savefig(fon + '/' + fn + '.png')
        plt.close(fig)
        
        fig,ax = plt.subplots()
        for i, Q in enumerate(Qlps):
            lb1= f'Q={int(Q)}lps_'
            lb2= 'ini' if i%2==1 else 'fin'
            ax.plot(d['Fr'][i], d['c'][i], mk[i], label=lb1+lb2)
            ax.set_xlabel('Froude number'); ax.set_ylabel('bedload volum. concentration')
            ax.set_xscale('log'); ax.set_yscale('log')
            ax.set_xlim(left=.1, right=2)
            ax.set_ylim(bottom=1e-7, top=1e-2)
            ax.legend(loc='upper left', ncol=1, fancybox=True, fontsize=8)
        fig.tight_layout()        
        fn = 'ej2Xcorrel' + case_n
        plt.savefig(fon + '/' + fn + '.png')
        plt.close(fig)
        
        fig,ax = plt.subplots()
        for i, Q in enumerate(Qlps):
            lb1= f'Q={int(Q)}lps_'
            lb2= 'ini' if i%2==1 else 'fin'
            ax.plot(d['grad_wgw'][i], d['gradNeg_hs'][i], mk[i], label=lb1+lb2)
            ax.set_xlabel('grad_wgw'); ax.set_ylabel('gradNeg_hs')
            ax.set_xlim(left=anom_mn, right=anom_mx); ax.set_ylim(bottom=anom_mn, top=anom_mx)
            ax.legend(loc='best', ncol=1, fancybox=True, fontsize=8)
        fig.tight_layout()        
        fn = 'ej3Xcorrel' + case_n
        plt.savefig(fon + '/' + fn + '.png')
        plt.close(fig)        
        
    def setThreshVals(self, v, tol, vc):
        f = np.array([1-tol, 1+tol])
        vlim = vc*f  
        print(f'1apr; f={f}; vlim[0]={vlim[0]}, v.shape={v.shape}')
        cond = np.logical_and(v>vlim[0], v<vlim[1]) 
        v[cond] = np.nan
        v_msk = np.ma.masked_invalid(v)
        return v_msk
            
    def plt_hc(self, which_flume, plot_whole, prt_wdw, dt_serie, Q, v, h, vf, w, c, D50, D84, D50c_ss, 
        D84c_ss, S, h_s, h_sf, w_c, w_b, w_v, ts, tc84, tsm,
        Fr, aspect, subm, Mgtm, gam, pi):
        
        #set last time NaN for all j, to avoid line fall at end of timeseries plots.
        Q[-1,:]=np.nan; v[-1,:]=np.nan; h[-1,:]=np.nan; vf[-1,:]=np.nan; w[-1,:]=np.nan; 
        c[-1,:]=np.nan; D50[-1,:]=np.nan; D84[-1,:]=np.nan; D50c_ss[-1,:]=np.nan; 
        D84c_ss[-1,:]=np.nan; 
        S[-1,:]=np.nan; h_s[-1,:]=np.nan; h_sf[-1,:]=np.nan; w_c[-1,:]=np.nan; 
        w_b[-1,:]=np.nan; ts[-1,:]=np.nan; tc84[-1,:]=np.nan; tsm[-1,:]=np.nan;
        Fr[-1,:]=np.nan; aspect[-1,:]=np.nan; subm[-1,:]=np.nan; Mgtm[-1,:]=np.nan; 
        gam[-1,:]=np.nan; pi[-1,:]=np.nan
        
        Lf = self.Lf
        dx = self.dx
        nR = len(w_v)
        nt = len(Q[:, 0])
        fowe = 'calib/' if which_flume == 1 else 'calib_flood/'
        fon1 = 'doct_fluvial/hg_model_flume/' + fowe 
        if not os.path.exists(fon1): os.makedirs(fon1)
        case_n = '_10wv' + str(int(round(10*w_v[0],0))) + \
            '_Lf' + str(int(Lf)) + '_10dx' + str(int(round(10*dx, 0))) + '_dt' + str(int(dt_serie[0]))
        
#        CORRELACION ESPACIAL ANCHO VS LECHO DE SHAWN CHARTRAND18jgr
        wXanom = np.zeros((nt, nR))
        FrXanom = np.zeros((nt, nR)) 
        grad_wgw  = np.zeros((nt, nR))
        gradNeg_hs  = np.zeros((nt, nR))
#        hsXanom = np.zeros((nt, nR))
        wmi= np.mean(w, axis=1) #mean along channel, given t.
#        hsmi = np.mean(h_s, axis=1)
        FrXanom = self.Xanomalizer(Fr, nR, nt)
        wXanom = self.Xanomalizer(w, nR, nt)
        hXanom = self.Xanomalizer(h, nR, nt)        
        vXanom = self.Xanomalizer(v, nR, nt)
        hsXanom = self.Xanomalizer(h_s, nR, nt)        
        for i in range(1, nt+1):
            wdif= np.diff(w[i-1, :])
            hsdif = -np.diff(h_s[i-1, :])
#            print(f'wdif={wdif}')
            gwl = [0]
            ghsl = [0]
            gwl.extend(wdif.tolist())
            ghsl.extend(hsdif.tolist())
#            print(f'gwl={gwl}')
            grad_wgw[i-1, :] = np.array(gwl)/dx
            gradNeg_hs[i-1,:] = np.array(ghsl)/dx        

        t = self.get_t_list(dt_serie)        
        t = np.array(t)
        t=t/3600

        Qprk = Q/(np.sqrt(16.5)*D50**2.5)
        bb = .44
        braidInd = S*Qprk**bb
        
        jam = w/D84
        ttc = ts/tc84
        ttm = ts/tsm

        cumOut = np.zeros((nt, nR))                
        for j in range(1, nR+1):
            cumOut[:-1, j-1] = np.cumsum(2650*c[:-1, j-1]*Q[:-1, j-1]*dt_serie[:-1]) #[:-1]
            
        if which_flume ==1: #trench            
            Qgoal=[.01, .01, .032, .032, .018, .018, .007, .007]
        else: #flood
            Qgoal=[.012, .012, .0242, .0242, .0318, .0318, .0422, .0422]
        indWhich = [nR+1, -1]    #if info on Q initial stages will be plotted, use tstep nR+1: when Q arrives to outlet.
            #m3/s as comes from f90. 18 is only in rising limb and 7 only in falling.
#        firstElseLast = [1,0,1,0]
        t_ej = []
        for k, qg in enumerate(Qgoal):
            absRelErrQ = np.absolute(Q[:,0] - qg)/qg
            ind,=np.where(absRelErrQ<.01)
#            print(f'30mar, len(ind)={len(ind)}, ind={ind}')
#            print(f'30mar, k{k}, firstElseLast[k]={firstElseLast[k]}')
            tejNew = ind[indWhich[k%2]]   # indWhich[ firstElseLast[k] ]  
            t_ej.append(tejNew)
#            print(f'30mar, t_ej={t_ej}')            
        xn = ['hsXanom', 'wXanom', 'gradNeg_hs', 'grad_wgw', 'vXanom', 'c', 'Fr']
        tup_ls = []
        for kk, varn in enumerate(xn):
            profs_choosen = np.zeros((len(Qgoal), nR))            
            for ii, ti in enumerate(t_ej):
                profs_choosen[ii, :] = eval(varn)[ti,:]
            tup2add = (varn, profs_choosen)
#            print(f'30mar, at varn{varn}, ti{ti}, profs_choosen[ii, :]={profs_choosen[ii, :]}')            
            tup_ls.append(tup2add)
            dic = dict(tup_ls)
        self.plotXcorr(xn, Qgoal, dic, case_n, fon1)
                
        if plot_whole==0:
            ii= int(prt_wdw[0])
            ifin=int(prt_wdw[1])
        else:
            ii= 1
            ifin= nt
        if self.plot_series:
            fon = fon1 + 't' + case_n
            if not os.path.exists(fon): os.makedirs(fon)
#            fc.cleanFolder(fon)        #no need, as every name is overwritten.
            for j in range(int(prt_wdw[2]), int(prt_wdw[3]) + 1):   #all channel reaches
                tr= t[ii:ifin]      
                D50j=D50[ii:ifin, j-1];   D50ss_j = D50c_ss[ii:ifin, j-1];  D84ss_j = D84c_ss[ii:ifin, j-1]
                cj=c[ii:ifin, j-1]; D84j=D84[ii:ifin, j-1]; jam_j = jam[ii:ifin, j-1]
                Qj= Q[ii:ifin, j-1]; hj=h[ii:ifin, j-1];
                vj=v[ii:ifin, j-1]; vfj=vf[ii:ifin, j-1];  
                wj=w[ii:ifin, j-1]; Sj=S[ii:ifin, j-1]; hsj= h_s[ii:ifin, j-1]; 
                hsfj=h_sf[ii:ifin, j-1]; wcj=w_c[ii:ifin, j-1]
                wbj=w_b[ii:ifin, j-1]; tsj= ts[ii:ifin, j-1]; tcj= tc84[ii:ifin, j-1]; tsmj= tsm[ii:ifin, j-1]
                Frj= Fr[ii:ifin, j-1]; asp_j= aspect[ii:ifin, j-1]; subm_j= subm[ii:ifin, j-1]
                pi_jk= pi[ii:ifin, j-1, :]; Mgtm_j = Mgtm[ii:ifin, j-1]; gam_j = gam[ii:ifin, j-1]
                #timeseries:
                self.series(fon, tr, cj, D50j, D84j, jam_j, D50ss_j, D84ss_j,
                    Qj, hj, vj, vfj, wj, Sj, hsj, hsfj, wcj, wbj, w_v[j-1], Frj, asp_j, subm_j, j)
                self.seriesSed(fon, dt_serie, tr, cumOut[:, j-1], Qj, braidInd[:, j-1], Sj, cj, tsj, tcj, tsmj, pi_jk,
                    Mgtm_j, gam_j, j)
        
                    
#        mn_an = -.1; mx_an = .1; tstep1 = 2000; tstep2 = 5000
#        fig,ax = plt.subplots(nrows = 2, ncols = 2)
#        ax[0,0].plot(wXanom[tstep1, :], hsXanom[tstep1, :], '.r') #formative flow
#        ax[0,0].plot(wXanom[tstep2, :], hsXanom[tstep2, :], '.k') #recesion flow
#        ax[0,0].plot(wXanom[3000, :], hsXanom[3000, :], '.m') #formative flow        
#        ax[0,0].plot(wXanom[4000, :], hsXanom[4000, :], '.',c='.4') #recesion flow                
#        ax[0,0].set_xlim(left=mn_an, right=mx_an)        
#        ax[0,0].set_ylim(bottom=mn_an, top=mx_an)
#        ax[0,0].set_xlabel('wXanom');         ax[0,0].set_ylabel('hsXanom')
#        ax[0,1].plot(grad_wgw[tstep1, :], gradNeg_hs[tstep1, :], '.r', label= f'Qfor, t#{tstep1}') #formative flow
#        ax[0,1].plot(grad_wgw[3000, :], gradNeg_hs[3000, :], '.m', label= f'Qfor, t#3000') #formative flow        
#        ax[0,1].plot(grad_wgw[tstep2, :], gradNeg_hs[tstep2, :], '.k', label = f'Qrec, t#{tstep2}') #recesion flow
#        ax[0,1].plot(grad_wgw[4000, :], gradNeg_hs[4000, :], '.',c='.4', label = f'Qrec, t#4000') #recesion flow        
#        ax[0,1].set_xlim(left=-.3, right=.5)
#        ax[0,1].set_ylim(bottom=-.06, top=.06)
#        ax[0,1].set_xlabel('grad_w');         ax[0,1].set_ylabel('gradNeg_hs')
#        ax[0,1].legend()        
#        ax[1,0].plot(wXanom[tstep1, :], vXanom[tstep1, :], '.r') #formative flow
#        ax[1,0].plot(wXanom[tstep2, :], vXanom[tstep2, :], '.k') #recesion flow
#        ax[1,0].plot(wXanom[3000, :], vXanom[3000, :], '.m') #formative flow        
#        ax[1,0].plot(wXanom[4000, :], vXanom[4000, :], '.',c='.4') #recesion flow                
#        ax[1,0].set_xlim(left=mn_an, right=mx_an)        
#        ax[1,0].set_ylim(bottom=mn_an, top=mx_an)
#        ax[1,0].set_xlabel('wXanom');         ax[1,0].set_ylabel('vXanom')
#        ax[1,1].plot(Fr[tstep1, :], c[tstep1, :], '.r') #formative flow
#        ax[1,1].plot(Fr[tstep2, :], c[tstep2, :], '.k') #recesion flow
#        ax[1,1].plot(Fr[3000, :], c[3000, :], '.m') #formative flow        
#        ax[1,1].plot(Fr[4000, :], c[4000, :], '.',c='.4') #recesion flow                
#        ax[1,1].set_yscale('log')
#        ax[1,1].set_xscale('log')        
#        ax[1,1].set_xlim(left=.1, right=2)
#        ax[1,1].set_ylim(bottom=1e-6, top=1e-3)        
#        ax[1,1].set_xlabel('Fr');         ax[1,1].set_ylabel('c')
#        fig.tight_layout()        
#        fn = 'ejXcorrel' + case_n
#        plt.savefig('doct_fluvial/hg_model_flume/tx/' + fn + '.png')
#        plt.close(fig)
                    
        mn_an = -.5; mx_an = .5
        tf = t[-1]
        xf = nR*dx
        
        #filters to highlight threshold values:
        filt_tol =.1
        Fr_msk = self.setThreshVals(Fr, filt_tol, 1)
        subm_msk = self.setThreshVals(subm, filt_tol, 4)
        jam_msk = self.setThreshVals(jam, filt_tol, 5)
        aspect_msk = self.setThreshVals(aspect, filt_tol, 50)
        braidInd_msk = self.setThreshVals(braidInd, filt_tol, .4)
        ttc_msk = self.setThreshVals(ttc, filt_tol, 1)
        ttm_msk = self.setThreshVals(ttm, filt_tol, 1)
        
        self.mapSaver(fon1, case_n, wXanom, 'wXanom', mn_an, mx_an, xf, tf)
        self.mapSaver(fon1, case_n, hXanom, 'hXanom', mn_an, mx_an, xf, tf)        
        self.mapSaver(fon1, case_n, FrXanom, 'FrXanom', mn_an, mx_an, xf, tf)
        self.mapSaver(fon1, case_n, vXanom, 'vXanom', mn_an, mx_an, xf, tf)        
        self.mapSaver(fon1, case_n, hsXanom, 'hsXanom', mn_an, mx_an, xf, tf)
        self.mapSaver(fon1, case_n, aspect_msk, 'aspect', 1, 100, xf, tf, log=True)
        self.mapSaver(fon1, case_n, subm_msk, 'subm', 1, 7, xf, tf)
        self.mapSaver(fon1, case_n, jam_msk, 'jam', 1e0, 1e2, xf, tf, log=True)
        self.mapSaver(fon1, case_n, ttc_msk, 'ttc', 1e-1, 1e1, xf, tf, log=True)
        self.mapSaver(fon1, case_n, ttm_msk, 'ttm', 1e-1, 1e1, xf, tf, log=True)
        outFin= cumOut[nt-2,nR-1] #-2 is feasible given timeSteps are thousands, 
                                  #so it is avoided the NaN set previously to print series
        titcum = f'VolSedOut = {int(round(outFin))}kg'
        self.mapSaver(fon1, case_n, Fr_msk, 'Fr', .5, 1.2, xf, tf, tit=titcum)
        cumOut_noZ = cumOut[cumOut<1e-1]=1e-1
        self.mapSaver(fon1, case_n, cumOut, 'cumOut', np.amin(cumOut), np.amax(cumOut), 
            xf, tf, log=True)
#        bmn = np.amin(braidInd[:-1,:]) 
#        bmx = np.amax(braidInd[:-1,:])
        bmn = 1e-1
        bmx = 1e0
        self.mapSaver(fon1, case_n, braidInd, 'braidInd', 
            bmn, bmx, xf, tf, log=False)   
            #avoids last row (t=tfin) which is NaN.
        self.mapSaver(fon1, case_n, w_b, 'wb', 0, .6, xf, tf)
        self.mapSaver(fon1, case_n, S, 'S', 0, .06, xf, tf)
        self.mapSaver(fon1, case_n, c, 'c', 1e-6, 1e-3, xf, tf, log=True)
        self.mapSaver(fon1, case_n, w, 'w', .2, 1, xf, tf)
        self.mapSaver(fon1, case_n, v, 'v', .3, 1.3, xf, tf)
        self.mapSaver(fon1, case_n, h, 'h', .03, .2, xf, tf)
        self.mapSaver(fon1, case_n, gradNeg_hs, 'gradNeg_hs', -.5, .5, xf, tf)
        self.mapSaver(fon1, case_n, grad_wgw, 'grad_wgw', -.5, .5, xf, tf)
        

#            #Hydraulic geometry of standard variables:
#            fig, ax = plt.subplots(nrows=2,ncols=1)
#            tmx = max(tr)
#            cl = [str(ti/(1.1*tmx)) for ti in tr]
#            lh=ax[0].scatter(Qj, hj, marker= 'x', c=cl)
#            lv=ax[0].scatter(Qj, vj, marker= '+', c=cl)
#            ax1 = ax[0].twinx()        
#            lw=ax1.scatter(Qj, wj, marker= '.', c=cl) #, label = '$w$'
#            ls = (lh, lv, lw)
#            labs = ('h', 'v', 'w')
#            leg = ax[0].legend(ls, labs, scatterpoints=1, loc='best', ncol=1, fancybox=True, fontsize=9)            
#            cj[cj<1e-8]=1e-8
#            D84j[D84j<1e-6]=1e-6
#            lc = ax[1].scatter(Qj, cj, marker= '.', c=cl)
#            ax2 = ax[1].twinx()
#            ld = ax2.scatter(Qj, D84j, marker= '+', c=cl)
#            ax[1].set_yscale('log')
#            ax2.set_yscale('log')            

#            ls = (lc, ld)
#            labs = ('c', '$D_{84}$')
#            leg = ax[1].legend(ls, labs, scatterpoints=1, loc='best', ncol=1, fancybox=True, fontsize=9)

#            ax[0].get_shared_x_axes().join(ax[0], ax[1])
#            ax[0].set_xticklabels([])           
#            ax[1].set_xlabel('$Q(m^3/s)$')
#            ax[0].set_ylabel('h(m) or v(m/s)')
#            ax[1].set_ylabel('c($m^3/m^3$)')
#            ax2.set_ylabel('$D_{84}(m)$')            
#            ax1.set_ylabel('w(m)')

##            ls2 = l12+l13
##            ls=l1+l4+l7
##            labs = [l.get_label() for l  in ls]; labs2 = [l.get_label() for l  in ls2]
##            leg = ax[0].legend(ls, labs, loc='best', ncol=1, fancybox=True, fontsize=9)
##            leg.get_frame().set_alpha(0.7)
##            leg = ax[1].legend(ls2, labs2, loc='best', ncol=1, fancybox=True, fontsize=9)
##            leg.get_frame().set_alpha(0.7)            
#            fn = 'HG_reach' + str(j)
#            plt.tight_layout()       
#            plt.savefig(f'doct_fluvial/hg_model_flume/{fn}.png')
#            plt.clf()
            
            
        
    def pct_dyn(self, dt, x):
        a = [(dt[i], x[i]) for i in np.argsort(x)]
        lenPerTuple = len(a[0]) #2 by definition of pct_dyn method args.
        lenTuple = len(a)
        for i in range(1, lenPerTuple+1):
            bi = [aj[i-1] for j, aj in enumerate(a)]
            b = bi if i==1 else np.stack((b, bi))
#        print(f'b{b}')            
        t = np.cumsum(b[0, :])
        tp = 1-t/t[-1]
        return tp, b[1, :]        
        
    def plt_ddt(self, flume_else_basin, ds, t, y, yd, dt, l, d2pi, d2pf, st, j, ymn):   
         #d2p: days to plot; ds: day-timeStep relation [day, ts].
        of = self.of         
        tp_y, ys = self.pct_dyn(dt, y)
        tp_dt, dts = self.pct_dyn(dt, dt)
#        ys = np.sort(y);  dts = np.sort(dt)
        yds = np.sort(yd)
        ny = len(y);  nyd = len(yd)        
#        imx = int(ny/20/365*d2p)   #jan8, 2020: 1mes
        imn_pair = ds[ds[:, 0]==d2pi]
        imx_pair = ds[ds[:, 0]==d2pf]
        imn = int(imn_pair[0][1]) - 1   #-1 for py indexes that start at 0
        imx = int(imx_pair[0][1]) - 1
        t_steps = np.arange(ny);  td_steps = np.arange(nyd)
        tpd = 1-td_steps/nyd    #tp = 1-t_steps/ny;  
        tr = t[imn:imx];  yr = y[imn:imx];  dtr = dt[imn:imx]
#        print(f'imn_pair, imx_pair ={imn_pair}, {imx_pair}')
#        print(f'imn_pair[0][1]={imn_pair[0][1]}')
#        print(f'imn{imn}, imx{imx}')
#        print(f'tr= {tr}')
#        print(f'yr= {yr}')
        td = np.arange(1, nyd+1);  tdr = td[d2pi-1 : d2pf-1+1];  ydr = yd[d2pi-1 : d2pf-1+1]
#        print(f'sizes in plot_cd: tpd{len(tpd)}, yds{len(yds)}, tp_y{len(tp_y)}, tp_dt{len(tp_dt)}, ys{len(ys)}, , dts{len(dts)}')
#        print(f'tp_y_mn{min(tp_y)}, tp_y_mx{max(tp_y)}')
#        print(f'tp_dt_mn{min(tp_dt)}, tp_dt_mx{max(tp_dt)}')

#        fig,ax=plt.subplots(nrows=2,ncols=1)
        fig = plt.figure()        
        grid = plt.GridSpec(2, 2) #, hspace=.4, wspace=.6)
                
        ax1 = fig.add_subplot(grid[0, :])
        l1 = ax1.plot(tr, yr, 'ok', label= l[1-1], markersize = 1.5);     ax1.set_ylabel(l[1-1] + ' (m3/s)'); 
        ax1.set_xlabel('time (s)')
#        ax1.plot(tdr, ydr, 'or', label= l[2-1], markersize = 3, markerfacecolor='none');
        axa = ax1.twinx()
        l2 = axa.plot(tr, dtr, '.r', label= l[3-1], markersize = 1.5);     axa.set_ylabel(l[3-1] + ' (s)')
        ls=l1+l2
        labs = [l.get_label() for l  in ls]        
        leg = ax1.legend(ls, labs, loc='best', ncol=1, fontsize=9) #ncol=2, , fancybox=True
        
        ax2 = fig.add_subplot(grid[1, 0]) if flume_else_basin == 0 else fig.add_subplot(grid[1, :])
        ax2.plot(tp_y, ys, 'ok', label = f'{l[1-1]}, Nsteps{ny}', markersize = 1.5)
#        ax2.plot(tpd, yds, '.r', label = f'{l[2-1]}', markersize = 1)
        ax2.set_xlabel('prob.exced');  ax2.set_ylabel(l[1-1] + ' [m3/s]')
        axb = ax2.twinx()
        axb.plot(tp_dt, dts, '.r', markersize = .3)
        if flume_else_basin == 1: axb.set_ylabel(l[3-1] + ' [s]')

        
        #as previos ax2 but no log scale in axis.
        if flume_else_basin == 0:
            ax3 = fig.add_subplot(grid[1, 1])
            ax3.plot(tp_y, ys, 'ok', label = f'{l[1-1]}, Nsteps{ny}', markersize = 1.5)
#            ax3.plot(tpd, yds, '.r', label = f'{l[2-1]}', markersize = 1)
            ax3.set_xlabel('prob.exced');  #ax3.set_ylabel(l[1-1] + ' [m3/s]')
            axb = ax3.twinx()
            axb.plot(tp_dt, dts, '.b', markersize = .3)
            axb.set_ylabel(l[3-1] + ' [s]')
        
        if flume_else_basin == 0:
            ax1.set_yscale('log')  
            axa.set_yscale('log')
            ax2.set_xscale('log')
            ax2.set_yscale('log')
            axb.set_yscale('log')        
            ax3.set_yscale('log')
        
#        if l[1-1] != 'Qhill':  
#            ax1.set_ylim(bottom = ymn)  #shows only sediment supply concentration >1e3ppm.
#            ax2.set_ylim(bottom = ymn);  ax3.set_ylim(bottom = ymn)

        path = f'{of}/{l[1-1]}_reach{j}.png'
        fig.suptitle(st);  
        fig.tight_layout()
        print(f'23mar, to save downscaled series and CDs, path: {path}')
        plt.savefig(path)
        plt.close()
