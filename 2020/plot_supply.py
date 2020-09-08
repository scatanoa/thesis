    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#dec27, 2019

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gsp
import sys
from scipy import optimize


class c:

    def __init__(self, n_days, output_folder, nDh, ncv):
        self.n_days = n_days;  self.output_folder =  output_folder;  self.nDh = nDh
        self.ncv = ncv

    def hgy_plotter(self, j, Pserie, Q, x1w, x2w, x3w, x4w, o1w):
        n_days = self.n_days; output_folder = self.output_folder
        t=np.arange(n_days)        
        fig,ax=plt.subplots(nrows=5,ncols=1)
        
        ax[0].set_title('mean(Q)='+str(round(Q.mean(),2))+', mean(x2w)='+str(x2w.mean())+\
        ', var(x2w)='+str(np.var(x2w)),fontsize=10)
        
    #        ax[0].plot(t,o1w,label='mean(o1w)='+str(round(o1w.mean(),1)))
    #        ax[0].set_ylabel("o1w[mm]")
        Pm = round(Pserie.mean(),1); CVP = round(np.std(Pserie)/Pm, 1)
        ax[0].plot(t, Pserie, label=f'mean{Pm}, CV{CVP}')
        ax[0].set_ylabel("P[mm]")
        
        ax[0].legend()
        ax[1].plot(t,x1w)
        ax[1].set_ylabel("x1w[mm]")
        ax[2].plot(t,x3w)
        ax[2].set_ylabel("x3w[mm]")
        for k in range(1,3+1):
            ax[k-1].get_xaxis().set_visible(False)
        ax[3].plot(t,x4w)
        ax[3].set_ylabel("x4w[mm]")
        ax[4].plot(t,o1w)
        ax[4].set_ylabel("o1w[mm]")
        ax[4].set_xlabel("t [days]",fontsize=10)        
        plt.savefig(f'{output_folder}/hill_water_stores_reach{j}.png')
        plt.close()


    def supply_x(self, Ptx, Q, Qsl, Qsl_jjAgr, hcr, hs, FS, slidsPerDay, coh_vect, ncv, i_cohMin, 
            nR, Abasin, slids_mx_perSubb):
        n_days = self.n_days;  output_folder = self.output_folder
        t=np.arange(n_days)
        fig,ax=plt.subplots(nrows=2,ncols=1)
        ii = np.arange(n_days-1 +1)        
        ii = ii/n_days
        for jj in range(ncv+1):
            ax[0].plot(t, slidsPerDay[:,jj-1])
            y4cd = np.sort(slidsPerDay[:,jj-1])[::-1]
            ax[1].plot(ii, y4cd)
        fig.tight_layout()
        plt.savefig(f'{output_folder}/slidsPerDay.png')
        plt.close()
        
        fig1,ax1=plt.subplots(nrows=2,ncols=1)
        cmn=np.amin(coh_vect)
        cmx=np.amax(coh_vect)
        for jj,coh in enumerate(coh_vect): #any subbasin j, p ej #7.
            y4cd = np.sort(FS[:, 7, jj]) #[::-1]
            y4cdr = y4cd[::-1]            
#            cohVal= coh*np.ones(n_days)
            r= coh/cmx
            col=(r,r,r)
#            print('23jul20; len x y z:', len(ii), len(y4cd), len(cohVal))
#            print('23jul20; jj, cohVal:', jj, cohVal)
#            im = ax1.scatter(ii, y4cd, c=cohVal, s=1, cmap='viridis', alpha=.7)
            ax1[1].plot(ii, y4cd, color=col, markersize=1, marker='o')
            ax1[1].set_xlabel('probability of non-exceedance')
            ax1[1].set_ylabel('Factor of Security')
            ax1[0].loglog(ii, y4cdr, color=col, markersize=1, marker='o')
            ax1[0].set_xlabel('probability of exceedance')
            ax1[0].set_ylabel('Factor of Security')            
#        ax1[0].set_ylim(bottom=.9, top=1.1)
#        ax1[0].set_xlim(left=1e-4, right=1e-1)
        ax1[1].set_ylim(bottom=.9, top=1.1)
        ax1[1].set_xlim(left=1e-4, right=1e-1)
        ax1[1].set_xscale('log');        
#        cbar = fig1.colorbar(im, ax=ax1)            
#        cbar.set_label(r'cohes $[N/m^2]$')
#            ax.plot(ii, y4cd)
        plt.grid(True, ls="--")                            
#        ax.set_yscale('log')
        fig1.tight_layout()       
        plt.savefig(f'{output_folder}/fg9_cdFS.png')
        plt.close()

        
#        fig1,ax1=plt.subplots(nrows=2,ncols=1) #25jul20
        fig = plt.figure()
        gs = gsp.GridSpec(2,2)
        ax0= fig.add_subplot(gs[0,0])
        ax1= fig.add_subplot(gs[1,0])
        ax2 = fig.add_subplot(gs[:,1])
        Q7 = Q[:, 7]
        for jj,coh in enumerate(coh_vect): #any subbasin j, p ej #7.
            r= coh/cmx
            col=(r,r,r)
            y4cd = np.sort( Qsl[:, 7, jj] ) [::-1]
            ax2.loglog(ii, y4cd, color=col, markersize=3, marker='o')
            y4cd = np.sort( Qsl[:, 7, jj] / Q[:, 7] ) [::-1] 
                #29jul20; slids_mx_perSubb so potential Qs to be weighted.
#            cohVal= coh*np.ones(n_days)

            ax0.loglog(ii, y4cd, color=col, markersize=3, marker='o')
        samplSubb = [7,17]                                
        for j in range(1, nR):
            Qslm_subb_j = Qsl[:, j-1, :] #time series for each cohesion.
            Qslm = np.mean(Qslm_subb_j, axis=1) #average along cohesion.
            c = Qslm / Q[:,j-1]
            y4cd = np.sort(c)[::-1]
            r= j/nR
            col=(r,r,r)            
            ax1.loglog(ii, y4cd, color=col, linewidth=.5)
        ax1.set_xlabel('probability of exceedance')
        ax0.set_ylabel('potential c[], per cohesion, j=7')
        ax2.set_ylabel('potential Qs [m3/s], per cohesion, j=7')
        ax2.set_xlabel('probability of exceedance')
        ax1.set_ylabel('actual c[], per subbasin j')
        Qslm_aver_coh = np.mean(Qsl, axis=2) #25jul20; transforms array from 3d to 2d.
        cont0=0
        ls0=[]
#        for i in range(1, n_days+1): #26jul20; ok, I discovered Q is not saved at last subbasin.
#            for j in range(1, nR+1):
#                if Q[i-1,j-1]==0: 
#                    cont0+=1                    
#                    ls0.append( (i,j) )
#        p0 = cont0/(n_days*nR)
#        print('ls0: ', ls0)
        c_slm_aver_j = np.mean(Qslm_aver_coh[:, :-1] / Q[:, :-1], axis=1)  #25jul20; transforms array from 2d to 1d.
        print(f'25jul20, Qslm_aver_coh.shape: {Qslm_aver_coh.shape}')
        print(f'25jul20, Q.shape: {Q.shape}')
        print(f'25jul20, np.amin(Q): {np.amin(Q)}')
#        Qslm = np.mean(Qslm_aver_j, axis=1)
        y4cd = np.sort(c_slm_aver_j)[::-1]
        Qslm_aver_j = np.mean(Qslm_aver_coh[:, :-1], axis=1)  #mean along space (subbasins)
        Qslmm = np.mean(Qslm_aver_j) #mean along time 
        rs = Qslmm*86400*365*2.65/Abasin #25jul20; t/km2/yr
        print('25jul20, int(round(Abasin, 0)):', int(round(Abasin, 0)) )
        lb = 'AggrSupply' + f'\nA={ int(round(Abasin, 0)) }km2' + f'\nsedYield={int(round(rs, 0))}t/km2/yr'
        ax1.loglog(ii, y4cd, color='r', markersize=3, marker='o', label=lb)
        ax1.legend(fontsize=8, loc='best')
        fig.tight_layout()
        plt.savefig(f'{output_folder}/fg10_cdQs.png')
        plt.close()
        
        fig,ax = plt.subplots() #25jul20
        linStyl = ['.k', '+r']
        jl = [7,17]
        for ind,j in enumerate(jl):
            lb = f'subbasin {j}'
            ax.plot(Qslm_aver_coh[-365:, j-1] / Q[-365:, j-1], linStyl[ind], 
                markersize=3, label=lb)
            ax.set_yscale('log')
        ax.set_xlabel('time (days)')
        ax.legend()
        ax.set_ylabel('c (m3/m3)')
        plt.savefig(f'{output_folder}/fg11_ctj.png')
        plt.close()        
        
        fig,ax=plt.subplots(nrows=1, ncols=2) #rating curve, so neither uncorrel (pure stoch slids), 
            #nor pure determinism (as rating curves for morphodyn inputs).
        ax[0].loglog(Q[:,7], Qsl_jjAgr[:,7], 'ok', markersize=1) #any subbasin j, p ej #7.
        ax[1].loglog(Q[:,17], Qsl_jjAgr[:,17], 'or', markersize=1) #any subbasin j, p ej #17.
        Q7m = np.mean(Q[:,7]); Q17m = np.mean(Q[:,17])
        Qm = [Q7m, Q17m]
        ax[0].set_ylabel('Qs (m3/s)')        
        for i in [0,1]: #concentration reference lines, to infer if flood of debris flow.
            ax[i].set_xlabel('Q (m3/s)')
            ax[i].loglog( [ Qm[i], Qm[i] ], [1e-5, 1e0], '--b')
            txt = f'Qm={ round(Qm[i],3) }'
            ax[i].text(Qm[i], 1e0, txt, ha='center', c='b')
            for c in [1e-4, 1e-3, 1e-2, 1e-1, 1e0]:
                qq = np.array([1e-1, 2e1])
                qqs = c*qq
                ax[i].loglog(qq, qqs, label= f'{c}')
        ax[0].legend(title='c:', ncol=2, fontsize=7) #, loc=4
        fig.tight_layout()
        plt.savefig(f'{output_folder}/fg8_ratsQQs.png')
        plt.close()        

        p = [.1,.25,.4,.55,.7,.85] #portion of whole basin suffering slides.
#        col = ['c','b','g','r','m','k']
        p = np.array(p)
        nsb = p*nR  #number of subbasins
        nsbi = nsb.astype(int)
#        print(f'nsbi: {nsbi}')
#        print(f'i_cohMin: {i_cohMin}')
#        i_cohInterm= int(ncv/2)
        t_malamud = []  #6jul20; p266r: sample days to get spatial variation of slide volumes
        tP=[]
        for y in np.nditer(nsbi):
            ind, = np.where( slidsPerDay[: , i_cohMin] == y) #9jul20; i_cohInterm: few events
            n=100 #whole bulk cloud of black lines to show scaling break.
            nP=4 #to show colored multized dots showing to trend with P.                        
            ii = ind[-n:]
            iiP = ind[-nP:]

#            iiP = np.zeros(n)
#            nnP = n-nP
#            iiP[nnP:] = 1 #last t cases are to be P (rain) analized.
            if len(ind)>0: 
                t_malamud.extend(ii)
                tP.extend(iiP)
                    #use last 2 time 'y' subbasins showed slide with weakest cohesion.
                            
#        fig,ax=plt.subplots(nrows=2,ncols=2)
        x=np.arange(nR)
        xx = x/nR
        y4cd_allCoh_ls = []; r4cd_allCoh_ls = []; x3_ls = []; nx_ls=[]
        for i,t in enumerate(t_malamud): #9jul20; each t for sample of %basin eroded p
            yu_all_coh = []; ru_all_coh = []
            for jj in range(ncv+1):
                y_unsort = Qsl[t,:,jj-1] #*86400 / hs[t,:,jj-1] #27jul20;
                r_unsort = Ptx[t,:]
                yu_all_coh.extend(y_unsort); ru_all_coh.extend(r_unsort)
#                y4cd = np.sort(y_unsort)[::-1]
#                ax[0,0].plot(x, y_unsort, c=col[i])                
#                ax[1,0].plot(xx, y4cd, c=col[i])            
#            print(f't_malamud #{i}: day {t}')
#            print(f' yu_all_coh:{yu_all_coh}')            
#            print(f' ru_all_coh:{ru_all_coh}')
            yr_u = [yu_all_coh, ru_all_coh]
                #9jul20; yr means sed_yield | rain; u means unsorted
            yr_u = np.array(yr_u)
            yr_u = yr_u.T         
#            print(f' yr_u:{yr_u}')
            yr_s = yr_u[yr_u[:,0].argsort()] #9jul20; sort by 1st column: sedYield, preserving asociated rain.
#            print(f' yr_s:{yr_s}')
            yr_s_noZ = yr_s[np.logical_not(yr_s[:,0]==0)] #survive rows where 1st column (yield) is not null.
#            print(f' yr_s_noZ:{yr_s_noZ}')            
#            sys.exit('rev yr')    #9jul20, 13:30; works only until here.

#            y_all_coh = np.array(y_all_coh)
#            y4cd = np.sort(y_all_coh)[::-1]            
#            y4cd = y4cd[y4cd != 0]
#            y4cd_allCoh_ls.append(y4cd)
            x= yr_s_noZ[:,0]
            y4cd_allCoh_ls.append( x.tolist() )
            r4cd_allCoh_ls.append( yr_s_noZ[:,1].tolist() )
#            nx = len(y4cd)
            nx = len(x)
            x3= np.arange(nx)
            nx_ls.append(nx)
            x3 = x3/nx #9jul20; duration curve become biased to left; notably when small sample.
            x3_ls.append(x3)
#            print(f'y_all_coh: {y_all_coh}')
#            print(f' noZeros SORT; y_all_coh: {y4cd}')
#        ax[0,0].set_yscale('log')
#        ax[1,0].set_xscale('log'); ax[1,0].set_yscale('log')
#        ax[1,1].set_xscale('log'); ax[1,1].set_yscale('log')
#        ax[1,1].set_xlim(left=1e-2, right=1e0)
#        ax[1,1].set_ylim(bottom=1e-5, top=1e-0)
#        fig.tight_layout()
#        plt.savefig(f'{output_folder}/malamud.png')
#        plt.close()
        
        xl = []; yl = []; zl = []; sl=[]
        fig,ax=plt.subplots(nrows=2, ncols=2)
#        print(f'len(t_malamud):{len(t_malamud)}; len(tP):{len(tP)}')
#        sys.exit('rev tP, t_malamud')    #9jul20, 15:30
        Qs_t_co1_ls = [];         Qs_t_co2_ls = [];         
        P_t_ls = [];
        Qs_t_co3_ls = [];         Qs_t_co4_ls = []; Qs_t_ls = []
        cn = [0, int(ncv/3), int(2*ncv/3), ncv-1]   # cohesion #
        cnc= ['oc', 'ob', 'og', 'or']
        nc = len(cnc)
        n_tev = len(t_malamud)
        Qs_t = np.zeros( (n_tev, nc) ); hcr_t = np.zeros( (n_tev, nc) )
        hs_t = np.zeros( (n_tev, nc) )
        Qs_xt = np.zeros( (n_tev, nR, nc) )
        hs_xt = np.zeros( (n_tev, nR, nc) )
        hcr_xt = np.zeros( (n_tev, nR, nc) )                
        P_xt_ev = np.zeros( (n_tev, nR) )
        xset = []; yset = []; ntm = len(t_malamud)
        a_ls = []; b_ls = []; R2_ls = []
        for i,t in enumerate(t_malamud):
            x = x3_ls[i]; 
            x = [1-xi for xi in x] #from Pr(no exc), as by argsort, to Pr(exc)
            x=np.array(x)
#            dx=-np.diff(x)            
            y = y4cd_allCoh_ls[i]
            y = np.array(y)*86400/slids_mx_perSubb #28jul20; INDIVIDUAL slids spectra.
            y = y.tolist()
#            Qsm_f=np.sum(y[1:]*dx) #12jul20; integrating area under pdf.
                #12jul20; freq domain deprecated as it considers only no0 values, while
                    #rain vs sedYield relation must show spatially averaged event yield.
            for jj,c in enumerate(cn):
                Qs_t[i, jj] = np.mean( Qsl[t,:,jj] ) #12jul20.
                hcr_t[i, jj] = np.mean( hcr[t,:,jj] ) #22jul2020.
                hs_t[i, jj] = np.mean( hs[t,:,jj] ) #22jul2020.
                Qs_xt[i, :, jj] = Qsl[t,:,jj] #22jul20.
                hs_xt[i, :, jj] = hs[t,:,jj] #23jul20.                                                
                hcr_xt[i, :, jj] = hcr[t,:,jj] #23jul20.                
#            Qsm=np.mean(Qsl[t,:,cn[1]]);             Qs_t_co2_ls.append(Qsm)
#            Qsm=np.mean(Qsl[t,:,cn[2]]);             Qs_t_co3_ls.append(Qsm)
#            Qsm=np.mean(Qsl[t,:,cn[3]]);             Qs_t_co4_ls.append(Qsm)
#            print(f'y[1:]: {y[1:]}; dx: {dx}')
#            print(f'Qsm_f: {Qsm_f}; Qsm_t: {Qsm_t}')
#            sys.exit('Qsm_freq=Qsm_tserie?')
            s = r4cd_allCoh_ls[i]; 
            s = [si+1 for si in s] #to avoid 0size marker in slides with 0 rain.
            P = np.mean(Ptx[t,:])
            P_t_ls.append(P)
            P_xt_ev[i, :] = Ptx[t,:]
            z = np.ones(nx_ls[i])*P
            z = z.tolist()
            ax[0,0].plot(x, y,'-k', linewidth=.3)
#            if len(y)>1 and np.std(y)>0: #28jul20; if not all y values are the same.
            cont_yValUniq = 0
            if len(y)>5: #fit to power regression only events with more than five slide sizes.
                yprev = y[0]
                for yj in y[1:]:
                    if yj>yprev:
                        cont_yValUniq+=1
                    yprev = yj
                if cont_yValUniq>5:
                    print( f'28jul20; event {i} of {ntm}')
                    print(' ny: ', cont_yValUniq)
                    print(' x: ', x)
                    print(' y: ', y)
                    pars, pars_cov = optimize.curve_fit(self.pow_fn, x, y, p0=[1e1, -.1]) 
                        #28jul20; exponent <0.
                    y_powEval = self.pow_fn(x, pars[0], pars[1])
                    R2 = self.get_r2(x, y, y_powEval, pars)
                    ax[0,1].plot(x, y_powEval)
                    a = pars[0]; b = pars[1]
                    print(' a,b,R2:', round(a, 1), round(b, 1), round(R2, 2) )
                    a_ls.append(a); b_ls.append(b); R2_ls.append(R2)            
            xset.extend(x)
            yset.extend(y)
#            print(f't_malamud #{i}: day {t}')
#            print(f' malamud x, len{len(x)}: {x}')
#            print(f' malamud y, len{len(y)}: {y}')
#            print(f' malamud z, len{len(z)}: {z}')
#            print(f' malamud s, len{len(s)}: {s}')
#            sys.exit('rev malam colorP')
            ind, = np.where( np.array(tP) == t )
            if len(ind)>0: #analize rain only for reduced sample per portion of subbasins failed.
                xl.extend( x )       
                yl.extend( y )     
                sl.extend( s )
                zl.extend(z)
#        ax[0,0].plot(xset, yset, 'or', markersize=1)
        pars, pars_cov = optimize.curve_fit(self.pow_fn, xset, yset, p0=[1e1, -.1]) #28jul20; exponent <0.
        y_powEval = self.pow_fn(xset, pars[0], pars[1])
        R2 = self.get_r2(xset, yset, y_powEval, pars)
#        lb = f'pars:{pars}, R2:{R2}'
#        ax[0,0].plot(xset, y_powEval) #, label=lb)
#        ax.plot([5e-2, 5e-1],[5e4, 1e3],'-r')
        lbVol = 'slide volume ' + r'$V_s$ $[m^3]$'
        ax[0,0].set_ylabel(lbVol)
        ax[0,0].set_xlabel('exceedance probability F')
        lb = 'regressions ' + r'$ F \sim V_s^{-\theta} $' + ' per event'
        ax[0,1].set_xlabel(lb)
        for i in [0,1]:
            ax[0,i].set_xscale('log');         ax[0,i].set_yscale('log')
#            ax[0,i].set_xlabel('probability of exceedance');         
            ax[0,i].set_xlim(left=1e-2, right=1e0);     
            ax[0,i].set_ylim(bottom=1e-3, top=1e3)
        expo = -1/np.array(b_ls)     
        ax[1,0].hist( -1/np.array(b_ls) ); #28jul20; p272ur: F=pss~Vol^-(1/beta) as F is x-axis.
        ax[1,0].set_xlim(left=0); 
        lb = 'exponent ' + r'$\theta$' + f', mean={ round(np.mean(expo), 1) }'
        ax[1,0].set_xlabel(lb); 
        ax[1,0].set_ylabel('count')
        ax[1,1].hist( np.array(R2_ls) ); ax[1,1].set_xlabel(r'$R^2$' + ' of power law regressions')
        fig.tight_layout()        
        plt.savefig(f'{output_folder}/fg6_malamud3.png')
        plt.close()
    
        fig,ax=plt.subplots(nrows=2,ncols=2)
        ny=0; nyx=0
        for jj,col in enumerate(cnc): #vertical hztal y diagonal collapse.
            den=hcr_t[:, jj]
            den3 = den**3
            x = hs_t[:, jj]/den
            y = Qs_t[:, jj]/den3            
            coh_value = coh_vect[ cn[jj] ]
            lb = str(int(coh_value)) + r' $N/m^2$'
            #22ju20; spatially averaged, to test marc18, which failured 
            ax[0,0].plot(P_t_ls, y, col, label=lb)
            ax[0,1].plot(x, y, col, label=lb)
            ny += len(y)
            yx = Qs_xt[:, :, jj]
            xx = hs_xt[:, :, jj]
            z= hcr_xt[:, :, jj]    
            #22ju20; spatially distributed, to test marc19.
            for j in range(1, nR+1):
                denx = z[:, j-1]
                denx3 = denx**3
                yp = yx[:, j-1] #/ denx3 #23jul20
                ax[1,0].plot(P_xt_ev[:, j-1], yp, col, label=lb, markersize=1)
                xxp = xx[:, j-1] #/ denx #23jul20
                ax[1,1].plot(xxp, yp, col, label=lb, markersize=1)
                nyx += len(yp)
#            ax[1,1].plot(x, yx, col, label=lb)            
        for i in [0,1]:
            for j in [0,1]:
                ax[i,j].set_xscale('log'); ax[i,j].set_yscale('log')
        ax[1,0].set_xlabel('P_day [mm]');         
        lb2= f'nt: {ny}'
        lb = r'$Q_s/h_{cr}^3$, ' + lb2  #23jul20; n data per cohesion value.
        ax[0,0].set_ylabel(lb) #$[m^3/s]$
        lb = r'$Q_s/h_{cr}^3$, ' + f'nxt: {nyx}'  #23jul20; n data per cohesion value.
        ax[1,0].set_ylabel(lb) #$[m^3/s]$        
        ax[0,1].set_xlabel('hs/hcr')      
        ax[0,0].legend()
#        ax[1].plot()
        fig.tight_layout()  
        plt.savefig(f'{output_folder}/fg7_rs_vs_P_Xaggr.png')
        plt.close()
        
        fig1,ax1=plt.subplots()
        im = ax1.scatter(xl, yl, s=sl, c=zl, cmap='viridis', alpha=.7) #c=col[i]
        ax1.legend( *im.legend_elements('sizes', num=5),
            loc="lower left", title=r'$Subbasin P_{day}$ [mm]')
            #https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/scatter_with_legend.html
        cbar = fig1.colorbar(im, ax=ax1)
        cbar.set_label(r'Basin $P_{day}$ [mm]')
        ax1.set_xscale('log');         ax1.set_yscale('log')
        ax1.set_xlabel('probability of excedence');         ax1.set_ylabel(lbVol)
        ax1.set_xlim(left=1e-2, right=1e0); ax1.set_ylim(bottom=1e-3, top=1e3)
        plt.savefig(f'{output_folder}/fg6_malamud2.png')
        plt.close()
        
    def pow_fn(self, x, a, b):
        return a*x**b       
        
    def get_r2(self, xdata, ydata, yEval, pars):
        res = ydata - yEval
        ss_res = np.sum(res**2)
        ss_tot = np.sum( ( ydata - np.mean(ydata) )**2 )
        rat = ss_res/ss_tot
        return 1-rat

    def supply_plotter(self, j, A_hs, FS, hcr, hmx, gm, pAslid, zsg, hw, hwg, o2w, o3w, ch, ch_wash, ch_slid, 
        P, Q, xsh, wash_sh, slid_sh, qs_mx_teor, q, asp, As, ncv, slidFlg):
        n_days = self.n_days;  output_folder = self.output_folder;  nDh = self.nDh; ncv = self.ncv
        #plot fluxes
        t=np.arange(n_days)     
        Qm = np.mean(Q)
        cn = [4, int(ncv/2), ncv-1]   # cohesion #
        col = ['r','b','c']
        mk = ['or','ob','oc']
        ctp = len(cn)   #number of cohesion cases to plot
        for jj in range(1,ctp+1):
            cnjj = cn[jj-1]; mkjj = mk[jj-1]
            FSjj = FS[:, cnjj]
            fig,ax = plt.subplots(nrows=2,ncols=2)
            ax[0,0].plot(P,hw[:, cnjj], mkjj, label = f'coh{cnjj}', markersize=1)
            ax[0,0].set_xscale('log');         ax[0,0].set_yscale('log'); 
            ax[0,0].set_xlabel('P_day [mm]');             ax[0,0].set_ylabel('hw [m]')
            ax[1,0].plot(hw[:, cnjj], FSjj, mkjj, label = f'coh{cnjj}', markersize=1)
            ax[1,0].set_ylabel('Fact.Secur');             ax[1,0].set_xlabel('hw [m]')                       
            ax[1,0].set_xscale('log');         ax[1,0].set_yscale('log');
            ax[1,0].set_ylim(top = 5)
            ax[0,1].plot(P, FSjj, mkjj, label = f'coh{jj}', markersize=1)
            ax[0,1].set_ylabel('Fact.Secur');             ax[0,1].set_xlabel('P_day [mm]')
            ax[0,1].set_xscale('log');         ax[0,1].set_yscale('log');
            ax[0,1].set_ylim(top = 5)
            ax[1,1].plot(zsg[:, cnjj], FSjj, mkjj, label = f'coh{cnjj}', markersize=1)
            ax[1,1].set_ylabel('Fact.Secur');             ax[1,1].set_xlabel('hs [m]')
            ax[1,1].set_yscale('log'); ax[1,1].set_xscale('log');
            ax[1,1].set_ylim(top = 5)
            fig.tight_layout()
            plt.savefig(f'{output_folder}/fg4_FSrats{j}_coh{cnjj}.png') #see ch2 phd fig definition in phdSCnotes 268.
            plt.close()
        ind_ini=[3,3,0]
        for jj in range(1,ctp+1):
            ii = ind_ini[jj-1]
            cnjj = cn[jj-1]; mkjj = mk[jj-1]; 
            coljj = col[jj-1]
            fig = plt.figure(tight_layout=True)
            gs = gsp.GridSpec(2,3)
            ax = fig.add_subplot(gs[0, :])
            tt = t[ii:]/365
            ax.plot(tt, zsg[ii:, cnjj], mkjj, markersize=1, label='zsg')
            FSjj = FS[:, cnjj]
            ax1=ax.twinx()
            ax1.plot(tt, FSjj[ii:], 'ok', markersize=1, label='FS')
            ax.set_xlabel('time [yr]'); 
            ax.set_ylabel('hs [m]', color=coljj)
            ax1.set_ylabel('Factor of security')
            ax.tick_params(axis = 'y', labelcolor=coljj)            
#            ax.legend()                        
            asp_jj = asp[:, cnjj]
            As_jj = As[:, cnjj]
            ind, = np.where(slidFlg[:, cnjj]>0)
            t_slid = np.array(ind)
            Tr_sl_jj = np.diff(t_slid)/365 #9jul20; day to [yr]
            vl = [Tr_sl_jj[ii:], As_jj[ii:], asp_jj[ii:]]
            nl = ['return period (yr)', 'scar area (m2)', 'aspect= L/b']
            for i in range(3):
                ax = fig.add_subplot(gs[1, i])
                vli = vl[i]
                vli_noNull = vli[vli>0]
                if i==1: #As
#                    if jj==ctp:
#                        print(f'vli_noNull: {vli_noNull}')
#                        sys.exit('9jul20; check As sample')
#                    As_mn = np.amin(vli_noNull);  As_mx = np.amax(vli_noNull)
#                    bn=np.logspace(np.log10(As_mn), np.log10(As_mx), 10)
#                    ax.hist(vli_noNull, bins=bn) #log=True,
                    ax.hist(vli_noNull) #log=True,
#                    ax.set_xscale('log')
                    ax.grid(which='minor', axis='x')
                else:                    
                    ax.hist(vli_noNull)          #, log=True
                ax.set_xlabel(nl[i])
                if i==0: ax.set_ylabel('count') #left subplot.            
            plt.savefig(f'{output_folder}/fg5_slidHists_rch{j}_coh{cnjj}.png')
            plt.close()

            
                
        
#        fig,ax=plt.subplots(nrows=3,ncols=2)
#        l1=ax[0,0].plot(t,FS, 'or', label = 'FS', markersize=1)
#        l2=ax[0,0].plot(t, gm/1e4, 'oc', label = 'gm [ton/m3]', markersize=1)
#        ax[0,0].set_ylabel("gm [ton/m3], FS")
#        ax[0,0].set_ylim(bottom = .9)
#    #            ax[0].set_yscale('log')
##        ax2=ax[0,0].twinx()
##        l3=ax2.plot(t, pAslid, 'ok', label = 'pAslid', markersize=1)
##        ax2.set_ylabel("pAslid")
#        ls = l1+l2#+l3
#        labs = [l.get_label() for l  in ls]        
##        leg = ax2.legend(ls, labs, loc='best', fancybox=True) #ncol=2, , ncol=4, , fontsize=8                
##        for jj in range(1,ncv+1):
##            ax[1].plot(t,zsg[:, jj-1],'k',label='zsg')
##            ax[1].plot(t,hcr[:, jj-1],'or',label='hcr', markersize=1)
##        for jj in range(1,ncv+1):
##            ax[1].plot(t,zsg[:, jj-1],'k',label='zsg')
##            ax[1].plot(t,hcr[:, jj-1],'or',label='hcr', markersize=1)
#        ax[1,0].plot(t,hw[:, 0],':b',label='hw')
#        ax[1,0].plot(t,zsg[:, 0],'k',label='zsg')
#        ax[1,0].plot(t,hcr[:, 0],'or',label='hcr', markersize=1)
#        ax[1,0].plot(t,hmx[:, 0],'oc',label='hcr', markersize=1)
##        ax[1].plot(t,hwg,':c',label='hwg')
##        leg = ax[1].legend(loc=3,fancybox=True, fontsize=10)
##        leg.get_frame().set_alpha(0.7)        
#        ax[1,0].set_ylabel("zsg: slide layer [m]") #,fontsize=6
#        ax[1,0].get_shared_x_axes().join(ax[0,0], ax[1,0])
#        ax[1,0].set_xticklabels([])
#        ax[1,0].set_xlabel('time (days)')       
#        
#        ax[2,0].plot(t,hw[:, 1],':b',label='hw')
#        ax[2,0].plot(t,zsg[:, 1],'k',label='zsg')
#        ax[2,0].plot(t,hcr[:, 1],'or',label='hcr', markersize=1)
#        ax[2,0].plot(t,hmx[:, 1],'oc',label='hcr', markersize=1)        
#        
#        ax[0,1].plot(t,hw[:, 2],':b',label='hw')
#        ax[0,1].plot(t,zsg[:, 2],'k',label='zsg')
#        ax[0,1].plot(t,hcr[:, 2],'or',label='hcr', markersize=1)
#        ax[0,1].plot(t,hmx[:, 2],'oc',label='hcr', markersize=1)        
#        
#        ax[1,1].plot(t,hw[:, 3],':b',label='hw')
#        ax[1,1].plot(t,zsg[:, 3],'k',label='zsg')
#        ax[1,1].plot(t,hcr[:, 3],'or',label='hcr', markersize=1)        
#        ax[1,1].plot(t,hmx[:, 3],'oc',label='hcr', markersize=1)                
#        
#        ax[2,1].plot(t,hw[:, 4],':b',label='hw')
#        ax[2,1].plot(t,zsg[:, 4],'k',label='zsg')
#        ax[2,1].plot(t,hcr[:, 4],'or',label='hcr', markersize=1)                                        
#        ax[2,1].plot(t,hmx[:, 4],'oc',label='hcr', markersize=1)
#        
#        fig.tight_layout()        
#        plt.savefig(f'{output_folder}/soilFailures_reach{j}.png')
#        plt.close()


        
        fig,ax=plt.subplots(nrows=3,ncols=1)
        ax[0].plot(t,ch,'k',label='total')
        ax[0].plot(t,ch_wash,'b',label='wash')
        ax[0].plot(t,ch_slid,'r',label='slide')        
        ax[0].set_ylabel("ch [m3/m3]")
        ax[0].set_yscale('log')
        ax[0].set_ylim(bottom=1e-2, top=1e1)
        ch_slid_mean='%1.1e'%np.mean(ch_slid)
        ch_wash_mean='%1.1e'%np.mean(ch_wash)
        ax[0].set_title(f'means: ch_slid= {ch_slid_mean}, ch_wash= {ch_wash_mean}')
        leg = ax[0].legend(loc=2,fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.7)
        cl = ['b', 'r', 'm', 'k', 'g']
        for kk in range(1, nDh + 1):
            ax[1].plot(t, slid_sh[:, kk-1], color = cl[kk-1], label = str(kk))
            ax[2].plot(t, wash_sh[:, kk-1], color = cl[kk-1], label = str(kk))
        ax[1].set_ylabel("slid_sh [m3/s]")            
        ax[1].set_yscale('log')
        ax[1].set_ylim(1e-4*Qm)
        leg = ax[1].legend(loc=2,fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.3)                        
        ax[2].set_ylabel("wash_sh [m3/s]")            
        ax[2].set_yscale('log')
        ax[2].set_ylim(1e-4*Qm)
#        ax[2].set_ylim(bottom = 1e-6)
        leg = ax[2].legend(loc=2,fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.3)            
        plt.savefig(f'{output_folder}/Qs_hill_reach{j}.png')
        plt.close()            
        
        #plot water flow
        Qspssp = (o2w + o3w)*A_hs/86400              
        fig,ax=plt.subplots(nrows=3,ncols=2)    #one image and file per run, no overlap        
        lb = ['P_day [mm]', 'Qhill_day [m3/s]']
        for k in [0,1]:
            lbk = lb[k]
            if k==0:
                Y = P
            else:
                Y = Q
            ax[0,k].plot(Y,'k')    #kind='bar'
            ax[0,k].set_ylabel(lbk)
            ii=np.arange(n_days-2 +1)
            ii=ii/n_days     #array between 0 and 1
        #            Q = Q[:-1]  #erase last value so size is compatible with ii for plotting.        
            Q4cd = np.sort(Y)[::-1]
            Q4cd = Q4cd[:-1]
            ax[1,k].plot(ii, Q4cd, 'k', label = 'Q')
            if k==1: 
                Qspssp4cd = np.sort(Qspssp)[::-1] 
                Qspssp4cd = Qspssp4cd[:-1]
                compl_pQbase = int(round(100*np.mean(Qspssp4cd)/Qm, 0))                
                ax[1,k].plot(ii, Qspssp4cd, '--k', label = f'QminusQb: {compl_pQbase}%')                
                ax[2,k].plot(ii, Qspssp4cd, '--k', label = f'QminusQb: {compl_pQbase}%')                
                ax[1,k].legend()  
                    #as Qm~40% exc.prob and Qbf around 1m3/s/km2 I am tranquil with the shape of the CDQ 
                        #(duration curve).                
            ax[1,k].set_ylabel(lbk)
#            ax[1,k].set_xlabel('prob_exced')
    #        ax[1].set_yscale('log')
    #        ax[1].set_xscale('log')
            ax[2,k].plot(ii, Q4cd, 'k', label = 'Q')
            ax[2,k].set_yscale('log')            
            if k==1: 
                ax[2,k].legend()
                ax[2,k].set_ylim(bottom=1e-2, top=1e1)
#                majors = [1e-2, 1e-1, 1e0, 1e1]
#                ax[2,k].yaxis.set_major_locator(ticker.FixedLocator(majors))
            ax[2,k].set_ylabel(lbk)
            ax[2,k].set_xlabel('prob_exced')
            ax[2,k].set_xscale('log')        
        CVQ=np.std(Q)/Qm; Qm=int(Qm*1e3); str_T = str(int(n_days))
        plt.tight_layout()
        plt.savefig(f'{output_folder}/Qh_reach{j}_Qm{Qm}lps_10CV{str(int(10*CVQ))}.png')
        plt.close()

        #plot states
        l = ["0.1mm", "1mm", "10mm", "100mm", "1000mm"]
        fig,ax=plt.subplots()
        for kk in range(1, nDh+1):
            ax.plot(t,xsh[:, kk-1], color = cl[kk-1], label = l[kk-1])            
        ax.set_ylabel("xsh [m]")
        leg = ax.legend(loc=2,fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.7)
    #        ax[0,1].legend()            
        plt.savefig(f'{output_folder}/Xs_hill_reach{j}.png')
        plt.close()            
        
#        #plot rating curve
#        fig,ax=plt.subplots()
##        ax.plot(washSlop**(1.5-2.)*qs_Rick01/.058, qs_wash*3600*2650, 'ok', markersize=.5)    #fig6_Rick01 = fig12_Zimm13.
##        ax.set_xscale('log')
##        ax.set_yscale('log')
##        ax.set_xlabel(" (q-qc)S^1.5 [m3/s/m]")
##        ax.set_ylabel("qs [kg/m/h]")        
#        ax.plot(q, ch_wash, 'ok', markersize=.5)    #fig6_Rick01 = fig12_Zimm13.
#        ax.set_xscale('log')
#        ax.set_yscale('log')
#        ax.set_xlabel(" q [m3/s/m]")
#        ax.set_ylabel("c wash")
#        ax.set_ylim(bottom=1e-4)
#        plt.savefig(f'{output_folder}/HillRatingCurve_reach{j}.png')
#        plt.close()
