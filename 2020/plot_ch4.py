    #!/usr/bin/env python
# -*- coding: UTF-8 -*-

#3sep, 2020. By: Santiago Cataño Alvarez.
#Attempt to compare num modeled patterns of both ubc and fld experiments. See p280dl notes.

import numpy as np
import matplotlib.pyplot as plt
import sys

class a:

    def __init__(self):
        x=1

    def plts(self):
        print('In plot_ch4.py/plts')
        of = 'doct_fluvial/hg_model_flume/ch4phd_bothFlumes'
        var_l = ['D84c', 'h_s', 'w', 'Q', 'c', 'ts_c', 'tc_c']
        nv = len(var_l)
        nx = 2

        flumeName = 'ubc'
        xl = [2, 16] #1m from upstream and downstream boundary condition.        
        nt_ubc = np.savetxt(f'{of}/{flumeName}_nt.txt')
        V_ubc = np.zeros((nt_ubc, nv, nx))        
        for ii, vr in enumerate(var_l):
            for jj, j in xl:
                V_ubc[:, ii, jj] = np.loadtxt(f'{of}/{flumeName}_{vr}_j{j}.txt')

        flumeName = 'fld'
        xl = [2, 23] #1m from upstream and downstream boundary condition.        
        nt_fld = np.savetxt(f'{of}/{flumeName}_nt.txt')
        V_fld = np.zeros((nt_fld, nv, nx))        
        for ii, vr in enumerate(var_l):
            for jj, j in xl:
                V_fld[:, ii, jj] = np.loadtxt(f'{of}/{flumeName}_{vr}_j{j}.txt')

        fig, ax = plt.subplots(nrows=2,ncols=2)
        ax.scatter(V_ubc[:, 5, 1], V_ubc[:, 4, 1]) #c vs tao*
        fnfull = fon + fn
        plt.savefig(fnfull)
        plt.close('all')        
        
        
    def dQ_serie_partit(self):
        
        
    def ratings(self, fon, Q, h, v, w, cs, vf, j):
        print(f'In plots2f.py/ratings; printing j={j}')
        fig, ax = plt.subplots(nrows=2,ncols=2)
        s = .1; m = 'o'
#        c = 'b'
        c = np.where(vf>0, 'r', 'k')

#        tmx = max(tr)
#        cl = [str(ti/(1.1*tmx)) for ti in tr]
        ax[0,0].scatter(Q, v, marker= m, s=s, c=c) #, c=cl
        ax[0,0].set_yscale('log');         ax[0,0].set_xscale('log')
        ax[0,0].set_ylabel('v')
        ax[0,0].set_xlim(  left=1e-2, right=1e2  )
        ax[0,0].set_ylim(bottom=1e-1, top=1e1)
        
        ax[0,1].scatter(Q, h, marker= m, s=s, c=c) #, c=cl
        ax[0,1].set_yscale('log');         ax[0,1].set_xscale('log')
        ax[0,1].set_ylabel('h')
        ax[0,1].set_xlim(  left=1e-2, right=1e2  )        
        ax[0,1].set_ylim(bottom=1e-2, top=1e1)
        
        ax[1,0].scatter(Q, w, marker= m, s=s, c=c) #, c=cl
        ax[1,0].set_yscale('log');         ax[1,0].set_xscale('log')
        ax[1,0].set_ylabel('w')
        ax[1,0].set_xlim(  left=1e-2, right=1e2  )        
        ax[1,0].set_ylim(bottom=1e-1, top=1e2)
        
        ax[1,1].scatter(Q, cs, marker= m, s=s, c=c) #, c=cl
        ax[1,1].set_yscale('log');         ax[1,1].set_xscale('log')
        ax[1,1].set_ylabel('c')
        ax[1,1].set_ylim(bottom=1e-7, top=1e-1)
        ax[1,1].set_xlim(  left=1e-2, right=1e2  )

        fn = '/ratings_x' + str(j) + '.png'
#        print(f'  17may20; just after namingFig')        
#        gs.tight_layout(fig) #17may20; deprec as just here comes coredump.
#        print(f'  17may20; just after tight_layout')
#        gs.update(hspace=.13)
#        print(f'  17may20; just after hspace')
        fnfull = fon + fn
#        print(f'  17may20; just before savefig')
        plt.savefig(fnfull)
#        print(f'  17may20; just after savefig')        
        plt.close('all')
        return c
        
    def dimless_ratings(self, fon, Qbf_dl, wbf_dl, hbf_dl, Qsbf_dl, S, cl, j):
        print(f'In plots2f.py/dimless_ratings; printing j={j}')
        fig, ax = plt.subplots(nrows=2,ncols=2)
        s = .1; m = 'o'
#        k1 = 1/np.sqrt(10*D50**5)
#        k2 = 10**.2/Qbf**.4
#        Qbf_dl = Qbf*k1  #dl: dimensionless
#        Qsbf_dl = Qsbf*k1

        ax[0,0].scatter(Qbf_dl, wbf_dl, marker= m, s=s, c=cl) #, c=cl
        ax[0,0].set_yscale('log');         ax[0,0].set_xscale('log')
        ax[0,0].set_ylabel('wbf*')
        
        ax[0,1].scatter(Qbf_dl, hbf_dl, marker= m, s=s, c=cl) #, c=cl
        ax[0,1].set_yscale('log');         ax[0,1].set_xscale('log')
        ax[0,1].set_ylabel('hbf*')
        
        ax[1,0].scatter(Qbf_dl, Qsbf_dl, marker= m, s=s, c=cl) #, c=cl
        ax[1,0].set_yscale('log');         ax[1,0].set_xscale('log')
        ax[1,0].set_ylabel('Qsbf*')
        ax[1,0].set_xlim(  left=1e2, right=1e8  ) #np.amin( Qbf_dl[:-1] )
        ax[1,0].set_ylim(  bottom=1e-3, top=1e6 )
        
        ax[1,1].scatter(Qbf_dl, S, marker= m, s=s, c=cl) #, c=cl
        ax[1,1].set_yscale('log');         ax[1,1].set_xscale('log')
        ax[1,1].set_ylabel('S')        

        fn = '/dimless_rats_x' + str(j) + '.png'
        fnfull = fon + fn
        plt.tight_layout()
        plt.savefig(fnfull)
        plt.close('all')
        
    def x_dimless_ratings(self, fon, Qbf_dl, S, wbf_dl, hbf_dl, Qsbf_dl, prt_wdw, A):
        varn = ['S', 'wbf_dl', 'hbf_dl', 'Qsbf_dl']
        print(f'In plots2f.py/x_dimless_ratings')
#        Amx= 1.1*np.amax(A)
        nt = len(Qbf_dl[:,0])
        xi = int(prt_wdw[2]); xf= int(prt_wdw[3])
#        nx = xf-xi+1 
#        n = nx*nt
        nv = len(varn)
        ons_t = np.ones(nt)
        x=[]; z1 = []; z2 = []
        for j in range( xi, xf + 1 ):
            xl = Qbf_dl[:, j-1].tolist()        
            z1l = S[:, j-1].tolist()
            z2l = ons_t*A[j-1]
            z2l = z2l.tolist()
            x.extend(xl)
            z1.extend(z1l)
            z2.extend(z2l)
        x = np.array(x); z1 = np.array(z1); z2 = np.array(z2);
        for k in range(1, nv+1):
            fig1, ax1 = plt.subplots()        
            vk = varn[k-1]
            yij = eval(vk)
            print(f'  vk:{vk}; yij.shape={yij.shape}')
            y = []
            for j in range( xi, xf + 1 ):
                yj = yij[:, j-1]
                yl = yj.tolist()
                y.extend(yl)
            si = np.where( z1>.03, .1, 5 )
            y = np.array(y)
            print(f'    preplot; x.shape:{x.shape}, y.shape:{y.shape}')
            im = ax1.scatter(x, y, c=z2, cmap = 'viridis', s=si, alpha=.3) # , , 
            cbar = fig1.colorbar(im, ax=ax1)
            cbar.set_label('A (km2)')
            ax1.set_yscale('log');         ax1.set_xscale('log')
            if vk=='Qsbf_dl': ax1.set_ylim(  bottom=1e-3, top=1e6 )
            fn = '/dimless_allX_' + vk + '.png'
            fnfull = fon + fn
            plt.savefig(fnfull)
            plt.close('all')    
    
    def basinSedBudget(self, fon, t, ci, ciw, cis, co, cdv):
        print(f'In plots2f.py/basinSedBudget')
#        print(f'  t.shape:{t.shape}; shapes: ci:{ci.shape}, co:{co.shape}, cdv:{cdv.shape}')
        col = 'tab:red'
        tr = t[:-1] #3jun20; to avoid last 0 y value.
        fig, ax = plt.subplots()
        l1 = ax.plot(tr, ci[:-1], '-c', markersize=1, label='cumIn')
        l3 = ax.plot(tr, ciw[:-1], ':c', markersize=1, label='cumIn_wash')        
        l4 = ax.plot(tr, cis[:-1], '-.c', markersize=1, label='cumIn_slid')        
        l2 = ax.plot(tr, co[:-1], '-k', markersize=1, label='cumOut')
        ax.set_ylabel('cumSedFlow(ton)')
        ax2 = ax.twinx()
        ax2.plot(tr, cdv[:-1], '-r', markersize=1)
        ax2.tick_params(axis = 'y', labelcolor=col)
        ax2.set_ylabel('cum'+r'$\Delta$alluv(ton)', color=col)
        ax.yaxis.set_major_formatter( ticker.FormatStrFormatter('%.1e') )
        ax2.yaxis.set_major_formatter( ticker.FormatStrFormatter('%.1e') )        
        ax.set_xlabel('time (hours)')
        ls = l1+l2+l3+l4
        labs = [l.get_label() for l  in ls]        
        leg = ax.legend(ls, labs, loc='best', fancybox=True) #ncol=2, , ncol=4, , fontsize=8
        fn = '/basinSedBudget.png'
        fnfull = fon + fn
        fig.tight_layout()        
        plt.savefig(fnfull)
        plt.close('all')
        
    def x_scaling(self, fon, xx, v_pct, varn, vmn2, dim, typ):
        print(f'In plots2f.py/x_scaling')
        s = 20
        m = ['*', '_', '+', 'o', '+', '_']
        c = ['r', 'r', 'r', 'k', 'b','b']
        print(f'  varn[3:]: {varn[3:]}')
#        typ = np.array(typ) #to support conditional
        for k1, var in enumerate( varn[3:] ): #1jun20; dont plot supply vars.
            k = k1+3
            typk = typ[k]
            if typk>0: #1jun20; do not plot Qbf_dl in y axis.
                xxk = xx[:,k]
                vmn1 = np.amin(v_pct[:,k,:])
                vmn = max(vmn2[k], vmn1);
                vmx = np.amax(v_pct[:,k,:])
                fig, ax = plt.subplots()
                vnk = varn[k]
                print(f'  k1={k1}, k={k}, vnk:{vnk}')
                print(f'  xxk:{xxk}')
        #        print(f'  1jun20; A{A}')
        #        ax.scatter(A, v_pct[:,4,-3], marker=m[-3], c=c[-3], s=s)
                for i,mp in enumerate(m):
                    ax.scatter(xxk, v_pct[:,k,i], marker=mp, c=c[i], s=s)
                    print(f'    pct#{i}, yi: { v_pct[:, k, i] }')
                ax.set_yscale('log');         ax.set_xscale('log')
                if typk==2: ax.set_xlim(left=1e2, right=1e8)
                xlab = 'A (km2)' if typk==1 else 'Qbf_dl_p50'
                ax.set_ylabel(f'{vnk} {dim[k]}');         ax.set_xlabel(xlab); 
                ax.set_ylim(bottom=vmn, top=vmx)
                fn = '/x_scaling_' + vnk + '.png'
                fnfull = fon + fn
                plt.savefig(fnfull)
                plt.close('all')
#        sys.exit('  just printed x_scaling')
        
    def seriesSed(self, flume_else_basin, fon, dt, t, cum_dAllu, cumIn_wash, cumIn_slid, cumIn, cumOut, Q, braidInd, S, c, ts, tc, 
        wproy,wv,wc,wb,colm,j): #, pi, Mgtm, gam....tsm, 
#        fig,ax=plt.subplots(nrows = 2, ncols = 1)
        print(f'In plots2f.py/seriesSed; printing j={j}')
        col = 'tab:red'        
        Dc = self.Dc
        fb = self.fb
        ms= 2
        gs = grsp.GridSpec(4, 1, height_ratios = [3, 3, 3, 3])
        fig = plt.figure()
        axd = fig.add_subplot(gs[0,:])
        axb = fig.add_subplot(gs[1,:])
        axa = fig.add_subplot(gs[2,:])                
        axc = fig.add_subplot(gs[3,:])
        axb.yaxis.set_major_locator(ticker.MultipleLocator(100))

        print(f'  17may20; before cumOut plot')
        tm1 = t[:-1]
        l1 = axb.plot(tm1, cumOut[:-1], 'k', label = 'cumOut', markersize=ms)
        l3 = axb.plot(tm1, cumIn[:-1], 'c', label = 'cumIn', markersize=ms)
        l4 = axb.plot(tm1, cumIn_wash[:-1], ':c', label = 'cumIn_wash', markersize=ms)
        l5 = axb.plot(tm1, cumIn_slid[:-1], '-.c', label = 'cumIn_slid', markersize=ms)
        print(f'  17may20; after cumOut plot')
        axb.set_ylabel('cumSedFlow(ton)', fontsize=9)
        cumOmx = np.amax(cumOut)
        cumImx = np.amax(cumIn)
        cummx = max(cumOmx, cumImx)
        majors = [0, .5*cummx, cummx]
        axb.yaxis.set_major_locator(ticker.FixedLocator(majors))
        axb.tick_params(axis = 'y', labelsize=8)
        axb.yaxis.set_major_formatter( ticker.FormatStrFormatter   ('%.0e') )
        ax1 = axb.twinx()
#        ax1.set_ylim(bottom=1e-7, top=1e-1)
#        ax1.set_yscale('log')        
#        majors = [1e-7, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
#        ax1.yaxis.set_major_locator(ticker.FixedLocator(majors)) #1apr20: LogLocator did not work for 2nd axis.
        l2 = ax1.plot(tm1, cum_dAllu[:-1], '-r', markersize=ms)
        print(f'  17may20; after alluv plot')        
        ax1.set_ylabel('cum'+r'$\Delta$alluv(ton)', color=col, fontsize=9)
        vmx = np.amax(cum_dAllu)
        vmn = np.amin(cum_dAllu)      
#        if vmn<0:
#            majors = [vmn, 0, vmx]
#        else: 
#            majors = [0, .5*vmx, vmx]
#        ax1.yaxis.set_major_locator(ticker.FixedLocator(majors))        
        ax1.yaxis.set_major_formatter( ticker.FormatStrFormatter('%.0e') )
        ax1.tick_params(axis = 'y', labelcolor=col, labelsize=8)  
#        ax1.yaxis.grid(which='both')        
#        ax2 = axb.twinx()
#        l3 = ax2.plot(t, D84, '.g', label = '$D_{84}$', markersize=ms)
#        ax2.set_ylabel('$D_{84}$(m) or S(m/m)')
#        ax2.set_yscale('log')
#        ax2.set_ylim(bottom=1e-3, top=1e-1)
#        axb.set_ylim(bottom=0, top=1e-3)
        ls= l1 + l3 + l4 + l5
        labs = [l.get_label() for l  in ls]        
        leg = axb.legend(ls, labs, loc='best', ncol=2, fontsize=8) #ncol=2, , fancybox=True

        lb1 = r'$\tau^*$'
        lb2 = r'$\tau^*_{c50}$'
        l1= axd.plot(t, ts, 'ob', markersize=1, label = lb1)
        l2= axd.plot(t, tc, ':k', markersize=ms, label = lb2)
#        lbtc = r'$\tau^*/\tau^*_{c50}$'        
#        l3 = axd.plot(t, ts/tc, '-.c', markersize=ms, label = lbtc)        
        yl = lb1 + '; ' + lb2
        axd.set_ylabel(yl)

        axd.yaxis.grid(which='both', color='0.8')

        if flume_else_basin==0: #29ago20; to better see exp decay by coarsening, especially in pp13 flood.    
            axd.set_yscale('log')
            axd.set_ylim(bottom=1e-2, top=1e-0)
        
        ax8=axd.twinx()
        print(f'  17may20; after taos plot')
        ax8.set_ylim(bottom=1e-5, top=1e-1)
        ax8.set_yscale('log')        
        majors = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
        ax8.yaxis.set_major_locator(ticker.FixedLocator(majors)) #1apr20: LogLocator did not work for 2nd axis.
        l4 = ax8.plot(t, c, '.r', label = 'c', markersize=1)
        print(f'  17may20; after c plot')
        ax8.set_ylabel('$c (m^3/m^3)$', color=col)        
#        lbtsm = r'$\tau^*/\tau^*_{m}$'
#        l4 = ax8.plot(t, ts/tsm, '--r', markersize=ms, label = lbtsm)
        ax8.tick_params(axis = 'y', labelcolor=col)
        ls = l1+l2 #+l4
        labs = [l.get_label() for l  in ls]        
        leg = axd.legend(ls, labs, loc='best', ncol=4, fancybox=True, fontsize=8) #ncol=2,        
#        ax8.set_yscale('log')
#        ax8.set_ylim(bottom=1e-1, top=1e1)
#        lbtc2 = lbtc # + r' ; $\tau^*/\tau^*_{m}$'
#        ax8.set_ylabel(lbtc2, color=col)
        
        l4=axa.plot(t, wproy, 'oc', label='wc_qsy', markersize=1.5)
        axa.set_ylabel('width(m)') #if infinite w_valley qsy would produce it.
        axa.set_yscale('log')
        l1 = axa.plot(t, wv, 'xb', markersize=1, label='wv' )
        majors = [1e-1, 1, 1e1, 1e2]
        axa.yaxis.set_major_locator(ticker.FixedLocator(majors)) #1apr20: LogLocator did not work for 2nd axis.                
        ax1 = axa.twinx()    
        ax1.plot(t, colm, '.r',markersize=3)
        ax1.set_ylabel('MphdynType', color=col)
        ax1.set_ylim(bottom=-5, top=5)
        majors = [-5,-2,0,2,5]
        ax1.yaxis.set_major_locator(ticker.FixedLocator(majors))
        ax1.tick_params(axis = 'y', labelcolor=col)        
        l2 = axa.plot(t, wc, 'og', label='wc', markersize=1.5) #2jun20; called w_proy in the f90 model.
        l3 = axa.plot(t, wb, '.k', label='wb', markersize=1)
        ls = l1+l2+l3+l4
        labs = [l.get_label() for l  in ls]        
        leg = axa.legend(ls, labs, loc='best', ncol=2, fancybox=True, fontsize=8) #ncol=2,                
        
#        nD = len(pi[0,:])
#        nD2 = nD+1
#        kl = np.arange(nD)
#        axa.set_yscale('log')
#        axa.set_ylim(bottom=1e-2, top=1e1)        
#        majors = [1e-2, 1e-1, 1, 10]
#        axa.yaxis.set_major_locator(ticker.FixedLocator(majors)) #1apr20: LogLocator did not work for 2nd axis.        
#        for k in kl:
#            cl = (kl[k]/nD2, kl[k]/nD2, kl[k]/nD2)
#            Dck = Dc[k]
#            lbn = np.around(Dck, decimals=2) if Dck<1 else int(Dck)
#            lb = str(lbn) + 'mm'
#            pik = pi[:,k]
#            rat = pik/fb[k]
#            rat[rat<1e-2] = np.nan
#            axa.plot(t, rat, '-', c=cl, markersize=ms, label = lb)
#        axa.legend(title='D(mm)',fontsize=6,ncol=6, fancybox=True, title_fontsize=7)
##        axa.set_ylim(bottom=0, top=pi_mx)
#        axa.set_ylabel('$p_i / f_{i,b}$')
#        ax1 = axa.twinx()
#        Mgtm[-1]=Mgtm[-2] #to avoid 0
#        gam[-1]=gam[-2] #to avoid 0
#        majors = [1e-1, 1, 1e1, 1e2]
#        ax1.set_yscale('log')
#        ax1.set_ylim(bottom=1e-1, top=1e2)
#        ax1.yaxis.set_major_locator(ticker.FixedLocator(majors))                                            
#        ax1.plot(t, Mgtm, '--r')
#        ax1.plot(t, gam, ':r')
#        ax1.set_ylabel(r'(--)$M_{gtm}$; (..)$\gamma_{gtm}$', color=col)
#        ax1.tick_params(axis = 'y', labelcolor=col)
        
        l1= axc.plot(t, Q, ':k', label = 'Q')
        axc.set_ylabel('$Q(m^3/s)$')
#        axc.set_xlabel('time (hours)') #dont work: , ha = 'left', va = 'top'
        axc.xaxis.set_label_coords(0,-.1)
        axc.set_xlabel('t(hours)', ha='right')
        ax1 = axc.twinx()
#        ax1.set_ylim(1e-3, 1e-1)
#        majors = [1e-3, 1e-2, 1e-1]
#        ax1.yaxis.set_major_locator(ticker.FixedLocator(majors))        
        lbi = 'i'
        ons = np.ones(len(t))   #reference braidInd for AB (anabranching) and B (braiding)
        mc=.4   #after Eaton10, threshold braidInd for multichannel mc is .28 or .4 for eqns 4 or 7 respect/. 
            #In this eq 7, threshold grows with bank resistance miu.
        ana_arr = mc*ons; bra_arr = 1.82*ana_arr
        l2 = ax1.plot(t, braidInd, 'or', label = lbi, markersize=1)
        l3 = ax1.plot(t, ana_arr, '-.r', label = 'i_ana')
        l4 = ax1.plot(t, bra_arr, ':r', label = 'i_bra')
        ax1.set_ylim(1e-1, 1e1)
        ax1.set_yscale('log')        
        print(f'  17may20; after all plots')
        lbi2 = lbi + r'$=SQ^{*0.44}$'
        ax1.set_ylabel(lbi2, color=col)
        ax1.tick_params(axis = 'y', labelcolor=col)        
        ls=l1+l2+l3+l4
        labs = [l.get_label() for l  in ls]        
        leg = axc.legend(ls, labs, loc='best', ncol=4, fancybox=True, fontsize=8) #ncol=2, loc=8

        axc.get_shared_x_axes().join(axb, axc, axd, axa)
        axa.set_xticklabels([])
        axb.set_xticklabels([])
        axd.set_xticklabels([])                
        print(f'  17may20; just shared axs')
        
        fn = '/sed_x' + str(j) + '.png'
        print(f'  17may20; just after namingFig')        
#        gs.tight_layout(fig) #17may20; deprec as just here comes coredump.
        fig.tight_layout()
        print(f'  17may20; just after tight_layout')
        gs.update(hspace=.13)
        print(f'  17may20; just after hspace')
        fnfull = fon + fn
        print(f'  17may20; just before savefig')
        plt.savefig(fnfull)
        print(f'  17may20; just after savefig')        
        plt.close('all')

        
    def series(self, flume_else_basin, fon, t, c, D50, D84, jam, D50ss, D84ss, Q, h, v, vf, w, S, h_s, h_sf, w_c, w_b, w_v,  
        Fr, asp, subm, j, nR):        
        col = 'tab:red'
#        fig,ax=plt.subplots(nrows = 2, ncols = 1)
#        print(f'28apr20; in plot_s2f_bin/series')        
#        print(f'  27apr20; before, check varType2; Q_mean={np.mean(Q[:-1])}, Q_median={np.median(Q[:-1])}') 
#            #ignore last time; set NaN.
#        print(f'  27apr20; before, check varType2; h={np.mean(h[:-1])}, h_median={np.median(h[:-1])}')
#        print(f'  27apr20; before, check varType2; v_mean={np.mean(v[:-1])}, v_median={np.median(v[:-1])}')
#        print(f'  27apr20; before, check varType2; vf={np.mean(vf[:-1])}, vf_median={np.median(vf[:-1])}')
#        print(f'  27apr20; before, check varType2; w={np.mean(w[:-1])}, w_median={np.median(w[:-1])}')
#        print(f'  27apr20; before, check varType2; S_mean={np.mean(S[:-1])}, S_median={np.median(S[:-1])}')
                
        par_wlog = self.par_wlog
        ms= 2
        gs = grsp.GridSpec(4, 1, height_ratios = [1, 1, 1, 1])
            #https://stackoverflow.com/questions/52560527/how-to-fully-customize-subplot-size-in-matplotlib
        fig = plt.figure()       
        axa = fig.add_subplot(gs[1,:])
        axaa = fig.add_subplot(gs[0,:])
        axb = fig.add_subplot(gs[2,:])
        axc = fig.add_subplot(gs[3,:])
        
        l1 = axa.plot(t, h_s + h, '.c', label = '$h+h_s$', markersize=1)
        l4 = axa.plot(t, h_sf, '-.', c='.6',  label = '$h_{sf}$', markersize=ms)
        l5 = axa.plot(t, h_s, '--k', label = '$h_s$', markersize=ms)
        axa.set_ylabel('h (m)')        
        ax12 = axa.twinx()        
        l11 = ax12.plot(t, subm, ':r', markersize=ms) #, label= '$R/D_{84}$'
        ax12.set_ylabel('subm$=R/D_{84}$', color = col)
        ax12.tick_params(axis = 'y', labelcolor=col)
#        ax12.set_ylim(bottom=0, top=10)
        l2 = axaa.plot(t, w, '.c', label = 'w', markersize = 1)
        l6 = axaa.plot(t, w_c, 'og', label = '$w_c$', markersize = 1.5) #, linestyle= 'None', marker='.'
        l7 = axaa.plot(t, w_b, '.k', label = '$w_b$', markersize=1)
#        l8 = axaa.plot(t, par_wlog[1]*Q**par_wlog[4], ':g', label = '$w_s$')
#        w_v_t = w_v*np.ones(len(t))
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
        if flume_else_basin==0: 
            axaa.set_yscale('log')        
            majors = [1e-2, 1e-1, 1e0, 1e1, 1e2]            
            axaa.yaxis.set_major_locator(ticker.FixedLocator(majors))        
        y2lb = lba + '. ' + lbj
        ax11.set_ylabel(y2lb, color=col)
        ax11.tick_params(axis = 'y', labelcolor=col)
        ls=l1+l5+l4 #+l11
        labs = [l.get_label() for l  in ls]
        leg = axa.legend(ls, labs, loc='best', ncol=3, fancybox=True, fontsize=9) #ncol=2, 
        leg.get_frame().set_alpha(0.7)
        ls = l2 + l6 + l7  + l10 + l12 #+ l9 + l8
        labs = [l.get_label() for l  in ls]
        leg = axaa.legend(ls, labs, loc='best', ncol=3, fancybox=True, fontsize=9) #ncol=2, 
        leg.get_frame().set_alpha(0.7)
        
#        print(f'28apr20; in plot_s2f_bin/series; after 1st plot code chunck')
        
        l2 = axb.plot(t, D84, '-k', label = '$D_{84}$', markersize=ms)
        l3 = axb.plot(t, D50, '-g', label = '$D_{50}$', markersize=ms)
        l5 = axb.plot(t, D84ss, ':k', label = '$D_{ss,84}$', markersize=ms)
        l6 = axb.plot(t, D50ss, ':g', label = '$D_{ss,50}$', markersize=ms)        
        axb.set_ylabel('$D$(m)')
        ax2 = axb.twinx()
        l4 = ax2.plot(t, S, '--r', label = 'S', markersize=ms)
        ax2.set_ylabel('S(m/m)', color=col)
        ax2.tick_params(axis = 'y', labelcolor=col)
        if flume_else_basin==0: #27ago20. 
            if j<nR: ax2.set_yscale('log')
            ax2.set_ylim(bottom=1e-2, top=1e-0) #bottom=0, top=.08        
            axb.set_ylim(bottom=1e-3, top=1e0)
            majors = [1e-3, 1e-2, 1e-1, 1e0]
            axb.yaxis.set_major_locator(ticker.FixedLocator(majors)) #1apr20: LogLocator did not work for 2nd axis.
            axb.set_yscale('log')
#        axb.set_ylim(bottom=.005, top=.035) #bottom=1e-3, top=1e-1        
        ls=l2+l3+l4+l5+l6
        labs = [l.get_label() for l  in ls]
        leg = axb.legend(ls, labs, loc='best', ncol=5, fancybox=True, fontsize=9) #ncol=2,


#        l1= axc.plot(t, Q, ':k', label = 'Q')
#        axc.set_ylabel('$Q(m^3/s)$')
#        axc.set_xlabel('time (hours)')
        axc.xaxis.set_label_coords(0,-.1)
        axc.set_xlabel('t(hours)', ha='right')
        l3 = axc.plot(t, v, 'ob', label = 'v', markersize = 1)
        l4 = axc.plot(t, vf, 'oc', label = '$v_f$', markersize = 1)
        axc.set_ylabel('v(m/s)')        
        ax3 = axc.twinx()
        l2 = ax3.plot(t, Fr, 'or', markersize = 1.5) #label = 'Fr',         
#        ax3.set_ylim(bottom=0, top=4)
        ax3.set_ylabel('Fr', color=col)
        ax3.tick_params(axis = 'y', labelcolor=col)
        ls=l3+l4
        labs = [l.get_label() for l  in ls]        
        leg = axc.legend(ls, labs, loc=1, ncol=2, fancybox=True, fontsize=8) #ncol=2,

        axa.get_shared_x_axes().join(axaa, axa, axb, axc)
        axa.set_xticklabels([])
        axaa.set_xticklabels([])        
        axb.set_xticklabels([])
              
        fn = '/hc_x' + str(j) + '.png'
        fnfull = fon + fn
#        print(f'28apr20; in plot_s2f_bin/series; fnfull: \n {fnfull}')
#        gs.tight_layout(fig)
        gs.update(hspace=.13)
#        print(f'28apr20; in plot_s2f_bin/series; after last plot code chunck, prev to save')
#        plt.show()
        plt.savefig(fnfull)
#        print('just after saving fig')
        plt.close('all')
                      
            
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
            
    def plt_hc(self, flume_else_basin, which_flume, case_n, plot_whole, flag_day_else_tStep, cdt,
        prt_wdw, dt_serie, Q, v, h, vf, w, wproy, c, D50, D84, 
        D50c_ss, D84c_ss, S, h_s, h_sf, w_c, w_b, w_v, ts, tc84,
        Fr, aspect, subm, colm, Qbf_dl, Qsbf_dl, hbf_dl, wc_dl, xx, v_pct, varn, vmn2, dim, typ, A,
        volIn, dvol, wash, slid): 
            #, Mgtm, gam, pi..., tsm
        print('27apr20: in plt_hc')
        #set last time NaN for all j, to avoid line fall at end of timeseries plots.
        Q[-1,:]=np.nan; 
        v[-1,:]=np.nan; 
        h[-1,:]=np.nan; 
        vf[-1,:]=np.nan; 
        w[-1,:]=np.nan;
        wproy[-1,:]=np.nan; colm[-1,:]=np.nan; volIn[-1,:]=np.nan; dvol[-1,:]=np.nan
        wash[-1,:]=np.nan; slid[-1,:]=np.nan
        Qbf_dl[-1,:]=np.nan; Qsbf_dl[-1,:]=np.nan; hbf_dl[-1,:]=np.nan; wc_dl[-1,:]=np.nan;
        c[-1,:]=np.nan; D50[-1,:]=np.nan; D84[-1,:]=np.nan; D50c_ss[-1,:]=np.nan; 
        D84c_ss[-1,:]=np.nan; 
        S[-1,:]=np.nan; h_s[-1,:]=np.nan; h_sf[-1,:]=np.nan; w_c[-1,:]=np.nan; 
        w_b[-1,:]=np.nan; ts[-1,:]=np.nan; tc84[-1,:]=np.nan; #tsm[-1,:]=np.nan;
        Fr[-1,:]=np.nan; aspect[-1,:]=np.nan; subm[-1,:]=np.nan; 
#        Mgtm[-1,:]=np.nan; gam[-1,:]=np.nan; pi[-1,:]=np.nan
#        print('28apr20: in plt_hc; after nan')
        
        Lf = self.Lf        
        nR = len(w_v)
        nt = len(Q[:, 0])

        if flume_else_basin==1:
            dx = self.dx
                #27apr20: as dx seems to be used in this code to calculate x,t map lims, 
                    #read only dx for flume cases. This map is not clear visualization for tree-like instead serial reaches.
            pwe = 'hg_model_flume/' #path for which environment        
            pwf = 'trench/' if which_flume == 1 else 'flood/'
        else:
            pwe = 'hg_model/'
            pwf = ''
        fon1 = 'doct_fluvial/' + pwe + pwf
#        if not os.path.exists(fon1): os.makedirs(fon1) #26apr20: folder has to exist, as it is the location of f90 calcs output.

#        case_n = '_10wv' + str(int(round(10*w_v[0],0))) + \
#            '_Lf' + str(int(Lf)) + '_10dx' + str(int(round(10*dx, 0))) + '_dt' + str(int(dt_serie[0]))
        
#        CORRELACION ESPACIAL ANCHO VS LECHO DE SHAWN CHARTRAND18jgr

#        print('28apr20: in plt_hc; after paths production')

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
        if flume_else_basin:            
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
#        ttm = ts/tsm

#        print('28apr20: in plt_hc; just before calc cumOut')
        cumOut = np.zeros((nt, nR)); cumIn = np.zeros((nt, nR)); cum_dAllu = np.zeros((nt, nR))
        cum_dAlluBasin = np.zeros(nt); cumInBasin = np.zeros(nt); 
        cumInBasin_wash = np.zeros(nt); cumInBasin_slid = np.zeros(nt);
        cumIn_wash = np.zeros((nt, nR)); cumIn_slid = np.zeros((nt, nR));
        for j in range(1, nR+1):
            cumOut[:-1, j-1] = 2.65*np.cumsum(c[:-1, j-1]*Q[:-1, j-1]*dt_serie[:-1]) #[:-1]
            cumIn[:-1, j-1] = 2.65*np.cumsum(volIn[:-1, j-1])
            cumIn_wash[:-1, j-1] = 2.65*np.cumsum(wash[:-1, j-1]*dt_serie[:-1])
            cumIn_slid[:-1, j-1] = 2.65*np.cumsum(slid[:-1, j-1]*dt_serie[:-1])
            cum_dAllu[:-1, j-1] = 2.65*np.cumsum(dvol[:-1, j-1])
        cum_dAlluBasin[:-1] = np.sum(cum_dAllu[:-1, :-1], axis=1) #2jun20; last reach is Fr1 bypass (matters if backwater),
            #so sed in=out so no dVolsed.
        slidBasin = np.sum(slid[:-1, :-1], axis=1); #Qsi = sum_j(Qsij)
        washBasin = np.sum(wash[:-1, :-1], axis=1);
        supplyHills = slidBasin + washBasin
        cumInBasin_wash[:-1] = 2.65*np.cumsum(washBasin*dt_serie[:-1]) #vol~ sum(Qsi*dti)
        cumInBasin_slid[:-1] = 2.65*np.cumsum(slidBasin*dt_serie[:-1])
        cumInBasin[:-1] = cumInBasin_wash[:-1] + cumInBasin_slid[:-1]
#        cumOut[-1,:]=np.nan; cumIn[-1,:]=np.nan; cum_dAllu[-1,:]=np.nan
            
#            print('cm,Qm,dt_serie',)

#        print('28apr20: in plt_hc; just before Xcorr')
#        #this Xcorr plotter np.where which gets unuseful for recalculated Q, which involves errors (<10%) from f90 numSolver.
#        if which_flume ==1: #trench            
#            Qgoal=[.01, .01, .032, .032, .018, .018, .007, .007]
#        else: #flood
#            Qgoal=[.012, .012, .0242, .0242, .0318, .0318, .0422, .0422]
#        indWhich = [nR+1, -1]    #if info on Q initial stages will be plotted, use tstep nR+1: when Q arrives to outlet.
#            #m3/s as comes from f90. 18 is only in rising limb and 7 only in falling.
##        firstElseLast = [1,0,1,0]
#        t_ej = []
#        for k, qg in enumerate(Qgoal):
#            print(f'28apr20: in Xcorr inputs production; qg: \n {qg}')        
#            print(f'28apr20: in Xcorr inputs production; Q[1:360,0]: \n {Q[1:360,0]}')
#            absRelErrQ = np.absolute(Q[:,0] - qg)/qg
#            ind,=np.where(absRelErrQ<.03)
##            print(f'28apr20: in Xcorr inputs production; ind: \n {ind}')
##            print(f'30mar, len(ind)={len(ind)}, ind={ind}')
##            print(f'30mar, k{k}, firstElseLast[k]={firstElseLast[k]}')
#            tejNew = ind[indWhich[k%2]]   # indWhich[ firstElseLast[k] ]  
#            t_ej.append(tejNew)
##            print(f'30mar, t_ej={t_ej}')            
#        xn = ['hsXanom', 'wXanom', 'gradNeg_hs', 'grad_wgw', 'vXanom', 'c', 'Fr']
#        tup_ls = []
#        for kk, varn in enumerate(xn):
#            profs_choosen = np.zeros((len(Qgoal), nR))            
#            for ii, ti in enumerate(t_ej):
#                profs_choosen[ii, :] = eval(varn)[ti,:]
#            tup2add = (varn, profs_choosen)
##            print(f'30mar, at varn{varn}, ti{ti}, profs_choosen[ii, :]={profs_choosen[ii, :]}')            
#            tup_ls.append(tup2add)
#            dic = dict(tup_ls)
#        self.plotXcorr(xn, Qgoal, dic, case_n, fon1)
#        print(f'prt_wdw={prt_wdw}')
#        sys.exit('17may20; checked ii,fin?')
        if plot_whole==0:
            if flag_day_else_tStep==0:
                ii= int(prt_wdw[0])
                ifin=int(prt_wdw[1])
            else:
                di= int(prt_wdw[0])
                dfin=int(prt_wdw[1])                
                imn_pair = cdt[cdt[:, 0]==di]
                imx_pair = cdt[cdt[:, 0]==dfin]
                ii = int(imn_pair[0][1]) - 1   #-1 for py indexes that start at 0
                ifin = int(imx_pair[0][1]) - 1                
        else:
            ii= 1
            ifin= nt
        if self.plot_series:
            fon = fon1 + 't' + case_n
            if not os.path.exists(fon): os.makedirs(fon)
#            fc.cleanFolder(fon)        #no need, as every name is overwritten.
            tr= t[ii:ifin]
            ntr = len(tr)
#            sys.exit(f'3jun20; ntr:{ntr}')            
            ons_t = np.ones(ntr)
            fld = np.zeros(ntr)
            if flume_else_basin==0: #27ago2020.
                self.x_scaling(fon, xx, v_pct, varn, vmn2, dim, typ)
                self.x_dimless_ratings(fon, Qbf_dl, S, wc_dl, hbf_dl, Qsbf_dl, prt_wdw, A)
            self.basinSedBudget(fon, tr, cumInBasin[ii:ifin], cumInBasin_wash[ii:ifin], cumInBasin_slid[ii:ifin], 
                cumOut[ii:ifin, nR-2], cum_dAlluBasin[ii:ifin])
#            sys.exit('3jun20; check basinSedBudget plot')
                #spatially lumped: In=hillslope supply, Out=fromXfin(j-1 as j=bypass).
            for j in range(int(prt_wdw[2]), int(prt_wdw[3]) + 1):   #all channel reaches
                D50j=D50[ii:ifin, j-1];   D50ss_j = D50c_ss[ii:ifin, j-1];  D84ss_j = D84c_ss[ii:ifin, j-1]
                cj=c[ii:ifin, j-1]; D84j=D84[ii:ifin, j-1]; jam_j = jam[ii:ifin, j-1]
                Qj= Q[ii:ifin, j-1]; hj=h[ii:ifin, j-1];
                vj=v[ii:ifin, j-1]; vfj=vf[ii:ifin, j-1];  
                wj=w[ii:ifin, j-1]; Sj=S[ii:ifin, j-1]; hsj= h_s[ii:ifin, j-1]; 
                wproy_j=wproy[ii:ifin, j-1]; colm_j=colm[ii:ifin, j-1];
                Qbf_dl_j = Qbf_dl[ii:ifin, j-1]; Qsbf_dl_j = Qsbf_dl[ii:ifin, j-1]; 
                hbf_dl_j = hbf_dl[ii:ifin, j-1]; wc_dl_j = wc_dl[ii:ifin, j-1]
                hsfj=h_sf[ii:ifin, j-1]; wcj=w_c[ii:ifin, j-1]
                wbj=w_b[ii:ifin, j-1]; tsj= ts[ii:ifin, j-1]; tcj= tc84[ii:ifin, j-1]; #tsmj= tsm[ii:ifin, j-1]
                Frj= Fr[ii:ifin, j-1]; asp_j= aspect[ii:ifin, j-1]; subm_j= subm[ii:ifin, j-1]
#                pi_jk= pi[ii:ifin, j-1, :]; Mgtm_j = Mgtm[ii:ifin, j-1]; gam_j = gam[ii:ifin, j-1]
                #timeseries:
                coj = cumOut[ii:ifin, j-1]
                cij = cumIn[ii:ifin, j-1]; ciwj = cumIn_wash[ii:ifin, j-1]; cisj = cumIn_slid[ii:ifin, j-1]
                cdaj = cum_dAllu[ii:ifin, j-1]                 
                bij = braidInd[ii:ifin, j-1]
                wvj = w_v[j-1]*ons_t
                print(f'28apr20: in plt_hc; prev to call series() for j={j}')
                self.series(flume_else_basin, fon, tr, cj, D50j, D84j, jam_j, D50ss_j, D84ss_j,
                    Qj, hj, vj, vfj, wj, Sj, hsj, hsfj, wcj, wbj, w_v, Frj, asp_j, subm_j, j, nR)
#                print(f'28apr20: in plt_hc; POS disq plot series() for j={j}')
                self.seriesSed(flume_else_basin, fon, dt_serie, tr, cdaj, ciwj, cisj, cij, coj, Qj, bij, Sj, cj, tsj, tcj, 
                    wproy_j,wvj, wcj, wbj, colm_j,j) 
                    #, tsmj                     #, pi_jk, Mgtm_j, gam_j
                if flume_else_basin==0:    
                    cl_fld = self.ratings(fon, Qj, hj, vj, wj, cj, vfj, j)
                    self.dimless_ratings(fon, Qbf_dl_j, wc_dl_j, hbf_dl_j, Qsbf_dl_j, Sj, cl_fld, j) 
                    #29may20; parker07 basin calibration, p233d
                            
                
        

                    
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
        if flume_else_basin:         
            xf = nR*dx        
            #filters to highlight threshold values:
            filt_tol =.1
            Fr_msk = self.setThreshVals(Fr, filt_tol, 1)
            subm_msk = self.setThreshVals(subm, filt_tol, 4)
            jam_msk = self.setThreshVals(jam, filt_tol, 5)
            aspect_msk = self.setThreshVals(aspect, filt_tol, 50)
            braidInd_msk = self.setThreshVals(braidInd, filt_tol, .4)
            ttc_msk = self.setThreshVals(ttc, filt_tol, 1)
    #        ttm_msk = self.setThreshVals(ttm, filt_tol, 1)
            
            self.mapSaver(fon1, case_n, wXanom, 'wXanom', mn_an, mx_an, xf, tf)
            self.mapSaver(fon1, case_n, hXanom, 'hXanom', mn_an, mx_an, xf, tf)        
            self.mapSaver(fon1, case_n, FrXanom, 'FrXanom', mn_an, mx_an, xf, tf)
            self.mapSaver(fon1, case_n, vXanom, 'vXanom', mn_an, mx_an, xf, tf)        
            self.mapSaver(fon1, case_n, hsXanom, 'hsXanom', mn_an, mx_an, xf, tf)
            self.mapSaver(fon1, case_n, aspect_msk, 'aspect', 1, 100, xf, tf, log=True)
            self.mapSaver(fon1, case_n, subm_msk, 'subm', 1, 7, xf, tf)
            self.mapSaver(fon1, case_n, jam_msk, 'jam', 1e0, 1e2, xf, tf, log=True)
            self.mapSaver(fon1, case_n, ttc_msk, 'ttc', 1e-1, 1e1, xf, tf, log=True)
    #        self.mapSaver(fon1, case_n, ttm_msk, 'ttm', 1e-1, 1e1, xf, tf, log=True)
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
            self.mapSaver(fon1, case_n, w_b, 'wb', 0, np.amax(w_b[:-1,:-1]), xf, tf)
            self.mapSaver(fon1, case_n, S, 'S', 0, np.amax(S[:-1,:-1]), xf, tf)
            self.mapSaver(fon1, case_n, c, 'c', 1e-6, 1e-3, xf, tf, log=True)
            self.mapSaver(fon1, case_n, w, 'w', .2, 2, xf, tf)
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
        
    def plt_ddt(self, flume_else_basin, plot_tseries, ds, t, y, dt, l, d2pi, d2pf, st, j, pct_req, 
        plog, dim, pad, yd=None):
         #d2p: days to plot; ds: day-timeStep relation [day, ts].
        print('en plt_ddt:')
        of = self.of
        tp_y, ys = self.pct_dyn(dt, y)

        ys = np.array(ys)
        ypl = []
        tol = .05; toli=1-tol; tolf=1+tol
        for i,p in enumerate(pct_req):
            y_pct = ys[ (tp_y > toli*p) & (tp_y < tolf*p) ] #17jun20; 1.1*p instead p+.01 to better catch pct .3%.
            if len(y_pct)<1: sys.exit(f'1jun20; pctile range with no value at p{p}')
            ypl.append(y_pct[0])
#        print('tp_y:', tp_y)
#        print('ys:', ys)
        print(f'  {l[0]} ypl:{ypl}')
#        sys.exit('31may, stop to check if ok pctils')

        tp_dt, dts = self.pct_dyn(dt, dt)
#        ys = np.sort(y);  dts = np.sort(dt)
        if pad==1: 
            yds = np.sort(yd)
            nyd = len(yd)
            td_steps = np.arange(nyd)
            tpd = 1-td_steps/nyd    #tp = 1-t_steps/ny;             
            td = np.arange(1, nyd+1);  tdr = td[d2pi-1 : d2pf-1+1];  ydr = yd[d2pi-1 : d2pf-1+1]            
        ny = len(y);  
#        imx = int(ny/20/365*d2p)   #jan8, 2020: 1mes
        imn_pair = ds[ds[:, 0]==d2pi]
        imx_pair = ds[ds[:, 0]==d2pf]
        imn = int(imn_pair[0][1]) - 1   #-1 for py indexes that start at 0
        imx = int(imx_pair[0][1]) - 1
        t_steps = np.arange(ny);   
        tr = t[imn:imx];  yr = y[imn:imx];  dtr = dt[imn:imx]
#        print(f'imn_pair, imx_pair ={imn_pair}, {imx_pair}')
#        print(f'imn_pair[0][1]={imn_pair[0][1]}')
#        print(f'imn{imn}, imx{imx}')
#        print(f'tr= {tr}')
#        print(f'yr= {yr}')

#        print(f'sizes in plot_cd: tpd{len(tpd)}, yds{len(yds)}, tp_y{len(tp_y)}, tp_dt{len(tp_dt)}, ys{len(ys)}, , dts{len(dts)}')
#        print(f'tp_y_mn{min(tp_y)}, tp_y_mx{max(tp_y)}')
#        print(f'tp_dt_mn{min(tp_dt)}, tp_dt_mx{max(tp_dt)}')

#        fig,ax=plt.subplots(nrows=2,ncols=1)        
        if plot_tseries==1:
            fig = plt.figure()        
            grid = plt.GridSpec(2, 2) #, hspace=.4, wspace=.6)
            ylab = f'{l[1-1]} {dim}'
                    
            ax1 = fig.add_subplot(grid[0, :])
            l1 = ax1.plot(tr, yr, 'ok', label= l[1-1], markersize = 1.5);     
            if l[0]=='c':
                ax1.set_ylim(bottom=1e-5, top=1e-1)
            ax1.set_ylabel(ylab); 
            ax1.set_xlabel('time (s)')
            if pad==1: ax1.plot(tdr, ydr, 'or', label= l[2-1], markersize = 3, markerfacecolor='none');
            axa = ax1.twinx()
            l2 = axa.plot(tr, dtr, '.r', label= l[3-1], markersize = 1.5);     axa.set_ylabel(l[3-1] + ' (s)')
            ls=l1+l2
            labs = [l.get_label() for l  in ls]        
            leg = ax1.legend(ls, labs, loc='best', ncol=1, fontsize=9) #ncol=2, , fancybox=True
            
            ax2 = fig.add_subplot(grid[1, 0]) if flume_else_basin == 0 else fig.add_subplot(grid[1, :])
            ax2.plot(tp_y, ys, 'ok', label = f'{l[1-1]}, Nsteps{ny}', markersize = 1.5)
            if l[0]=='c':
                ax2.set_ylim(bottom=1e-5, top=1e-1)        
            if pad==1: ax2.plot(tpd, yds, '.r', label = f'{l[2-1]}', markersize = 1)
            ax2.set_xlabel('prob.exced');  ax2.set_ylabel(ylab)
            axb = ax2.twinx()
            axb.plot(tp_dt, dts, '.r', markersize = .3)
            if flume_else_basin == 1: axb.set_ylabel(l[3-1] + ' (s)')
            
            #as previos ax2 but no log scale in axis.
            if flume_else_basin == 0:
                ax3 = fig.add_subplot(grid[1, 1])
                ax3.plot(tp_y, ys, 'ok', label = f'{l[1-1]}, Nsteps{ny}', markersize = 1.5)
                if l[0]=='c':
                    ax3.set_ylim(bottom=1e-5, top=1e-1)                    
                if pad==1: ax3.plot(tpd, yds, '.r', label = f'{l[2-1]}', markersize = 1)
                ax3.set_xlabel('prob.exced');  #ax3.set_ylabel(l[1-1] + ' [m3/s]')
                axb = ax3.twinx()
                axb.plot(tp_dt, dts, '.b', markersize = .3)
                axb.set_ylabel(l[3-1] + ' [s]')
            
            if flume_else_basin == 0:
                if plog==1:
                    ax1.set_yscale('log')            
                    ax2.set_yscale('log')
                    ax3.set_yscale('log')  
                axa.set_yscale('log')
                ax2.set_xscale('log')
    #            ax2.set_ylim(bottom=1e-2, top=1e2)
                axb.set_yscale('log')        
    #            ax3.set_ylim(bottom=1e-2, top=1e2)
            
    #        if l[1-1] != 'Qhill':  
    #            ax1.set_ylim(bottom = ymn)  #shows only sediment supply concentration >1e3ppm.
    #            ax2.set_ylim(bottom = ymn);  ax3.set_ylim(bottom = ymn)

            path = f'{of}/{l[1-1]}_reach{j}.png'
            fig.suptitle(st);  
            fig.tight_layout()
            print(f'23mar, to save downscaled series and CDs, path: {path}')
            plt.savefig(path)
            plt.close()
        return ypl
