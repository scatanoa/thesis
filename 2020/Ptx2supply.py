#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#mar17, 2018
#updated 22dic2019  #even if fluvial GHmodel is not grainFractional, sedupply here is 1 textural & 2aggregated, to recicle code.

import numpy as np
import matplotlib.pyplot as plt
#import random as rn

import mapPlot as mp
import plot_supply as ps
import sys



class lumped: 
    def __init__(self, Ptx, n_days, A_hs, elong_ratio, hgy, geotech, geotech2, weath_coefPerD, weath_xrefPerD, pr, flags):
        self.i_mn_show = pr[0];  self.i_mx_show = pr[1];  self.j_mn_show = pr[2];  self.j_mx_show = pr[3]
        self.flag_printDetails = flags[5]                                
        drain_dens_km = geotech2[0]; slidSlop = geotech2[1]; cwash = geotech2[2]
        ms = geotech2[3]; mw = geotech2[4]  
            #wat&sed 'catena': porcentual comparisons gully/rill, i.e. mw amplifies water gradient related to sed one.
        self.nDh = int(geotech2[5])
        self.ncv = 20 #6jul20; # of cohesion values inside subbasin.        
            # #diameters in hills (1e-4, -3, -2, -1, 1e0 m). Get fract supply even to lump for fluvHGmodel.
        X = np.loadtxt('doct_topol/Topol_X.txt'); n_reaches = len(X[:, 0])
        Abasin = X[-1, 4]        #25jul20; km2
        t_ev = np.loadtxt('doct_rain/t_ev.txt'); t_ev = t_ev.astype(int)
        self.output_folder = 'doct_supply/supply_wat_sed'
        self.n = int(n_days); self.elong_ratio = elong_ratio; self.A_hs = A_hs       #'hillslope' area [m2]
        self.psc = ps.c(self.n, self.output_folder, self.nDh, self.ncv)
        Qhill = np.zeros((self.n, n_reaches)); wash = np.zeros((self.n, n_reaches)); 
        slid = np.zeros((self.n, n_reaches, self.ncv)); hcr = np.zeros((self.n, n_reaches, self.ncv))
        zsg = np.zeros((self.n, n_reaches, self.ncv)); FS = np.zeros((self.n, n_reaches, self.ncv))
        slidFlg = np.zeros((self.n, n_reaches, self.ncv))
        dayOfEv = np.zeros(self.n)
        for ii in t_ev:
            dayOfEv[ii-1] = 1
        coh = np.linspace(1e3, 13e3, self.ncv)     #9jul20;     
#        coh = np.random.default_rng().uniform(1e3, 13e3, self.ncv) #6jul20;
            #dietrich95 p 17 range from grass to forest (sc: 1 OM!). proportional to hcr odorico03. 
            #odorico03. pr13: zs(t) study so coh>0.
            #pr29; coh ~ 1/zs by less effect of roots (see odorico veget slids)                    
        for j in range(1, (n_reaches-1) + 1):    #last topol reach denotes final basin point. #points = #subbasins - 1.
            print()
            print()
            print(f'reach{j}')
            Pt = Ptx[:, j-1]
            wash_j, slid_j, Qhill_j, slidFlg_j, hcr_j, zsg_j, FS_j, slids_mx_perSubb = \
                self.hill_s_t(j, dayOfEv, Pt, hgy, drain_dens_km, slidSlop, cwash, ms, mw, \
                    geotech, coh, weath_coefPerD, weath_xrefPerD)
            Qhill[:, j-1] = Qhill_j;  wash[:, j-1] = self.cumD(wash_j);             
            for jj in range(1, self.ncv+1):
                slid[:, j-1, jj-1] = self.cumD(slid_j[:,:,jj-1])
                slidFlg[:, j-1, jj-1] = slidFlg_j[:, jj-1]
                hcr[:, j-1, jj-1] = hcr_j[:, jj-1]
                zsg[:, j-1, jj-1] = zsg_j[:, jj-1]
                FS[:, j-1, jj-1] = FS_j[:, jj-1]                                                
        slidsPerDay = np.sum(slidFlg, axis = 1)
        coh = np.array(coh)        
        print(f'6jul20; coh={coh}')
        i_cohMin, = np.where(coh==np.amin(coh))
        slid_cohAve = np.mean(slid, axis=2)
        self.psc.supply_x(Ptx, Qhill, slid, slid_cohAve, hcr, zsg, FS, slidsPerDay, coh, self.ncv, i_cohMin[0], 
            n_reaches, Abasin, slids_mx_perSubb)
        f = ['Qhill_tx', 'wash_tx', 'slid_tx']
        of = self.output_folder
        mp.mapPlotter(Qhill, of, f[0], 'Hillslope Q(m3/s)', '# subbasin', 'time (days)', log_cbar=True);  
#        mp.mapPlotter(wash, of, f[1], 'surface erosion i.e. wash flow (m3/s)', '# subbasin', 'time (days)', 
#            mn=1e-3, log_cbar=True);  
#        mp.mapPlotter(slid[:,:,0], of, f[2], 'landslide flow(m3/s)', '# subbasin', 'time (days)') #, 
##            mn=1e-5, log_cbar=True) #23jul20; discarded as no clear visual. So, show Qmap & rating cloud Qs-Q.
        np.savetxt(f'{of}/{f[0]}.txt', Qhill); np.savetxt(f'{of}/{f[1]}.txt', wash); 
        np.savetxt(f'{of}/{f[2]}.txt', slid_cohAve) #12jul20; mean of all cohesion values.
            #23dic2019: saving only aggr. supply for GH fluvial model that does not account for fractional sed balance.            
        self.output_folder = of
        
    def cumD(self, x):        #cumD: cumulate all grain sizes D.
        return np.sum(x, axis = 1)

        
#    def hill_w_plotter
       
        
    def hill_w_t(self, j, dayOfEv, P, hgy): #parameter units [L=mm, T=dia]    
    #tropicalVals(Marinilla): 1-4, 20-600(100), 1-100(20), .01-10(4), 1-10(1), 1-10(3), 50-200(80)
        etp = hgy[0]; hu=hgy[1]; kssp=hgy[2]; kb=hgy[3]; tsp=hgy[4]; tssp=hgy[5]; t_aqf=hgy[6]
        output_folder=self.output_folder
    
        x1w=np.zeros((self.n)); x2w=np.zeros((self.n)); x3w=np.zeros((self.n)); x4w=np.zeros((self.n))
        o1w=np.zeros((self.n)); o2w=np.zeros((self.n)); o3w=np.zeros((self.n))
        Pserie = np.zeros((self.n)); Q=np.zeros((self.n)) # mm and m3/s
        x1w[0]=0.  # water at residual soil [0-.2]*1e3 mm
        x2w[0]=0.   # water for fast runoff   
        x3w[0]=0.   # water at saprolite[.1-1] *1e3 mm
        x4w[0]=.1e3     # water at fractured rock(.1-2) *1e3 mm
        ev = 0
        balance = 0
        for i in range(1,self.n):        #hillslope hidrology
            if dayOfEv[i-1] == 1:
                Pserie[i-1] = P[ev]
                ev += 1                    
            else:
                Pserie[i-1] = 0            
            Ptoday = Pserie[i-1]                 
#            print(f'in hill_w_t, tot{self.n}days, i{i}')
            i1w=min(Ptoday*(1-x1w[i-1]/hu)**2,hu-x1w[i-1])
            x1w[i]=x1w[i-1]+i1w
            o1w[i]=min(x1w[i],etp*(x1w[i]/hu)**.6)
            x1w[i]=x1w[i]-o1w[i]
            i2w=max(0,Ptoday-i1w-kssp)
            x2w[i]=x2w[i-1]+i2w
            o2w[i]=x2w[i]/tsp   
            x2w[i]=x2w[i]-o2w[i]
            i3w=max(0,Ptoday-i1w-i2w-kb)
            x3w[i]=x3w[i-1]+i3w
            o3w[i]=x3w[i]/tssp    
            x3w[i]=x3w[i]-o3w[i]
            i4w=max(0,Ptoday-i1w-i2w-i3w)
            x4w[i]=x4w[i-1]+i4w
            o4w=x4w[i]/t_aqf
            x4w[i]=x4w[i]-o4w
            Q[i]=o2w[i]+o3w[i]+o4w    #mm/d
#            if j < self.j_mx_show: print(f'in hill_w_t, Q[mm/d]={Q[i]}')
            balance += x1w[i]+x2w[i]+x3w[i]+x4w[i]-(x1w[i-1]+x2w[i-1]+x3w[i-1]+x4w[i-1])+ \
                (o1w[i]+o2w[i]+o3w[i]+o4w)-Ptoday   #must be 0                
        Q=Q/1e3*self.A_hs/86400    #m3/s
        Q[0]=.05*(self.A_hs/1.e6)    #average yield as artificial starting water discharge for sediment simulation.
        if (j >= self.j_mn_show and j <= self.j_mx_show):
            self.psc.hgy_plotter(j, Pserie, Q, x1w, x2w, x3w, x4w, o1w)
            if i<=self.i_mx_show and i>=self.i_mn_show:
                print(f'balance={balance}')
    #            print ("in hill, P=", P)
                print ("daily mean ETR=", o1w.mean())
    #            print ("daily timeseries of ETR=", o1w)
    #            print ("from fictBasin py, Qhill= ", Q)
        results_hillW = (Q, o2w, o3w, x3w)
        return results_hillW
        
        
#    def hill_s_plotter        23dec2019: no time to organize.
        
    
    def hill_s_t(self, j, dayOfEv, Pt, hgy, drain_dens_km, slidSlop, cwash, ms, mw, \
        geotech, coh, weath_coefPerD, weath_xrefPerD):
        #produces sediment supply from hillslope, triggered by surface runoff wash and saturation planar lansdlide.
        
        
        #general parameters:        
        ntd = self.n
        A_hs = self.A_hs; A_hs_km2 = A_hs/1e6
        output_folder=self.output_folder
        nDh = self.nDh
        ncv = self.ncv
        L_hill = .5/drain_dens_km*1e3
        tssp = hgy[5]
        ks = L_hill/(tssp*86400) #7jul20; [m3/s]
        #p267 notes and eq5 in https://pubs.usgs.gov/pp/0422c/report.pdf

        #slide parameters:        
        coefMalamud = .024; expMalamud = 1/.368 #from eq 18 of Malamud04 [km2, km3].
        mmYr = 1e-3/365; gw=1.e4; gs = 2.65*gw
        soilPor=geotech[0] #; coh=geotech[1]; coh_cv =geotech[2] #coeficient of variation (std_desv/mean) of cohesion 
#        print(f'coh={coh}')
#        sys.exit('check coh unif 10 samples')
        tanFi = geotech[1];  zsIniPerSize = geotech[2] 
#        tanFi = np.tan(frictAngle*3.14/180)
        cwsh = mmYr*np.array(weath_coefPerD);  xsh_ref = weath_xrefPerD
        solidPort = 1-soilPor; voidRat = soilPor / (1-soilPor); f_gm = gs / (voidRat + 1)
        zsSolidIniPerSize = zsIniPerSize*solidPort #(1.3 - np.random.rand())*.. desv<=50% for each subbasin j, given diffHistory.
        slidSlop_angle = np.arctan(slidSlop); sinS = np.sin(slidSlop_angle); cosS = np.cos(slidSlop_angle)
        gwtf = gw*tanFi; S_exc = slidSlop - tanFi
        dif = 100*1e-5 #7jul20; p266d 3e-3 dietrich95 [m2/yr] (1e-5m2/day)
        tb_ta = .8 #assumed
        Adr_slid = 3.7e3     #3700m2 as pr27 of odorico03 in case steep suplim 
            #(vary from 7e2 to 4e3 among cases).
        slids_mx_perSubb = A_hs/Adr_slid  #  # max of slides in subbasin.
            #as marcHovius slides amalgamation may occur.
            #Debe corregirse pa restar las ya deslizadas con cohesiones anteriores.
        hollowSidSlop = slidSlop/tb_ta
        seedLayer = .01; Aseed = seedLayer**2/hollowSidSlop  #1cm; see geom in fig3 odorico03.
        kdif = 2*dif*cosS*(hollowSidSlop**2 - slidSlop**2)    #eq 9 odorico03
#        print(f'kdif: {kdif}')
#        sys.exit('check kdif~dif')
        
        #wash parameters:
        fWh = 2*drain_dens_km*A_hs_km2*1e3

#        #bedload wash parameters: (the famous jarret84 eqn is for slopes from .2 to 3%: milder than steppools)
#        #31dec2019 I concluded that these search for qc (~3 so Q=qW=20!! is 10 times Q2.33)
#            #to apply Rickenmann01 equation is not needed,
#            #as order1 reach already computes fluvial mass balance. Then kilincRicharson wash was recovered.
##        e = 1.2 #        #to estimate S0 without macroroughness for Chiari n0/ntot equation.
###        abs_exp_nPart = .36; coef_nPart = .102  # Recking11hc: vpn and nPart from eqn19b.
##        b2= e_vpn_h*(1.-2./e);  b1 = .17 + 1./e - b2;  b3= .67 - b1;  b4= .5 - 1./e;  b5= 1./(1.-2./e)
##        eR = b3*b5; eD = -b2*b5; eS = b4*b5; eg = eS; ec = 2*b4*b5
###        b1 = 1.5 + exp_vpn;  b2 = .5 - abs_exp_nPart
#        c_vpn_l = 2.5; e_vpn_l = 1.0; c_vpn_h = 6.5; e_vpn_h = .167    #params Ferguson07 eq 16..
#            #.. that Recking11hc suggest, for for Low and High submerg relative to  R/D ~4, which Recking11hc eqn 19b shows.
#        W = 7.6*A_hs_km2**.2   #1.8*A_hs_km2**.42....consistente con dispersion W-A en paper 'continental bankfull width'.
#            #fig4.14 in Jimenez15 dissertation, supplyLim to get higher W to avoid overstimate dc (and S>2.5%).
#        D84 = 5*washSlop #~fig3, Recking2013.
#            #critical discharge, see deduction in pg2019_37:
#        tao_star_cr_84 = .56*washSlop + .021 #eq12 recking2016gtm.        
#        Rc = 1.65*D84*tao_star_cr_84/washSlop;  vfr = np.sqrt(10*Rc*washSlop); rRD = Rc/D84;  r_v_vfr = c_vpn_l*rRD**e_vpn_l
#        vc = r_v_vfr*vfr;   dc = Rc*W/(W-2*Rc); qc = vc*dc
##        vc = c_vpn_h**(2*b4)*10**b4*D84**eD*Rc**eR*washSlop**eS  #pg 2019_48        
##        qc = coef_nPart*coef_vpn*hc**b1*washSlop**b2/D84**abs_exp_nPart 
#            #1 OM closer qc to table1 Rickenmann01 than Zimm13 qc below.
###        qc = 11*D84**1.5*washSlop**(-2)*tao_star_cr_84**2.5*np.sqrt(10) #alternative deduction from v in Zimm13 (my pg2019_45).
#            #last figs in Ferguson12: d/D ~2 for Fr=1 where So~.1. f is around 1 according to Zimm13 and fig2-15 CatañoMSc15
#        if j==1: 
#            f = 8/r_v_vfr**2;  n = (f/8*10)**.5/Rc**.17;  Frc = vc/np.sqrt(10*dc);
#            v0 = vfr*6.5*(Rc/D84)**.167   #assuming high submergence.
#            f0 = 8*(vfr/v0)**2;  r_S0_S = f0/8*vc**2/(10*Rc)/washSlop; n0 = (f0/8*10)**.5*Rc**.17            
#            Wchr = round(W, 3); Rcr = round(Rc, 3); dcr = round(dc, 3); rWd= round(W/dc, 3 ); D84r = round(D84, 3)
#            print(f'D84{D84r}, Rc/dc={Rcr}/{dcr}, W{Wchr}, W/d{rWd}, R/D{round(rRD, 3)}, thcr{round(tao_star_cr_84, 3)}, \
#                vc{round(vc, 2)}, Frc{round(Frc, 2)}, f0/f={round(f0, 2)}/{round(f, 2)}, n0/n={round(n0, 2)}/{round(n, 2)} \
#                S0/S={round(r_S0_S, 2)}, qc{round(qc, 3)}')


        #variable initialization and initial conditions:
#        self.n=int(self.n)
        xsh=np.zeros((ntd, nDh, ncv)); wash_sh=np.zeros((ntd, nDh)); slid_sh=np.zeros((ntd, nDh, ncv)); 
        slidFlg = np.zeros( (ntd,ncv) )
        ch=np.zeros((ntd)); ch_wash=np.zeros((ntd)); ch_slid=np.zeros((ntd))
        FS=np.zeros((ntd, ncv)); pAslid = np.zeros((ntd)); hw=np.zeros((ntd, ncv)); hwg=np.zeros((ntd)); 
        asp=np.zeros((ntd, ncv));         As=np.zeros((ntd, ncv))
        zsg=np.zeros((ntd,ncv))
        hcr=np.zeros((ntd, ncv)); hmx = np.zeros((ntd, ncv))
        q=np.zeros((ntd)); qs_wash=np.zeros((ntd)); qs_mx_teor=np.zeros((ntd)); gm=np.zeros((ntd))
        #initial storages of sediment in hillslope, by grain size (can be up to 5, increasing 1 OM from .1mm)
        for jj in range(1, ncv+1):
            xsh[0,:, jj-1] = np.ones(nDh)*zsSolidIniPerSize                
        pSat=0.; ii= []


        #produce hillslope hidrology input:
        results_hillW = self.hill_w_t(j, dayOfEv, Pt, hgy)
        Q=results_hillW[0]; o2w=results_hillW[1]/1e3; o3w=results_hillW[2]/1e3; x3w=results_hillW[3]/1e3
            #[o3w] got from mm to m, so now it is [m/day]

        #get sediment from wash and slides:
#        p = 1e-3 #initial, cos each t, it is updated by Malamud law.
#        cps= ms/(1-p);  cpw = mw*cps #5jun20; cps is c_prime for seds as aux coef defined in p41ul.
#        cp = [cps, cpw];  cp = np.array(cp)  #1st component for sed and 2nd for wat.
#        yg = cp/(1+cp*p)    # pg 2019_42 in notes.        
        for i in range(1,ntd):
            ii.append(i-1)
            for jj,co in enumerate(coh): #6jul20; each cohesion value (each subsubbasin class)                            
                #get p, which is the portion of hillslope area likely to slide (gully, instead of rill):
                tot_soil = sum(xsh[i-1,:,jj]);  zsg[i,jj] = tot_soil/solidPort    #soil thickness in gully. 
                    #7jul20; deprec: = yg[0]*
#                p = (zsg[i,jj]/1e3/coefMalamud)**expMalamud / A_hs_km2 #deduced at p2019_42d
                    #does pAslid allows study pdfs as Malamud04? (no: #>1 per subbasin). 
                    #pAslid derived from Malamud04's eqn 18(NGuinea), which is tropical instead of eqn 19(NZ).            
#                p = max(1e-8, p); p = min(.99, p)
#                pAslid[i] = p
#                cps= ms/(1-p);  cpw = mw*cps
#                cp = [cps, cpw];  cp = np.array(cp)  #1st component for sed and 2nd for wat.
#                yg = cp/(1+cp*p)    # pg 2019_42 in notes.
#                hwg = yg[1]*x3w[i]/soilPor   #water column in soil in gully.
                
                Qdh = o3w[i]/86400*Adr_slid  #[o3w]= m/day
                    #Q via soil darcy permeability though crosssection in hollow. p267.
                hwg[i] = np.sqrt(hollowSidSlop*Qdh/ks)/soilPor
                    #ks is hztal, while hgy kssp is vertical (~dm/day). if k is darcy perm, ks=k*dh/dx.
#                if Q[i]>0.1:
#                    print( 'Q[i],o3w[i],ks,Qdh,hwg:', Q[i],o3w[i],ks,Qdh,hwg[i] )
#                    sys.exit('check hw')
                
                hw[i,jj] = min(zsg[i,jj], hwg[i])
    
                #get bedload 'wash' load. Assume no transient sediment store inside drainage network of subbasin.
                if j <= self.j_mx_show and j >= self.j_mn_show and i<=self.i_mx_show and i>=self.i_mn_show:
                    print()           
                    print(f'reach{j}, day{i}')
                if hwg[i] > zsg[i,jj]: #6jul20; wash was deprect to better analyze slids. p266.
                    pSat += 1
#                    Wh = fWh*p
#                    q[i] = Q[i]/Wh    #[m2/s] #pg2019_44.
                    
    #                p_bed_wash = xsh[i-1,:]/tot_soil;  #q_eff = max(0,q[i]-qc)
    #    #            qs_Rick01[i] = 5.8/100*q_eff*washSlop**2.     #[m2/s]

    #                qs_mx_teor[i] = cwash*8e3*q[i]**2.0*slidSlop**1.7  #coef deduced for qs and q in m2/s at pg 2019_50. 
    #                    #cwash repair dt: q daily might be at least 25% of instantaneous, so qs might be 1/4^2~ 6% (1/.06= 17)!!
    #                washmx = p_bed_wash*qs_mx_teor[i]*86400    #m2/d    #/100 due to 2 OM below flume eqn (Rickenmann01).
    #                    #18jun20; rick01 eqn 9 was likely discarded due to the uncertain qc.
    #                washmx_denud = washmx*Wh/A_hs  #[m/d]
    #                wash_sh[i, :] = np.minimum(washmx_denud, xsh[i-1,:])            
                    if j <= self.j_mx_show and j >= self.j_mn_show and i<=self.i_mx_show and i>=self.i_mn_show:
#                        print(f'Wh{Wh}, slidSlop{slidSlop}')
                        qr= round(q[i],2)  #; qcr= round(qc, 2)
                        qsr= round(qs_mx_teor[i],6)
                        print(f'surface&subsurf runoff q = {qr}m2/s')  #, qc{qcr}')
                        print(f'qs_mx = {qsr}m2/s')  #, qc{qcr}')                                        
    #            if j <= self.j_mx_show and j >= self.j_mn_show and i<=self.i_mx_show and i>=self.i_mn_show:
    #                xp = ['%1.3e'%a for a in xsh_prelim]
    #                print(f'after wash, xsh_prelim{xp}')
                            
                #get supply from lansdlide material. Assume no store, as for wash. 
                    #Take care: BendaDunne97 have buffer (store).
                    #yg is ratio of sediment[0] or water[1] depth in gully compared to the hillslope (
                        #gullies + rills) mean.
    #            xsh_prelim = np.zeros((ncv))
                xsh_prelim = xsh[i-1,:,jj] - wash_sh[i,:]        #6jul20; wash was turned off.
                hW_net = hw[i, jj]*soilPor;  hS_net = zsg[i,jj]*solidPort
                moist = gw*hW_net/(gs*hS_net)
                gm[i] = f_gm*(moist + 1)
#                print('jj,co[jj]: ' jj,co[jj])
                FS[i,jj] = (co + (gm[i]*zsg[i,jj] - gw*hw[i,jj])*cosS*tanFi) / (gm[i]*zsg[i,jj]*sinS) #no cos2 but cos, 
                hcr[i,jj] = co / ( cosS*( gwtf + gm[i]*S_exc ) ) #eq3 odorico03
                hmx[i,jj] = co / ( gm[i]*cosS*S_exc )   #eq4 odorico03
                    #as h,H are normal to slope. see Odorico03 (this allows rectangVolCalcs)
                if FS[i,jj]<1: #lansdlide is triggered
#                    slid_sh[i,:,jj] = p*(yg[0]*xsh_prelim) #all gully material falls.
#                    yr = 1/(1-p)*(1-yg*p) #p41cr
                    As_sx = zsg[i,jj]**2/hollowSidSlop #7jul20; p266dl. triangle prism in hollow.
                    b_ind = 2*zsg[i,jj]/hollowSidSlop #7jul20; ind mean individual slide.
                    L_ind = 50*zsg[i,jj]**2
                    As_pl = L_ind*b_ind; As[i,jj] = As_pl
                    asp[i,jj] = L_ind/b_ind
                        #individual slide; expon from vol~ A**1.3 ~ LarsenMontgomery10.
                        #define coef to accomplinsh aspect l/b ratio e.g. taylorMalamud
                    L_slid = slids_mx_perSubb*L_ind
                        #7jul20 (at this time at this cohesion) (p265r).
                    slid_sh[i,:,jj] = L_slid*(As_sx-Aseed)/86400 #7jul20 [m3/s]
#                    print( 'zsg[i,jj], b_ind, L_ind, As_pl, As_sx, Aseed, slid_sh[i,:,jj]: ',
#                        zsg[i,jj], b_ind, L_ind, As_pl, As_sx, Aseed, slid_sh[i,:,jj] )
#                    sys.exit('check slid vol')

#                    xsh_prelim = yr[0]*xsh_prelim*(1-p)     #rill material remains in hillslope.
#                    print(f'xsh_prelim: {xsh_prelim}')
                    xsh_prelim = seedLayer*np.ones(nDh)    #7jul20; 1cm remains as seed for weathering.
                    slidFlg[i,jj] = 1
                else:
                    slid_sh[i,:,jj] = 0                            
                if j <= self.j_mx_show and j >= self.j_mn_show and i<=self.i_mx_show and i>=self.i_mn_show:
                    moist_r = round(moist, 2); 
#                    zsr= round(zsg[i], 2); hwir= round(hwg, 3)
#                    FSr= round(FS[i], 3); pasr= round(p, 6)
                    print(f'moist{moist_r}, gm{gm}')
#                    print(f'zsg{zsr}, yg{yg}, hwg{hwir}, FS{FSr}, pAslid{pasr}')
                
#            #slid suply ensemble due to uniform random nabla coh inside subbasin, to widen Aslid spectra per event.    
#            hW_net = hw[i]*soilPor; hS_net = zsg[i]*solidPort
#            moist = gw*hW_net/(gs*hS_net)
#            gm[i] = f_gm*(moist + 1)
#            FS[i] = (coh[0] + (gm[i]*zsg[i] - gw*hw[i])*cosS*tanFi) / (gm[i]*zsg[i]*sinS) #no cos2 but cos, 
#            hcr[i] = coh[0] / ( cosS*( gwtf + gm[i]*S_exc ) ) #eq3 odorico03
#                #as h,H are normal to slope. see Odorico03 (this allows rectangVolCalcs)
#            if FS[i]<1: #lansdlide is triggered
#                slid_sh[i,:] = p*(yg[0]*xsh_prelim) #all gully material falls.
#                yr = 1/(1-p)*(1-yg*p) #p41cr
#                xsh_prelim = yr[0]*xsh_prelim*(1-p)     #rill material remains in hillslope.
#            else:
#                slid_sh[i,:] = 0                            
#            if j <= self.j_mx_show and j >= self.j_mn_show and i<=self.i_mx_show and i>=self.i_mn_show:
#                moist_r = round(moist, 2); 
#                zsr= round(zsg[i], 2); hwir= round(hwg, 3)
#                FSr= round(FS[i], 3); pasr= round(p, 6)
#                print(f'moist{moist_r}, gm{gm}')
#                print(f'zsg{zsr}, yg{yg}, hwg{hwir}, FS{FSr}, pAslid{pasr}')

            
            
            
                #sediment replenishment:                
#                recov_sh = cwsh*np.exp(-xsh_prelim/xsh_ref)   #pelletier16 globalDz.asc says xsh_ref~.5m
                recov_sh = .5*kdif/xsh_prelim   #7jul20; dif is the unique parameter. *dt=1day
                
                
                
                #mass balance of hillslope sediment stores:
                xsh[i, :, jj] = xsh_prelim + recov_sh            
                if j <= self.j_mx_show and j >= self.j_mn_show and i<=self.i_mx_show and i>=self.i_mn_show: 
#                    xsh_now = [round(x, 2) for x in xsh[i, :]]
                    recov_sh = [np.format_float_scientific(x, precision=2, exp_digits=1) for x in recov_sh]            
                    print(f'recov_sh{recov_sh}')
#                    print(f'xsh_now{xsh[i, :]}')
    #                print(f'wash_sh{wash_sh[i, :]}')
    #                print(f'slid_sh{slid_sh[i, :]}')                
            
        #save supply as sediment discharge in m3/s for fluvial model.                
        Q[Q == 0] = 1e-6    #just for avoid ch = #/0= NaN.
        f = A_hs/86400;  fQ = f/Q
#        wash_sh = f*wash_sh;  slid_sh = slid_sh*f   #from [m/day] to [m3/s]
#        wash_sh[wash_sh < 1e-10] = 1e-10;  slid_sh[slid_sh < 1e-10] = 1e-10
        wash_aggrD = np.sum(wash_sh, axis=1); slid_aggrD = np.sum(slid_sh[:,:,0], axis=1)
        ch_wash = wash_aggrD*fQ;  qs_wash = ch_wash*q
        ch_slid = slid_aggrD*fQ    #sum for all grain sizes.        
        ch = ch_wash + ch_slid;  ch [ch < 1e-6] = 1e-6
#        wash_sh[0, :] = 1e-10;  slid_sh[0, :] = 1e-10


        #main report to review coherence:
        if j <= self.j_mx_show and j >= self.j_mn_show:
            Qm = '%1.3e'%np.mean(Q)
            print ("pSat= ", pSat/ntd)
            print(f'Qmean[m3/s]{Qm}')                    
            print ("daily wash [m3/s]=", wash_aggrD.mean())
            print ("daily slid [m3/s]=", slid_aggrD.mean())
            print ("min slid_mperday=", min(slid_aggrD))
            print ("max slid_mperday=", max(slid_aggrD))
            print ("daily c=", ch.mean())
#            print ("zsg mean", zsg.mean())
            print('next reach j -------------------------------------------------')
            print()
            print()
        print(f'reach{j} was finished')            
#        np.savetxt(f'{output_folder}/results_hillS_'+str(self.n)+'steps.txt',results_hillS,fmt='%d')        
        
        if (j >= self.j_mn_show and j <= self.j_mx_show):
            self.psc.supply_plotter(j, A_hs, FS, hcr, hmx, gm, pAslid, zsg, hw, hwg, o2w, o3w, ch, 
                ch_wash, ch_slid, Pt, Q, xsh[:,:,0], wash_sh, slid_sh, qs_mx_teor, q, asp, As, ncv,
                slidFlg)
        
        return wash_sh, slid_sh, Q, slidFlg, hcr, zsg, FS, slids_mx_perSubb
