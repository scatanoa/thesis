#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#mar17, 2018

import numpy as np
from scipy.stats import genpareto
import pandas as pd
import matplotlib.pyplot as plt
import math
import scipy.stats as st


class ideal: 
    def __init__(self, dT, dt, A_hs, output_folder, stoch_rain, rain_fileName):
        self.rain_fileName=rain_fileName
        self.stoch_rain=stoch_rain    
        self.output_folder= output_folder
        self.dT=dT    #[hours]
        self.dt=dt    #[hours]
        self.n=dT/dt
        self.A_hs=A_hs    #'hillslope' area
        print ("n= ", self.n)
        print ("instantiation from Fluv..py to fictBasin.py, output_folder= ", output_folder)
        
        
    def t_daily_Salas(self):    #copied from rain.py module. Old comments were erased.        
    
        output_folder=self.output_folder
    
        clim_change=1
        if clim_change==0:
            params=[.1,5.5,1,4000]  #RFG2018 presentation, default: 1st 3 params must be 1 according to Salas2015 Msc
#             cclim: var 6 instead 5.5
            days_per_yr_P=1.3*(params[3])**.65
        else:
            params=[.1,7,1,4000]  #RFG2018 presentation, default: all 1
            days_per_yr_P=1.0*(params[3])**.65            
        #https://www.climatestotravel.com/climate/colombia
        #now: days_per_yr_P ~ 1.3*(P[mm/y])**.65; P<4e3mm/y. Hypoth: 1.3 change to 1 with increased variance due to climate change.
        eventT=365/days_per_yr_P
        k=params[0]*.4    #.4, dt day   #(-.02*self.dt+.914)    #shape parameter.
        sigma=params[1]*3.      #3, dt day.    015*self.dt**1.706#/2    #variance parameter
        theta=params[2]*(-.7)     #-.7, dt day  (-.037*self.dt+.099)    Salas2015 Msc
        Pt=[]
        t_ev=[1]
        n_ev=0    #number of events. It starts in 0, as index
        t=0
        while t_ev[n_ev]<=self.n:
            addRain=genpareto.rvs(k,loc=theta,scale=sigma,size=1)    #[mm]
            addRain=int(round(addRain[0]))
            addRain=max(1,addRain)    #trick: to record notable rain[mm] in integer series            
            Pt.append(addRain)
            t_arriv=np.random.poisson(lam=eventT)        
            t_arriv=max(1,t_arriv)
            t_ev.append(t_ev[n_ev]+t_arriv)    #time of next event to be simulated            
            if t_arriv>1:
                nZeros=min(t_arriv-1,self.n-t_ev[n_ev])
                [Pt.append(0) for i in range(0,nZeros)]    #do not exceed required series lenght        
            n_ev=n_ev+1    #increase number of events
        if len(Pt)<self.n:
            [Pt.append(0) for i in range(len(Pt)+1,self.n+1)]
        if t_ev[len(t_ev)-1]>self.n:
            del t_ev[-1]    #event times until series lenght
        np.savetxt(f'{output_folder}/t_ev.txt',np.transpose(t_ev),fmt='%d')
        #plot
        self.Pt=np.array(Pt)
        Pserie=pd.Series(self.Pt, index=pd.date_range('10/1/1990', periods=self.n)) #note capital S to call function [month, d, y]
        fig,ax=plt.subplots(nrows=2,ncols=1)    #one image and file per run, no overlap        
        ax[0].plot(Pserie,'k')    #kind='bar'
        ax[0].set_ylabel('P[mm/day]')
#        ax.set_yscale('log')
#        fig.autofmt_xdate()        
        ii=[]
        for i,p in enumerate(self.Pt):
            ii.append(i)
        ii=np.array(ii)
        ii=ii/self.n     #array between 0 and 1
        ax[1].plot(np.sort(self.Pt)[::-1],ii,'k')
        ax[1].set_xlabel('P[mm/day]')
        ax[1].set_ylabel('prob_exced')
        ax[1].set_yscale('log')
#        ax[1].set_xscale('log')   #TO COMPARE: https://www.hydrol-earth-syst-sci.net/14/2559/2010/hess-14-2559-2010.pdf
        Pm=self.Pt.mean()*24*365/self.dt
        CV_P=np.std(self.Pt)/self.Pt.mean()
        Pm=int(Pm)
        plt.savefig(f'{output_folder}/P'+str(self.n)+'days_Pm'+str(Pm)+'_10CV'+str(int(10*CV_P))+'.png')
        plt.clf
        np.savetxt(f'{output_folder}/P'+str(self.n)+'steps.txt',np.transpose(self.Pt),fmt='%d')
        print ("P[mm/y]= ", Pm)
        print ("CV_P= stdP/meanP= ", CV_P)
        return self.Pt        
        
    def hill_w_t(self,etp,hu,kssp,kb,tsp,tssp,tb): #parameter units [L=mm, T=dia]    
    #tropicalVals(Marinilla): 1-4, 20-600(100), 1-100(20), .01-10(4), 1-10(1), 1-10(3), 50-200(80)
    
        output_folder=self.output_folder
        rain_fileName=self.rain_fileName
    
        x1w=np.zeros((self.n))
        x2w=np.zeros((self.n))
        x3w=np.zeros((self.n))
        x4w=np.zeros((self.n))
        o1w=np.zeros((self.n))
        o2w=np.zeros((self.n))
        o3w=np.zeros((self.n))
        Q=np.zeros((self.n)) # m3/s
        x1w[0]=0.  # water at residual soil [0-.2]*1e3 mm
        x2w[0]=0.   # water for fast runoff   
        x3w[0]=0.   # water at saprolite[.1-1] *1e3 mm
        x4w[0]=.1e3     # water at fractured rock(.1-2) *1e3 mm
        #produce rain input
        if self.stoch_rain==1:
            P=self.t_daily_Salas()        
        else:
            P=np.loadtxt(f'{output_folder}/{rain_fileName}')  
            #saved timeseries, 6 years (1990-1996) to calibrated San Carlos(Col) sed wave
            
#        P.astype(float)    #to produce more variable output
        for i in range(1,self.n):        #hillslope hidrology
            i1w=min(P[i]*(1-x1w[i-1]/hu)**2,hu-x1w[i-1])
            x1w[i]=x1w[i-1]+i1w
            o1w[i]=min(x1w[i],etp*(x1w[i]/hu)**.6)
            x1w[i]=x1w[i]-o1w[i]
            i2w=max(0,P[i]-i1w-kssp)
            x2w[i]=x2w[i-1]+i2w
            o2w[i]=x2w[i]/tsp   
            x2w[i]=x2w[i]-o2w[i]
            i3w=max(0,P[i]-i1w-i2w-kb)
            x3w[i]=x3w[i-1]+i3w
            o3w[i]=x3w[i]/tssp    
            x3w[i]=x3w[i]-o3w[i]
            i4w=max(0,P[i]-i1w-i2w-i3w)
            x4w[i]=x4w[i-1]+i4w
            o4w=x4w[i]/tb
            x4w[i]=x4w[i]-o4w
            Q[i]=o2w[i]+o3w[i]+o4w    #mm/d
            balance=x1w[i]+x2w[i]+x3w[i]+x4w[i]-(x1w[i-1]+x2w[i-1]+x3w[i-1]+x4w[i-1])+(o1w[i]+o2w[i]+o3w[i]+o4w)-P[i]    #must be 0
        Q=Q/1e3*self.A_hs/(3600*self.dt)    #m3/s
        Q[0]=.05*(self.A_hs/1.e6)    #average yield as artificial starting water discharge for sediment simulation.
        print ("in hill, P=", P)
        print ("daily mean ETR=", o1w.mean())
        print ("daily timeseries of ETR=", o1w)
        print ("from fictBasin py, Qhill= ", Q)
        #Q duration curve
        
        #plot        
        t=np.arange(self.n)        
        fig,ax=plt.subplots(nrows=4,ncols=1)
        
        ax[0].set_title('mean(Q)='+str(round(Q.mean(),2))+', mean(x2w)='+str(x2w.mean())+\
        ', var(x2w)='+str(np.var(x2w)),fontsize=10)
        
        ax[0].plot(t,o1w,label='mean(o1w)='+str(round(o1w.mean(),1)))
        ax[0].set_ylabel("o1w[mm]")        
        ax[0].legend()
        ax[1].plot(t,x1w)
        ax[1].set_ylabel("x1w[mm]")
        ax[2].plot(t,x3w)
        ax[2].set_ylabel("x3w[mm]")
        for k in range(1,3+1):
            ax[k-1].get_xaxis().set_visible(False)
        ax[3].plot(t,x4w)
        ax[3].set_ylabel("x4w[mm]")
        ax[3].set_xlabel("t [days]",fontsize=10)
        plt.savefig(f'{output_folder}/hill_water_stores.png')            
#        self.Q=Q
#        self.o2w=o2w
#        self.o3w=o3w
#        self.x3w=x3w
        results_hillW=(Q,o2w,o3w,x3w)        
        return results_hillW
    
    def hill_s_t(self,washSlop,slidSlop,cwash1,qc,pAwash,geotech):
        #fb=fictBasin.ideal(24*365*100,24,1.e5)
        #fb.hill_s_t(.2,.5,1,1e-5,.003,[.3,1.8e4,1.e4,.6,1e-3/365/3,1e-3/365/3,1e-3/365/3,1e-3/365/3,1.,.9,.8,.7])
        #cwash1 and pAwash are correction factors due to aggregation
        #correction factor 1 due to aggregation
        
        output_folder=self.output_folder
        
        soilPor=geotech[0]
        gm=geotech[1]
        coh=geotech[2]
        tanFi=geotech[3]
        cw1sh=geotech[4]
        cw2sh=geotech[5]
        cw3sh=geotech[6]
        cw4sh=geotech[7]
        x1sh_ref=geotech[8]
        x2sh_ref=geotech[9]
        x3sh_ref=geotech[10]
        x4sh_ref=geotech[11]
        pAslid=pAwash
        gw=1.e4
        gs=2.65*gw
        slidSlop=np.arctan(slidSlop)
        sinS=np.sin(slidSlop)
        cosS=np.cos(slidSlop)
        cwash2=cwash1    #/10 due to size, but *10 due to exposure. Non-cohesive wash load sizes are only .1mm and 1mm 
        Lhill=math.sqrt(self.A_hs)
        Whill=Lhill    #assume aspect ratio=1 (basin is a square)
        #------------------------------------------
        #initial conditions
        self.n=int(self.n)
        x1sh=np.zeros((self.n))
        x2sh=np.zeros((self.n))
        x3sh=np.zeros((self.n))
        x4sh=np.zeros((self.n))
        hw=np.zeros((self.n))
        zsh=np.zeros((self.n))
        ch=np.zeros((self.n))
        ch_wash=np.zeros((self.n))
        ch_slid=np.zeros((self.n))
        FS=np.zeros((self.n))
        wash_1sh=np.zeros((self.n))
        wash_2sh=np.zeros((self.n))
        slid_1sh=np.zeros((self.n))
        slid_2sh=np.zeros((self.n))
        slid_3sh=np.zeros((self.n))
        slid_4sh=np.zeros((self.n))
        q=np.zeros((self.n))
        qs_wash=np.zeros((self.n))
        fsq=np.zeros((self.n))
        #initial storages of sediment in hillslope, by grain size (can be up to 5, increasing 1 OM from .1mm)        
        x1sh[0]=.01*(1-soilPor)    #initial total soil depth= .4m
        x2sh[0]=.01*(1-soilPor) 
        x3sh[0]=.01*(1-soilPor)
        x4sh[0]=.01*(1-soilPor)
        #produce hillslope hidrology input
        results_hillW=self.hill_w_t(4.,100.,20.,4.,1.,3.,80.)
        Q=results_hillW[0]
        o2w=results_hillW[1]/1e3    #mm to m
        o3w=results_hillW[2]/1e3
        x3w=results_hillW[3]/1e3
        pSat=0.
        secs_per_dt=3600*self.dt
        for i in range(1,self.n):
#            print ("i= ", i)
            q[i]=(o2w[i]+o3w[i])*self.A_hs/pAwash/Whill/(3600*self.dt)    #m2/s. o3w 'emerges' to surface
            if q[i]<qc:
                q[i]=qc
#            fsq[i]=25500*washSlop**1.67*q[i]**2.    #eqn KilincRich, cited by Montoya2008 in dissertation
            fsq[i]=5.8/100*max(0,q[i]-qc)*washSlop**2.    #/100 due to 2 OM below flume eqn (Rickenmann01). 
            #Take care, eqn 9 there is better in steep
            
            p1_bed_wash=x1sh[i-1]/(x1sh[i-1]+x2sh[i-1])
            p2_bed_wash=1-p1_bed_wash
            washmx1=p1_bed_wash*cwash1*fsq[i]    #m2/s if cw~2e4~KilincRich
            washmx2=p2_bed_wash*cwash2*fsq[i]
            washmx1=washmx1/Lhill*3600*self.dt    #m/d
            washmx2=washmx2/Lhill*3600*self.dt
            wash_1sh[i]=max(0,min(x1sh[i-1],washmx1))
            wash_2sh[i]=max(0,min(x2sh[i-1],washmx2))
#            if washmx1>0:
#                wash_1sh[i]=washmx1/(3.14/2)*np.arctan(x1sh[i-1]/washmx1)
#            else:
#                wash_1sh[i]=0
#            if washmx2>0:
#                wash_2sh[i]=washmx2/(3.14/2)*np.arctan(x2sh[i-1]/washmx2)
#            else:
#                wash_2sh[i]=0
            zsh[i]=(x1sh[i-1]+x2sh[i-1]+x3sh[i-1]+x4sh[i-1])/(1-soilPor)    #soil evolves slower
            pAslid=pAwash*zsh[i]**2     #depth of soil to be in slide increases lenght and width of surficial area of slide
            hw_int=x3w[i]/pAslid/soilPor
            hw[i]=zsh[i]/(3.14/2)*np.arctan(hw_int/zsh[i])    #'slidable(~high_permeab?)' soil thickness acts as maximum water column
            if hw_int>zsh[i]:
                pSat+=1
            FS[i]=(coh+(gm*zsh[i]-gw*hw[i])*cosS**2*tanFi)/(gm*zsh[i]*cosS*sinS)
            if FS[i]<1:
                slid_1sh[i]=.95*(x1sh[i-1]-wash_1sh[i])                       
                slid_2sh[i]=.95*(x2sh[i-1]-wash_2sh[i])   
                slid_3sh[i]=.95*x3sh[i-1]
                slid_4sh[i]=.95*x4sh[i-1]
            recov_1sh=1./pAwash*cw1sh*np.exp(-x1sh[i-1]/x1sh_ref)
            recov_2sh=1./pAslid*cw2sh*np.exp(-x2sh[i-1]/x2sh_ref)
            recov_3sh=1./pAslid*cw3sh*np.exp(-x3sh[i-1]/x3sh_ref)
            recov_4sh=1./pAslid*cw4sh*np.exp(-x4sh[i-1]/x4sh_ref)
            x1sh[i]=x1sh[i-1]-wash_1sh[i]-slid_1sh[i]+recov_1sh
            x2sh[i]=x2sh[i-1]-wash_2sh[i]-slid_2sh[i]+recov_2sh
            x3sh[i]=x3sh[i-1]-slid_3sh[i]+recov_3sh
            x4sh[i]=x4sh[i-1]-slid_4sh[i]+recov_4sh
            if Q[i]==0:
                Q[i]=1e-10    #just to make sure ch element is not #/0=NaN
            ch_wash[i]=(wash_1sh[i]+wash_2sh[i])*pAwash*self.A_hs/(Q[i]*3600*self.dt)
            ch_slid[i]=(slid_1sh[i]+slid_2sh[i]+slid_3sh[i]+slid_4sh[i])*pAslid*self.A_hs/(Q[i]*3600*self.dt)
            ch[i]= ch_wash[i] + ch_slid[i]
            if ch[i]<1e-6:
                ch[i]=1e-6
            qs_wash[i]=ch_wash[i]*q[i]                
            wash_1sh[i]=max(1e-10,wash_1sh[i]*pAwash*self.A_hs/secs_per_dt)  # from [m/day] to [m3/s]
            wash_2sh[i]=max(1e-10,wash_2sh[i]*pAwash*self.A_hs/secs_per_dt)
            slid_1sh[i]=max(1e-10,slid_1sh[i]*pAslid*self.A_hs/secs_per_dt)
            slid_2sh[i]=max(1e-10,slid_2sh[i]*pAslid*self.A_hs/secs_per_dt)
            slid_3sh[i]=max(1e-10,slid_3sh[i]*pAslid*self.A_hs/secs_per_dt)
            slid_4sh[i]=max(1e-10,slid_4sh[i]*pAslid*self.A_hs/secs_per_dt)
        wash_1sh[0]=1e-15
        wash_2sh[0]=1e-15
        slid_1sh[0]=1e-15
        slid_2sh[0]=1e-15
        slid_3sh[0]=1e-15    
        slid_4sh[0]=1e-15
        #plot fluxes
        t=np.arange(self.n)        
        fig,ax=plt.subplots(nrows=2,ncols=1)
        ax[0].plot(t,FS)
        ax[0].set_ylabel("FS")
        ax[0].set_yscale('log')
        ax[1].plot(t,zsh,'k',label='zsh')
        ax[1].plot(t,hw,':b',label='hw')
        leg = ax[1].legend(loc=2,fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.7)        
        ax[1].set_ylabel("slide layer [m]",fontsize=6)
        plt.savefig(f'{output_folder}/soilFailures.png')            
        fig,ax=plt.subplots(nrows=3,ncols=1)
        ax[0].plot(t,ch,'k',label='total')
        ax[0].plot(t,ch_wash,'b',label='wash')
        ax[0].plot(t,ch_slid,'r',label='slide')        
        ax[0].set_ylabel("ch [m3/m3]")
        ax[0].set_yscale('log')
        ax[0].set_ylim(1e-4)
        ch_slid_mean='%1.1e'%np.mean(ch_slid)
        ch_wash_mean='%1.1e'%np.mean(ch_wash)
        ax[0].set_title(f'means: ch_slid= {ch_slid_mean}, ch_wash= {ch_wash_mean}')
        leg = ax[0].legend(loc=2,fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.7)
        
        ax[1].plot(t,slid_1sh,'b',label='1')
        ax[1].plot(t,slid_2sh,'r',label='2')
        ax[1].plot(t,slid_3sh,'m',label='3')
        ax[1].plot(t,slid_4sh,'k',label='4')
        ax[1].set_ylabel("slid_sh [m3/s]")            
        ax[1].set_yscale('log')
        leg = ax[1].legend(loc=2,fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.3)
        ax[2].plot(t,wash_1sh,'b',label='1')
        ax[2].plot(t,wash_2sh,'r',label='2')
        ax[2].set_ylabel("wash_sh [m3/s]")            
        ax[2].set_yscale('log')
        leg = ax[2].legend(loc=2,fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.3)
        plt.savefig(f'{output_folder}/Qs_hill.png')            
        
        #plot water flow
        fig,ax=plt.subplots()
        ax.plot(t,Q)
        ax.set_ylabel("Qhill_day [m3/s]")
#        ax[0,1].legend()            
        plt.savefig(f'{output_folder}/Qhill_day.png')
                
        #plot states
        fig,ax=plt.subplots()
        ax.plot(t,x1sh,'b',label="0.1mm")
        ax.plot(t,x2sh,'r',label="1mm")
        ax.plot(t,x3sh,'m',label="10mm")
        ax.plot(t,x4sh,'k',label="100mm")
        ax.set_ylabel("xsh [m]")
        leg = ax.legend(loc=2,fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.7)
#        ax[0,1].legend()            
        plt.savefig(f'{output_folder}/Xs_hill.png')
        #plot rating curve
        fig,ax=plt.subplots(nrows=2,ncols=1)
        ax[0].plot((q-qc)*washSlop**1.5,qs_wash*3600*2650,'ok',markersize=.5)    #to compare with figure 12 in Zimmermann2013, e.g. supp.lim.
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        ax[0].set_title("xlabel: (q-qc)S^1.5 [m3/s/m]")
        ax[0].set_ylabel("qs wash [kg/m/h]")
        ax[1].plot(q,qs_wash,'ok',markersize=.5)    #to compare with figure 12 in Zimmermann2013, e.g. supp.lim.
        ax[1].set_xscale('log')
        ax[1].set_yscale('log')
        ax[1].set_xlabel("q [m2/s]")
        ax[1].set_ylabel("qs_wash [m2/s]")        
        plt.savefig(f'{output_folder}/HillRatingCurve.png')
        #main report to review coherence
        print ("pSat= ", pSat/self.n)
        print ("daily c=", ch.mean())
        print ("zsh mean", zsh.mean())
        print ("Lhill= ", Lhill)            
#        np.savetxt(f'{output_folder}/results_hillS_'+str(self.n)+'steps.txt',results_hillS,fmt='%d')
        return wash_1sh,wash_2sh,slid_1sh,slid_2sh,slid_3sh,slid_4sh,Q
