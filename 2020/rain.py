    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#jan30, 2018

import numpy as np
import random as rn
import scipy.ndimage as im
import matplotlib.pylab as plt
from scipy.stats import genpareto
import pandas as pd
import matplotlib.dates as mdates
import math
import time
import sys

import mapPlot as mp


class P:  #required to call one function in other
    def __init__(self,dT,dt,eventT,p,N,M,smoothFact, R0flag, Pm_mmyr, a_zcit):
        #global variables, shared in several functions in this class:
        self.roSm=[]
        self.dT=dT	#[hours]
        self.dt=dt	#[hours]
        self.eventT=eventT
        self.Pm_mmyr = Pm_mmyr	
        self.a_zcit = a_zcit
            #must be integer?, 
            #http://www.eltiempo.com/colombia/otras-ciudades/las-ciudades-en-las-que-mas-llueve-en-colombia-129436
        self.p=p
        self.N=N
        self.M=M
        self.smoothFact=smoothFact	#factor of aggregation of field. Power of 2 is recommended.
        self.n=dT/dt       
        self.output_folder = 'doct_rain'
        self.R0flag = R0flag

	
    def xt_genP(self):
        #notes by pen, 2018, pg12, area reduction factor calibrated:
        
#        ns=self.N - math.log(self.smoothFact,2)        
#        kk=-math.log(float(self.dt)/24,2)	#increase nt if dt finer than 1 day        
#        nt_day=9	#9<nt<17, after calibration by Catano, with daily parameterization
#        nt=nt_day+kk
#        print('relative spatial averaging degree: nt-ns=', nt-ns)
#        x=nt-ns	#comparison between temporal and spatial resolution of binary cascade (Veneziano et al, 2006)
#        y=min(0,-.53*x-.26)	#from fig 5b, (Veneziano et al, 2006)
#        ARF=2**y
#        print('ARF=', round(ARF, 2))
        ARF = 1

        #time amplifier of dimensionless rain fields:
        if self.R0flag == 0:
            nd = int(self.n)
            self.Pt = self.Pm_mmyr/365*np.ones(nd)
        else:
            print('va pa salas')
            self.Pt=self.t_daily_Salas(self.Pm_mmyr)
            #or analytical sinusoidal, which Over95 solves for event (intraday) rain fields. 
#        #en sgte loop llamarás fortran pa construir bin files de rain fields. 22DEC2019: DEPRECATED COS F90 HAS NOT IM.ZOOM().
#        fields.Pxt()
        
        last_days = 0
        last_day = 0
        t1 = time.time()
        cum_t = 0
        xml = []
        for k in self.t_ev: #produce field only if rain occurs in t
            if k == self.t_ev[-1:][0]: last_day = 1
            if k > self.t_ev[-5:][0]: last_days = 1   #22dec2019: plot fields4only (effic) last 5 events, to review spatial peaks.
            intMean=self.Pt[k-1]	#[mm] accumulated in dt
            print('t stepping[days] ', self.dt, ', field at t=', k)
            print(f'last_day = {last_day}')
            print(f'last_days = {last_days}')
            print("Pmean[mm]= ", intMean)
            Px=self.x_casc(last_day)
            while Px.mean()==0:	#repeat spatial cascade if null field is produced
                print('null field was produced')
                Px=self.x_casc(last_day)
            xm = Px.mean()
            Px = Px/xm #5jun20; star at p242u; forced normalization to 1 so rand P source in (t,x) is (Salas, p_casc).
            xm = Px.mean()            
            print( "5jun20; cascBalancCheck; must be 1; Px_mean[adim]= ", xm)
            xml.append(xm)
            Px=ARF*intMean*Px
            sides=Px.shape[0]
            for i in range(sides):
	            for j in range(sides):
		            Px[i,j]=int(round(Px[i,j]))
            np.savetxt('doct_rain/Px_day'+str(int(k*self.dt/24))+'.txt',Px,fmt='%d')
#            np.save('doct_rain/Px_day' + str(int(k*self.dt/24)) + '.npy', Px)
            t2 = time.time()    #for 1 year daily N=10, M=3, smoothFact=2**3; save/savetxts take 110/120s. 
                #Run also shortened by 1s if only needed programs are open: no word nor firefox: bah!
            dt = t2-t1
            print(f'run_dt[s] = {dt}')
            cum_t += dt
            print(f'run cum_t = {cum_t} segs.')
            print()
            t1 = t2
            if last_days==1: mp.mapPlotter(Px, 'doct_rain', 'Px_day' + str(int(k*self.dt/24)), 'P(mm)',
                'easting', 'northing')
                #redundant, as Pxt se ve en folder supply.
        xml = np.array(xml)
        xmlm = xml.mean()
#        print(xml)
        print( "5jun20; cascBalancCheck; must be 1 with ~0 desv: mean(xm), std(cv): ", xmlm, np.std(xml) )
            
    def t_daily_Salas(self, P_mmyr):    #copied from rain.py module. Old comments were erased.        
        #17dic2019: esta es la version mas actual a la q me refiero en la version deprecated de este archivo.
    
        output_folder=self.output_folder
    
        clim_change=0
        if clim_change==0:
            params=[-3, 5, 1, 0]  #RFG2018 presentation, default: 1st 3 params must be 1 according to Salas2015 Msc            
#             cclim: var 6 instead 5.5
            Pm = params[3]  #mm/yr
#            days_per_yr_P=1.3*Pm**.65 #30dic2019: deprecated if casc_fields are to be produced, where p accounts for 0days.
            days_per_yr_P = 365   #goal, after applying casc_field, is to have 150-200 daysP/year, i.e. rain falls each 2 days.
        else:
            params=[.1,7,1,4000]  #RFG2018 presentation, default: all 1
            days_per_yr_P=1.0*(params[3])**.65            
        #https://www.climatestotravel.com/climate/colombia
        #now: days_per_yr_P ~ 1.3*(P[mm/y])**.65; P<4e3mm/y. Hypoth: 1.3 change to 1 with increased variance due to climate change.
        eventT=365/days_per_yr_P
        k=params[0]*.4    #.4, dt day   #(-.02*self.dt+.914)    #shape parameter.
        sigma=params[1]*3.      #3, dt day.    015*self.dt**1.706 #/2    #variance parameter
        theta=params[2]*(-.8)     #-.8, dt day  (-.037*self.dt+.099)    Salas2015 Msc
        Pt=[]
        t_ev=[1]
        n_ev=0    #number of events. It starts in 0, as index
        t=0
        Pmd = P_mmyr/365
        Pmdr = Pmd*eventT   #6jul20; mean rain in rainy days; grow is more t-spaced.
        print(f'Pmd={Pmd}')
#        a_zcit = self.a_zcit #p260r #depr 6jul20; zcit sinus set exp pdf shape of Pt of area average.
#                #(i) that, and given..
#                    #(ii)isagen notes hgy saying bimodal is x-t sensible, by geogr & ENSO, 
#                    #(iii) stefa tucker odorico refs simplifying rain as memoryless exp pdf.
#                    #SO I SET EXP PDF. no for time between events, as these are created by x-0's of casc
#                        #casc better than odorico00 x-poisson to set easier tested correl x-0's vs P x_mean.
#                            #casc is set for amazon by salas.
##        sys.exit(f'stop to rev a_zcit={a_zcit}')
#        T_zcit = 183
#        T_zc_sin = 2*np.pi/T_zcit
        while t_ev[n_ev]<=self.n:
#            td = t_ev[n_ev]        
##            addRain=genpareto.rvs(k,loc=theta,scale=sigma,size=1)    #[mm]
##            addRain=int(round(addRain[0]))
##            addRain=max(1,addRain)    #trick: to record notable rain[mm] in integer series            
#            zcit = Pmd*( 1 - a_zcit*np.sin(T_zc_sin*td) ) #27jun20; may improve via seasonal p_casc?
#            print(f'td={td}')
#            print(f'zcit={zcit}')

#            addRain = rn.uniform(1 , zcit*2) #27jun20; uniform instead random, to include rain min = 1mm (mean bias).

            addRain = rn.expovariate(1/Pmdr)
                    
            addRain = int(addRain) #27jun20; round avoided; who cares if rain was 3 or 2mm/day?
                #27jun20; p260ul; random instead pareto for eficiency, as this value just amplifies map of whole basin.
            

            Pt.append(addRain)
            if days_per_yr_P == 365: 
                t_arriv = 1
#                sys.exit('so no poisson')
            else:
                t_arriv=np.random.poisson(lam=eventT)        
                t_arriv=max(1,t_arriv)
            print(f't_arriv={t_arriv}')                
            t_ev.append(t_ev[n_ev]+t_arriv)    #time of next event to be simulated            
            if t_arriv>1:
                nZeros = int(min(t_arriv-1, self.n-t_ev[n_ev]))
                [Pt.append(0) for i in range(0,nZeros)]    #do not exceed required series lenght        
            n_ev=n_ev+1    #increase number of events
        if len(Pt)<self.n:
            [Pt.append(0) for i in range(len(Pt)+1,self.n+1)]
        if t_ev[len(t_ev)-1]>self.n:
            del t_ev[-1]    #event times until series lenght
        np.savetxt(f'{output_folder}/t_ev.txt',np.transpose(t_ev),fmt='%d')
        self.t_ev = t_ev
        #plot
        self.Pt=np.array(Pt)
        
        Pserie=pd.Series(self.Pt, index=pd.date_range('10/1/1990', periods=self.n)) #note capital S to call function [month, d, y]
        
        fig,ax=plt.subplots(nrows=2, ncols=1)    #one image and file per run, no overlap        
        ax[0].plot(Pserie,'k')    #kind='bar'
        ax[0].set_ylabel('P[mm/day]')
#        ax.set_yscale('log')
#        fig.autofmt_xdate()        
        ii=[]
        for i,p in enumerate(self.Pt):
            ii.append(i)
        ii=np.array(ii)
        ii=ii/self.n     #array between 0 and 1
        ax[1].plot(ii, np.sort(self.Pt)[::-1], 'k')
#        ax[1].set_xlabel('P[mm/day]')
        ax[1].set_xlabel('prob_exced')
        ax[1].set_ylabel('P[mm/day]')
#        ax[1].set_xscale('log')
#        ax[1].set_yscale('log')
#        ax[1].set_xscale('log')   #TO COMPARE: https://www.hydrol-earth-syst-sci.net/14/2559/2010/hess-14-2559-2010.pdf
        Pm=self.Pt.mean()*24*365/self.dt
        CV_P=np.std(self.Pt)/self.Pt.mean()
        Pm=int(Pm)
        str_P = str(int(self.n))
        plt.savefig(f'{output_folder}/P'+ str_P +'days_Pm'+str(Pm)+'_10timesCVeq'+str(int(10*CV_P))+'.png')
        plt.clf
        
        np.savetxt(f'{output_folder}/P'+ str_P + 'steps.txt',np.transpose(self.Pt),fmt='%d')
        print ("P[mm/y]= ", Pm)
        print ("CV_P= stdP/meanP= ", round(CV_P, 2))
        return self.Pt

	
    def x_casc(self, last_day):
#		R0,L0 are not used arguments yet

        '''	goal: multiplicative binary cascade, Over y Gupta(1994)
        and smoothing to get more than two values of rain intensity.
        by: Santiago Cataño Alvarez.
        Parameters:
        p: controls binary probability distribution.
        b: branches, e.g. use 2 for rainfall field
        d: dimension, e.g. use 2 for rainfall field
        R0: spatially averaged rainfall.
        L0: map side lenght.
        N, float: number of iterations -'ni'- to get finer resolution.
        M, float (<N, suggested: M=N-3): 'ni' to get coarser resolution.
        agrFact: factor of aggregation of field. 
        Use power of 2 because sides number of cascade field is too'''

        #		zoomFact=b**(1./d)
        #note1: loss generality without free b and d, but gain robustness for using code as rain field generator.
        zoomFact=2
        #frequency: N=high, M: low
        x=[self.N,self.M]
        for k in [0,1]:	#each spatial frequency: N=fine, M=coarse
            ro=np.ones((2,2))
            for k2 in np.arange(1,x[k]+1):	#each iter to refine field resolution 
#				sides=b**(k2/d)
				#see note 1 in this function
                sides = int(4**(k2/2))
                W=np.zeros((sides,sides))
                #weigths W, E[W]=1, bynary
                for i in np.arange(0,sides): #matrix starts in 				i=0&j=0
                    for j in np.arange(0,sides):
                        if rn.random()<self.p:
                            W[i,j]=0	
                        else:
                            W[i,j]=1./(1-self.p)
				#apply weigths to last -coarser- density 				map	
                ro=np.multiply(ro,W)
                if k2<x[k]: #resample density map to increase resolution
                    ro=im.zoom(ro,zoomFact,order=0)
            if k==0: 
                roN = ro
                if last_day==1: mp.mapPlotter(roN, 'doct_rain', 'roN','roN', 'easting', 'northing')
            else: #given M<N
                roM=ro
                if last_day==1: mp.mapPlotter(roM, 'doct_rain', 'roM','roM', 'easting', 'northing')
                #				print "mean_roM=",roM.mean()			
                roM_zoom=im.zoom(ro,zoomFact**(self.N-self.M),order=0)
                ro=np.multiply(roN,roM_zoom) #integrate maps of coarse and fine resolution
                if last_day==1: mp.mapPlotter(ro, 'doct_rain', 'ro','ro','easting', 'northing')
                
#        #get smoothed rain field and save as image
#        self.roSm=self.x_smooth(ro)
#        #		print "roSm=",self.roSm
#        #check mean value of smoothed field:
#        meanSm=self.roSm.mean()
#        #rescale roSm to get mean=1 in field
#        if meanSm>0:
#            self.roSm=self.roSm/meanSm
##        meanSm=self.roSm.mean()
##        print ("mean_roSm=",meanSm)	#check if =1
#        if last_day==1: mp.mapPlotter(self.roSm, 'doct_rain', 'roSm') 		#plot rescaled, smoothed field and save as image.
#        #Closes smooth. Leaves roSm available for other function		

        return ro	
		
			
    def x_smooth(self,field):
        '''goal: get more diverse values in field through upscaling.  			proposed by Catano, S.
        Parameters: 
        field=2D array'''

        #number of sides of the field
        dim=field.shape
        dimSm = int(dim[0]/self.smoothFact)
        self.roSm=np.zeros((dimSm,dimSm))  #smoothed field
        #fill smoothed field. u=up, d=down, l=left, r=right
        for i in np.arange(1,dimSm+1): #rows
            iu=(i-1)*self.smoothFact+1
            ido=i*self.smoothFact
            for j in np.arange(1,dimSm+1): #cols
                jl=(j-1)*self.smoothFact+1
                jr=j*self.smoothFact
                #following indexes are translated from 'real' to 				python indexing standard
                self.roSm[i-1,j-1]=field[iu-1:ido,jl-1:jr].mean()
        return self.roSm	#Closes smooth. Leaves roSm available for other function	
