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


class P:  #required to call one function in other
#functions contained:
#xt_genP;t_hourly_Huff, t_daily_Salas; x_casc, x_smooth; x_zoom; basin_vector

#to fix/review: 
#P>0
#coherence of rain(dt) pdf and #zeros
#coherence of area reduction factor
#save fields with date in xt_genP


	def __init__(self,dT,dt,eventT,p,N,M,smoothFact):
		#global variables, shared in several functions in this class:
		self.roSm=[]
		self.dT=dT	#[hours]
		self.dt=dt	#[hours]
		self.eventT=eventT	#must be integer?, http://www.eltiempo.com/colombia/otras-ciudades/las-ciudades-en-las-que-mas-llueve-en-colombia-129436
		self.p=p
		self.N=N
		self.M=M
		self.smoothFact=smoothFact	#factor of aggregation of field. Power of 2 is recommended.
		self.n=dT/dt

#		print "n=",self.n
		
		
	def xt_genP(self):
		#notes by pen, 2018, pg12, area reduction factor calibrated:
		ns=self.N-math.log(self.smoothFact,2)
		nt_day=9	#9<nt<17, after calibration by Catano, with daily parameterization
		k=-math.log(float(self.dt)/24,2)	#increase nt if dt finer than 1 day
		nt=nt_day+k
		print 'relative spatial averaging degree: nt-ns=', nt-ns
		x=nt-ns	#comparison between temporal and spatial resolution of binary cascade (Veneziano et al, 2006)
		y=min(0,-.53*x-.26)	#from fig 5b, (Veneziano et al, 2006)
		ARF=2**y
		print 'ARF=', ARF
		#time amplifier of dimensionless rain fields:
		self.Pt=self.t_daily_Salas()
#		for k in range (0,self.n):
		for k in self.t_ev: #produce field only if rain occurs in t
			intMean=self.Pt[k-1]	#[mm] accumulated in dt
			print "field at t=",k			
			print "intMean=",intMean
			Px=self.x_casc()
			while Px.mean()==0:	#repeat spatial cascade if null field is produced
				Px=self.x_casc()
			Px=ARF*intMean*Px
			sides=Px.shape[0]
			print "sides_field_xtProcess=",sides
			for i in range(sides):
				for j in range(sides):
					Px[i,j]=int(round(Px[i,j]))
			print "Px_tGiven_integer=",Px.astype(int)
			np.savetxt('doct_rain/Px_day'+str(k*self.dt/24)+'.txt',Px,fmt='%d')	
				
					
	def t_daily_Salas(self):
		'''started by Santiago Cataño. April 25th, 2017. From MatLab
		#goal: P[m] by seasonal parameters + random Pareto
		#dt=time step in hours, 1<dt<24
		#Dt=time window in hours
		#Dt/dt must be integer
		#ref: fig21 Salas2013MSc'''

		k=-.02*self.dt+.914	#shape parameter
		sigma=.015*self.dt**1.706	#variance parameter
		theta=-.037*self.dt+.099	#(miu) mean parameter

		#seasonality:
		#P=max(0,a*(1+b*sin(w.*t))+normrnd(0,a*cvn,1,n));

		#following 2 lines are example from
		#https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.genpareto.html#scipy.stats.genpareto
		##r = genpareto.rvs(c, size=1000)
		#“frozen” continuous RV object: rv = genpareto(c, loc=0, scale=1)
		
		#generate GPrandom per time step
#		P=[]
#		for i in range(0,n):		
#			P=001*max(0,gp(k,size=n)); #[m]
		Pt=[]
		self.t_ev=[1]
		n_ev=0	#number of events. It starts in 0, as index
		t=0
		while self.t_ev[n_ev]<=self.n:	#poisson arrival of event and random genPareto for rain value (Salas2013MSc)
			addRain=genpareto.rvs(k,loc=theta,scale=sigma,size=1)	#[mm], to consider efficient integer values for files
			addRain=int(round(addRain[0]))
			addRain=max(1,addRain)	#trick: to record notable rain[mm] in integer series			
			Pt.append(addRain)
			t_arriv=np.random.poisson(lam=self.eventT)		
			t_arriv=max(1,t_arriv)
			self.t_ev.append(self.t_ev[n_ev]+t_arriv)	#time of next event to be simulated			
#			print "event #=",n_ev
#			if self.t_ev[n_ev]>self.n:
#				t_arriv=self.n-self.t_ev[n_ev-1]	#trick: close series with rain
#			print "t_arriv=",t_arriv 
#			addZeros=[0 for i in range(0,t_arriv-1)]
#			print "addZeros=",addZeros
			if t_arriv>1:
				nZeros=min(t_arriv-1,self.n-self.t_ev[n_ev])
#				print "nZeros=",nZeros
				[Pt.append(0) for i in range(0,nZeros)]	#do not exceed required series lenght
#			print "Pt=",Pt
#			print "lenght Pt=",len(Pt)
#			print "t_event_plusNext=",self.t_ev
#			t=t+1
			n_ev=n_ev+1	#increase number of events
#			print "--------------------------"			
		if len(Pt)<self.n:
			[Pt.append(0) for i in range(len(Pt)+1,self.n+1)]
#			print "filled with 0s at end of array"
		if self.t_ev[len(self.t_ev)-1]>self.n:
			del self.t_ev[-1]	#event times until series lenght
#		print "t_event",self.t_ev
		np.savetxt('doct_rain/t_ev.txt',np.transpose(self.t_ev),fmt='%d')
		self.Pt=np.array(Pt)
#		print "Pt=",self.Pt
#		print "mean_Pt=",self.Pt.mean()
#		self.Pt=theta+sigma*self.Pt
		#plot Pt frequencies
#		fig,ax=plt.subplots(1,1)
#		ax.hist(self.Pt)
#		plt.show()	
		#plot series, longer steps to auto format date axis
#		#Use PANDAS to translate 1D array to series, assigning date in daily format:
		#assume reasonable # of timesteps
		if self.dt>12:	#e.g. daily during year(s) 
			fr=str(self.dt/24)+'D'	#will bring bugs when dt!=24, because pandas allow only certain integer frequencies
			res='days'
		else:
			fr=str(self.dt)+'H'
			res='hours'
#		print "fr=",fr
		dates=pd.date_range('20/2/2018',periods=self.n,freq=fr)
#		print dates
#		df=pd.DataFrame({res:dates,'values':self.Pt})
#		print df
#		df[res]=pd.to_datetime(df[res])
#		print df
#		df=df.set_index(res)

		Pserie=pd.Series(self.Pt,index=dates)	#note capital S to call function
#		print "Pserie",Pserie	
#		print "dates=",dates	
		#plot series, longer steps to auto format date axis
		
		fig,ax=plt.subplots()	#one image and file per run, no overlap		
		Pserie.plot(kind='bar')
#		ax.plot(df.index,df['values'],kind='bar')
#		loc=mdates.AutoDateLocator()
#		ax.xaxis.set_major_locator(loc)	#exist MonthLocator(interval=number)...
#		ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(loc))		
		
#		ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
		
#		ax.xaxis_date()		

		fig.autofmt_xdate()
#		plt.show()
		plt.savefig('doct_rain/P'+str(self.dt)+'h.png')
		plt.clf
		return self.Pt		
		
		
#		figure
#		set(gcf,'color','w');
#		subplot(2,2,[1 2]),plot(P),xlabel(['time steps=' num2str(dt) 'h']),ylabel('P(m)');
#		x=flipud(sort(P));
#		p=1/n:1/n:1;p=p';
#		subplot(2,2,3),loglog(p,x),xlabel('exced.probability'),ylabel('P(m)');
#		subplot(2,2,4),loglog(p,x),xlim([1e-5 1e-1]),ylim([1e-3 1e-0]),grid on,xlabel('exced.probability'),ylabel('P(m)');

#		#mean annual P [m]:
#		Pm_t=mean(P)
#		Pma=Pm_t*8760/dt
#		CV=Pm_t/std(P)
#		p_zeros=(n-nnz(P))/n
	
	
	def x_casc(self):
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
#			print "k=",k 
			for k2 in np.arange(1,x[k]+1):	#each iter to refine field resolution 
#				sides=b**(k2/d)
				#see note 1 in this function
				sides=4**(k2/2)	
				W=np.zeros((sides,sides))
				#weigths W, E[W]=1, bynary
				for i in np.arange(0,sides): #matrix starts in 				i=0&j=0
					for j in np.arange(0,sides):
						if rn.random()<self.p:
							W[i,j]=0	
						else:
							W[i,j]=1./(1-self.p)
				#apply weigths to last -coarser- density 				map	
#				print "k2=",k2
#				print "dims ro=",ro.shape
#				print "dims W=",W.shape
				ro=np.multiply(ro,W)

				#resample density map to increase resolution
				if k2<x[k]:
					ro=im.zoom(ro,zoomFact,order=0)
			if k==0:
				roN=ro
#				print "mean_roN=",roN.mean()
				#plot rain field and save as image
				fig=plt.figure()	
				plt.imshow(roN,interpolation='nearest',cmap='gray')
				plt.colorbar()
				plt.savefig('doct_rain/roN.png')
				plt.clf
			else: #given M<N
				roM=ro
#				print "mean_roM=",roM.mean()			
				fig=plt.figure()	
				plt.imshow(roM,interpolation='nearest',cmap='gray')
				plt.colorbar()
				plt.savefig('doct_rain/roM.png')
				plt.clf
				roM_zoom=im.zoom(ro,zoomFact**(self.N-self.M),order=0)
				#integrate maps of coarse and fine resolution
				ro=np.multiply(roN,roM_zoom)
				fig=plt.figure()	
				plt.imshow(ro,interpolation='nearest',cmap='gray')
				plt.colorbar()
				plt.savefig('doct_rain/ro.png')	
#				print "mean_ro=",ro.mean()
#				print "ro=",ro
		#get smoothed rain field and save as image
		self.roSm=self.x_smooth(ro)
#		print "roSm=",self.roSm
		#check mean value of smoothed field:
		meanSm=self.roSm.mean()
		#rescale roSm to get mean=1 in field
		if meanSm>0:
			self.roSm=self.roSm/meanSm
		meanSm=self.roSm.mean()
		print "mean_roSm=",meanSm	#check if =1
		#plot rescaled, smoothed field and save as image
		fig=plt.figure()	
		plt.imshow(self.roSm,interpolation='nearest',cmap='gray')
		plt.colorbar()
		plt.savefig('doct_rain/roSm.png')
		plt.clf
		#Closes smooth. Leaves roSm available for other function		
		return self.roSm	
		
			
	def x_smooth(self,field):
		'''goal: get more diverse values in field through upscaling.  			proposed by Catano, S.
		Parameters: 
		field=2D array'''
		
		#number of sides of the field
		dim=field.shape
#		print "bin_field_dim=",dim[0]
		dimSm=dim[0]/self.smoothFact
#		print "smooth_field_dim=",dimSm
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
#		print "roSm",self.roSm
		return self.roSm	#Closes smooth. Leaves roSm available for other function
		




					
		#save values density maps
	
		#np.savetxt(ro) 
		
		
	#feb 1, 2018
	#resample to finer resolution, as required by basin

	#def resample()

	#mask basin		 		
		
