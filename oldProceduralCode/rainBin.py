#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#feb 6, 2018

import numpy as np
import scipy.ndimage as im
import matplotlib.pylab as plt

from cu import *
#from wmf import wmf


class Pread:

	def __init__(self,dT,sideCoarse,dxM,nodata):
		#dxM: spatial resolution of Map
		self.dT=dT	#[hours]
		self.sideCoarse=sideCoarse	#side lenght, #cells, of coarse rain field
		#load basin file [files,cols,drain2] and get size required for rain fields to be used in hydrology model:
		self.basin=np.loadtxt('doct_DEM/ArenBasin.txt',dtype='int')
		Nrows=max(self.basin[:,2])-min(self.basin[:,2])
		Ncols=max(self.basin[:,1])-min(self.basin[:,1])
		Ncells=self.basin.shape[0]
		print 'Nrows=', Nrows
		print 'Ncols=', Ncols
		print self.basin		
		self.Nrows=Nrows
		self.Ncols=Ncols
		self.Ncells=Ncells
		
	def Pxt_fine(self):
		#basin mask to be applied to rain field in each time:
		inBasin=self.mask()	#ones just in pixels inside basin
		#downscaling
		Lmx=max(self.Nrows,self.Ncols)
		zoomFact=Lmx/self.sideCoarse
		print 'zoomFact=', zoomFact
		t_ev=np.loadtxt('doct_rain/t_ev.txt',dtype='int')
		t_ev=t_ev.tolist()
		for i,t in enumerate(t_ev):
			Pcoarse=x=np.loadtxt('doct_rain/Px_day' + str(t) + '.txt',dtype='int')
			#downscaling
			Pfine=im.zoom(Pcoarse,zoomFact,order=0)
			Lx=inBasin.shape[1]
			Ly=inBasin.shape[0]
			#clip squared field to become hztal or vertical rectangle:
			if Lx<Lmx:	#less columns than rows in basin mask
				clip=Lmx-Lx
				Pfine=Pfine[:,:-clip]	#clip unnecessary rain columns
			else:	#viceversa, then clip unnecessary rain rows
				clip=Lmx-Ly
				Pfine=Pfine[:-clip,:]
			#apply mask
			Pfine=Pfine*inBasin
			np.savetxt('doct_rainBin/PxDownscal_day'+str(t)+'.txt',Pfine,fmt='%d')	
			fig=plt.figure()	
			plt.imshow(Pfine,interpolation='nearest',cmap='gray')
			plt.colorbar()
			plt.savefig('doct_rainBin/PxDownscal_day'+str(t)+'.png')
			plt.clf
#			print 't=', t
#			print Pcoarse
			#map to basin, basin is like vector of topology, which is compatible with SHIA:
#			wmf.Transform_Map2Basin(Pfine,[self.Ncols,self.Nrows,,,])
			vecPfine=cu.basin_map2basin(self.basin,self.Ncells,Pfine,1,self.Ncols,self.Ncols,self.Nrows,
			dxM,nodata)			

	def mask(self):
		X=np.zeros((self.Nrows,self.Ncols))
		print 'shape X =', X.shape
		print 'Ncells=', self.Ncells
		for k in range(0,self.Ncells):
			print 'cell#', k
			#cell coordinates vary as 0<=x<=Ncols-1 and 0<=y<=Nrows-1
			i=self.basin[k,2]-1			
			j=self.basin[k,1]-1
			print 'i', i
			print 'j', j
			X[i,j]=1
		fig=plt.figure()	
		plt.imshow(X,interpolation='nearest',cmap='gray')
		plt.colorbar()
		plt.savefig('doct_rainBin/basinMask.png')
		plt.clf					
		return X  
		
	
	
		
		
		
		
	def	vectBasin2bin(self,vec=None,ruta_out=None,fecha=None,dt=None,status='update',umbral = 0.01):	
	#modified from method rain_radar2basin_from_array in module wmf.py, developed by Nicolas Velazquez. 
	#https://github.com/scatanoa/WMF/blob/master/wmf/wmf.py
	
		'Descripcion: Genera campos de lluvia a partir de archivos array\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : .\n'\
		'vec: Array en forma de la cuenca con la informacion.\n'\
		'ruta_out: Ruta donde escribe el binario con la lluvia.\n'\
		'fecha: Fecha del registro actual.\n'\
		'dt: Intervalo de tiempo entre registros.\n'\
		'status: Estado en el cual se encuentra el binario que se va a pasar a campo.\n'\
		'	update: (Defecto) Con este estado se genera un binario y se agregan nuevos campos.\n'\
		'	old: Estado para abrir y tomar las propiedades de self.radar.. para la generacion de un binario.\n'\
		'	close: Cierra un binario que se ha generado mediante update.\n'\
		'	reset: Reinicia las condiciones de self.radar... para la creacion de un campo nuevo.\n'\
		'Retornos\n'\
		'----------\n'\
		'Guarda el binario, no hay retorno\n'\
		'meanRain :  La serie de lluvia promedio.\n'\
		'\n'\
		'Mirar Tambien\n'\
		'----------\n'\
		'rain_interpolate_idw: interpola campos mediante la metodologia idw.\n'\
		'rain_radar2basin_from_asc: Mete campos de lluvia mediante multiples arrays.\n'\
		#Edita la ruta de salida 
		if ruta_out <> None:
			if ruta_out.endswith('.hdr') or ruta_out.endswith('.bin'):
				ruta_bin = ruta_out[:-3]+'.bin'
				ruta_hdr = ruta_out[:-3]+'.hdr'
			else:
				ruta_bin = ruta_out+'.bin'
				ruta_hdr = ruta_out+'.hdr'
		#Establece la cantidad de elementos de acuerdo al tipo de cuenca
#		if self.modelType[0] is 'c':
#			N = self.ncells
#		elif self.modelType[0] is 'h':
#			N = self.nhills
#			try:
#				if vec.shape[0]  == self.ncells:
#					vec = self.Transform_Basin2Hills(vec,sumORmean=1)		
#			except:
#				pass
		# De acuerdo al estado actualiza las variables o guarda el binario final 
		actualizo = 1
		if status == 'update':
			#Entrada 1 es la entrada de campos sin lluvia 
			if len(self.radarDates) == 0:
				models.write_int_basin(ruta_bin,np.zeros((1,N)),1,N,1)
			if vec.mean() > umbral:
				#Actualiza contador, lluvia media y pocisiones 
				self.radarCont +=1
				self.radarMeanRain.append(vec.mean())
				self.radarPos.append(self.radarCont)
				#Guarda el vector 
				vec = vec*1000; vec = vec.astype(int)
				models.write_int_basin(ruta_bin,np.zeros((1,N))+vec,
					self.radarCont,N,1)
				actualizo = 0
			else:
				#lluvia media y pocisiones 
				self.radarMeanRain.append(0.0)
				self.radarPos.append(1)
			self.radarDates.append(fecha)
		#Si ya no va a agregar nada, no agrega mas campos y genera el .hdr 
		elif status == 'close':
			self.radarMeanRain = np.array(self.radarMeanRain)
			self.radarPos = np.array(self.radarPos)
			#Guarda un archivo con informacion de la lluvia 
			f=open(ruta_hdr[:-3]+'hdr','w')
			f.write('Numero de celdas: %d \n' % self.ncells)
			f.write('Numero de laderas: %d \n' % self.nhills)
			f.write('Numero de registros: %d \n' % self.radarMeanRain.shape[0])
			f.write('Numero de campos no cero: %d \n' % self.radarPos.max())
			f.write('Tipo de interpolacion: radar \n')
			f.write('IDfecha, Record, Lluvia, Fecha \n')
			c = 1
			for d,pos,m in zip(self.radarDates,
				self.radarPos,self.radarMeanRain):
				f.write('%d, \t %d, \t %.2f, %s \n' % (c,pos,m,d.strftime('%Y-%m-%d-%H:%M')))
				c+=1
			f.close()
			#Vuelve las variables listas de nuevo 
			self.radarMeanRain = self.radarMeanRain.tolist()
			self.radarPos = self.radarPos.tolist()
		elif status == 'reset':
			#Variables de radar
			self.radarDates = []
			self.radarPos = []
			self.radarMeanRain = []
			self.radarCont = 1
		elif status == 'old':
			#si es un archivo viejo, lo abre para tomar las variables y continuar en ese punto 
			f=open(ruta_hdr[:-3]+'hdr','r')
			Lista = f.readlines()
			self.radarCont = int(Lista[3].split()[-1])
			f.close()
			#Abre con numpy para simplificar las cosas 
			a = np.loadtxt(ruta_hdr,skiprows=6,dtype='str').T
			self.radarPos = [int(i.split(',')[0]) for i in a[1]]
			self.radarMeanRain = [float(i.split(',')[0]) for i in a[2]]
			for i in a[3]:
				d = datetime.datetime.strptime(i,'%Y-%m-%d-%H:%M')
				self.radarDates.append(d)
		return actualizo 

