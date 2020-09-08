    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#mar19, 2018

#pre-existing modules
import numpy as np
import matplotlib.pyplot as plt
#'home made' Python modules
import fictBasin
#'home made' Fortran modules
from mod_model_sc import *
import read_model_results as rmr
from scipy import stats

#from byte2strarray import *


class idealStreams: 
    def __init__(self,dT,dt_minu_ev,dt_hours_postEv,flume_else_basin=0):
        #dt_minu_ev [minutes] is fine resolution during event. It must be divisor of 60; e.g.: 1,2,3,5,10,12,20,30,60
        
        #dt_hours_postEv [hours] is semi-fine resolution in day after event. It must be divisor of 'coarse chunker' dt,
        #e.g. if dt=24 then dt_hours_postEv can be: 1,2,3,4,6,8,12.
        
        #PARAMETER DEFINITION
        if flume_else_basin==1:
            self.output_folder='trench_mw2' #must be in home folder, so you avoid to use '/'
        else:
            self.output_folder='doct_fictBasin_cclim'        
        self.dT=dT    #[hours], integer
        self.dt_dQ_coarse=2    #[hours], integer, default=24. It denotes time period in which Q changes.
        self.n_coarse_dt=int(dT/self.dt_dQ_coarse)
        self.nW=9   #default nW is 5, if basin. Can vary from 1 to 9 (from June8, 2018) if no flume. If jamming, recommended nW8
        self.nD=8  #Move between 2 and 8 if flume_else_basin=1, else leave 5
        trenchFlag=1
        self.stochGSD=1 #turn ON to produce percentile values to be manually written as GSD array. Then, turn OFF to run model.
        self.Dmean=np.log(7.3e-3)   #statistical parameters of Shawn's dissertation. Coarser bed than Maria&Tobias
        self.Dstdg=np.log(2.1) #equivalent to sigma. If GSD is lognormal then D84=sigma*D50. 15jul2019, edited from 2.5 to get maxD as ltgreen~20mm
        #call spatial topology
        if flume_else_basin==1:
            Sof=.02    #maria: .022
            self.amplif_D = 1    #for sensitivity analysis in 3rd paper with Marwan, comparing flume and numeric model of lateral input channel.
            L_flume=1*9.  #maria: 12
            repos_angle_ini = 30.
            base_level = .2
            self.Wflume = 1
            self.dx_flume=1    #maria: 1
            triangl_base = self.Wflume/2.
            rectang_store = base_level*self.Wflume*L_flume
            triang_store = .5*triangl_base**2*np.tan(3.14/180*repos_angle_ini)*L_flume
            self.ini_store_m3 = 2*triang_store + rectang_store
            self.ini_store = 2650*(1-.3)*self.ini_store_m3
#            print('rectang_store[m3]', rectang_store, 'triang_store', triang_store, 'ini_store[kg_sed]', self.ini_store)
            n_nodes=int(L_flume/self.dx_flume) #L must be multiple of dx             
            self.X=self.serial_topol(n_nodes,self.dx_flume,Sof)
        else:
            stoch_rain=1
            rain_fileName='P36500steps.txt'  
            wash_flag=1
            self.A_hs=4e6    #'hillslope' (subbasin) area
            self.Omega=3    #order of drainage network at outlet
            self.X=self.fictNetwork(1.,2.)
            nR=len(self.X[:,0])     #number of reaches                       
            
        #CALL geometric variables
        self.flume_else_basin=flume_else_basin                    
        Wfeas, Wflow, Dcentr, Dpct, pctD = self.assumed_vars(self.nW, self.Dmean, self.Dstdg, self.nD, self.stochGSD, trenchFlag)
        Dcentr = self.amplif_D*Dcentr
        Dpct = self.amplif_D*Dpct        
        
        self.pctD=pctD
        nD=len(Dcentr)     
        nW=self.nW  

        path_list = self.path_list_builder(nW, nD, flume_else_basin)
        self.path_list=path_list            
                
        #CALL water and sediment supply
        if flume_else_basin==1:
            self.Qhill, self.wash_1sh, self.wash_2sh, self.slid_sh = \
            self.supply_calc_flume(self.n_coarse_dt, n_nodes, trenchFlag, nD, pctD)
        else:
            fb=fictBasin.ideal(self.dT,self.dt_dQ_coarse,self.A_hs,self.output_folder, stoch_rain, rain_fileName)
            wash_1sh_j,wash_2sh_j,slid_1sh,slid_2sh,slid_3sh,slid_4sh,Qhill_j = fb.hill_s_t(.15,\
            np.tan(35*3.14/180),wash_flag,2e-1,.005,\
            [.3,1.8e4,6.5e3,np.tan(34*3.14/180),1e-3/365*1,\
            1e-3/365*2.3,1e-3/365*.8,1e-3/365*.03,.1,.15,.09,.08])
            slid_sh_j=np.stack((slid_1sh,slid_2sh,slid_3sh,slid_4sh))
            slid_sh_j=slid_sh_j.T             
            self.Qhill=np.zeros((self.n_coarse_dt,nR))
            self.slid_sh=np.zeros((self.n_coarse_dt,nR,nD_hill)) #5th=greatest size, boulder, is supplied in f90 fluvial model
            self.wash_1sh=np.zeros((self.n_coarse_dt,nR))
            self.wash_2sh=np.zeros((self.n_coarse_dt,nR))
            for j in range(1,nR+1): #June3, 2018: short time before RFG18 implies to assume no rain gradient
                self.Qhill[:,j-1]=Qhill_j
                self.wash_1sh[:,j-1]=wash_1sh_j
                self.wash_2sh[:,j-1]=wash_2sh_j
                for k in range(1,nD_hill+1):   
                    self.slid_sh[:,j-1,k-1]=slid_sh_j[:,k-1]        
        #CALL dt dynamic             
        self.Wfeas=Wfeas
        self.Wflow=Wflow
        self.Dcentr=Dcentr
        self.Dpct=Dpct
        slid_sh_sum=np.sum(self.slid_sh, axis=2)   #aggregatte load of all grain sizes
        self.nt, self.dt_fine, self.coarse_dt = self.refine_dt(self.Qhill,slid_sh_sum,dt_minu_ev,dt_hours_postEv,self.dt_dQ_coarse)
        self.dt_minu_ev=dt_minu_ev
        self.dt_hours_postEv=dt_hours_postEv
#        print ("slid_sh= ", self.slid_sh)
        np.set_printoptions(precision=3, suppress=True)    
#        print("nD_hill= ", nD_hill)   
        print ("X= ", self.X)        
        print ("Qhill= ", self.Qhill)   
#        print ("flume_feed(n_coarse_dt,nD)= ", self.slid_sh[:,0,:])
#        print ("total feed rate", sum(self.slid_sh[1,0,:]))    
        print("dt[hours]= ", self.dt_fine)
        print("n_coarse_dt=", self.n_coarse_dt) 
        print("nt=", self.nt)
        self.nR = len(self.Wfeas[:,0])
        self.trenchFlag = trenchFlag
        
    def run_model_sc(self,iref_ini, iref_fin, jref_ini, jref_fin):
        #neither wash nor .1mm grain size: 0 initial store, 0 rate of recover    
        np.set_printoptions(precision=3, suppress=False)    
#        print "py, inside setup inputs, in Mphdyn, Q= ", Q                    
        #initialize variables previous to run time steps
        nW=len(self.Wfeas[0,:])
        nR=self.nR
        nD=len(self.Dcentr)
#        np.savetxt('doct_fictBasin/Qhill.txt',self.Qhill)    #,fmt='%d' #saved to be used in plot method
#        print "dt="
#        for    i in range(1,nt+1):
#            print dt[i-1]            
        np.set_printoptions(precision=3, suppress=True)    
#        print "py, inside setup inputs, in Mphdyn, D= ", D                    
        #call topology
#        print "py, inside setup inputs, in Mphdyn, X= ", self.X                    
#        print "nW= ", nW                                    
#        print "nR= ", nR
#        Wfeas=np.asfortranarray(Wfeas)
#        Wflow=np.asfortranarray(Wflow)
#        D=np.asfortranarray(D)        
#        print "X as numpy= ", self.X
#        X=np.asfortranarray(self.X)
#        print "X as f90= ", X
#        dt_fine=1
#        dt_fine.astype(int)

        if self.flume_else_basin==1:
            nD_hillslope=nD
        else:
            nD_hillslope=nD-1
            
        mod_model_sc.model_sc(iref_ini,iref_fin,jref_ini,jref_fin, self.output_folder,\
        self.path_list2f90,self.pos_ini_var,self.flume_else_basin, self.dt_dQ_coarse,self.dt_minu_ev,self.dt_hours_postEv,\
        self.dt_fine,self.wash_1sh,self.wash_2sh,self.slid_sh,self.Qhill, self.Wfeas, self.X, self.Wflow,\
        self.Dcentr,self.Dpct,self.pctD,self.ini_store_m3,self.nt,self.n_coarse_dt, nW, nR, nD, nD_hillslope,\
        self.num_paths,self.nvars,self.num_char)
            #h_kk1, v_kk1, Fr_kk1, Q,E,rel_submg,W_h=
        print("review if effectively results are stored in binary files now ")        

    def supply_calc_flume(self, n_coarse_dt, n_nodes, trenchFlag, nD, pctD):
        Qhill=np.zeros((self.n_coarse_dt,n_nodes))
        if trenchFlag==0: 
            Q2reach=np.ones((self.n_coarse_dt))*5e-3  #maria: 65e-3
        else:
            Q2reach=[4, 6, 8, 10, 14, 18, 24, 28, 32, 32, 14, 8, 7, 6, 4]            
            Q2reach = np.array(Q2reach)
            Q2reach = Q2reach*1e-3            
        #test for 2-fold increase in flume discharge after 2nd day
#            for i in range(self.n_coarse_dt+1-3,self.n_coarse_dt+1):
#                Q2reach[i-1]=Q2reach[i-1]            
        for j in range(1,n_nodes+1):
            if j==1:    #only water input in upstream boundary of flume
                Qhill[:,j-1]=Q2reach
            else:
                Qhill[:,j-1]=Q2reach/1e3 #to initialize water transit in flume, <1%err if 10 spatial nodes
        #daily timesteps imply to pass only constant fractional feed rate: 
        #Intraday pulses and zeros are assigned inside fortran:
        slid_sh=np.zeros((n_coarse_dt,n_nodes,nD))  
        feed2reach=0/2.65                   
        slid_sh[:,0,0]=pctD[0]/100.*feed2reach                        
        for k in range(2, nD+1):
            slid_sh[:,0,k-1]=(pctD[k-1]-pctD[k-2])/100.*feed2reach
        wash_1sh=np.zeros((n_coarse_dt,n_nodes))
        wash_2sh=np.zeros((n_coarse_dt,n_nodes))    
        if (self.flume_else_basin==1 and trenchFlag==1): self.Q2reach = Q2reach  #to plot hydraulic geometry, method plot_hg
        return Qhill, wash_1sh, wash_2sh, slid_sh
        
    def serial_topol(self,n_nodes,dx_flume,Sof):       
        #fill columns 4, 6, 7 and 8 (fortran format) with drain2, Sor, dx and zr respectively
        #flume has 'serial' instead of branched topology        
        X=np.zeros((n_nodes,8))
        X[:,6-1]=Sof  #flume slope as bedrock slope
        X[:,7-1]=dx_flume #dx        
        #cumulate altitude from downstream  
        for i in range(1,n_nodes):
            ii=n_nodes-i
            X[ii-1,8-1]=X[ii,8-1]+X[ii,6-1]*X[ii,7-1]       
        #drain to..?                    
        for i in range(1,n_nodes):  #node in downstream boundary drains to '0' reach
            X[i-1,4-1]=i+1
        return X
                                    
    def refine_dt(self, Qhill, slid_sh_sum, dt_minu_ev,dt_hours_postEv, dt_dQ_coarse):
        #fine timesteps according to Courant criteria when high event (exced.prob~2%)        
        n=int(self.n_coarse_dt)        
        if (self.flume_else_basin==1 and n>1):
            Qref=.9*np.amin(Qhill)  #trick: set Qref such that code decide to produce finest dt always 
        elif (self.flume_else_basin==0):   
        #note!!: code block is to be refined as Qhill vary spatially. Get goal t's when event happens
            refine_pctile=1
            Qref=np.percentile(Qhill,100-refine_pctile)
        t_ev_spatial=[]
        Nreaches=len(self.X[:,0])
        np.set_printoptions(precision=3, suppress=False)    
        print("Qhill= ", Qhill)
        print("Qref= ", Qref)
        slid_sh_sum_LocalMean=np.mean(slid_sh_sum, axis=0)
        print("slid_sh_sum_LocalMean= ", slid_sh_sum_LocalMean)
        for j in range(1,Nreaches+1): 
            t_ev,=np.where(Qhill[:,j-1] >= Qref)
            t_ev2,=np.where(slid_sh_sum[:,j-1] > slid_sh_sum_LocalMean[j-1])
#            print("j=", j)
#            print("t_ev", t_ev)
#            k=1    !next while was deprecated in Jun4, 2018 given huge rain can last >1 day (e.g. Salgar2015, Col)
#            for i,value in enumerate(t_ev):
#                while (i<len(t_ev)-1 and t_ev[i+k]<t_ev[i+k-1]+1):                                
#                    k=k+1                                       
            t_ev_spatial.extend(t_ev)
            t_ev_spatial.extend(t_ev2)                   
#            print("t_ev_noConsecu", t_ev)
        t_ev_spatial = np.array(t_ev_spatial)
        t_ev_spatial = np.unique(t_ev_spatial)
        print("t_ev_spatial", t_ev_spatial)                        
#        print ("Qhill= ", Qhill)
#        print ("Qref= ", Qref)
        dt_fine=[]        
        steps_per_hour=int(60/dt_minu_ev) #[minutes]
        dtfine_list=[dt_minu_ev/60.]*dt_dQ_coarse*steps_per_hour    #units=[h], as required by fortran model. Fine graining events to dt=dt_minu_ev[min] in coarse_dt (e.g. day) which contains the event.
        dtinterm_list=[dt_hours_postEv]*int(dt_dQ_coarse/dt_hours_postEv)  
        i_last_ev=0
        
        coarse_dt=np.zeros((n,2))  #to have relation [day,time_step*] *(at end of day) for plotting series with variable dt
        len_dtfine_list=len(dtfine_list)
        len_dtinterm_list=len(dtinterm_list)
        
        for i in range(1,n+1):    #here i means coarse time chunks (days except for trench experiment -hours)
            coarse_dt[i-1,0]=i
#            if (n==1 or Qhill[i-1]>Qref):    #if nablaP then Qhill must be array and 'if' question must be for max(Qhill, given i, in all X)
            any_reach_ev,=np.where(t_ev_spatial==i-1)
#            print("any_reach_ev=", any_reach_ev)                        
            if (i==1 or len(any_reach_ev>0)):   #any reach with event? 
                i_last_ev=i
                dt_fine.extend(dtfine_list)
                coarse_dt[i-1,1]=1*(0 if i==1 else coarse_dt[i-2,1])+len_dtfine_list
#            elif (i>1 and Qhill[i-2]>Qref):    #Use hourly resolution after event, before going again to daily dt.
            elif (i>1 and i==i_last_ev+1):
                dt_fine.extend(dtinterm_list) 
                coarse_dt[i-1,1]=1*(0 if i==1 else coarse_dt[i-2,1])+len_dtinterm_list
            else: #current Qhill > Q
                dt_fine.append(dt_dQ_coarse)
                coarse_dt[i-1,1]=1*(0 if i==1 else coarse_dt[i-2,1])+1
        print('coarse_dt= ', coarse_dt)
        dt_fine=np.array(dt_fine)
#        dt=dt.T    !no considered to test if fortran read this new array as argument (12Apr2018)
        nt=len(dt_fine)
        return nt, dt_fine, coarse_dt                

    def path_list_builder(self, nW, nD, flume_else_basin):
        #DYNAMIC PATHS TO SAVE AND PLOT KEY TIME SERIES:
        #dictionary. var_name: number_of_copies. Copies are needed when 3rd dimension needs to be plotted, because I only know\
        #method to save 2D arrays in binary format from fortran code.        
        
        #variables 1 to 14. # in ':#' means 3rd dimension of array, additional to t=i and x=j
        #new var to plot: (1) py: add to this dict. (2) f90: add to save bin with its nvar. (3) py: declare, read and plot
        nD_hill= 1*(nD if flume_else_basin==1 else nD-1)  #if basin then boulder is not supplied in hillslope module but in river
        self.nD_hill=nD_hill
        dict_vars={'Q':'unique','h':'nW',\
        'v':'nW',\
        'z':'nW',\
        'Da84':'nW',\
        'Da50':'nW',\
        'Dss50':'nW',\
        'So':'unique',\
        'indW':'unique',\
        'th50a':'unique',\
        'thm84B':'unique',\
        'thratio84':'unique',\
        'c':'nD',\
        'Qhill':'unique',\
        'Qs_feed':'nD_hill'}
        dict_copies={'unique':1, 'nW':nW, 'nD':nD, 'nD_hill':nD_hill}
        print('dict_vars= ', dict_vars)
        vars2plot=[]        #get list from dictionary to iterate with numeric values
        n_copies_name=[]
        for i,j in dict_vars.items():
            vars2plot.append(i)
            n_copies_name.append(j)
        nvars=len(vars2plot)
        n_copies=[]
        for i in range(1,nvars+1):
            n_copies.append(dict_copies[n_copies_name[i-1]])
        print('n_copies= ', n_copies)                    
        len_vars=[len(var) for var in vars2plot]
        len_max_var=max(len_vars)+6 #4: maybe double subindex letter and index greater than 10         
        path_list=[]
        pos_ini_var=[1] #contains as many elements as variables. 
        #Each element says position (f90 format i.e. start=1) of first path related to each variable.        
        fill='o'*(len_max_var-len_vars[0])
        path_list.append(f'{fill}{vars2plot[0]}')        
        for i in range(2,nvars+1):  #each var                        
            pos_ini_var.append(pos_ini_var[i-2]+n_copies[i-2])
            if n_copies_name[i-1]=='unique':            
                fill='o'*(len_max_var-len_vars[i-1])     #Needed given f90 requires fixed lenght of array string elements.
                path_list.append(f'{fill}{vars2plot[i-1]}')
            elif n_copies_name[i-1]=='nW':
                for j in range(1,nW+1):    #each copy of var
                    fill='o'*(len_max_var-len_vars[i-1]-3-(1 if j<10 else 2))   #3 char in '_kk'                    
                    path_list.append(f'{fill}{vars2plot[i-1]}_kk{j}')
            elif n_copies_name[i-1]=='nD':
                for j in range(1,nD+1):    #each copy of var                        
                    fill='o'*(len_max_var-len_vars[i-1]-2-(1 if j<10 else 2))   #2char in '_k'                       
                    path_list.append(f'{fill}{vars2plot[i-1]}_k{j}')                                   
            elif n_copies_name[i-1]=='nD_hill':
                for j in range(1,nD_hill+1):    #each copy of var                        
                    fill='o'*(len_max_var-len_vars[i-1]-3-(1 if j<10 else 2))   #2char in '_k'                       
                    path_list.append(f'{fill}{vars2plot[i-1]}_kh{j}')
        self.num_paths=len(path_list)
        self.num_char=len(path_list[0])
#        path_list2f90=np.empty((self.num_paths,self.num_char),dtype='c') #1st arg is tuple of (#elem,#charPerElem)
        print('just pre pass string array to f90 format, path_list= ', path_list)
        path_list2f90=np.array(path_list,dtype='c').T   #trick: https://stackoverflow.com/questions/22293180/passing-numpy-string-format-arrays-to-fortran-using-f2py/23204241
#        for i in range(1,nvars+1):        
#            path_list2f90[i-1]=path_list[i-1]            
        self.path_list2f90=path_list2f90                        
        self.pos_ini_var=pos_ini_var
        self.nvars=nvars                                                          
        for i in range(1,self.num_paths+1):  #path list in py format is filled with folder, to be used in plotting binaries
            path_list[i-1]= self.output_folder + '/' + path_list[i-1]                
        return path_list
        
    def assumed_vars(self, nW, Dmean, Dstdg, nD, stochGSD, trenchFlag):    
    #GSD is produced in certain key pctiles. Then please choose nD such that 3<nD<6 predefined 
        X=self.X
#        print "in funct 'assumed...'"
        nR=len(X[:,0])
        if self.flume_else_basin==0:    
            Akm2=X[:,4]
    #        Sor(1:nR)=.7e-3;
    #        cw4sc=.1*.5*ones(nR,1);    #Mar23/2018: grain size #4 =100mm must be provided by slides, at 1st experiment, as Hassan said in March2018.
    #        cw5sc=.1*.1*ones(nR,1);x4sca_ref=.1*ones(nR,1);x5sca_ref=.1*ones(m,1);    #Mar23/2018: no boulders are supplied
            Qm=.05*Akm2  #water yield .05 m3/s/km2
            Qbf=2*Akm2**.75   #rbf~2m3/s/km2 while rm~.04. Then: Qbf~50Qm (in fact, 90Qm by A scaling)  
            Wv=6.4*Akm2**.5  #eqn in fig 6.19 Takahashi(2014), coef=6.4, exp=.5
            Wbf=1.8*Akm2**.42  #transport limited, Jimenez2013
                        
            #composite cross section of valley: 1= low flow, 2= bankfull flow (expected), 3= 2*bf, 4= 5*bf, 5= Wv
            Wfeas_menu=np.vstack([.5*np.minimum(Wv,Wbf)*(Qm/Qbf)**.2,\
            np.minimum(Wv,Wbf)*(Qm/Qbf)**.2,\
            np.minimum(Wv,Wbf),\
            np.minimum(Wv,Wbf),\
            np.minimum(Wv,1.5*Wbf),\
            np.minimum(Wv,2*Wbf),\
            np.minimum(Wv,3*Wbf),\
            np.minimum(Wv,5*Wbf),\
            Wv])           
            
            Wfeas = self.Wfeas_calc(Wfeas_menu)
                            
#            Wfeas=Wfeas.T    #feasible flow widths
            nW=len(Wfeas[0,:])

            self.plot_W_feas(self, nW, Akm2, Wfeas)
            Wflow = self.Wflow_calc(Wfeas, nW, nR)
                        
            #grain sizes                           
            D=np.array([1e-4,1e-3,1e-2,1e-1,1e0])    #lower2greater grain sizes. Mar23, 2018: only 3(4) classes after Hassan
            pctD=np.array([5,50,80,95,100])   #uniform mixture of scales just after big event of deposition (say San Carlos, Colombia, 1990)
        else: #flume
#            stochGSD=1 #deterministic GSD (0) given high sensitivity of model avoids to locate bug when I try to reproduce it
            Wflume=self.Wflume
            if trenchFlag==1:
                Wfeas=np.zeros((nR,nW))
#                Wmx_freq = .95   #4jul2019: refine width discretization in values you know occur more often.      
#                Xsec=np.arange(1,nW-1+1)/nW*Wmx_freq
#                Xsec = np.concatenate((Xsec, Wflume*np.ones(1)))
                Xsec=np.arange(1,nW+1)/nW*Wflume
                for j in range(1,nR+1):
                    Wfeas[j-1,:] = Xsec
                Wflow = self.Wflow_calc(Wfeas, nW, nR)
            else:
                Wfeas=np.ones((nR,1))*Wflume                                    
                Wflow=Wfeas
                    #grain sizes
            if nD==2:
                pctD=[50,100]        
            elif nD==3:
                pctD=[50,84,100]
            elif nD==4:
                pctD=[16,50,84,100]    
            elif nD==5:
                pctD=[16,50,65,84,100]    
            elif nD==6:
                pctD=[16,50,65,84,90,100]    
            elif nD==7:
                pctD=[16,30,50,65,84,90,100]
            elif nD==8:
                pctD=[10,16,30,50,65,84,90,100]    
#            elif nD==9:
#                pctD=[10,16,30,50,65,84,90,95,100]    
#            elif nD==10:
#                pctD=[10,16,30,50,65,84,90,95,98,100]
                
            Dsamples=np.random.lognormal(Dmean, Dstdg, 10000000)  #virtual Wolman (#=100) count in logNormal pool of grains                        
            if stochGSD==1:
                D=np.zeros((nD))                         
                for k in range(1,nD+1):
                    if pctD[k-1]<100:
                        D[k-1]=np.percentile(Dsamples,pctD[k-1])            
                    else:
                        D[k-1]=np.percentile(Dsamples,pctD[k-1]-1)  #not really 100, but 98%. Avoid bias compared to flume while preserving total probability to pass to fortran            
            else:   #generated when stochGSD flag was turned on. Parameters miu and sigma of Shawn's dissertation
                #1st list: Shawn Experiment, D50=7.3mm
                #2st list: MariaExperiment, D50=5mm
                #3rd list: AlexExp.  #3sep2019
                if nD==2:
#                    D=[0.00729017, 0.04851702]
#                    D=[0.00564174, 0.04322329]                                          
                    D=[0.008,0.032]
                elif nD==3:
#                    D=[0.00751987, 0.01856935, 0.04761879]
#                    D=[0.00564174, 0.01793813, 0.04322329]
                    D=[0.008,0.017,0.032]
                elif nD==4:
#                    D=[0.00292903, 0.00736181, 0.01810215, 0.04537524]
#                    D=[0.00177475, 0.00564174, 0.01793813, 0.04322329]                                        
                    D=[0.0035,0.008,0.017,0.032]
                elif nD==5:
#                    D=[0.00291723, 0.00720183, 0.01024248, 0.01786403, 0.0477301 ]
#                    D=[0.00177475, 0.00564174, 0.00883529, 0.01793813, 0.04322329]
                    D=[0.0035,0.008,0.017,0.020,0.032]
                elif nD==6:
#                    D=[0.00300337, 0.00737362, 0.01030656, 0.01812114, 0.02342613, 0.0487648 ]
#                    D=[0.00177475, 0.00564174, 0.00883529, 0.01793813, 0.02504979, 0.04322329]
                    D=[0.0035,0.008,0.011,0.017,0.020,0.032]
                elif nD==7:
#                    D=[0.00295069,0.00454375,0.00730662,0.0104271,0.01829968,0.02363937,0.0499629]
#                    D=[0.00177475, 0.00306532, 0.00564174, 0.00883529, 0.01793813, 0.02504979, 0.04322329]
                    D=[0.0035,0.005,0.008,0.011,0.017,0.020,0.032]
                elif nD==8:
#                    D=[0.00227894,0.00297331,0.0046032,0.00731512,0.01024277,0.01849975,0.02376357,0.04848185]
#                    D=[0.0012712, 0.00177475, 0.00306532, 0.00564174, 0.00883529, 0.01793813, 0.02504979, 0.04322329]
                    D=[0.0028,0.0035,0.005,0.008,0.011,0.017,0.020,0.032]   #3sep2019

                D=np.array(D)    
#            D=np.array([2e-3,7e-3,2e-2,6e-2])    #lower2greater grain sizes. Mar23, 2018: only 3 classes after Hassan
        #texture for calculations by size fraction must be based in central value of each histogram bin of GSD
        Dpct=D
        Dcentr=np.zeros((nD))                         
        Dcentr[0]=np.sqrt(1e-3*D[0])   #size such that <= 1% sample is below (GSD Alex, UBC)
        for k in range(2,nD+1):
            Dcentr[k-1]=np.sqrt(D[k-1]*D[k-2])
        print("Wfeas= ", Wfeas)        
        print("Wflow= ", Wflow)
#        print('Dsamples', Dsamples)
        print('D_pct= ', Dpct)
        print('D_centr= ', Dcentr)
        print('pctD= ',pctD)                
        return Wfeas, Wflow, Dcentr, Dpct, pctD        
        
    def Wfeas_calc(self, Wfeas_menu):            
        Wfeas_menu[2]=Wfeas_menu[1]+.5*(Wfeas_menu[3]-Wfeas_menu[1])            
        if nW==1:
            Wfeas_chosen=[4]
        elif nW==2:
            Wfeas_chosen=[4,9]
        elif nW==3:
            Wfeas_chosen=[2,4,9]
        elif nW==4:
            Wfeas_chosen=[2,4,6,9]
        elif nW==5:
            Wfeas_chosen=[2,4,6,7,9]
        elif nW==6:
            Wfeas_chosen=[2,4,6,7,8,9]
        elif nW==7:
            Wfeas_chosen=[2,4,5,6,7,8,9]
        elif nW==8:
            Wfeas_chosen=[2,3,4,5,6,7,8,9]
        elif nW==9:            
            Wfeas_chosen=[1,2,3,4,5,6,7,8,9]
        Wfeas=np.zeros((nR,nW))            
        for kk in range(1,nW+1):
            Wfeas[:,kk-1]=Wfeas_menu[Wfeas_chosen[kk-1]-1]                                    
        return Wfeas
        
    def Wflow_calc(self, Wfeas, nW, nR):
        #width occuped effectively by the flow in each portion of valley
        Wflow=np.zeros((nR,nW))
        for j in range(1,nR+1):
            for kk in range(1,nW+1):
                if kk==1:  
                    Wflow[j-1,kk-1]=Wfeas[j-1,kk-1]
                else:
                    Wflow[j-1,kk-1]=Wfeas[j-1,kk-1]-Wfeas[j-1,kk-2]
        return Wflow
                
    def plot_W_feas(self, nW, Akm2, Wfeas):
        col_kk=list()
        for kk in range(1,nW+1):
            tuplRGB=((kk-1)/(1.2*nW),(kk-1)/(1.2*nW),(kk-1)/(1.2*nW))     #R,G,B. R=G=B means gray              
            col_kk.append(tuplRGB)                
        
        fig,ax=plt.subplots()
        
        for kk in range(1,nW+1):
            ax.plot(Akm2,Wfeas[:,kk-1],color=col_kk[kk-1])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('drainage area (km2)')
        ax.set_ylabel('feasible" channel widths (m)')
        plt.savefig(f'{self.output_folder}/Widths_order'+ str(self.Omega) + '.png')      
        
    def fictNetwork(self,a,c):
        Rl=2 #Rl~2 (http://seismo.berkeley.edu/~kirchner/reprints/1993_18_Hortons_laws.pdf)
        l_mn=1
        Omega=self.Omega
        #--------------------------------------------------------------------
        #step 1: stream network
        T=np.zeros((Omega-1))
        for k in range(Omega-1):
            T[k]=a*c**k
#        print "T= ", T
        #Tokunaga matrix 'N', to review coherent number of reaches per w
        N=np.zeros((Omega,Omega))
        N[Omega-1,Omega-1]=1 #stream order i arriving to order j 
        for i in range(1,Omega):    #1:Omega-1
            w=Omega-i    #-1?
            for j in range(1,Omega-w+1):       #=1:Omega-w #always arrive to greater order         
                y=np.sum(N,1)    #sum rows, i.e. vertically
                N[w-1,w+j-1]=T[j-1]*y[w+j-1]
                if j==1:
                    N[w-1,w+j-1]=N[w-1,w+j-1]+2*y[w+j-1]    #add upstream branches                    
#        print "N= ", N
        #indicate # of internal nodes to preserve in each stream, given its order
        nnod=np.zeros((Omega))
        l=np.zeros((Omega))
        for w in range(1,Omega+1):
            nnod[w-1]=np.sum(T[:w-1])
            l[w-1]=l_mn*Rl**(w-1)
#        print "nnod= ", nnod         
#        print "l= ", l
        fig,ax=plt.subplots()
        pi=np.pi
        #plot branches
        X=[[0, 0, Omega, 3*pi/2, 0, 1],[0, l[Omega-1], Omega, 3*pi/2, 1, 1]]    #list: [x,y,w,dir_azimutRads,tail=1(ifNot:Head=0)orInternalFree=2orInternalOccupied=3,StreamID_assocToTail]
        X=np.array(X)
        ax.plot(X[:,0],X[:,1],'b')
        a=[1,-1]  #to alternate dir of each bifurcation
        ID=1
        deg45=np.pi/4
        for i in range(1,Omega):
            w=Omega-i
#            print "w= ", w
            #bifurc branches
            ind,=np.where((X[:,2]==w+1) & (X[:,4]==1.))  #tail points of headwater branches
#            print "ind= ", ind
            for kk in range(1,len(ind)+1):  #each headwater branch
                nptos=len(X[:,0])
                newZeros=np.zeros((2,len(X[0,:])))
                X=np.vstack([X,newZeros])    #create rows for new pairs of banches
                for k in [1,2]:   #2 branches
                    X[nptos+k-1,3]=X[ind[kk-1],3]+a[k-1]*deg45
                    X[nptos+k-1,2]=w
                    X[nptos+k-1,4]=1.
                    ID=ID+1
                    X[nptos+k-1,5]=ID
                    X[nptos+k-1,0]=X[ind[kk-1],0]-l[w-1]*np.cos(X[nptos+k-1,3])
                    X[nptos+k-1,1]=X[ind[kk-1],1]-l[w-1]*np.sin(X[nptos+k-1,3])            
                    ax.plot([X[ind[kk-1],0], X[nptos+k-1,0]],[X[ind[kk-1],1], X[nptos+k-1,1]],'b')
            #internal branches
#            nptos=len(X[:,0])
            np.set_printoptions(precision=1, suppress=True)
#            print "X= ", X
            ind1,=np.where((X[:,2]==w+1) & (X[:,4]==1.))  #tail points of headwater branches
#            print "ind1= ", ind1
            for kk in range(1,len(ind1)+1):  #each stream with order w+1
#                print "kk=", kk
                nptos=len(X[:,0])
#                print "nptos= ", nptos
                newZeros=np.zeros((int(nnod[w+1-1]),len(X[0,:])))
                X=np.vstack([X,newZeros])
                np.set_printoptions(precision=1, suppress=True)
#                print "with stackZeros for kkth stream order w+1, X= ", X
                for k in range(1,int(nnod[w+1-1])+1):
                    #set intermediate points for confluences in streams with order w+j            
#                    print "k=", k
                    r=(k/(nnod[w+1-1]+1))*l[w+1-1]
#                    print "for k=", k, " , r= ", r
                    X[nptos+k-1,0]=X[ind1[kk-1],0]+r*np.cos(X[ind1[kk-1],3])
                    X[nptos+k-1,1]=X[ind1[kk-1],1]+r*np.sin(X[ind1[kk-1],3])
                    X[nptos+k-1,2]=w+1
#                    print "X[",nptos+k-1,",2]= ", X[nptos+k-1,2] 
                    X[nptos+k-1,3]=X[ind1[kk-1],3]
                    X[nptos+k-1,4]=2 #report that this internal node is available
                    X[nptos+k-1,5]=X[ind1[kk-1],5]   #preserve streamID
            np.set_printoptions(precision=1, suppress=True)
#            print "internal point were set in X= ", X
            for j in range(1,Omega-w+1):  #each order difference                   
                IDj=X[X[:,2]==w+j,5]    #ID of streams w+j
                IDj=np.unique(IDj)
#                print "starting to build internal branches. We are in w= ", w
#                print "IDj= ", IDj
                for k in range(1,len(IDj)+1):  #each stream w+j                
                    nptos=len(X[:,0])
                    newZeros=np.zeros((int(T[j-1]),len(X[0,:])))
                    X=np.vstack([X,newZeros])
                    for jj in range(1,int(T[j-1])+1): 
                        ind2,=np.where((X[:,4]==2) & (X[:,5]==IDj[k-1]))                     
                        if j==1:
                            kk=int(round(len(ind2)/2.))    #internal nodes must be taken first by higher order affluents
                        elif j>1: 
                            if jj==1:
                                kk=1
                            else:
                                kk=len(ind2)                                                    
                        aIndex=int(round(np.random.rand()))
                        randSide=a[aIndex]    #choose randomly left or right side for confluence
                        X[nptos+jj-1,3]=X[ind2[kk-1],3]+randSide*np.pi*(1./2-1./8)  
                        X[nptos+jj-1,2]=w
                        X[nptos+jj-1,0]=X[ind2[kk-1],0]-l[w-1]*np.cos(X[nptos+jj-1,3])
                        X[nptos+jj-1,1]=X[ind2[kk-1],1]-l[w-1]*np.sin(X[nptos+jj-1,3])
                        X[nptos+jj-1,4]=1
                        ID=ID+1
                        X[nptos+jj-1,5]=ID
                        X[ind2[kk-1],4]=3 #report that this internal node is already occupied
                        ax.plot([X[ind2[kk-1],0],X[nptos+jj-1,0]],[X[ind2[kk-1],1],X[nptos+jj-1,1]],'b')
#                    print "for this stream, BRANCHES to internal points were set in X= ", X    
#        print "#reaches= ", len(X[:,0])-1                
#        plt.savefig('doct_fictBasin/FictTopol_Order'+ str(Omega) + '.png')    
        X=np.delete(X,0,axis=0)    #remove outlet. Topology is reaches in X defined by their upstream point.
        #--------------------------------------------------------------------
        #step2: topologic matrix of reaches and drainage areas
        ind2pi,=np.where(X[:,3]>2*3.14159)
        X[ind2pi,3]=X[ind2pi,3]-2*3.14159    #needed in following lines to rank topology for flow accum
        n=len(X)
        cheqID=np.zeros((n,2))    ##[ID, exist?]
        cheqID[:,0]=np.arange(1,n+1)
#        print "cheqID= ", cheqID
#        print "pre unique IDs, X= ", X
        contRepet=0
        for i in range(1,n+1): #get unique IDs
            if cheqID[int(X[i-1,5])-1,1]==1:    #available?
                ID=ID+1
                X[i-1,5]=ID    #new
                contRepet+=1
#            print "i= ", i    
#            print "contRepet= ", contRepet                         
            cheqID[int(X[i-1,5])-1,1]=cheqID[int(X[i-1,5])-1,1]+1 #reserve ID to be unique (cheq<=1)
#            print cheqID
#            print "---------------"
#        print "every node ID is different, X= ", X
        dist=np.zeros((n,n))
        X=np.concatenate((X,np.zeros((n,5))),axis=1)
        X[:,7]=X[:,7]+self.A_hs/1.e6    #every reach has at least its local drainage area
#        print "including A, X= ", X
        L=l[0]
        #first solve low order reaches to simplify procedure to get cumulative area
        X=X[X[:,2].argsort()]
#        print "sorted by w, X= ", X
        for i in range(1,n+1):
#            print "------------------------------------------"
#            print "X_row= ", i-1
            x_down=X[i-1,0]+L*np.cos(X[i-1,3])
            y_down=X[i-1,1]+L*np.sin(X[i-1,3])
            #search nearest point        
            for j in range(1,n+1):
                if i!=j:    #downstream neighbor must be a point different to the current one
                    dist[i-1,j-1]=abs(x_down-X[j-1,0])+abs(y_down-X[j-1,1])            
                else:
                    dist[i-1,j-1]=100 #some value such that minimum dist will be in other index surely
            goalDist=min(dist[i-1,:])        
            ind,=np.where(dist[i-1,:]==goalDist)
#            print "nearest ID to downstream calculated point= ", ind        
            if np.sqrt(x_down**2+y_down**2)>0.01*l_mn:    
                X[i-1,6]=X[ind,5]
            else:    #drain to x,y=origin
                X[i-1,6]=0
#        print "PRE double sort X= ", X
        w=1    #create nw_pts 
        nw_pts=np.zeros((Omega))
        for i in range(1,n+1):
#            print "creating nw, step i=", i
            if X[i-1,2]>w:
                w=w+1
            nw_pts[w-1]=nw_pts[w-1]+1
#        print "nw_pts= ", nw_pts
#        nw_pts_cum=np.zeros((Omega))
#        ok_drain=0    #start to account reaches arranged in series (e.g. several w=1 arriving to same w=2 while flowing downstream)
#        for w in range(1,Omega+1):
#            index=range(0,w)
#            nw_pts_cum[w-1]=nw_pts[index].sum()
#            print "nw_pts_cum= ", nw_pts_cum
        ind=np.lexsort((X[:,3],X[:,2]))    #1st sort by w at first, 2nd by flowDir
#        print "index list to perform double sort= ", ind    
        temp=np.zeros((n,len(X[0,:])))
        for i in range(1,n+1):
            temp[i-1,:]=X[ind[i-1],:]
        X=temp
#        print "POS double sort X= ", X
        i=int(nw_pts[0])+1    #only streams w>2 have chance to be re-ranked according to flow direction
        totPts=len(X[:,0])
#        print "totPts= ", totPts
        ip=0        
        while i<totPts:    #whole X
#            print "i_iniStream= ", i-1
            cont=0
            iniDir=i
            yN=0
            xE=0
            if X[i-1,3]<3.1416:    #drain to north
                yN=1
            if (X[i-1,3]<3.1416*1./2) or (X[i-1,3]>3.1416*3./2):
                xE=1    
            while (X[i,2]==X[i-1,2]) & (X[i,3]==X[i-1,3]):    #same stream preserve w and direction
                i+=1
                if i==totPts: 
#                    print "in break due to max i"
                    break
            S=X[iniDir-1:i,:]
#            print "S= ", S
            S=S[S[:,1].argsort()]    #default sort by increasing north
            if yN==0:
                S=S[::-1]    #reverse rank
#            print "i= ", i
            if i==totPts:
                ip=i-1
            if (X[ip,1]==X[ip-1,1]) & (X[ip,2]==X[ip-1,2]):    #flow direction, in w, with no northing component
                S=S[S[:,0].argsort()]    ##default sort by increasing east
                if xE==0:
                    S=S[::-1]    #reverse rank                                            
#            print "flow dir OK. S= ", S
            X[iniDir-1:i,:]=S    #reset points ranked according to flow dir
            i=i+1
#        print "sorted by flow dir, X= ", X
        for i in range(1,n+1):        
#            print "-----------------------------------------------------------"
            ind,=np.where(X[:,5]==X[i-1,6])    #unique IDs, as stated lines above, imply unique index
#            print "reach downstr, with ID ", X[ind,5], ", has previous Aacum= ", X[ind,7]
#            print "A_bringing [row= ", i-1, "]=",  X[i-1,7]
            X[ind,7]=X[ind,7]+X[i-1,7]    #add drainage area of new tributary        
        X[:,8]=.04*X[:,7]**-.3    #local bedrock slope in reaches (default coef 08 and exp -.3 from Colombian regression in ideam2017)
        X[n-1,8]=1e-10  #downstream base level fixes alluvium thickness (June8, 2018)
        X[:,9]=X[:,9]+np.sqrt(self.A_hs/1.e6)*1e3    #lenght of reaches [m] (then x y coor)
#        print "full X= ", X
        X=X[X[:,7].argsort()]    #sort by drainage area
        #before trim X, come back to IDs defined by rows. It is easier to deal with this in fluvial model
#        print "sorted by A, X= ", X
        indVect=np.zeros((totPts))
        for i in range(1,n+1-1):     #-1 because outlet tail point do not have downstream node
#            print "restating IDs, i=", i
            X[i-1,5]=0
            ind,=np.where(X[:,5]==X[i-1,6])
            indVect[i-1]=ind[0]
#            print "ind=", ind
        for i in range(1,n+1):
            X[:,6]=indVect    #defined in Python format i.e. arrays starting in index 0 
#        print "IDs as rows, X= ", X
        X=np.vstack([X[:,0],X[:,1],X[:,2],X[:,6],X[:,7],X[:,8],X[:,9]])
        X=X.T
#        X=[X[:,2], X[:,6], X[:,7], X[:,8], X[:,9]]    #useful columns for morphodynamics: w, downstrID, A(km2), S and L(km)  
        #upstream altitude of bedrock in reaches
        X=np.concatenate((X,np.zeros((totPts,1))),axis=1)    #column to store bedrock elevation [m]
        X[totPts-1,7]=0
        for jj in range(1,totPts+1-1):    #avoid to find drainage dir to "0", which is outlet          
            j=totPts-jj+1    #cumulate elevation in upstream direction
            ind,=np.where(X[:,3]==j-1)    #find flow direction
#            print "ind=", ind            
            nn=len(ind)    #1 or 2, asuming no trifurcation upstream
            for k in range(1,nn+1):
                X[ind[k-1],7]=X[j-1,7]+X[ind[k-1],5]*X[ind[k-1],6]
#            print "jj=", jj
#            print "with zs, X= ", X                


#        fig,ax=plt.subplots()
        maxZ=max(X[:,7])
        color = [str(item/maxZ) for item in X[:,7]]
#        print "color=", color
        plt.scatter(X[:,0],X[:,1],c=color)
        plt.colorbar()
        plt.savefig('doct_fictBasin/TopolZcolor_order'+ str(self.Omega) + '.png')
        plt.clf()

#        print "elevations added in last column, X= ", X
        for j in range(1,totPts+1):
            if X[j-1,3]>0:
                X[j-1,3]=X[j-1,3]+1                
#        print "downstream ID + 1 for Fortran compatibility, X= ", X
        return X                              
        
    def simul_plots(self,t_ref_ini_raw=1, t_ref_fin_raw=0, t_scale_log=0, coarse_dt_flag=1):
        #goal: plot binary files of results
        #[t_ref]=days if (coarse_dt_flag==1) else [t_ref]=Tsteps
         
        flume_else_basin=self.flume_else_basin
        X=self.X
        coarse_dt=self.coarse_dt
        trenchFlag = self.trenchFlag
        
        event=1 
        if t_ref_fin_raw==0:  #default argument means user wants to plot the whole simulation time
            if coarse_dt_flag==1:
                t_ref_fin_raw=self.n_coarse_dt
            else:
                t_ref_fin_raw=self.nt
            event=0 
        n=self.nt        
        
        if coarse_dt_flag==1:         
            ii1,=np.where(coarse_dt[:,0]==t_ref_ini_raw)   #find step i related to the day indicated as argument
            ii2,=np.where(coarse_dt[:,0]==t_ref_fin_raw)
            t_ref_ini=int(coarse_dt[ii1,1]) if t_ref_ini_raw>1 else 1  #include fine initial process if plotting whole simuTime
            t_ref_fin=int(coarse_dt[ii2,1])
            print('t_ref_ini= ', t_ref_ini)
            print('t_ref_fin= ', t_ref_fin)
        else:
            t_ref_ini=t_ref_ini_raw
            t_ref_fin=t_ref_fin_raw
            
            
        #Prepare evenly spaced* time array to plot longitudinal profiles:
        #*in days if raw arguments are so and flag is on    

        t_prof_interv=(t_ref_fin_raw-t_ref_ini_raw)/10            
        t_profile=np.arange(t_ref_ini_raw,t_ref_fin_raw+t_prof_interv,t_prof_interv)    
        #t[units] according to time units of raw args
        
        t_profile=t_profile.astype(int)
        nt_profile=len(t_profile)
        print('raw t_profile= ', t_profile)
        if coarse_dt_flag==1:  #translate t_profile from [days] to [Tsteps]    
            t_profile_days=t_profile
            t_profile_ok=np.zeros((nt_profile))
            for i in range(1,nt_profile+1): 
                ii,=np.where(coarse_dt[:,0]==t_profile[i-1])
                t_profile_ok[i-1]=coarse_dt[ii,1]
            t_profile=t_profile_ok
            t_profile=t_profile.astype(int)    
            print('saved[days] t_profile= ', t_profile_days)        
        print('[Tsteps] t_profile= ', t_profile)        
        
#        print(str('days' if flag_coarse_dt_flag==1 else 'Tsteps')+' for t_profile= ', t_profile)             
            
            
            
        Wfeas, Wflow, Dcentr, Dpct, pctD = self.assumed_vars(self.nW, self.Dmean, self.Dstdg,\
        self.nD, self.stochGSD, self.trenchFlag)
#        Qhill=np.loadtxt('doct_fictBasin/Qhill.txt')
        Dcentr = self.amplif_D*Dcentr
        Dpct = self.amplif_D*Dpct 
        nW=len(Wfeas[0,:])
        nR=len(Wfeas[:,0])
        nD=len(Dcentr)
        jref=[2, int(.3*nR), int(.7*nR), int(1.*nR)-1]   #[1,3,5,8]
        njref=len(jref)
        
        dtf=self.dt_fine*60  #from hours to minutes
        t=[]
        t.append(0)
        for i in range(2,n+1):
            t.append(t[i-2]+dtf[i-2]) #[t]=minutes. t is estimated at the end of time step
        t=np.array(t)
        #        t=np.arange(1,self.n+1)
        if flume_else_basin==1:
            t=t/60. #[hours]
        else:
            t=t/1440.   #[days]                           
        indW_max=np.zeros((nR))
        for j in range(1,nR+1):
            kk=1
            while (kk<=nW and Wflow[j-1,kk-1]>0):         
                kk=kk+1
                indW_max[j-1]=indW_max[j-1]+1                            
        path_list=self.path_list
#        Q,res=mod_model_sc.read_FloatArrIJ(path_Q,1,self.n,nR)  #(var_bin,posFecha+1,N,5)
        Q=np.zeros((n,nR)); h=np.zeros((n,nR,nW)); v=np.zeros((n,nR,nW))
        z=np.zeros((n,nR,nW)); Da84=np.zeros((n,nR,nW)); Da50=np.zeros((n,nR,nW)) 
        zw=np.zeros((n,nR)); So=np.zeros((n,nR)); W=np.zeros((n,nR))
        indW=np.zeros((n,nR)); th50a=np.zeros((n,nR)); thm84B=np.zeros((n,nR));
        Dss50=np.zeros((n,nR,nW)); thratio84=np.zeros((n,nR)) 
        ck=np.zeros((n,nR,nD)); Qhill=np.zeros((n,nR)); slid_sh=np.zeros((n,nR,nD))
        fileunit=10
        for i in range(t_ref_ini, t_ref_fin+1):
            Q[i-1,:],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[1-1]-1],i,nR)  #(var_bin,posFecha+1,N,5)
            So[i-1,:],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[8-1]-1],i,nR)
            indW[i-1,:],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[9-1]-1],i,nR)
            th50a[i-1,:],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[10-1]-1],i,nR)
            thm84B[i-1,:],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[11-1]-1],i,nR)
            thratio84[i-1,:],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[12-1]-1],i,nR)            
            Qhill[i-1,:],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[14-1]-1],i,nR)            
            for kk in range(1,nW+1):
                h[i-1,:,kk-1],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[2-1]-1+kk-1],i,nR)            
                v[i-1,:,kk-1],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[3-1]-1+kk-1],i,nR)
#                v[i-1,:,1],res=rmr.read_float_arr(fileunit,path_v_kk2,i,nR)
                z[i-1,:,kk-1],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[4-1]-1+kk-1],i,nR)
                Da84[i-1,:,kk-1],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[5-1]-1+kk-1],i,nR)                
                Da50[i-1,:,kk-1],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[6-1]-1+kk-1],i,nR)
                Dss50[i-1,:,kk-1],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[7-1]-1+kk-1],i,nR)
#                z[i-1,:,1],res=rmr.read_float_arr(fileunit,path_z_kk2,i,nR)
#                z[i-1,:,2],res=rmr.read_float_arr(fileunit,path_z_kk3,i,nR)
#                z[i-1,:,3],res=rmr.read_float_arr(fileunit,path_z_kk4,i,nR)
            for k in range(1,nD+1):
                ck[i-1,:,k-1],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[13-1]-1+k-1],i,nR)
            for k in range(1,self.nD_hill+1):
                slid_sh[i-1,:,k-1],res=rmr.read_float_arr(fileunit,path_list[self.pos_ini_var[15-1]-1+k-1],i,nR)

            for j in range(1,nR+1):        
                W[i-1,j-1]=Wfeas[j-1,int(indW[i-1,j-1]-1)]
        Da84=1000*Da84; Da50=1000*Da50; Dss50=1000*Dss50                
        zw=z[:,:,0]+h[:,:,0]
        np.set_printoptions(precision=4, suppress=False)    
        print ("finally, in py, after read binary, Q=", Q)        
            
        t=t[t_ref_ini-1:t_ref_fin]  
        
        ck_sum=np.zeros((n,nR))
        for j in range(1,nR+1):
            ck_sum[t_ref_ini-1:t_ref_fin,j-1]=np.sum(ck[t_ref_ini-1:t_ref_fin,j-1,:], axis=1)

        #plot waves (f(t), curve per x) of water and sediment        
        fig,ax=plt.subplots(nrows=2,ncols=1)
        yw=np.zeros((n,nR))                
        ys=np.zeros((n,nR))                
        col_x=[]
        if flume_else_basin==0:
            Amx=self.X[nR-1,4]        
            for j,A in enumerate(self.X[:,4]):
                tuplRGB=(A/(1.2*Amx), A/(1.2*Amx), A/(1.2*Amx))    #R,G,B. R=G=B means gray. Degrade color by drainA
                col_x.append(tuplRGB)
        else: 
            for j in range(1,nR+1):
                tuplRGB=(j/(1.2*(nR+1)), j/(1.2*(nR+1)), j/(1.2*(nR+1))) #Degrade color downstream
                col_x.append(tuplRGB)
        convFact=2600*3600                
        if flume_else_basin==0:
            for j in range(1,nR+1):
                yw[t_ref_ini-1:t_ref_fin,j-1]= 3.6*Q[t_ref_ini-1:t_ref_fin,j-1]/self.X[j-1,4]
                ax[0].plot(t,yw[t_ref_ini-1:t_ref_fin,j-1],color=col_x[j-1],\
                marker='o',markersize=.5,ls='None',label='j'+str(j-1))            
                ys[t_ref_ini-1:t_ref_fin,j-1]=ck_sum[t_ref_ini-1:t_ref_fin,j-1]*Q[t_ref_ini-1:t_ref_fin,j-1]/\
                self.X[j-1,4]*2.65*31.536e6            
                ax[1].plot(t,ys[t_ref_ini-1:t_ref_fin,j-1],color=col_x[j-1],\
                marker='o',markersize=.5,ls='None',label='j'+str(j-1))
        else:
            for j in range(1,nR+1):            
                yw[t_ref_ini-1:t_ref_fin,j-1]= Q[t_ref_ini-1:t_ref_fin,j-1]*1e3 #l/s
                ax[0].plot(t,yw[t_ref_ini-1:t_ref_fin,j-1],color=col_x[j-1],\
                marker='o',markersize=.5,ls='None',label='j'+str(j-1))            
                ys[t_ref_ini-1:t_ref_fin,j-1]=ck_sum[t_ref_ini-1:t_ref_fin,j-1]*Q[t_ref_ini-1:t_ref_fin,j-1]*convFact
                ax[1].plot(t,ys[t_ref_ini-1:t_ref_fin,j-1],color=col_x[j-1],\
                marker='o',markersize=.5,ls='None',label='j'+str(j-1))        
            
#            ax.plot(t,Q[:,12-1],'r',label='j12')
#            ax.plot(t,Q[:,18-1],'g',label='j18')
#            ax.plot(t,Q[:,21-1],'k',label='j21')
#        leg = ax[0].legend(loc=1,fancybox=True, fontsize=6)
#        leg.get_frame().set_alpha(0.3)
#        leg = ax[1].legend(loc=1,fancybox=True, fontsize=6)
#        leg.get_frame().set_alpha(0.3) 
        
        if t_scale_log==1:
            ax[0].set_xscale('log')        
        if flume_else_basin==0:
            t_equil=int(t_ref_fin/5.) #initial 10% of time is warming period of system
            ax[0].set_title('Lighter color with A. yw_equil='+\
            str(round(np.mean(yw[t_equil:t_ref_fin,nR-1])*24*365, 0))+'mm/yr. ys_eq='+\
            str(round(np.mean(ys[t_equil:t_ref_fin,nR-2]), 0))+'ton/km2/yr', fontsize=10)        
            ax[0].set_ylabel('water yield (mm/h)',fontsize=10)
            ax[0].set_ylim(ymin=1e-2)
            ax[1].set_ylabel('sedim yield (ton/km2/yr)',fontsize=10)
            if t_scale_log==1:  ax[1].set_xscale('log')
            ax[1].set_ylim(ymin=1e0)                
            ax[0].set_yscale('log')
            ax[1].set_yscale('log')
        else:            
            t_equil=1
            ax[0].set_title('Lighter color downstream. Qmean='+\
            str(round(np.mean(yw[t_equil:t_ref_fin,nR-1]), 0))+'l/s. Qs_mean='+\
            str(round(np.mean(ys[t_equil:t_ref_fin,nR-2]), 0))+'kg/h', fontsize=10)        
            ax[0].set_ylabel('Q (l/s)',fontsize=10)
#            ax[0].set_ylim(ymin=1e-2)
            ax[1].set_ylabel('Qs (kg/h)',fontsize=10)
            if t_scale_log==1:  ax[1].set_xscale('log')
            ax[1].set_ylim(ymin=1e-1)
#            ax[0].set_yscale('log')
            ax[1].set_yscale('log')            
#        ax.set_ylim(0,120)
        if t_scale_log==1:  ax[1].set_xscale('log')
        plt.savefig(f'{self.output_folder}/Q_Qs_xt.png')
        plt.clf()
        
        #plot downstream incision:
        col_kk=list()
        for kk in range(1,nW+1):
            tuplRGB=((kk-1)/(1.2*nW),(kk-1)/(1.2*nW),(kk-1)/(1.2*nW))                   
            col_kk.append(tuplRGB)        
        if coarse_dt_flag==1:            
            jpath=[]
            j=1
            jpath.append(j)        
            abscissa=[0,X[j-1,6]]   #1st calc node of path has x=0. Append position of next node adding local length
            j=int(X[j-1,3])
            jpath.append(j)
            while j>0:  #start in reach j1 finding downstream node (while existing) and appending cumulative distance
    #            print('j= ', j) 
                abscissa.append(abscissa[-1]+X[j-1,6])            
                j=int(X[j-1,3])  #j_downstream
                jpath.append(j)
    #            print('new abscissa= ', abscissa[-1])
    #            print('jnext= ', j)
            jpath=jpath[:-1]
            print('abscissa= ', abscissa) 
            print('jpath= ', jpath)
            nR_profile=len(jpath)         
            
            zw_abs=np.zeros((nt_profile,nR_profile))
            zs_abs=np.zeros((nt_profile,nR_profile,nW))
            Da84_downstr=np.zeros((nt_profile,nR_profile,nW))
            Da50_downstr=np.zeros((nt_profile,nR_profile,nW))
            Dss50_downstr=np.zeros((nt_profile,nR_profile,nW))
            
            for ii,iprofi in enumerate(t_profile): #absolute elev: add local bedrock elevation to w and s thickness
                for jj,jprofi in enumerate(jpath):  #jj produced by 'enumerate' starts at '0', then ii index does not need: -1
                    zw_abs[ii,jj]= X[jprofi-1,7] + zw[iprofi-1,jprofi-1]
                    zs_abs[ii,jj,:]= X[jprofi-1,7] + z[iprofi-1,jprofi-1,:]
                    Da84_downstr[ii,jj,:]= Da84[iprofi-1,jprofi-1,:]
                    Da50_downstr[ii,jj,:]= Da50[iprofi-1,jprofi-1,:]
                    Dss50_downstr[ii,jj,:]= Dss50[iprofi-1,jprofi-1,:] 
                        
            fig,ax=plt.subplots(nrows=2,ncols=1)
            abscissa=abscissa[:-1]           
            for ii in range(1,nt_profile+1):    
                ax[0].clear()
                ax[1].clear()
                ax[0].set_title('after ' + str(round(t_profile_days[ii-1]/365,1)) + 'years')
                ax[0].plot(abscissa,zw_abs[ii-1,:],color='b')        #,label='zw_abs'
                ax[0].plot(abscissa,[X[jprofi-1,7] for jj,jprofi in enumerate(jpath)],':k',marker='o',markersize=3)
#                ax[0].set_ylim(ymax=600)
                for kk in range(1,int(indW_max[jprofi-1])+1):
                    ax[0].plot(abscissa,zs_abs[ii-1,:,kk-1],color=col_kk[kk-1])
                    ax[1].plot(abscissa,Da84_downstr[ii-1,:,kk-1],color=col_kk[kk-1],ls='-')
                    ax[1].plot(abscissa,Da50_downstr[ii-1,:,kk-1],color=col_kk[kk-1],ls=':')
                    ax[1].plot(abscissa,Dss50_downstr[ii-1,:,kk-1],color=col_kk[kk-1],ls='--')            
                ax[0].set_ylim(ymin=0, ymax=.6)
                ax[1].set_yscale('log')
                max_file_number=int(np.log10(self.n_coarse_dt))
                fill_digits=max_file_number-int(np.log10(t_profile_days[ii-1]))
                fill=str(0)*int(fill_digits)
                print('max_file_number= ', max_file_number)
                print('fill_digits= ', fill_digits)
                print('fill= ', fill)                        
                
                    
                #to leave images files ready to create ordered gif from Bash shell
                #https://stackoverflow.com/questions/5417979/batch-rename-sequential-files-by-padding-with-zeroes
                plt.savefig(f'{self.output_folder}/incis_downstr_day'+fill+str(t_profile_days[ii-1])+'.png')

        #RUN FROM BASH IN IMAGES FOLDER: 
#        cd folder
#        convert -delay 100 incis_downstr_day*.png incis_downstr.gif       
#           bah: BUG: does not restart cycle with right time sequence
        
        
            
        #plot event downscaled WATER supply
        
        fig,ax=plt.subplots(nrows=2, ncols=1)
        j_suppl=[0,5]
        for j in range(1,2+1):
            ax[j-1].plot(t,Qhill[t_ref_ini-1:t_ref_fin,j_suppl[j-1]-1],'ok',markersize=.5)   #,'.')  
            #python reads array slice with interval like this: [ , )                        
            ax[j-1].set_ylabel('Qhill_'+str(j_suppl[j-1])+'[m3/s]',fontsize=10)
            if flume_else_basin==0:
                ax[j-1].set_yscale('log')
            if t_scale_log==1:
                ax[j-1].set_xscale('log') 
        print('in py to plot, j=0, Qhill= ', Qhill[:,0])      
        plt.savefig(f'{self.output_folder}/Qhill_eventDownscaled.png')
        plt.clf()  
        
        
        
        #plot event downscaled SEDIMENT supply
        
        fig,ax=plt.subplots(nrows=2,ncols=1)
        slid_sh=slid_sh*2.65e6  #[g/s]
        for j in range(1,2+1):
            for k in range(1,nD+1):    #plot series of each grain size, only in initial reach, by now(29may2018), looking at flume 
                col=((k-1)/nD,(k-1)/nD,(k-1)/nD)    #R,G,B. R=G=B means gray
                ax[j-1].plot(t,slid_sh[t_ref_ini-1:t_ref_fin,j-1,k-1],color=col,label=str(round(Dcentr[k-1]*1e3,2))+'mm')
            ax[j-1].plot(t,np.sum(slid_sh[t_ref_ini-1:t_ref_fin,j-1,:], axis=1),'--k',label= 'total') 
            leg = ax[j-1].legend(loc=1,fancybox=True, fontsize=6)
            leg.get_frame().set_alpha(0.3)                             
            ax[j-1].set_ylabel('feed_j'+str(j)+'[g/s]',fontsize=10)
            ax[j-1].set_yscale('log')
            if t_scale_log==1:
                ax[j-1].set_xscale('log')
        plt.savefig(f'{self.output_folder}/feed_eventDownscaled.png')
        plt.clf()                 
                         
#        fig,ax=plt.subplots(nrows=2,ncols=1)
#        ax[0].plot(t,h[t_ref_ini-1:t_ref_fin,jref[0]-1,0],'b',t,h[t_ref_ini-1:t_ref_fin,jref[1]-1,0],'r',t,\
#        h[t_ref_ini-1:t_ref_fin,jref[2]-1,0],'g',t,h[t_ref_ini-1:t_ref_fin,jref[3]-1,0],'k')
#        
#        ax[0].plot(t,v[t_ref_ini-1:t_ref_fin,jref[0]-1,0],'--b',t,v[t_ref_ini-1:t_ref_fin,jref[1]-1,0],'--r',t,\
#        v[t_ref_ini-1:t_ref_fin,jref[2]-1,0],'--g',t,v[t_ref_ini-1:t_ref_fin,jref[3]-1,0],'--k')
#        
#        if t_scale_log==1:
#            ax[0].set_xscale('log')    
#        if flume_else_basin==0:
#            ax[0].set_yscale('log')                            
#        ax[0].set_ylabel('[-]h(m), [--]v(m/s); in kk1',fontsize=10)
##        if self.flume_else_basin==1:       
##            h_exp=[.073, .08, .083, .072, .075, .073]    #value at the end of each run (ElguetaHassan2017)
##        ax[0].set_yscale('log')
#            
##        if nW>1:             
##            ax[1].plot(t,h[:,jref[0]-1,1],'b',t,h[:,jref[1]-1,1],'r',t,h[:,jref[2]-1,1],'g',t,h[:,jref[3]-1,1],'k')
##            ax[1].plot(t,v[:,jref[0]-1,1],'--b',t,v[:,jref[1]-1,1],'--r',t,v[:,jref[2]-1,1],'--g',t,v[:,jref[3]-1,1],'--k')    
##            ax[1].set_ylabel('[-]h(m), [--]v(m/s); in kk2',fontsize=10)
##    #        ax[1].set_yscale('log')
##            ax[1].set_ylim(ymin=1e-2)
##            if t_scale_log==1:
##                ax[1].set_xscale('log')            
#        plt.savefig(f'{self.output_folder}/v_h.png')
#        plt.clf()
        
        fig,ax=plt.subplots(nrows=2,ncols=1)
#        ax[0].plot(t,thm84a[:,jref[0]-1],'b',t,thm84a[:,jref[1]-1],'r',t,thm84a[:,jref[2]-1],'g',t,thm84a[:,jref[3]-1],'k')
#        ax[0].plot(t,thm84B[:,jref[0]-1],'--b',t,thm84B[:,jref[1]-1],'--r',t,thm84B[:,jref[2]-1],'--g',t,thm84B[:,jref[3]-1],'--k')
#        ax[0].set_ylabel('thm84, [-]act, [--]Budget', fontsize=10)
        if t_scale_log==1:
            ax[0].set_xscale('log')        
        ax[0].grid(which='major',axis='y',linewidth=.75)
        ax[0].grid(which='minor',axis='y',linewidth=.25)        
#        ax[0].set_yscale('log')
        for jj in range(1,njref+1):
            col=((jj-1)/njref,(jj-1)/njref,(jj-1)/njref)    #R,G,B. R=G=B means gray
            ax[0].plot(t,th50a[t_ref_ini-1:t_ref_fin,jref[jj-1]-1],color=col,label='j'+str(jref[jj-1]))
            ax[1].plot(t,thratio84[t_ref_ini-1:t_ref_fin,jref[jj-1]-1],color=col,label='j'+str(jref[jj-1]))        
        ax[0].set_ylabel('th50a', fontsize=10)            
        if flume_else_basin==0:
            ax[0].set_yscale('log')        
#        ax[1].plot(t,thratio84[:,12-1],'r',label='j12')        
#        ax[1].plot(t,thratio84[:,18-1],'g',label='j18')        
#        ax[1].plot(t,thratio84[:,21-1],'k',label='j21')
        ax[1].grid(which='major',axis='y',linewidth=.75)
        ax[1].grid(which='minor',axis='y',linewidth=.25)                
        leg = ax[1].legend(loc=1,fancybox=True, fontsize=6)
        leg.get_frame().set_alpha(0.3)                                                        
        ax[1].set_ylabel('thratio84', fontsize=10)
        if t_scale_log==1:
            ax[1].set_xscale('log')        
        ax[1].set_yscale('log')                
#        ax[1].set_ylim(ymin=1e-2)
        plt.savefig(f'{self.output_folder}/mobility.png')
        plt.clf()
        
        #plot phase states. W vs So
        

        

#        #plot bedload and hydraulics

        fig,ax=plt.subplots(nrows=2,ncols=1)
        nD=len(Dcentr)
        
#        col_W=[]        
#        for jj in range(1,nW+1):
#            tuplRGB=((jj-1)/nD,(jj-1)/nD,(jj-1)/nD)    #R,G,B. R=G=B means gray
#            col_W.append(tuplRGB)
        for j in range(1,nR+1):
            ax[0].clear()
            for jj in range(1,nD+1):    #plot series of each grain size
                col=((jj-1)/nD,(jj-1)/nD,(jj-1)/nD)    #R,G,B. R=G=B means gray
                ax[0].plot(t,ck[t_ref_ini-1:t_ref_fin,j-1,jj-1],color=col,label= str(int(Dcentr[jj-1]*1e3)) + 'mm')
            ax[0].plot(t,ck_sum[t_ref_ini-1:t_ref_fin,j-1],'--k',label= 'total')    
            leg = ax[0].legend(loc=1,fancybox=True, fontsize=6)
            leg.get_frame().set_alpha(0.3)                                                                
            ax[0].set_yscale('log')                
            if t_scale_log==1:  ax[0].set_xscale('log')            
            ax[0].set_ylim(ymin=1e-7, ymax=1e0)   #if j<nR:   #june10, 2018: to review basin outlet to make sediment balance
            ax[0].grid(which='major',axis='y',linewidth=.75)
            ax[0].grid(which='minor',axis='y',linewidth=.25)
            ax[0].set_ylabel('fractional cb [m3/m3]',fontsize=10)
            ax[0].set_title('reach j= '+str(j))

            ax[1].clear()
            
            #Froude only in main channel:
            ax[1].plot(t,v[t_ref_ini-1:t_ref_fin,j-1,0]/\
            np.sqrt(10*h[t_ref_ini-1:t_ref_fin,j-1,0]),'k')
            
            ii=0                        
            if event==1: #i loop to find flood out of main channel only if event, for time eficciency
                for i in range(t_ref_ini,t_ref_fin+1):
                    ii=ii+1
                    for jj in range(1,int(indW[i-1,j-1])+1):    #plot series of hydraulics in wet crosssection
                        ax[1].plot(t[ii-1],h[i-1,j-1,jj-1],color=col_kk[jj-1],ls='None',marker='+',markersize=1)
                        ax[1].plot(t[ii-1],v[i-1,j-1,jj-1],color=col_kk[jj-1],ls='None',marker='_',markersize=1)
                ax[1].set_ylabel('Fr|(r)h|(g)v|',fontsize=10)        
            else: #plot only hydraulics in main channel
                ax[1].plot(t,h[t_ref_ini-1:t_ref_fin,j-1,0],'r',ls='None',marker='+',markersize=1)
                ax[1].plot(t,v[t_ref_ini-1:t_ref_fin,j-1,0],'g',ls='None',marker='_',markersize=1)
                ax[1].set_ylabel('Fr|+h|_v|, darker deeper',fontsize=10)    
                
#            leg = ax[1].legend(loc=1,fancybox=True, fontsize=6)
#            leg.get_frame().set_alpha(0.3)                
            ax[1].grid(which='major',axis='y',linewidth=.75)
            ax[1].grid(which='minor',axis='y',linewidth=.25)
            
            ax[1].set_yscale('log')
            if t_scale_log==1: ax[1].set_xscale('log')            
            plt.savefig(f'{self.output_folder}/ck_hydrau_j'+str(j)+'.png')        
                                    
        #if trench experiment, plot spatially averaged hydraulic geometry and wet (channel) & dry (banks) store balance:
        if (flume_else_basin==1 and trenchFlag==1):    #  and t_ref_fin==0 #t_ref_fin==0 means you are not plotting certain time window but the whole simulation time.
            print('va a plotear lump vars of exper trench')
            Qtrf = self.Q2reach*1e3
            d_xmean, W_xmean, v_xmean, t_last_meas_rising, ind_t_list = self.xmean_hg(t, self.Q2reach,\
            h, W, indW, Wflow, v, nR)
            
            self.plot_hg(Qtrf, d_xmean, W_xmean, v_xmean, self.output_folder)
            ck_out = ck[:, nR-1, :]
            c_out_allk = np.sum(ck_out, axis = 1)
            n_tlmr = len(t_last_meas_rising)
            
            z_abs2 = [z[:,j-1,0] + X[j-1,7] for j in jpath]
            z_abs2 = np.array(z_abs2)
            z_abs2 = z_abs2.T
            print('dims z_abs2= ', z_abs2.shape)
            So_agr, R2_zx = self.gradients_during_exp(self.dx_flume, nR, z_abs2, ind_t_list, n_tlmr)
            self.plot_gradients(Qtrf, So_agr, R2_zx, n_tlmr, self.output_folder, 'elev')
            grad, R2_cx = self.gradients_during_exp(self.dx_flume, nR, np.sum(ck, axis = 2), ind_t_list, n_tlmr)
            self.plot_gradients(Qtrf, grad, R2_cx, n_tlmr, self.output_folder, 'concentr')
            
            pi, fi = self.fract_tpt_get(pctD, ck_out, c_out_allk, nD, ind_t_list, n_tlmr)
            self.plot_fract_tpt(pi, fi, Dcentr, Qtrf, self.output_folder, n_tlmr)
            Q_out = Q[:, nR-1]
            t_meas_intraQ, Qs_mean, c_mean, n_tm = self.c4trapCalib_get(t, c_out_allk, Q_out)
            self.plot_c4trapCalib(t_meas_intraQ, Qs_mean, Qtrf, c_mean, n_tm, self.output_folder)            
            store, cum_out = self.balance_flume_trench(t, z, Q_out, c_out_allk, indW, Wflow, nR, nW,\
            t_last_meas_rising, self.dx_flume, dtf[0])  #trench experiment has constant dx and dt
            
            #about last line: 2 rows; for channel and bank respectively
            self.plot_stores(Qtrf, store, cum_out, self.output_folder)
            
            sing_store = self.singl_balance_flume_trench(t, z, Wflow, nR, t_last_meas_rising, \
            self.dx_flume)  #trench experiment has constant dx and dt            
            self.plot_single_store(Qtrf, sing_store, cum_out, self.output_folder)
        else:
            print('eh, no plotea lump vars of exper trench')            
            
        
        #plot frequencies of water and sediment yield
#        ii=[]
#        r=np.zeros()
#        for i,p in enumerate(self.Pt):
#            ii.append(i)
#        ii=np.array(ii)
#        ii=ii/self.n     #array between 0 and 1
#        ax[1].plot(np.sort(self.Pt)[::-1],ii,'k')
#        ax[1].set_xlabel('P[mm/day]')
#        ax[1].set_ylabel('prob_exced')
#        ax[1].set_yscale('log')
        
        
        
        #plot hydraulic geometry (1st for all j to debug, then only for jref to show spatial scales)
        
                                               
#        #plot validation (may29, 2018: Maria's data) 
#        if self.flume_else_basin==1:       
#            t_Qs=np.arange(.25,40*7+.25,.25)   #[hours], each quarter
#            Qs_real=np.loadtxt(f'{self.output_folder}/data/Qbtotal_15min.txt')
#            t_bed_sizes=np.loadtxt(f'{self.output_folder}/data/time4ph.txt')
#            D90_bed_real=np.loadtxt(f'{self.output_folder}/data/FD90_4ph.txt')
#            Dg_bed_real=np.loadtxt(f'{self.output_folder}/data/FDg_4ph.txt')
#            t_slope=np.loadtxt(f'{self.output_folder}/data/time_slope.txt')        
#            So_real=np.loadtxt(f'{self.output_folder}/data/slope.txt')
#            fig,ax=plt.subplots(nrows=4,ncols=1)
#            ax[0].plot(t_Qs,Qs_real,'--k',label= 'exp_flume')
#            ax[0].plot(t,np.sum(ck[:,nR-1,:], axis=1)*Q[:,nR-1]*2.65e6,label= 'num_model')            
#            ax[0].set_yscale('log')
#            ax[0].set_ylim(ymin=5e-2)                               
#            ax[0].set_ylabel('Qs (g/s)',fontsize=8)  
#            leg = ax[0].legend(loc=2,fancybox=True, fontsize=6)
#            leg.get_frame().set_alpha(0.3)             
#            ax[1].plot(t,np.mean(Da50[:,4-1:6-1,0],axis=1))
#            ax[1].plot(t_bed_sizes,Dg_bed_real,'--k')
#            ax[2].plot(t,np.mean(Da84[:,4-1:6-1,0],axis=1))        
#            ax[2].plot(t_bed_sizes,D90_bed_real,'--k')                
#            ax[3].plot(t,np.mean(So[:,1:nR-2],axis=1))        #excluding boundary nodes as suggested by Tobias (May30/18)
#            ax[3].plot(t_slope,So_real,'--k')        
#            ax[3].set_ylabel('slope (m/m)',fontsize=8)  
#            plt.savefig(f'{self.output_folder}/calibration.png')
#            



#        #plot balance (if data to calibrate)
##        if self.flume_else_basin==1:
#        if coarse_dt_flag==1:
#            cum_in=np.zeros((n))
#            cum_out=np.zeros((n))
#            if self.flume_else_basin==1: cum_out_exp=np.zeros((len(t_Qs)))
#            store=np.zeros((n))
#            slid_sh=slid_sh/2.65e6  #[m3/s]      
#            i=1
#            cum_in[i-1]=sum(slid_sh[i-1,0,:])*dtf[i-1]*60
#            cum_out[i-1]=sum(ck[i-1,nR-1,:])*Q[i-1,nR-1]*dtf[i-1]*60  
#            fact_Qs_real=15*60/2.65e6   #t_Qs has constant dt=15min
#            #replace NaN by average neighbor value to get cumulative volume from Qs_real
#            inan=1
#            if self.flume_else_basin==1:
#                Qs_real[0]=Qs_real[1] #Nan replacement exception to following case-specific rule:
#                for i in range(1+10,len(t_Qs)+1): 
#                    #I reviewed manually that forward mobile averaging window with 10 is feasible. Last 10 values are not NaN in data
#                    if np.isnan(Qs_real[i-1]):                
#                        if i==inan+1:    #is not first NaN occuring consecutively 
#                            Qs_real[i-1]=Qs_real[i-2]   #smooth fill of NaNs with same value, to make visible the fill in plot
#                        else:
#                            Qs_real[i-1]=np.mean(Qs_real[i-2:i-11])
#                inan=i
#            if self.flume_else_basin==1: cum_out_exp[i-1]=Qs_real[i-1]*fact_Qs_real
#            for i in range(2,n+1):    
#                #cumulative sediment output in downstream reach            
#                cum_in[i-1]=cum_in[i-2]+sum(slid_sh[i-1,0,:])*dtf[i-1]*60   #all grain sizes            
#                cum_out[i-1]=cum_out[i-2]+sum(ck[i-1,nR-1,:])*Q[i-1,nR-1]*dtf[i-1]*60    #all grain sizes                        
#            if self.flume_else_basin==1:
#                for i in range(2,len(t_Qs)+1):
#                    cum_out_exp[i-1]=cum_out_exp[i-2]+Qs_real[i-1]*fact_Qs_real
#    #        if nW>1:
#            dx=1*(self.dx_flume*np.ones((nR)) if self.flume_else_basin==1 else X[:,6])
#            for i in range(1,n+1):    
#                #sediment storage in system
#                for j in range(1,nR+1):    #integrate longitudinal            
#                    for kk in range(1,int(indW_max[j-1])+1):   #integrate transversal
#                        store[i-1]=store[i-1]+Wflow[j-1,kk-1]*z[i-1,j-1,kk-1]*dx[j-1]
#            fig,ax=plt.subplots(nrows=2,ncols=1)
#            ax[0].clear()                    
#            ax[0].plot(t,store*(1-.3),label='storage')  #.3= assumed porosity
#    #        ax[0].plot(t,store*(1-.3)+)
#            ax[0].set_ylabel('sedim.store[m3]',fontsize=10)        
#            leg = ax[0].legend(loc=1,fancybox=True, fontsize=6)
#            leg.get_frame().set_alpha(0.3)                                                          
#            if t_scale_log==1:
#                ax[0].set_xscale('log')        
#            ax[1].clear()
#            ax[1].plot(t,cum_in,'--k',label='Input')        
#            ax[1].plot(t,cum_out,'k',label='Output')
#            if self.flume_else_basin==1: ax[1].plot(t_Qs,cum_out_exp,'.k',label='Output_exp')
#            ax[1].set_ylabel('cumul.sedim.flux [m3]',fontsize=10)
#            leg = ax[1].legend(loc=1,fancybox=True, fontsize=6)
#            leg.get_frame().set_alpha(0.3)                                                  
#            if t_scale_log==1:
#                ax[1].set_xscale('log')              
#            plt.savefig(f'{self.output_folder}/sed_mass_balanc.png')                                           
                    
        #plot local incision                          
        fig,ax=plt.subplots(nrows=3,ncols=1)
        ymn = np.min(z)
        ymx = max(np.max(z), np.max(zw)) 
        for j in range(1,nR+1):
            ax[0].clear()
            ax[0].plot(t,zw[t_ref_ini-1:t_ref_fin,j-1],'b',label='zw',marker='',markersize=1.,linewidth=.5)
            for kk in range(1,int(indW_max[j-1])+1):
                ax[0].plot(t,z[t_ref_ini-1:t_ref_fin,j-1,kk-1],color=col_kk[kk-1],\
                label= 'z_kk' + str(kk), marker='', markersize=1., linewidth=.5)
                
            leg = ax[0].legend(loc='best',fancybox=True, fontsize=6)
            leg.get_frame().set_alpha(0.3) 
            ax[0].set_ylim(ymin=ymn, ymax=ymx)                               
            ax[0].set_ylabel('alluvium thickness (m)',fontsize=10)
            ax[0].set_title(str(j) + ', darker deeper')  
            if t_scale_log==1:
                ax[0].set_xscale('log')                              
            ax[1].clear()
            for kk in range(1,int(indW_max[j-1])+1):                                     
                ax[1].plot(t,Da84[t_ref_ini-1:t_ref_fin,j-1,kk-1],'-',color=col_kk[kk-1],label='pct84, kk'+str(kk))
                ax[1].plot(t,Da50[t_ref_ini-1:t_ref_fin,j-1,kk-1],':',color=col_kk[kk-1],label='pct50, kk'+str(kk))            
#            leg = ax[1].legend(loc='best',fancybox=True, fontsize=6)
#            leg.get_frame().set_alpha(0.3)                                
            ax[1].set_ylabel('Da|.50|84 (mm)',fontsize=10)
            if t_scale_log==1:
                ax[1].set_xscale('log')            
#            ax[1].set_ylim(ymin=1e-1)            
            if flume_else_basin==0:
                ax[1].set_yscale('log')            
            ax[2].clear()
            for kk in range(1,int(indW_max[j-1])+1):                                     
                ax[2].plot(t,Dss50[t_ref_ini-1:t_ref_fin,j-1,kk-1],color=col_kk[kk-1],label='kk'+str(kk))
#            leg = ax[2].legend(loc='best',fancybox=True, fontsize=6)
#            leg.get_frame().set_alpha(0.3)    
            ax[2].set_ylabel('Dss50(mm)',fontsize=10)
#            if flume_else_basin==0:
#                ax[2].set_yscale('log')            
            if t_scale_log==1:
                ax[2].set_xscale('log')            
#            ax[2].set_ylim(ymin=1e-2)
#            ax[1].set_yscale('log')                                                            
            plt.savefig(f'{self.output_folder}/incis_j'+str(j)+'.png')
                      
        #plot long term geomorphic traits: W and So            
        fig,ax=plt.subplots(nrows=2,ncols=1)    
        for j in range(1,nR+1):
            ax[0].clear()                                     
            ax[0].plot(t,So[t_ref_ini-1:t_ref_fin,j-1],'ok',markersize=.5)
            ax[0].set_ylabel('So[m/m]',fontsize=10)
            ax[0].set_title('j= ' + str(j))  
            if t_scale_log==1:
                ax[0].set_xscale('log')            
#            ax[0].set_ylim(ymin=1e-4,ymax=1e-1)  
#            ax[0].set_yscale('log')
            ax[1].clear()                                     
            ax[1].plot(t,W[t_ref_ini-1:t_ref_fin,j-1],'ok',markersize=.5)
            ax[1].set_ylabel('W[m]',fontsize=10)
            if t_scale_log==1:
                ax[1].set_xscale('log')            
            plt.savefig(f'{self.output_folder}/midTerm_evol_j'+str(j)+'.png')
        print ('Wfeas= ', Wfeas)
#        print ('cum_out_exp= ',cum_out_exp)
        print ('nD= ', nD)

    def xmean_hg(self, t, Q, d, W, indW, Wflow, v, nR):
        t_last_meas_rising = .9 + np.arange(0,10)
        t_last_meas_falling = .9 + np.arange(10,15)
        t_last_measure = np.concatenate([t_last_meas_rising, t_last_meas_falling])
        n_tlm = len(t_last_measure)
        d_gauge = np.zeros((nR))
        v_gauge = np.zeros((nR))
        d_xmean = np.zeros((n_tlm))
        W_xmean = np.zeros((n_tlm))
        v_xmean = np.zeros((n_tlm))
        ind_t_list = []
        for i in range(1, n_tlm + 1):
            ind, = np.where(t == t_last_measure[i-1])
            ind_t_list.append(ind[0])
            W_xmean[i-1] = np.mean(W[ind,:])
            for j in range(1, nR+1):
                kk_flow= int(indW[ind, j-1])
                d_gauge[j-1] = np.mean(d[ind, j-1, 0:kk_flow-1])
                cross_area = np.sum(np.multiply(d[ind, j-1, 0:kk_flow-1], Wflow[j-1, 0:kk_flow-1]))
                v_gauge[j-1] = Q[i-1]/cross_area
            d_xmean[i-1] = np.mean(d_gauge)
            v_xmean[i-1] = np.mean(v_gauge)   
        return d_xmean, W_xmean, v_xmean, t_last_meas_rising, ind_t_list        
        
    def fract_tpt_get(self, pctD, ck_out, c_out_allk, nD, ind_t_list, n_tlmr):
        #bed grainSize distribution
        fi = np.zeros((nD))
        fi[0] = pctD[0]
        for i in range(2, nD+1):
            fi[i-1] = pctD[i-1] - pctD[i-2]
        fi = fi/100
        #bed grainSize distribution            
        pi = np.zeros((n_tlmr, nD))
        for i in range(1, n_tlmr+1):
            ind_inf = 1*(0 if i==1 else ind_t_list[i-2]+1)
            vol_k = np.sum(ck_out[ind_inf:ind_t_list[i-1], :], axis = 0)
            vol = np.sum(c_out_allk[ind_inf:ind_t_list[i-1]])
            pi[i-1, :] = vol_k/vol
        return pi, fi        
        
    def c4trapCalib_get(self, t, c_out, Q_out):
        #array of each measurement time
        add_t = np.array([15, 30, 55])/60
        t_hours = np.arange(0, 9)
        t_meas_intraQ = []
        for i, ta in enumerate(t_hours):
            t_new = ta + add_t
            t_meas_intraQ.append(t_new.tolist())
        t_meas_intraQ = np.ravel(t_meas_intraQ)
        n_tm = len(t_meas_intraQ)
        #aggregation of transport according to t_meas_intraQ        
        c_mean = np.zeros((n_tm, 3))
        Qs_mean = np.zeros((n_tm))
        ind_prev = 0
        ii_meas_given_Q=0
        t_meas_intraQ = np.array(t_meas_intraQ)
        t_meas_intraQ = np.around(t_meas_intraQ, decimals = 3)
        t = np.around(t, decimals = 3)
        for i in range(1, n_tm + 1):
            t_inf=t_meas_intraQ[i-1]-1e-3  #dealing with rounding
            t_sup=t_meas_intraQ[i-1]+1e-3
            ind, = np.where((t>t_inf) & (t<t_sup))
            c = np.mean(c_out[ind_prev:ind[0]])
            Qs_mean[i-1] = c*Q_out[ind[0]]*2.65e6    #g/s)
            ii_meas_given_Q = 1*(0 if ii_meas_given_Q==2 else ii_meas_given_Q + 1)
            c_mean[i-1, ii_meas_given_Q] = c
            ind_prev = ind[0] + 1
        return t_meas_intraQ, Qs_mean, c_mean, n_tm       

    def balance_flume_trench(self, t, z, Q_out, c_out_allk, indW, Wflow, nR, nW, t_last_meas_rising, dx, dt):
        n_tlmr = len(t_last_meas_rising)
        store = np.zeros((n_tlmr, 2))  #rows 1 and 2 are channel and bank respectively respectively*
        cum_out = np.zeros((n_tlmr))
        ind_prev = 0
        for i in range(1, n_tlmr + 1):
            ind, = np.where(t == t_last_meas_rising[i-1])
            output_increms = c_out_allk[ind_prev:ind[0]]*Q_out[ind_prev:ind[0]]*dt*60
            cum_out[i-1] = np.sum(output_increms)
            if i>1: cum_out[i-1] = cum_out[i-2] + cum_out[i-1]
            ind_prev = ind[0] + 1
            for j in range(1, nR+1):
                kk_flow= int(indW[ind, j-1])
                kk_wetdry = [[0, kk_flow-1], [kk_flow, nW-1]]     #*
                for kkwd in range(1, 2+1):  #both wet and dry regions of crosssection
                    kki = kk_wetdry[kkwd-1][0]
                    kkf = kk_wetdry[kkwd-1][1]
                    sed2add = 1*(0 if (kkwd==2 and kki>kkf) else np.sum(np.multiply(z[ind, j-1, kki:kkf], Wflow[j-1, kki:kkf])))
                    #past conditional do not add sediment if no bank exist locally due to flume width occupying whole flume W.
                    store[i-1, kkwd-1] = store[i-1, kkwd-1] + sed2add*dx
        store = store*2650*(1-.3)  #m3 to kg, porosity (.3?) discouded
        for kkwd in range(1, 2+1):          
            store[:, kkwd-1] = store[:, kkwd-1] - store[0, kkwd-1]  #assume store starts at 0.
        cum_out = cum_out*2650
        return store, cum_out
        
    def singl_balance_flume_trench(self, t, z, Wflow, nR, t_last_meas_rising, dx):
        n_tlmr = len(t_last_meas_rising)
        sing_store = np.zeros((n_tlmr))  #rows 1 and 2 are channel and bank respectively respectively*
        ind_prev = 0
        for i in range(1, n_tlmr + 1):
            ind, = np.where(t == t_last_meas_rising[i-1])
            ind_prev = ind[0] + 1
            for j in range(1, nR+1):
                sed2add = np.sum(np.multiply(z[ind, j-1, :], Wflow[j-1, :]))
                sing_store[i-1] = sing_store[i-1] + sed2add*dx
        print('sing_store', sing_store)
        sing_store = sing_store*2650*(1-.3)  #m3 to kg, porosity (.3?) discouded
        sing_store = sing_store - self.ini_store  #assume store starts at 0.
        return sing_store           
        
    def gradients_during_exp(self, dx, nR, y_raw, ind_t_list, n_tlmr):
        #raw variables are spatio-temporal. Here, spatial long profiles are picked at target times.
        R2 = []
        grad = []
        x = np.arange(0, nR)*dx        
        print('x=', x)
        for i in range(1, n_tlmr+1):
            y = y_raw[ind_t_list[i-1], :]   #get longitudinal profile in specific time
            print('y=', y)
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y) 
            R2.append(r_value**2)
            grad.append(slope)
        return grad, R2   
        
    def plot_gradients(self, Q, grad, R2, n_tlmr, output_folder, name_GradVar):
        fig,ax = plt.subplots(nrows=2)
        q_ok= Q[0:n_tlmr]
        ax[0].plot(q_ok, grad)
        ax[1].plot(q_ok, R2)
        plt.savefig(f'{output_folder}/' + 'grad_' + f'{name_GradVar}' + '.png')
        plt.clf()

    def plot_hg(self, Q, d, W, v, output_folder):
        fig,ax=plt.subplots(nrows=3,ncols=1)
        ax[0].plot(Q[0:10] , d[0:10], 'ok', Q[10:15] , d[10:15], '+k')       #diff plotChar if raising or falling limb
        ax[0].set_ylabel('d(m)')
        ax[1].plot(Q[0:10] , W[0:10], 'ok', Q[10:15] , W[10:15], '+k')
        ax[1].set_ylabel('W(m)')
        ax[2].plot(Q[0:10] , v[0:10], 'ok', Q[10:15] , v[10:15], '+k')
        ax[2].set_ylabel('v(m/s)')
        ax[2].set_xlabel('Q (l/s)')        
        plt.savefig(f'{output_folder}/hidra_geom.png')
        plt.clf()    

    def plot_fract_tpt(self, pi, fi, Dcentr, Qtrf, output_folder, n_tlmr):
        fig,ax = plt.subplots()
        col_ref = 1.1*max(Qtrf)
        D = Dcentr*1e3
        for i in range(1, n_tlmr+1):
            q = Qtrf[i-1]
            col = (q/col_ref, q/col_ref, q/col_ref) 
            print('i, q', i, q)
            ax.plot(D, pi[i-1,:]/fi, color = col, label = str(int(q)) + 'l/s')
        ax.set_xlabel('D (mm)')        
        ax.set_ylabel('pi/fi')            
        ax.set_xscale('log')        
        ax.set_yscale('log')
        leg = ax.legend(loc='best', fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.3)        
        plt.savefig(f'{output_folder}/fract_tpt.png')
        plt.clf()                        
        
    def plot_stores(self, Q, store, cum_out, output_folder):
        fig,ax = plt.subplots()
        ls = ['-k', '--k', '+k']
        lbl = ['channel', 'hillslope', 'output']
        for kkwd in range(1, 2+1):
            st = store[0:10, kkwd-1]*1*(-1 if kkwd==2 else 1)    #hillslope supply is negative of bank store change
            ax.plot(Q[0:10] , st, ls[kkwd-1], label = lbl[kkwd-1])
        ax.plot(Q[0:10] , cum_out[0:10], ls[2], label = lbl[2])
        ax.plot(Q[0:10], -store[0:10, 0] - store[0:10, 1] - cum_out[0:10], label= 'errBalance')    #hillSupply - chStore - output == 0
        ax.set_xlabel('Q (l/s)')
        ax.set_ylabel('kg')
        leg = ax.legend(loc=2  ,fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.3)
        plt.savefig(f'{output_folder}/storeVsQ.png')
        plt.clf()
        
    def plot_single_store(self, Q, st, cum_out, output_folder):
        fig,ax = plt.subplots()
        st = st[0:10]
        ax.plot(Q[0:10] , -st, '-k', label = 'alluv.erosion')
        ax.plot(Q[0:10] , cum_out[0:10], '+k', label = 'output')
        ax.plot(Q[0:10], -st - cum_out[0:10], label= 'errBalance')    #hillSupply - chStore - output == 0
        ax.set_xlabel('Q (l/s)')
        ax.set_ylabel('kg')
        leg = ax.legend(loc=2  ,fancybox=True, fontsize=10)
        leg.get_frame().set_alpha(0.3)
        plt.savefig(f'{output_folder}/singleStoreVsQ.png')
        plt.clf()             
        
    def plot_c4trapCalib(self, t, Qs, Q, c, n_tm, output_folder):
        fig, ax = plt.subplots(nrows=2,ncols=1)
        ax[0].plot(t, Qs, 'ok')
        ax[0].set_xlabel('time (hours)')
        ax[0].set_ylabel('Qs (g/s)')
        ax[0].set_ylim(ymin=1e-2, ymax=1e2)
        ls = ['.k', '+k', '*k']
        lbl = ['15min', '30min', '55min']
        indx_Q = t.astype(int)  #get hour of measurement
        Q_allMeasures = np.zeros((n_tm))
        for i in range(1, n_tm+1):        
            Q_allMeasures[i-1] = Q[indx_Q[i-1]]
        for i in range(1, 3+1): #each measure given Q
            ax[1].plot(Q_allMeasures, c[:,i-1], ls[i-1], label = lbl[i-1])
        ax[1].set_xlabel('Q (l/s)')
        ax[1].set_ylabel('c (m3/m3)')
        ax[1].set_ylim(ymin=1e-6, ymax=1e-2)        
        leg = ax[1].legend(loc='best'  ,fancybox=True, fontsize=6)
        leg.get_frame().set_alpha(0.3)
        for b in [0, 1]:
            ax[b].set_yscale('log')            
        plt.savefig(f'{output_folder}/Qs_aggreg.png')
        plt.clf()
