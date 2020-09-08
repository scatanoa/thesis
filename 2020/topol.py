    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#mar19, 2018

#pre-existing modules
import numpy as np
import matplotlib.pyplot as plt
import readParamFromDict as rpfd
import list_str2float as lsf
import join_pp13exp_database as pp13
import sys

#from byte2strarray import *


class Streams: 
    def __init__(self, realBasin_else_fict, flume_else_basin, which_flume):
        pd = rpfd.readParam('param_dict.txt')
        #call spatial topology
        if flume_else_basin==1:
            Sof=.02    #maria: .022
            L_flume=1*9.  #maria: 12
            Sof=float(pd['Sof'])
            L_flume=float(pd['L_flume'])
            self.Wflume = float(pd['Wflume'])
            self.dx_flume=float(pd['dx_flume'])

            self.amplif_D = 1    #for sensitivity analysis in 3rd paper with Marwan, comparing flume and numeric model of lateral input channel.

            repos_angle_ini = 30.
            base_level = .2
            triangl_base = self.Wflume/2.
            rectang_store = base_level*self.Wflume*L_flume
            triang_store = .5*triangl_base**2*np.tan(3.14/180*repos_angle_ini)*L_flume
            self.ini_store_m3 = 2*triang_store + rectang_store
            self.ini_store = 2650*(1-.3)*self.ini_store_m3
#            print('rectang_store[m3]', rectang_store, 'triang_store', triang_store, 'ini_store[kg_sed]', self.ini_store)

            n_nodes=int(L_flume/self.dx_flume) #L must be multiple of dx             
            self.serial_topol(which_flume, n_nodes,self.dx_flume,Sof) #5mar20: self.X=
            print('5mar20, guardó?')
        elif realBasin_else_fict == 1:
            yyy = 1
        else:   #toy basin
#        
#            stoch_rain=1
#            rain_fileName='P36500steps.txt'  
#            wash_flag=1
#            print('param_dict= ', pd)        
                
            self.A_hs= float(pd['A_hs'])  #4e6    #'hillslope' (subbasin) area
            print('27jul20; A_hs: ', self.A_hs)
            self.elong_ratio = float(pd['elong_ratio'])
#            print('from dict, elong_ratio= ', self.elong_ratio)
            self.Omega = int(pd['Omega'])
#            self.Omega=3    #order of drainage network at outlet
            offset_angle_degrees = 0
            
            self.pars_wv_A = lsf.f(pd['fluv_wv_A'])
            self.pars_wc_A = lsf.f(pd['fluv_wc_A'])
            self.Sor_A = lsf.f(pd['fluv_rock_S_A'])  #25may20            
            
            self.X=self.fictNetwork(1.,2., offset_angle_degrees)      #15jun20; a,c values  below eq 5 in mesaGupta14.      
            
            nR=len(self.X[:,0])     #number of reaches                       
            print('X= ', self.X)
                  
        
    def serial_topol(self, which_flume, n_nodes, dx_flume, Sof):       
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
        ons = np.ones((n_nodes, 1))     
        Wv = self.Wflume*ons
        if which_flume == 1:
            Wc = Wv
        else:
            df = pp13.get_df()
            wciv = df['W_half_cm'].to_numpy()
            wci = wciv[0]*2/100
#            print(f'initial channel width, wciv= \n {wci}m')
#            sys.exit('stop here to check Wcini')
            Wc = wci*ons
            Wv[-1] = 1.2    #in pitlick13 photo outlet width seems ~2/5 of flume width.
            Wc[-1] = Wv[-1] #16apr20: to assume simple rigid rectangular critical flow section at outlet.
        X = np.hstack((X, Wv))            
        X = np.hstack((X, Wc))
        np.savetxt('doct_topol/Topol_X_flume.txt', X, \
            fmt= ' '.join(['%1.1f']*2 + ['%i']*3 + 1*['%1.3f'] + 4*['%1.2f']))
                                            
        
    def fictNetwork(self,a,c, offset_angle_degrees):
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
        deg45= pi/4 + offset_angle_degrees * pi/180
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
        
        noisA = np.zeros(n)
        ampl = 1
        nAmn = 1 - .5*ampl
        for i in range(1, n+1):
            noisA[i-1] = nAmn + ampl*np.random.rand()
            #1jun20; to spread var vs A samples in x_scaling validation plots.
#            print(nAmn + np.random.rand()*ampl)
        print('1jun20; noisA: ', noisA)
#        sys.exit('check if noisA OK')
        
        X[:,7] = X[:,7] + noisA*self.A_hs/1.e6    #every reach has at least its local drainage area
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
        X[:,8] = self.Sor_A[0]*X[:,7]**self.Sor_A[1]
            #local bedrock slope in reaches (default coef 08 and exp -.3 from Colombian regression in ideam2017)
        X[n-1,8]=1e-10  #downstream base level fixes alluvium thickness (June8, 2018)
        X[:,9] = X[:,9] + 2*np.sqrt(self.A_hs/np.pi)/self.elong_ratio    #lenght of reaches [m] (then x y coor)
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

#        print "elevations added in last column, X= ", X
        for j in range(1,totPts+1):
            if X[j-1,3]>0:
                X[j-1,3]=X[j-1,3]+1                
#        print "downstream ID + 1 for Fortran compatibility, X= ", X
#        print ("X= ", X)

#        fig,ax=plt.subplots()

        #w_v: valley width
        Wv = self.pars_wv_A[0]*X[:,4]**self.pars_wv_A[1]
        ons = np.ones((n, 1))        
        Wv_mtx = [oi*Wv[i] for i, oi in enumerate(ons)]
        X = np.hstack((X, Wv_mtx))        

        #w_c: initial channel top width
        Wc = []
        for j in range(1, n+1):
            #25may20; w~A scaling params depend on supply or tpt lim, which is S threshold. 
                #See flores06 in jimenez15phd.
            Srj = X[j-1,5] #local bedrock slope
            a = self.pars_wc_A[0] if Srj<.025 else self.pars_wc_A[2] 
            b = self.pars_wc_A[1] if Srj<.025 else self.pars_wc_A[3]
            Wc.append(a*X[j-1,4]**b)
        Wc = np.array(Wc)
        Wc = np.minimum(Wv,Wc) #4jun20; account for rockwall canyons.
        Wc_mtx = [oi*Wc[i] for i, oi in enumerate(ons)]
        X = np.hstack((X, Wc_mtx))

        maxZ=max(X[:,7])
        color = [item/maxZ for item in X[:,7]]
#        print "color=", color
        plt.scatter(X[:,0],X[:,1],c=color)
        plt.colorbar()
        plt.savefig('doct_topol/TopolZcolor_order'+ str(self.Omega) + '.png')
        plt.clf()

        np.savetxt(  'doct_topol/Topol_X.txt',X,fmt= ' '.join( ['%1.1f']*2 + ['%i']*2 + ['%1.2f'] + ['%1.3f'] + 
            ['%i']*2 + ['%1.1f']*2 )  )
#        sys.exit('1jun20; check saved topol')
        
        #X contains these 8 columns (5integer 1real 2integer): 0&1: x,y; 2: hortonStrahlerOrder; 3: IDrowDownstreamReach; 
        # 4:DrainArea_km2; 5: ReachSlope; #6: reachLenght_m; 7: elevationReachIniPoint_masl.        
        return X
