    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#dec18, 2019

#pre-existing modules
import numpy as np
import matplotlib.pyplot as plt
import sys

import mapPlot as mp
import rasterizeStraighLine as rsl
import readParamFromDict as rpfd

class Intercepter: 
#for spatial matrixes, (x, y) denotes (down, right) positions. See method of Bresenham's line algorithm in wiki article.
#SubbPerim: SubbasinPerimeter
    def __init__(self):    #elong_ratio increases with tectonic activity (To explore with SergioRest?).
        X = np.loadtxt('doct_topol/Topol_X.txt') #, dtype = {'formats': } !!deprec:fmt= ' '.join(['%1.1f']*2 + ['%i']*3 + ['%1.3f'] + ['%i']*2))
        ampl_X = 3
        X[:, 0: 1+1] = ampl_X*X[:, 0:1+1]    #as Topol_X.txt saved initial reach points (1st 2 cols of X matrix) with precision: float with 1 decimal.
            #also bresenham algorithm demands integer inputs (extreme points of straight line).
#        print('X= ', X)            
        self.X = X
        param_dict = rpfd.readParam('param_dict.txt')
        A_hs= float(param_dict['A_hs'])
        L_reach = X[:, 6]
        self.L_hs = np.round(A_hs/(2*L_reach))
        self.lhlr = self.L_hs/L_reach
#        print('L_reach= ', L_reach)
#        print('L_hs= ', self.L_hs)
#        print('lhlr= ', self.lhlr)        
        self.t_ev = np.loadtxt('doct_rain/t_ev.txt')
        self.Px0 = np.loadtxt('doct_rain/Px_day' + str(int(self.t_ev[0])) + '.txt')
        self.Px0_sh = self.Px0.shape #1st and 2nd components say rows/y/northing and columns/x/easting lenght respectively.
        
        
    def reachExtremPts(self, reachID):
        X= self.X
        point_upstream = (X[reachID-1, 0], X[reachID-1, 1])
        downstrRow = int(X[reachID-1, 3])
        point_downstream = (X[downstrRow-1, 0], X[downstrRow-1, 1]) 
        print('in reachExtremPts, reachID= ', reachID)
        print('point_upstream= ', point_upstream)
        print('point_downstream= ', point_downstream)
        print('downstrRow= ', downstrRow)        
        return point_upstream, point_downstream
        
    def SubbPerim_rectang_vector(self, Lhill_Lreach_ratio, pt_u, pt_d): #pt_u, pt_d: point_upstream, point_downstream.   
    #For rectangular subbasin perimeter, get corner points (x,y) from left hill upstream and rotating clockwise.
        dx_reach = pt_d[0] - pt_u[0]
        dy_reach = pt_d[1] - pt_u[1]                
        #get rectangle corners
        dx_wdh_l = -dy_reach*Lhill_Lreach_ratio   #wdh_l = waterDivideHillsl_left
        dy_wdh_l = dx_reach*Lhill_Lreach_ratio
        dx_wdh_r = dy_reach*Lhill_Lreach_ratio
        dy_wdh_r = -dx_reach*Lhill_Lreach_ratio        
        pt_ul = [pt_u[0] + dx_wdh_l, pt_u[1] + dy_wdh_l]
        pt_dl = [pt_d[0] + dx_wdh_l, pt_d[1] + dy_wdh_l]
        pt_dr = [pt_d[0] + dx_wdh_r, pt_d[1] + dy_wdh_r]
        pt_ur = [pt_u[0] + dx_wdh_r, pt_u[1] + dy_wdh_r]        
        SubbCornerPts = [pt_ul, pt_dl, pt_dr, pt_ur]
#        print('in SubbPerim_rectang_vector, REAL SubbCornerPts', SubbCornerPts) 
        return SubbCornerPts

    def rasterize(self, reachID, offset_x, y_mx, rescaling_factor, SubbCornerPts):
        #get intermediate points between corners, via Bresenham's line algorithm.
        perimPts = []
        rot_path = [2, 3, 4, 1] #path of arrow head locations of clockwise rotation through rectangle corners.        
        print('in rasterize, rescaling_factor = ', rescaling_factor)
        SubbCornerPts = rescaling_factor*np.array(SubbCornerPts)
#        print('in rasterize, REAL Rescaled SubbCornerPts= ', SubbCornerPts)        
        SubbCornerPts = [[int(round(x)) for x in y] for y in SubbCornerPts]        
#        print('in rasterize, ROUND INTEGER rescaled SubbCornerPts= ', SubbCornerPts)
        pt_ini = SubbCornerPts[0]        
        print()        
        for i, ptNum in enumerate(rot_path):            
            pt_fin = SubbCornerPts[ptNum-1]
#            print('in rasterize, i= ', i, 'pt_ini = ', pt_ini, 'pt_fin = ', pt_fin)            
            line = rsl.bresenham(pt_ini, pt_fin)
            line_pts = line.path
            for item in line_pts: perimPts.append(item)            
            pt_ini = pt_fin
        perimPts = np.array(perimPts)
#        print('perimPts= ', perimPts)
#        print('offset_x= ', offset_x)
        perimPts[:, 0] = perimPts[:, 0] + offset_x
            #as already rescaled, traslate perim_raster to fit with Px in same origin.
#        print('traslated perimPts= ', perimPts)
        perim_raster = np.zeros((self.Px0_sh[1], self.Px0_sh[0]))
        perim_raster[y_mx - perimPts[:,1], perimPts[:,0]] = 1
#        print('perim_raster= ', perim_raster)
#        self.mapPlotter(perim_raster, 'perim_raster_reachID' + str(reachID))
        return perim_raster 
               

    def addSubbMask(self, reachID,  Mask, rescaling_factor, offset_x, y_mx):
        print('in addSubbMask, reachID= ', reachID)
        pt_u, pt_d = self.reachExtremPts(reachID)
        SubbMask = self.rasterize(reachID, offset_x, y_mx, rescaling_factor, \
            self.SubbPerim_rectang_vector(self.lhlr[reachID-1], pt_u, pt_d))
        #FillCellsInSubbPerim:
        SubbMask_aggrAlongX = np.sum(SubbMask, axis = 1) #sum along columns "x".
        ind_y, = np.where(SubbMask_aggrAlongX > 0)
        print('ind_y= ', ind_y)
        ind_y_mn = min(ind_y)
        ind_y_mx = max(ind_y)
        print('ind_y_mn= ', ind_y_mn)
        print('ind_y_mx= ', ind_y_mx)
        for y in range(ind_y_mn, ind_y_mx+1):
            ind_x, = np.where(SubbMask[y, :] > 0)     #as subbasin is rectangular, 2 indexes must be found.
            print('y= ', y, ', ind_x= ', ind_x)
            ind_x_mn = ind_x[:1][0]
            ind_x_mx = ind_x[-1:][0]
            print('y= ', y, ', ind_x_mn= ', ind_x_mn, ', ind_x_mx= ', ind_x_mx)            
            SubbMask[y, ind_x_mn:ind_x_mx+1] = 1
#        mp.mapPlotter(SubbMask, 'doct_supply', 'SubbMask_reachID' + str(reachID))
        #overlap with previos subbasins. Unfortunately, this reduces whole basin area, even more if subbasin elong_ratio is low:
        UnavailCells = np.zeros((self.Px0_sh[1], self.Px0_sh[0])); AvailCells = np.zeros((self.Px0_sh[1], self.Px0_sh[0]))
        UnavailCells = SubbMask*Mask    #to get cells that are already draining to another subbasin.
        AvailCells[UnavailCells == 0] = 1
        SubbMask = SubbMask*AvailCells
        Mask = Mask + SubbMask*reachID          
        return Mask        
    
    def getRasterSubbasins(self):
        reachID = 0     #must start in 1, instead of 0, to avoid confusion with mask cells not assigned to any subbasin.
        Px0_sh = self.Px0_sh
        X = self.X       
        Mask = np.zeros((Px0_sh[1], Px0_sh[0]))    #each Mask cell has reachID where rain drains to
        n_reaches = len(X[: , 0])
        max_lhlr = max(self.lhlr)
        max_X_sideSize = max(max(X[:,0]) - min(X[:,0]), max(X[:,1]) - min(X[:,1]))
        rescaling_factor = .9*Px0_sh[0]/(max_X_sideSize + 4*max_lhlr) #max_lhlr to avoid boundary troubles.
            #lhlr instead Lh because Lreach=1 in topol X.
        if rescaling_factor<.1: sys.exit('.1 > rescaling_factor = ' + str(round(rescaling_factor, 1)) + '.\
            Improve spatial res parameter N of rain field, to get rescaling_factor > .1 to reduce subb rasterization mistakes')
        offset_x = rescaling_factor*abs(min(X[:,0])) + 2*2*rescaling_factor*max_lhlr
            #topology only has negative values in x coordinate.
        y_mx = int(rescaling_factor*(max(X[:,1]) + 2*rescaling_factor*max_lhlr))
        print('in getRasterSubbasins, rescaling_factor = ', rescaling_factor)
        print('in getRasterSubbasins, offset_x = ', offset_x)
        print('in getRasterSubbasins, y_mx = ', y_mx)
        print('in getRasterSubbasins, size_P_raster = ', Px0_sh[0])
        print('in getRasterSubbasins, size_topol = ', max_X_sideSize + 4*max_lhlr)
            #basin is squared like rain field. Factor 0.9 is for topol2fitRainField.
        for j in range(1, (n_reaches-1)+1): #because last X row is not upstream point of a reach.
            reachID+=1
            print()
            print()
            print()
            print()
            print('-----------------------------------')
            print('in getRasterSubbasins, reachID= ', reachID)            
            Mask = self.addSubbMask(reachID, Mask, rescaling_factor, offset_x, y_mx)
#            print('Mask= ', Mask)
            Mask_bin_currentSubb = np.zeros((Px0_sh[1], Px0_sh[0]))
            Mask_bin_currentSubb[Mask == reachID] = 1
            mp.mapPlotter(Mask_bin_currentSubb, 'doct_supply/geom', 'Mask_bin_currentSubb_reachID' + str(reachID),
                'isSubb','easting','northing')
            np.savetxt('doct_supply/geom/Mask_bin_currentSubb_reachID' + str(reachID)  + '.txt', Mask_bin_currentSubb, fmt='%d')
        np.savetxt('doct_supply/geom/subbasins.txt', Mask, fmt='%d')
        mp.mapPlotter(Mask, 'doct_supply/geom', 'subbasins', 'Subbasin #','easting','northing')
