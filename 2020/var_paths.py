    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#jan1, 2020

import numpy as np
import sys


class a:
    def __init__(self, HGelseFract):
        self.HGelseFract = HGelseFract
    
        
    def get_vars_dict(self, nD):
        #dict_copies is necessary to traslate copies name to copies amount, as nD might coincide with nW.
        if  self.HGelseFract == 1:
            dict_vars={'Qh': 'unique',\
            'wash_sh': 'unique',\
            'slid_sh': 'unique',\
#            'Mgtm': 'unique',\
#            'gam': 'unique',\
            'v': 'unique',\
            'h': 'unique',\
            'vf': 'unique',\
            'w': 'unique',\
            'wproy': 'unique',\
            'c': 'unique',\
            'h_s': 'unique',\
            'h_sf': 'unique',\
            'w_c': 'unique',\
#            'tsm_c': 'unique',\
            'D50c': 'unique',\
            'D84c': 'unique',\
            'D50c_ss': 'unique',\
            'D84c_ss': 'unique',\
            'Qbf': 'unique',\
            'Qsbf': 'unique',\
            'Q': 'unique',\
            'So': 'unique',\
            'w_b': 'unique',\
            'ts_c': 'unique',\
            'tc_c': 'unique',\
            'Fr': 'unique',\
            'aspect': 'unique',\
            'subm': 'unique',\
            'h_bf': 'unique',\
            'volIn': 'unique',\
            'dvol': 'unique',\
            'colm': 'unique'}
#            'pi': 'nD'}
            dict_copies={'unique':1, 'nD':nD} #unique means only one copie of 2D array is needed; 
                #i.e. fortran output variable has rank<3.
        else:
#            nD_hill= 1*(nD if flume_else_basin==1 else nD-1)  #if basin then boulder is not supplied from subbasin but in river.
            nD_hill = nD
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
        return dict_vars, dict_copies
        

    def path_list_builder(self, nW, nD, of):    #, flume_else_basin
        #DYNAMIC PATHS TO SAVE AND PLOT KEY TIME SERIES:
        #dictionary. var_name: number_of_copies. Copies are needed when 3rd dimension needs to be plotted, because I only know\
        #method to save 2D arrays in binary format from fortran code.        
        
        #variables 1 to 14. # in ':#' means 3rd dimension of array, additional to t=i and x=j
        #new var to plot?: (1) py: add to this dict. (2) f90: add to save bin with its nvar. (3) py: declare, read and plot
             
        dict_vars, dict_copies = self.get_vars_dict(nD)
        print('25apr20: dict_vars= ', dict_vars)
        vars2plot=[]        #get list from dictionary to iterate with numeric values
        n_copies_name=[]
        for i,j in dict_vars.items():
            vars2plot.append(i)
            n_copies_name.append(j)
        nvars=len(vars2plot)
        n_copies=[]
        for i in range(1,nvars+1):
            n_copies.append(dict_copies[n_copies_name[i-1]])
#        print('n_copies= ', n_copies)                    
        len_vars=[len(var) for var in vars2plot]
        len_max_var=max(len_vars)+6 #4: maybe double subindex letter and index greater than 10.  
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
#        path_list2f90=np.empty((self.num_paths,self.num_char),dtype='c') #1st arg is tuple of (#elem,#charPerElem)
        print('25apr20: just PRE pass string array to f90 format, path_list= ', path_list)
        path_list2f90=np.array(path_list,dtype='c').T
        print('25apr20: just POS pass string array to f90 format, path_list2f90= ', path_list2f90)
            #trick: https://stackoverflow.com/questions/22293180/passing-numpy-string-format-arrays-to-fortran-using-f2py/23204241
        num_paths = len(path_list2f90[0,:])
#        for i in range(1, num_paths+1):  #path list in py format is filled with folder, to be used in plotting binaries
#            path_list[i-1]= of + '/' + path_list[i-1]
        return path_list, path_list2f90, pos_ini_var, nvars, len(path_list2f90[:,0]), num_paths #, nvars
