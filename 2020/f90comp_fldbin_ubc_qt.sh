#!/bin/bash

#updated with newton nonlinear solver by mar17, 2020.
gfortran -c ddt_class_fld.f90   balance_fldbin_ubc_qt.f90
f2py -c  ddt_class_fld.f90   balance_fldbin_ubc_qt.f90   hg_model_fldbin_ubc_qt.f90 -m hg_model_fldbin_ubc_qt

#17apr: *** Error in `/home/scatanoa/venv_2020phd/bin/python3.6': double free or corruption (out): 0x0000000003ca8990 ***
#    flags (options, after '-c') --debug --verbose did not help.     





##try3: no uppercase in filenames: (updated with balance_class by jan11, 2020)
#gfortran -c ddt_class.f90 balance_class.f90
#f2py -c ddt_class.f90 balance_class.f90 hg_model_class.f90 -m hg_model_class


#jan3, 2020. I think this code worked in pcu by dec2020:
#try0:
#gfortran -c ddt_class.f90
#f2py -c ddt_class.f90 HGmodel_class.f90 -m HGmodel_class


#jan3, 2020:
#try1:
#f2py -c ddt_class.f90 -m ddt_class
#this code produces same bug, so it is not related to multicompile to have HGmodel_class using module ddt_class.
#Maybe is an issue related to numpy version, because the following name of extension module shows 35 instead 36, 
#even with python version being 3.6:
    #ddt_class.cpython-35m-x86_64-linux-gnu
    
#try2:
#"import 'module'" instead of "from 'module' import *": Did not work.
