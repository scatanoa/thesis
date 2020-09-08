import numpy as np

import readParamFromDict as rpfd
import Ptx2supply as s
import FolderCleaner as fc
import list_str2float as lsf


fc.cleanFolder('doct_supply/supply_wat_sed') #dedeprec by 27jun20, to clean old Qh figs as they are titled with Qm so no updated.

Ptx = np.loadtxt('doct_supply/PperSubb/P_tx.txt') #produced by topolPlusRain_PperSubb.py
pd = rpfd.readParam('param_dict.txt')
n = ['dT_days', 'A_hs', 'elong_ratio', 'hgy', 'geotech', 'geotech2', 'weath_coefPerD', 'weath_xrefPerD', 'print_range', 'flags']
s.lumped(Ptx, int(pd[n[0]]), float(pd[n[1]]), float(pd[n[2]]), lsf.f(pd[n[3]]), lsf.f(pd[n[4]]), \
    lsf.f(pd[n[5]]), lsf.f(pd[n[6]]), lsf.f(pd[n[7]]), lsf.f(pd[n[8]]), lsf.f(pd[n[9]]))
