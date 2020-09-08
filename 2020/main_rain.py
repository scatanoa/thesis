import rain
import readParamFromDict as rpfd
import FolderCleaner as fc
import list_str2float as lsf
#import rainBin

fc.cleanFolder('doct_rain')
pd = rpfd.readParam('param_dict.txt')
N = float(pd['N']); M = float(pd['M']); p = float(pd['p']); two2theEXPth_FieldSmoothing = int(pd['two2theEXPth_FieldSmoothing'])
dT_days = int(pd['dT_days']); dt = 24 #always daily, at least for dissertation, as required by aggregated hgy model SHIA (tanks).
daysAverBetweenEvents = int(pd['daysAverBetweenEvents'])

rainPars = lsf.f(pd['rainPars'])
flagR0var = int(rainPars[0]); Pm_mmyr = int(rainPars[1]); a_zcit = rainPars[2]

r=rain.P(dT_days*dt, dt, daysAverBetweenEvents, p, N, M, 2**two2theEXPth_FieldSmoothing, flagR0var, Pm_mmyr, a_zcit)
r.xt_genP()


#r=rainBin.Pread(24*30,16.)
#r.Pxt_fine()
