import topolPlusRain_geom as tpr_g
import FolderCleaner as fc

fc.cleanFolder('doct_supply/geom')

P_tx = tpr_g.Intercepter()
P_tx.getRasterSubbasins()
