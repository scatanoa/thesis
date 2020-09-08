import topol
import FolderCleaner as fc
import readParamFromDict as rpfd
import list_str2float as lsf

def fn(flume_else_basin, which_flume):
    top = topol.Streams(0, flume_else_basin, which_flume)
