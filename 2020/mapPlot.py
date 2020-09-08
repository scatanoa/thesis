    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#dec20, 2019
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
            
def mapPlotter(matrix, relative_folder_path, filename, cbarLabel, 
    xl, yl, mn=None, mx=None, log_cbar=False):
    cm = plt.get_cmap('magma')
    fig, ax1 =plt.subplots()
    if log_cbar:  
        im = ax1.imshow( matrix, interpolation='nearest', cmap=cm, aspect = 'auto',
            norm = LogNorm(vmin = mn) )
    else:
        im = ax1.imshow(matrix, interpolation='nearest', cmap=cm, aspect = 'auto',
            vmin=mn, vmax=mx)
    cbar= fig.colorbar(im, ax=ax1)
    cbar.set_label(cbarLabel)
    ax1.set_xlabel(xl)
#    ax1.set_ylabel(yl, rotation=270)
    ax1.set(ylabel=yl) #, rotation=270)    
    fig.tight_layout()        
    plt.savefig(relative_folder_path + '/' + filename + '.png')
    plt.clf
