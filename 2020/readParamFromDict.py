    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#dec18, 2019

#pre-existing modules
import numpy as np
def readParam(filename):
    d = {}
    with open(filename) as f:
        for line in f:
            lsp = line.split() #17mar2020: splits where space by default.
            if len(lsp) > 1:    #read neither title (1 word) nor white (0 words) lines.
                (var_name, val) = lsp 
                d[str(var_name)] = val
    return d
