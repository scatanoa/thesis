    #!/usr/bin/env python
# -*- coding: UTF-8 -*-
#jan11, 2020.

import numpy as np
import matplotlib.pyplot as plt

A_km2= [13,23,59,119];  tc_min= [120,219,350,940]
fig,ax=plt.subplots()
ax.plot(A_km2, tc_min, '.k')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(" A [km2]")
ax.set_ylabel("tc [min]")

#ax.set_ylim(bottom=1e-4)

#plt.savefig(f'{output_folder}/HillRatingCurve_reach{j}.png')
#plt.close()

plt.tight_layout()
plt.show()
