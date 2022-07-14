# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 10:56:27 2021

@author: Mirko
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 10:34:08 2021

@author: Mirko
"""

### ERT Inversion Code ### 
### NB before runnig the code, be sure that you selected the right path (where you saved this code and the file.ohm) in the top right of the screen (e.g. C:\Users\Mirko\Desktop\Schafberg_4PM\cond_inv_inv\ert)
# Import packages and functions
import numpy as np
import pygimli as pg
from pygimli.physics.ert import ERTManager, createGeometricFactors
import matplotlib.pyplot as plt

# Insert the dataset
data = pg.load("ert_filtered2.ohm")
print(data)

# NB The data file does not contain geometric factors (no ‘k’ field is present)
# We can create 'k' based on the given topography
data['k'] = createGeometricFactors(data, numerical=True)

# Initialize the ERTManager for further steps and inversion
ert = ERTManager(sr=False, useBert=True, verbose=True, debug=False)

# The data container has no apparent resistivities (no ‘rhoa’ field data is present)
# We can let the Manager fix and add all the columns (err, i, ip, iperr, k, r, rhoa, u, valid) to perform the inversion:
ert.checkData(data)

# define the inversion error based on the reciprocity check (here for example is 10%)
data['err'] = ert.estimateError(data, absoluteError=0.001, relativeError=0.10)
# remove the negative apparent resistivity values
data.remove(data["rhoa"] <= 0)
print(data)

nv = np.linspace(1,11,8)
li = 5
s = 0
for i in range(len(nv)):
    li = li*2 
    mod = ert.invert(data, lam=li, zweight = 1,
                     paraDX=0.5, paraMaxCellSize=2, quality=33.6)

    # Plot the resistivity section
    # plt.rcParams['figure.figsize'] = [8, 6]
    # ax , cbar = ert.showResult(cMap='RdYlBu', cMin = 5000, cMax = 200000)
    # # Plot electrodes position (use "diam" to increase/decrease the size of the electrodes in the figure)
    # pg.viewer.mpl.drawSensors(ax, data.sensors(), diam=1,
    #                          facecolor='white', edgecolor='black')
    # ax.set_xlim(-1, 95)
    # ax.set_ylim(2735, 2775 )
    # plt.savefig('ressec_%s.png' %s, dpi=600)

    # # Plot the quality of the obtained model
    # plt.rcParams['figure.figsize'] = [8, 6]
    # fig = ert.showFit(cMap='jet')
    # plt.savefig('resfit_%s.png' %s, dpi=600)
    
    file = open("parameters_%s.txt" % s, "w")
    lambdaa = repr(li)
    chii = ert.inv.chi2()
    chiii = repr(chii)
    relrmss = ert.inv.relrms()
    relrmsss = repr(relrmss)
    file.write("\n" + "lamdav = " + lambdaa + "\n" + "chi2 = " + chiii + "\n" + "rrms = " + relrmsss)
    file.close()
    
    s = s+1
    

# run the inversion defining the parameters for the inversion and the mesh 
# lam = Global regularization parameter lambda
# zweight = smoothing in lateral direction (e.g. 0.25 = four times more in horizontal direction)
# paraDX = relative distance for refinement nodes between two sensors (1=none), e.g., 0.5 means 1 additional node between two neighboring sensors e.g., 0.33 means 2 additional equidistant nodes between two sensors
# paraMaxCellSize = maximum size for parametricsize in m*m
# quality = 2D triangle quality sets a minimum angle constraint, be careful with values above 34 degrees
