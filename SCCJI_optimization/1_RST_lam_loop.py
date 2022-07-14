# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:35:12 2021

@author: Mirko Pavoni
"""

#### Seismic refraction Inversion ###
### NB before runnig the code, be sure that you selected the right path (where you saved this code and the file.sgt) in the top right of the screen (e.g. C:\Users\Mirko\Desktop\Schafberg_4PM\cond_inv_inv\rst)
# Import the packages needed to realize the inversion modelling
import pygimli as pg
from pygimli.physics import TravelTimeManager
import matplotlib.pyplot as plt
import numpy as np
import pygimli.meshtools as mt

# Load the picking file (.sgt format)
data = pg.load("rst_filtered2.sgt")
# remove the measurements which have the same shoot and geophone position 
data.remove(data["g"] == data["s"])
# Define the error in the seismic dataset in sec (e.g. 0.0003 means 0.3 msec)
data.set("err", pg.Vector(data.size(), 0.0006))
print(data)

# see traveltime curves  
fig, ax = pg.plt.subplots()
pg.physics.traveltime.drawFirstPicks(ax, data)

# Initialize the refraction manager
mgr = TravelTimeManager()

# matrix of apparent velocities 
ax, cbar = mgr.showData(data)

nv = np.linspace(1,11,8)
li = 5
s = 0
for i in range(len(nv)):
    li = li*2 
    # run the inversion, the mesh is created based on the sensors and shots positions
    inverted = mgr.invert(data, lam = li, secNodes=2, paraMaxCellSize=2.0, maxIter=20, zWeight=1, vTop=500, vBottom=5000,
               verbose=True)

    # plt.rcParams['figure.figsize'] = [8, 6]
    # vp , cbar = mgr.showResult(# color map settings, Cmin = Vp minimum of the colorbar, Cmax = Vp maximum value of the colorbar
    #                       cMap='seismic_r',cMin=350,cMax=3500,
    #                       sens=True, contour=True,
    #                       # colorbar settings
    #                       logScale = False, orientation='horizontal',
    #                       # show ray paths
    #                       rays = False)
    # vp.set_xlim(-1, 95)
    # vp.set_ylim(2735, 2775 )
    # # Sensors and shot (use "diam" to increase/decrease the size of sensors/shot)
    # pg.viewer.mpl.drawSensors(vp, data.sensors(), diam=1,
    #                          facecolor='white', edgecolor='black')
    # plt.savefig('vpsec_%s.png' %s, dpi=600)
    
    # plt.rcParams['figure.figsize'] = [8, 6]
    # fig = mgr.showFit() 
    # plt.savefig('vpfit_%s.png' %s, dpi=600)
    
    file = open("parameters_%s.txt" % s, "w")
    lambdaa = repr(li)
    chii = mgr.inv.chi2()
    chiii = repr(chii)
    relrmss = mgr.inv.relrms()
    relrmsss = repr(relrmss)
    file.write("\n" + "lamdav = " + lambdaa + "\n" + "chi2 = " + chiii + "\n" + "rrms = " + relrmsss)
    file.close()
    
    s = s+1
    
