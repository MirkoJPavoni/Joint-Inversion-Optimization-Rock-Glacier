# SCRIPT TO PERFORM CONVENTIONAL INDEPENDENT INVERSION AND SCC JOINT INVERSION

import os
import numpy as np
import pygimli as pg
import matplotlib.pyplot as plt
# from pygimli import meshtools as mt
# from pygimli.physics.petro import transFwdArchieS as ArchieTrans
# from pygimli.physics.petro import transFwdWyllieS as WyllieTrans
from pygimli.physics import ert
from pygimli.physics import traveltime as tt
from scci import SCCI
from plotting import drawCWeight

# %% define mesh and datasets
rMesh = pg.load("mesh_1.bms")
vMesh = pg.load("paraDomain_1.bms")
rKW = dict(logScale=True, cMin=5000, cMax=200000, cMap="Spectral_r")
vKW = dict(logScale=True, cMin=3500, cMax=3500, cMap="Spectral_r")

ertData = ert.load("ert_filtered2.ohm")
ERT = ert.ERTManager(ertData, verbose=True, sr=False)
ERT.setMesh(rMesh)

ttData = tt.load("rst_filtered2.sgt")
TT = tt.TravelTimeManager(ttData, verbose=True)
TT.errIsAbsolute = True
TT.setMesh(vMesh)

# %% ert inversion
ERT.invert(zWeight=1, lam=80, cType=1)

# %% plot resistivity section
plt.rcParams['figure.figsize'] = [8, 6]
ax , cbar = ERT.showResult(cMap='RdYlBu', cMin = 5000, cMax = 200000)
# Plot electrodes position (use "diam" to increase/decrease the size of the electrodes in the figure)
pg.viewer.mpl.drawSensors(ax, ertData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
ax.set_xlim(-1, 95)
ax.set_ylim(2735, 2775)
plt.savefig('res.png', dpi=600)

# Plot the quality of the obtained model
plt.rcParams['figure.figsize'] = [8, 6]
fig = ERT.showFit(cMap='jet')
plt.savefig('res_fit.png', dpi=600)

# %% rst inversion
TT.invert(zWeight=1, lam=180, vTop=500, vBottom=5000,
           verbose=True)

# %% plot vp section
plt.rcParams['figure.figsize'] = [8, 6]
vp , cbar = TT.showResult(# color map settings, Cmin = Vp minimum of the colorbar, Cmax = Vp maximum value of the colorbar
                      cMap='ocean',cMin=500,cMax=5000,
                      sens=True, contour=True,
                      # colorbar settings
                      logScale = False, orientation='horizontal',
                      # show ray paths
                      rays = False)
vp.set_xlim(-1, 95)
vp.set_ylim(2735, 2775 )
# Sensors and shot (use "diam" to increase/decrease the size of sensors/shot)
pg.viewer.mpl.drawSensors(vp, ttData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
plt.savefig('vp.png', dpi=600)
# plot the quality of inversion
plt.rcParams['figure.figsize'] = [8, 6]
fig = TT.showFit(cmap = 'jet' ) 
plt.savefig('vp_fit.png', dpi=600) 

# %% setup scci class and give manager
# NB a/b must be higher than 0 but lower than 1

scci = SCCI([ERT, TT], names=["ERT", "TT"])
scci.a = 0.05
scci.b = 0.05
scci.c = 2
scci.cmin = scci.b
# scci.cmax = 3.0
cw = scci.singleCWeights()
print(min(cw[0]), min(cw[1]), np.mean(cw[0]), np.mean(cw[1]))
# ax, cb = ERT.showResult(**rKW)
# drawCWeight(ax, ERT.paraDomain, cw[0])
# plt.savefig("ERT.png")
# ax, cb = TT.showResult(**vKW)
# drawCWeight(ax, TT.paraDomain, cw[1])
# plt.savefig("TT.png")
scci.runCoupled(maxIter=20)  # save=True)
print(min(ERT.inv.inv.cWeight()), min(TT.inv.inv.cWeight()))

# %% SCC Joint Inversion
ERT.inv.model = ERT.inv.inv.model()
TT.inv.model = 1./TT.inv.inv.model()

# %% Plot the new sections obtained with SCCI
plt.rcParams['figure.figsize'] = [8, 6]
ax, cb = ERT.showResult(cMap='RdYlBu', cMin = 5000, cMax = 200000)
# Plot electrodes position (use "diam" to increase/decrease the size of the electrodes in the figure)
pg.viewer.mpl.drawSensors(ax, ertData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
ax.set_xlim(-1, 95)
ax.set_ylim(2735, 2775)
plt.savefig('res_scci.png', dpi=600)
# Plot the quality of the obtained model
plt.rcParams['figure.figsize'] = [8, 6]
fig = ERT.showFit(cMap='jet')
plt.savefig('ert_scci_fit.png', dpi=600)

# drawCWeight(ax, ERT.paraDomain, ERT.inv.inv.cWeight())
# plt.savefig("ERT_coupled.png")
plt.rcParams['figure.figsize'] = [8, 6]
ax, cb = TT.showResult(# color map settings, Cmin = Vp minimum of the colorbar, Cmax = Vp maximum value of the colorbar
                      cMap='ocean',cMin=500,cMax=5000,
                      sens=True, contour=True,
                      # colorbar settings
                      logScale = False, orientation='horizontal',
                      # show ray paths
                      rays = False)
ax.set_xlim(-1, 95)
ax.set_ylim(2735, 2775 )
# Sensors and shot (use "diam" to increase/decrease the size of sensors/shot)
pg.viewer.mpl.drawSensors(ax, ttData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
plt.savefig('vp_scci.png', dpi=600)
# plot the quality of inversion
plt.rcParams['figure.figsize'] = [8, 6]
fig = TT.showFit()
plt.savefig('vp_scci_fit.png', dpi=600)  

# drawCWeight(ax, TT.paraDomain, TT.inv.inv.cWeight())
# plt.savefig("TT_coupled.png")
# ax, cbc =  TT.showResult(**vKW)
# TT.drawRayPaths(ax)
