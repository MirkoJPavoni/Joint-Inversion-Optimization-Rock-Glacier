# SCRIPT TO PERFORM CONVENTIONAL INDEPENDENT INVERSION AND SCC JOINT INVERSION
# WITH A LOOP TO SET DIFFERENT VALUES FOR THE PARAMETERS a,b, AND c

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

# %% define mesh and dataset
rMesh = pg.load("mesh_1.bms")
vMesh = pg.load("paraDomain_1.bms")
rKW = dict(logScale=True, cMin=2000, cMax=200000, cMap="Spectral_r")
vKW = dict(logScale=True, cMin=400, cMax=4000, cMap="Spectral_r")

ertData = ert.load("ert_filtered.ohm")
ERT = ert.ERTManager(ertData, verbose=True, sr=False)
ERT.setMesh(rMesh)

ttData = tt.load("rst_filtered.sgt")
TT = tt.TravelTimeManager(ttData, verbose=True)
TT.errIsAbsolute = True
TT.setMesh(vMesh)

# %%ert inversion
ERT.invert(zWeight=0.25, lam=100, cType=1)

# %% plot resistivity section
plt.rcParams['figure.figsize'] = [20, 6]
ax , cbar = ERT.showResult(cMap='jet_r', cMin = 2000, cMax = 200000)
# Plot electrodes position (use "diam" to increase/decrease the size of the electrodes in the figure)
pg.viewer.mpl.drawSensors(ax, ertData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
ax.set_xlim(-1, 95)
ax.set_ylim(2730, 2775)

# Plot the quality of the obtained model
plt.rcParams['figure.figsize'] = [20, 6]
fig = ERT.showFit(cMap='jet')

# %% rst inversion
TT.invert(zWeight=0.25, lam=250, vTop=500, vBottom=5000,
           verbose=True)

# %% plot vp section
plt.rcParams['figure.figsize'] = [20, 6]
vp , cbar = TT.showResult(# color map settings, Cmin = Vp minimum of the colorbar, Cmax = Vp maximum value of the colorbar
                      cMap='jet_r',cMin=350,cMax=3500,
                      sens=True, contour=True,
                      # colorbar settings
                      logScale = False, orientation='horizontal',
                      # show ray paths
                      rays = False)
vp.set_xlim(-1, 95)
vp.set_ylim(2730, 2775 )
# Sensors and shot (use "diam" to increase/decrease the size of sensors/shot)
pg.viewer.mpl.drawSensors(vp, ttData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
# plot the quality of inversion
plt.rcParams['figure.figsize'] = [12, 6]
fig = TT.showFit() 

# %% define the range for parameters a,b and c
# NB a/b must be higher than 0 and lower than 1

scci = SCCI([ERT, TT], names=["ERT", "TT"])
# a interval
ai = np.linspace(0.05,0.4, 4)
# b interval 
bi= np.linspace(0.05, 0.4, 4)
# c interval
ci = np.linspace(0.5, 2, 4)

# loop counter
s = 0
for i in range(len(ai)):
    scci.a = ai[i]    
    for j in range(len(bi)):
        scci.b = bi[j]
        for k in range(len(ci)):
            scci.c = ci[k]
            cw = scci.singleCWeights()
            print(min(cw[0]), min(cw[1]), np.mean(cw[0]), np.mean(cw[1]))
            # ax, cb = ERT.showResult(**rKW)
            # drawCWeight(ax, ERT.paraDomain, cw[0])
            # plt.savefig("ERT.png")
            # ax, cb = TT.showResult(**vKW)
            # drawCWeight(ax, TT.paraDomain, cw[1])
            # plt.savefig("TT.png")
            
            scci.runCoupled(maxIter=10)  # save=True)
            print(min(ERT.inv.inv.cWeight()), min(TT.inv.inv.cWeight()))
            
            ERT.inv.model = ERT.inv.inv.model()
            TT.inv.model = 1./TT.inv.inv.model()
            
            plt.rcParams['figure.figsize'] = [20, 6]
            ax, cb = ERT.showResult(cMap='jet_r', cMin = 1000, cMax = 200000)
            # Plot electrodes position (use "diam" to increase/decrease the size of the electrodes in the figure)
            pg.viewer.mpl.drawSensors(ax, ertData.sensors(), diam=1,
                                     facecolor='white', edgecolor='black')
            ax.set_xlim(-1, 95)
            ax.set_ylim(2730, 2775)
            plt.savefig('resec_%s.png' %s, dpi=600)
            # Plot the quality of the obtained model
            plt.rcParams['figure.figsize'] = [20, 6]
            fig = ERT.showFit(cMap='jet')
            plt.savefig('resfit_%s.png' %s, dpi=600)

            # drawCWeight(ax, ERT.paraDomain, ERT.inv.inv.cWeight())
            # plt.savefig("ERT_coupled.png")
            plt.rcParams['figure.figsize'] = [20, 6]
            ax, cb = TT.showResult(# color map settings, Cmin = Vp minimum of the colorbar, Cmax = Vp maximum value of the colorbar
                                  cMap='jet_r',cMin=300,cMax=3000,
                                  sens=True, contour=True,
                                  # colorbar settings
                                  logScale = False, orientation='horizontal',
                                  # show ray paths
                                  rays = False)
            ax.set_xlim(-1, 95)
            ax.set_ylim(2730, 2775 )
            # Sensors and shot (use "diam" to increase/decrease the size of sensors/shot)
            pg.viewer.mpl.drawSensors(ax, ttData.sensors(), diam=1,
                                     facecolor='white', edgecolor='black')
            plt.savefig('vpsec_%s.png' %s, dpi=600)

            plt.rcParams['figure.figsize'] = [12, 6]
            fig = TT.showFit()
            plt.savefig('vpfit_%s.png' %s, dpi=600)
            
            file = open("parameters_%s.txt" % s, "w")
            a_is = repr(scci.a)
            b_is = repr(scci.b)
            c_is = repr(scci.c)
            file.write("a=" + a_is + "\n" + "b=" + b_is + "\n" + "c=" + c_is + "\n")
            file.close()

            s = s+1

            # drawCWeight(ax, TT.paraDomain, TT.inv.inv.cWeight())
            # plt.savefig("TT_coupled.png")
            # ax, cbc =  TT.showResult(**vKW)
            # TT.drawRayPaths(ax)



