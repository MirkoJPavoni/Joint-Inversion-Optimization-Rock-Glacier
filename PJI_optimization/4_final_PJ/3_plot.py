# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 08:25:56 2021

@author: Mirko
"""

from string import ascii_uppercase
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid
import pygimli as pg
from fpinv import add_inner_title, logFormat, rst_cov, set_style
from pygimli.viewer.mpl import drawModel
import matplotlib.pyplot as plt


fs = 5.5
set_style(fs, style="seaborn-dark")   
# Load data
mesh1 = pg.load("paraDomain_1.bms")
joint = np.load("joint_inversion_0.npz")
sensors = np.loadtxt("sensors.npy")
def to_sat(fw, fi, fa, fr):
    phi = 1 - fr
    return fw / phi, fi / phi, fa / phi

veljoint1, rhojoint1, faj1, fij1, fwj1, frj1, maskj1 = joint["vel"], joint[
    "rho"], joint["fa"], joint["fi"], joint["fw"], joint["fr"], joint["mask"]

phij1 = 1-frj1
fwj1, fij1, faj1 = to_sat(fwj1, fij1, faj1, frj1)
 
# for the sensors
ertData = pg.DataContainerERT("ert_filtered.ohm")
plt.rcParams['figure.figsize'] = [8,6]

resi, cbar = pg.show(mesh1, data=rhojoint1, cMap='RdYlBu', cMin=5000, cMax=200000)
pg.viewer.mpl.drawSensors(resi, ertData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
resi.set_xlim(-1, 95)
resi.set_ylim(2735, 2775)
plt.savefig('res.png', dpi=600)

veli, cbar = pg.show(mesh1, data=veljoint1, cMap='ocean', cMin=500, cMax=5000)
pg.viewer.mpl.drawSensors(veli, ertData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
veli.set_xlim(-1, 95)
veli.set_ylim(2735, 2775)
plt.savefig('vp.png', dpi=600)

wci, cbar = pg.show(mesh1, data=fwj1, cMap='YlGnBu', cMin=0.06, cMax=0.21)
pg.viewer.mpl.drawSensors(wci, ertData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
wci.set_xlim(-1, 95)
wci.set_ylim(2735, 2775)
print('fw min is', min(fwj1))
print('fw max is', max(fwj1))
plt.savefig('water.png', dpi=600)

aci, cbar = pg.show(mesh1, data=faj1, cMap='Greens', cMin=0.04, cMax=0.81)
pg.viewer.mpl.drawSensors(aci, ertData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
aci.set_xlim(-1, 95)
aci.set_ylim(2735, 2775)
print('fa min is', min(faj1))
print('fa max is', max(faj1))
plt.savefig('air.png', dpi=600)


ici, cbar = pg.show(mesh1, data=fij1, cMap='turbo_r', cMin=0.13, cMax=0.88)
pg.viewer.mpl.drawSensors(ici, ertData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
ici.set_xlim(-1, 95)
ici.set_ylim(2735, 2775)
print('fi min is', min(fij1))
print('fi max is', max(fij1))
plt.savefig('ice.png', dpi=600)

poro, cbar = pg.show(mesh1, data=phij1, cMap='Oranges', cMin=0.47, cMax=0.80)
pg.viewer.mpl.drawSensors(poro, ertData.sensors(), diam=1,
                         facecolor='white', edgecolor='black')
poro.set_xlim(-1, 95)
poro.set_ylim(2735, 2775)
print('phi min is', min(phij1))
print('phi max is', max(phij1))
plt.savefig('phi.png', dpi=600)


