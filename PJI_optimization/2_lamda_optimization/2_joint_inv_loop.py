# SCRIPT FOr THE PETROPHISICAL JOINT INVERSION
# HERE YOU CAN CHOOSE DIFFERENT VALUES FOR THE PARAMETRERS AND PERFORM THE INVERSIONS WITH A LOOP

import numpy as np
import pygimli as pg
from fpinv import FourPhaseModel, JointInv, JointMod
from pygimli.physics import Refraction, ERTManager
from pygimli.physics.traveltime import createGradientModel2D
from settings import *

# m interval 
mi= np.linspace(1.9, 1.9, 1)
# n interval
ni = np.linspace(2.4, 2.4, 1)
# vr interval
vri = np.linspace(5000, 5000, 1)
# rhow interval
rhowi = np.linspace(100,100,1)

nv = np.linspace(1,10,10)
li = 5
s = 0

for i in range(len(nv)):
    li = li*2   
    for j in range(len(mi)):
        ml = mi[j]
        for k in range(len(ni)):
            nl = ni[k]
            for w in range(len(vri)):
                vrl = vri[w]
                for h in range(len(rhowi)):
                    rhowl = rhowi[h]
                    
                    fpm = FourPhaseModel(phi=poro, va=300., vi=3500., vw=1500, m=ml, n=nl,
                                         rhow=rhowl, vr=vrl)

                    args = sys.argv
                    lam = li
                    case = 1
                    weighting = False

                    if case == 2:
                        case = 2
                        constrained = True
                        mesh = pg.load("mesh_2.bms")
                        paraDomain = pg.load("paraDomain_2.bms")
                    else:
                        case = 1
                        constrained = False
                        mesh = pg.load("mesh_1.bms")
                        paraDomain = pg.load("paraDomain_1.bms")

                    pg.boxprint("Calculating case %s" % case)

                    # Load meshes and data
                    ertScheme = pg.DataContainerERT("ert_filtered.ohm")

                    fr_min = 0.1
                    fr_max = 0.9
                    phi = np.ones(paraDomain.cellCount()) * poro

                    # Setup managers and equip with meshes
                    ert = ERTManager()
                    ert.setMesh(mesh)
                    ert.setData(ertScheme)
                    ert.fop.createRefinedForwardMesh()

                    ttData = pg.DataContainer("rst_filtered.sgt", "s g")
                    rst = Refraction()
                    rst.setMesh(paraDomain)
                    rst.setData(ttData)
                    rst.fop.createRefinedForwardMesh()

                    # Set errors
                    ttData.set("err", np.ones(ttData.size()) * rste)
                    ertScheme.set("err", np.ones(ertScheme.size()) * erte)

                    if constrained:
                        # Find cells around boreholes to fix ice content to zero
                        fixcells = []
                        for cell in paraDomain.cells():
                            x, y, _ = cell.center()
                            if (x > 9) and (x < 11) and (y > -depth_5198):
                                fixcells.append(cell.id())
                            elif (x > 25) and (x < 27) and (y > -depth_5000):
                                fixcells.append(cell.id())
                        fixcells = np.array(fixcells)
                    else:
                        # Do not fix ice
                        fixcells = False

                    # Setup joint modeling and inverse operators
                    JM = JointMod(paraDomain, ert, rst, fpm, fix_poro=False, zWeight=zWeight,
                                  fix_ice=fixcells)

                    data = pg.cat(ttData("t"), ertScheme("rhoa"))

                    if weighting:
                        n_rst = ttData.size()
                        n_ert = ertScheme.size()
                        avg = (n_rst + n_ert) / 2
                        weight_rst = avg / n_rst
                        weight_ert = avg / n_ert
                    else:
                        weight_rst = 1
                        weight_ert = 1

                    error = pg.cat(
                        ttData("err") / ttData("t") / weight_rst,
                        ertScheme("err") / weight_ert)

                    minvel = 1000
                    maxvel = 5000
                    velstart = 1 / createGradientModel2D(ttData, paraDomain, minvel, maxvel)
                    rhostart = np.ones_like(velstart) * np.mean(ertScheme("rhoa"))
                    fas, fis, fws, _ = fpm.all(rhostart, velstart)
                    frs = np.ones_like(fas) - fpm.phi
                    frs[frs <= fr_min] = fr_min + 0.01
                    frs[frs >= fr_max] = fr_max - 0.01
                    if fixcells is not False:
                        fis[fixcells] = 0.0
                    startmodel = np.concatenate((fws, fis, fas, frs))

                    # Fix small values to avoid problems in first iteration
                    startmodel[startmodel <= 0.001] = 0.001

                    inv = JointInv(JM, data, error, startmodel, frmin=fr_min, frmax=fr_max,
                                   lam=lam, maxIter=maxIter)

                    # Run inversion
                    model = inv.run()
                    print(("Chi squared fit:", inv.getChi2()))
                    print("RST rms:", inv.relrms())
                    RMS = repr(inv.relrms())
                    CHI2 = repr(inv.getChi2())
    
                    # Save results
                    fwe, fie, fae, fre = JM.fractions(model)
                    fsum = fwe + fie + fae + fre

                    print("Min/Max sum:", min(fsum), max(fsum))

                    rhoest = JM.fpm.rho(fwe, fie, fae, fre)
                    velest = 1. / JM.fpm.slowness(fwe, fie, fae, fre)

                    array_mask = np.array(((fae < 0) | (fae > 1 - fre))
                                          | ((fie < 0) | (fie > 1 - fre))
                                          | ((fwe < 0) | (fwe > 1 - fre))
                                          | ((fre < 0) | (fre > 1))
                                          | (fsum > 1.01))
                    
                
                    np.savez("joint_inversion_%s.npz" % s, vel=np.array(velest),
                             rho=np.array(rhoest), fa=fae, fi=fie, fw=fwe,fr=fre, lambdaa = lam, vrock = vrl, wateres = rhowl, mexp = ml, nexp = nl, mask=array_mask)
                    
                    ertchi = JM.ERTchi2(model, error)
                    rstchi = JM.RSTchi2(model, error, ttData("t"))
                    rmsrho = ertchi[1]
                    rmsrhop = repr(rmsrho)
                    chi2rho = ertchi[0]
                    chi2rhop = repr(chi2rho)
                    rmsrst = rstchi[1]
                    rmsrstp = repr(rmsrst)
                    chi2rst = rstchi[0]
                    chi2rstp = repr(chi2rst)
                    
                    file = open("parameters_%s.txt" % s, "w")
                    lambdaa = repr(lam)
                    vrock = repr(vrl)
                    wateres = repr(rhowl)
                    mexp = repr(ml)
                    nexp = repr(nl)
                    file.write("\n" + "lamda = " + lambdaa + "\n" + "rhow = " + wateres + "\n" + "vr = " + vrock + "\n" + "m = " + mexp + "\n" + "n = " + nexp + "\n" + "rms = " + RMS + "\n" + "chi2 = " + CHI2 + "\n"  + "ertRMS = " + rmsrhop + "\n" + "ertCHI2 = " + chi2rhop + "\n"  + "rstRMS = " + rmsrstp + "\n" + "rstCHI2 = " + chi2rstp)
                    file.close()

                    s = s+1