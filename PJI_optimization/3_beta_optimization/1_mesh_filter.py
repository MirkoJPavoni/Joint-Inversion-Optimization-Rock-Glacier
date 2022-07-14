# SCRIPT TO CREATE THE MESH 

import numpy as np
import pygimli as pg
pg.verbose = print # temporary
import pygimli.meshtools as mt
from pygimli.physics.ert import ERTManager, createGeometricFactors

# Load the ERT input file (must be in the same folder of the code)
ertData = pg.DataContainerERT("ert_filtered.ohm")
print(ertData)
# We can create 'k' based on the given topography
ertData['k'] = createGeometricFactors(ertData, numerical=True)
# Initialize the ERTManager for further steps and inversion
ert = ERTManager(sr=False, useBert=True, verbose=True, debug=False)
# The data container has no apparent resistivities (no ‘rhoa’ field data is present)
# We can let the Manager fix and add all the columns (err, i, ip, iperr, k, r, rhoa, u, valid) to perform the inversion:
ert.checkData(ertData)
#  Define the error inversion (e.g here 10%)
ertData['err'] = ert.estimateError(ertData, absoluteError=0.001, relativeError=0.10)
# Remove the negative values off apparent resistivity
ertData.remove(ertData["rhoa"] <= 0)
print(ertData)

# Load the seismic input file (must be in the same folder of the code)
rstData = pg.DataContainer("rst_filtered.sgt", "s g")
# Define the error in the seismic dataset in sec (e.g. 0.0003 means 0.3 msec)
rstData.set("err", pg.Vector(rstData.size(), 0.0006))
# Remove the shoots in the same position of the geophones
rstData.remove(rstData["g"] == rstData["s"])
print(rstData)

# %% run this module only if the investigation lines don't have the same lenght 
# In Schafberg L2 the seismic line is 96 m and also the ert line, so you don't need to run this module
# See Schafberg L1 for different investigation lines lenght
# idx = []
# for i, sensor in enumerate(ertData.sensors()):
#    if sensor[0] >= 96.0:
#       idx.append(i)        
# ertData.removeSensorIdx(idx)

# %% run always this module (even if the investigation lines have the same lenght)
ertData.removeInvalid()
ertData.removeUnusedSensors()
ertData.save("ert_filtered.ohm")
rstData.save("rst_filtered.sgt")
ertData = pg.load("ert_filtered.ohm")
rstData = pg.load("rst_filtered.sgt")

# %% create a common file the position of the electrodes and geophones/shoots
ertData = pg.DataContainerERT("ert_filtered.ohm")
print(ertData)

print(len(ertData.sensorPositions()))
for pos in ertData.sensorPositions():
    print(pos)

def is_close(pos, data, tolerance=0.1):
    for posi in data.sensorPositions():
        dist = pos.dist(posi)
        if dist <= tolerance:
            return True
    return False

combinedSensors = pg.DataContainer()
for pos in ertData.sensorPositions():
    combinedSensors.createSensor(pos)

for pos in rstData.sensorPositions():
    if is_close(pos, ertData):
        print("Not adding", pos)
    else:
        combinedSensors.createSensor(pos)

combinedSensors.sortSensorsX()
x = pg.x(combinedSensors.sensorPositions()).array()
y = pg.y(combinedSensors.sensorPositions()).array()

np.savetxt("sensors.npy", np.column_stack((x, y)))

print("Number of combined positions:", combinedSensors.sensorCount())

# %%  Create a common mesh for the ERT and RST inversions
# combinedSensor = nodes for the position of the electrodes and geophones/shoots
# paraDepth: maximum depth for parametric domain, 0 (default) means 0.4 * maximum sensor range
# paraDX: relative distance for refinement nodes between two sensors (1=none), e.g., 0.5 means 1 additional node between two neighboring sensors e.g., 0.33 means 2 additional equidistant nodes between two sensors
# paraMaxCellSize = maximum size for parametric cells (size in m*m)
# paraBoundary: margin/offset for parameter domain in absolute sensor distances
# boundary: boundary width to be appended for domain prolongation in absolute para domain width. Values <=0 force the boundary to be 4 times para domain width.

plc = mt.createParaMeshPLC(combinedSensors, paraDX=0.3, boundary=4,
                               paraDepth=20, paraBoundary=3,
                               paraMaxCellSize=10)

mesh = mt.createMesh(plc, quality=33.8)
mesh.save("mesh_1.bms")
# Extract the inner domain where parameters should be estimated. 
# Outer domain is only needed for ERT forward simulation, not for seismic traveltime calculations.
paraDomain = pg.Mesh(2)
paraDomain.createMeshByMarker(mesh, 2)
paraDomain.save("paraDomain_1.bms")