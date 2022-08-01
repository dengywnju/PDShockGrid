import numpy as np
import os
from basicFunc import *
from gridUtils import *


masterDir = '/Users/dengyw/Workspace/Shock_1.1_rev122/'
workDir = '/Users/dengyw/Workspace/Shock_1.1_rev122/PDShockGrid'

vel = GridVar('Vs_km', 20, 40, 10)
dens = GridVar('nH_init', 3, 5, 1, True)
fixed_variables = inputVars(nameFixedVars = ['shock_type', 'Bbeta', 'Zeta', 'timeJ', 'duration_max'], valueFixedVars = ['C', 3, 3.7E-15, 5000, 10000])

runGrid2D(masterDir, workDir, 'TestGrid', vel,dens,fixed_variables)