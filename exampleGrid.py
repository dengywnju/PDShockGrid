from basicFunc import *
from gridUtils import *

#define directories
masterDir = '/Users/dengyw/Workspace/Shock_1.1_rev122/'
workDir = '/Users/dengyw/Workspace/Shock_1.1_rev122/PDShockGrid'

#example Grid: Non-radiative, non-steady C-type shock:
#Vs from 20 to 40 km/s with steps of 10 km/s, 
#density from 10^3 to 10^5 /cc with steps of 0.5 dex
vel = GridVar('Vs_km', 20, 40, 10)
dens = GridVar('nH_init', 3, 5, 1, True)
fixed_variables = inputVars(nameFixedVars = ['shock_type', 'Bbeta', 'Zeta', 'timeJ', 'duration_max'], valueFixedVars = ['C', 3, 3.7E-15, 5000, 10000])

#run the Grid
runGrid2D(masterDir, workDir, 'TestGrid', vel,dens,fixed_variables)
