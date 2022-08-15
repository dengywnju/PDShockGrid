# PDShockGrid
Grid scripts for Paris-Durham shock code
============
This is a python-based grid tool for running the Paris-Durham shock code
in a high-dimensional parameter space.

Author: Yunwei Deng

Requirement
===========
python3
numpy
scipy
astropy
matplotlib

Quick start
===========
0. $vim exampleGrid.py
1. Set the grid variables
2. Set the fixed parameters
3. $:wq
4. $python exampleGrid.py
6. See AnalyseExample.ipynb for more examples

To run with personlized parameters, eidt the input values of 

Func: runGrid2D(masterDir, Workdir, GridName, GridVar1, GridVar2, FixedVars, UV_field_on = False, fix_TJ_Tmax = False)

masterDir: where you install Shock_1.1_rev122

Workdir: your workdir

GridName: name your grid

GridVar1 & 2: should be two 'GridVar' objects, including name, upper+lower limits, sizes of steps, and scales of the grid varibles 

FixedVars: a inputVars object

UV_field_on: turn on UV field, must mannually set the UV parameters as FixedVars

fix_TJ_Tmax: make young shocks, if true, max_duration = 2 * timeJ 
