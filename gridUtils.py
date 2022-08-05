import numpy as np
import os
from basicFunc import *

#class: grid varaible
class GridVar:
    def __init__(self, nameVar, minVar, maxVar, stepVar, logVar = False):
        self.name = nameVar
        self.min = minVar
        self.max = maxVar
        self.step = stepVar
        self.logFlag = logVar
        

#class 2D-grid       
class Grid2D:
    def __init__(self, GridVar1, GridVar2):
        self.GridVar1 = GridVar1
        self.GridVar2 = GridVar2
        
    def getGrid2D(self):
        var1list = np.arange(self.GridVar1.min, self.GridVar1.max+self.GridVar1.step, self.GridVar1.step)
        if self.GridVar1.logFlag == True:
            var1list = 10**var1list

        var2list = np.arange(self.GridVar2.min, self.GridVar2.max+self.GridVar2.step, self.GridVar2.step)
        if self.GridVar2.logFlag == True:
            var2list = 10**var2list

        gridShape = var1list.shape[0]*var2list.shape[0]
        grid = np.zeros([gridShape,2])
        xx, yy = np.meshgrid(var1list, var2list)
        grid[:,0] = xx.reshape(gridShape)
        grid[:,1] = yy.reshape(gridShape)
        return grid
    
    def dimGrid(self):
        return (int((self.GridVar1.max+self.GridVar1.step-self.GridVar1.min)/self.GridVar1.step)*int((self.GridVar2.max+self.GridVar2.step-self.GridVar2.min)/self.GridVar2.step),2)
    
#class dictionary for Vars
class inputVars:
    Vars = {}
    
    def __init__(self, nameFixedVars, valueFixedVars):
        for i in enumerate(nameFixedVars):
            self.Vars[i[1]] = valueFixedVars[i[0]]
            
#run 2D grid
def runGrid2D(masterDir, Workdir, GridName, GridVar1, GridVar2, FixedVars, UV_field_on = False, fix_TJ_Tmax = False):
    #initialize log
    lenStr = 20
    modellog = open(GridName+'-models.txt','w', buffering = 1)
    modellog.write(GridName+'\n')
    modellog.write(50*'='+'\n')
    modellog.write('GridVar1:'+GridVar1.name+(lenStr-len(GridVar1.name))*' '+str(GridVar1.min)+(lenStr-len(str(GridVar1.min)))*' '
                   +'----'+str(GridVar1.max)+(lenStr-len(str(GridVar1.max)))*' '+str(GridVar1.step)+(lenStr-len(str(GridVar1.step)))*' '+'log:'+str(GridVar1.logFlag)+'\n')
    modellog.write('GridVar2:'+GridVar2.name+(lenStr-len(GridVar2.name))*' '+str(GridVar2.min)+(lenStr-len(str(GridVar2.min)))*' '
                   +'----'+str(GridVar2.max)+(lenStr-len(str(GridVar2.max)))*' '+str(GridVar2.step)+(lenStr-len(str(GridVar2.step)))*' '+'log:'+str(GridVar2.logFlag)+'\n')
    modellog.write(50*'='+'\n')

    lineFixedName = 'FixedVars:'+(lenStr-len('FixedVars'))*' '

    for key in list(FixedVars.Vars):
        lineFixedName += (key+': '+str(FixedVars.Vars[key])+(lenStr-len(str(FixedVars.Vars[key])))*' ')
    lineFixedName+='\n'
    
    modellog.write(lineFixedName)
    modellog.write(50*'='+'\n')
    
    modellog.write('model'+(lenStr-len('model'))*' '+GridVar1.name+(lenStr-len(GridVar1.name))*' '+GridVar2.name+'\n')

    #initialize vars
    Params = FixedVars.Vars    
    
    Params[GridVar1.name] = GridVar1.min
    Params[GridVar2.name] = GridVar2.min
    
    Grid = Grid2D(GridVar1, GridVar2)
    iterGrid = Grid.getGrid2D()
    
    
    for i, comb in enumerate(iterGrid):
        Params[GridVar1.name] = comb[0]
        Params[GridVar2.name] = comb[1]
        if fix_TJ_Tmax == True:
            Params['duration_max'] = 2*Params['timeJ']
        modelName = GridName + '-' + format(i, '04d' )
        staticParams, runParams = setModel(modelName, Params, UV_field_on = UV_field_on)
        modellog.write(modelName+(lenStr-len(modelName))*' '+format(comb[0],'.2E')+(lenStr-len(format(comb[0],'.2E')))*' '+format(comb[1],'.2E')+'\n')
        runModel(staticParams, runParams, masterDir, Workdir)
    modellog.close()
    print('Complete! See model logs in '+ GridName+'-models.txt')
    
