from scipy import array
import numpy as np
import os
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
from PartFunc import *


uVel = u.cm/u.s
uEmi = u.erg/u.s/u.cm**3
uLen = u.cm
uT = u.K

def get_runList(modellog):
    with open(modellog,'r') as log:
        argue0 = log.readline().split()[0]
        while argue0 != 'model':
            line = log.readline().split()
            argue0 = line[0]
            
        nameGridVar1 = line[1]
        nameGridVar2 = line[2]
        lines = log.readlines()
        models = {}
        models['model'] = np.full(len(lines),40*' ')
        models[nameGridVar1] = np.zeros(len(lines))
        models[nameGridVar2] = np.zeros(len(lines))
        models['Vars'] = [nameGridVar1,nameGridVar2]
        for (line_count, line) in enumerate(lines):
            items = line.split()
            for (col_count, col_name) in enumerate(['model',nameGridVar1,nameGridVar2]):
                value = items[col_count]
                models[col_name][line_count] = value
    return models


def plot_H2_emissivity(model, lines = ['0-0S(1)','1-0S(1)','2-1S(1)'], logX = False, logY = False, save_PDF = False, formatVar1 = '.0f', formatVar2 = '.0f'):
    numlines = len(lines)

    fig = plt.figure(figsize = (6.6,3))
    
    rows = math.ceil(numlines/3)
    cols = math.ceil(numlines/rows)
    for i,line in enumerate(lines):
        ax = fig.add_subplot(rows,cols,i+1)
        ax.set_title(line[:3]+' '+line[3:], weight = 'bold')
        for j,run in enumerate(model['model']):
            xmin, xmax = 1e-12,1e-12
            resultDir = masterDir+'output/'+run
            os.chdir(resultDir)
            H2_lines = get_H2line()
            label = model['Vars'][0] + '=' + format(model[model['Vars'][0]][j], formatVar1) + ', '+ model['Vars'][1] + '=' + format(model[model['Vars'][1]][j],formatVar2)
            xaxis = ((H2_lines['distance']*u.cm).to(u.pc)).value
            ax.plot(xaxis, H2_lines[lines[i]], label = label, linewidth = 1)
            if np.max(xaxis) > xmax:
                xmax = np.max(xaxis)
        
        if i == len(lines) - 1:
            ax.legend(framealpha = 0)
        ax.set_xlabel('z [pc]') 
        ax.set_xlim(xmin,xmax)
    
        if logX == True:
            ax.set_xscale('log')
        
        if logY == True:
            ax.set_yscale('log')
        ax.set_ylabel(r'Emissivity [erg s$^{-1}$ cm$^{-3}$]') 
    if save_PDF == True: 
        plt.savefig('H2_emissivity.pdf',bbox_inches = 'tight')

def checkValue(value):
    # Check if value should be a float or flagged as missing
    if value == "---":
        value = np.ma.masked
    else:
        value = float(value)
    return value



    
def get_H2line():
    fH2line = open('H2_line.out','r')
    H2line_names=fH2line.readline().split()
    H2line_lines = fH2line.readlines()
    H2_data = {}
    for H2line_name in H2line_names:
        H2_data[H2line_name] = np.ma.zeros(len(H2line_lines), 'f', fill_value = -999.999)
    for (line_count, line) in enumerate(H2line_lines):
        items = line.split()
        for (col_count, col_name) in enumerate(H2line_names):
            value = items[col_count]
            H2_data[col_name][line_count] = checkValue(value)
    return H2_data

def get_excit():
    fexcit = open('excit.out','r')
    excit_names=fexcit.readline().split()
    excit_lines = fexcit.readlines()
    excit_data = {}
    for excit_name in excit_names:
        excit_data[excit_name] = np.ma.zeros(len(excit_lines), 'f', fill_value = -999.999)
    for (line_count, line) in enumerate(excit_lines):
        items = line.split()
        for (col_count, col_name) in enumerate(excit_names):
            value = items[col_count]
            excit_data[col_name][line_count] = checkValue(value)
    return excit_data

def get_mhd_phys():
    fmhd_phys = open('mhd_phys.out','r')
    mhd_phys_names = fmhd_phys.readline().split() 
    mhd_phys_lines = fmhd_phys.readlines()
    mhd_phy_data = {}
    for mhd_phys_name in mhd_phys_names:
        mhd_phy_data[mhd_phys_name] = np.ma.zeros(len(mhd_phys_lines), 'f', fill_value = -999.999)
    for (line_count, line) in enumerate(mhd_phys_lines):
        items = line.split()
        for (col_count, col_name) in enumerate(mhd_phys_names):
            value = items[col_count]
            mhd_phy_data[col_name][line_count] = checkValue(value)
    return mhd_phy_data

def get_mhd_speci():
    fmhd_speci = open('mhd_speci.out','r')
    dum=fmhd_speci.readline()
    mhd_speci_names = fmhd_speci.readline().split() 
    mhd_speci_lines = fmhd_speci.readlines()
    mhd_speci_data = {}
    for mhd_speci_name in mhd_speci_names:
        mhd_speci_data[mhd_speci_name] = np.ma.zeros(len(mhd_speci_lines), 'f', fill_value = -999.999)
    for (line_count, line) in enumerate(mhd_speci_lines):
        items = line.split()
        for (col_count, col_name) in enumerate(mhd_speci_names):
            value = items[col_count]
            mhd_speci_data[col_name][line_count] = checkValue(value)
    return mhd_speci_data

def get_AtomicIntensity():
    fintensity = open('intensity.out', 'r')
    intensity_names=fintensity.readline().split()
    intensity_lines = fintensity.readlines()
    intensity_data = {}
    for intensity_name in intensity_names:
        intensity_data[intensity_name] = np.ma.zeros(len(intensity_lines), 'f', fill_value = -999.999)
    for (line_count, line) in enumerate(intensity_lines):
        items = line.split()
        for (col_count, col_name) in enumerate(intensity_names):
            value = items[col_count]
            intensity_data[col_name][line_count] = checkValue(value)
    return intensity_data

def get_Intensity(H2_line, mhd_phys, Vs, i=0):
    cosi = np.cos(i)
    E = (H2_line[1:]+H2_line[:-1])*uEmi/2
    #read
    Vn = (mhd_phys['Vn'][1:]+mhd_phys['Vn'][:-1])*uVel/2
    z = mhd_phys['distance']*uLen
    Tn = (mhd_phys['Tn'][1:]+mhd_phys['Tn'][:-1])*uT/2
    X_H2 = (mhd_phys['x(H2)'][1:]+mhd_phys['x(H2)'][:-1])/2
    X_H = (mhd_phys['x(H)'][1:]+mhd_phys['x(H)'][:-1])/2
    X_Hii = (mhd_phys['x(H+)'][1:]+mhd_phys['x(H+)'][:-1])/2
    wavel = 2.1218*u.um
    dz = z[1:]-z[:-1]
    mu = 1/(X_H2+X_H+X_Hii)
    #thermal dispersion
    #sigma = np.sqrt(2*c.k_B*Tn/(mu)/c.m_p+(3.5*u.km/u.s)**2).to(uVel)
    sigma = np.sqrt(2*c.k_B*Tn/(mu)/c.m_p).to(uVel)
    #initialise
    Vmax = 2*np.max(Vn)
    V = np.arange(-Vmax.value,Vmax.value,10000)*uVel
    Iv = np.zeros(V.shape)
    #integration
    for i in range(len(V)):
        Iv[i] = (np.sum(E/(np.sqrt(2*np.pi)*sigma)*np.exp(-((V[i]-(Vn-Vs)*cosi))**2/2/sigma**2)*dz).to(u.erg/u.s/u.cm**2/u.km*u.s)).value
    return V,Iv

class singleModel():
    def __init__(self, outDir, Var1 = None, Var2 = None, Var1value = None, Var2value = None, name = 'dirName'):
        self.outDir = outDir
        self.Var1 = Var1
        self.Var2 = Var2
        self.Var1value = Var1value
        self.Var2value = Var2value
        if name != 'sdirName':
            self.name = name
        else:
            self.name = os.getcwd().split('/')[-1]
    
    def H2_line(self):
        os.chdir(self.outDir)
        return get_H2line()
    
    def H2_population(self):
        os.chdir(self.outDir)
        return get_excit()
    
    def mhd_phys(self):
        os.chdir(self.outDir)
        return get_mhd_phys()
    
    def mhd_speci(self):
        os.chdir(self.outDir)
        return get_mhd_speci()
    
    def AtomicIntensity(self):
        os.chdir(self.outDir)
        return get_AtomicIntensity
    
    def H2_intensity_velocity(self, Vs, i):
        os.chdir(self.outDir)
        return get_Intensity(H2_line, mhd_phys, Vs, i=0)
    
    def plot_H2_emissivity_vs_Tn(self, lines = ['0-0S(1)','1-0S(1)','2-1S(1)'], magFac = [1,1,1], xmin = None, xmax = None ,ymin = None,ymax = None,logX = False, logY = False, Normalize = False ,save_PDF = False):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        H2_lines = self.H2_line()
        if self.Var1:
            title = self.Var1+'$=$'+format(self.Var1value,'.0f')+', '+self.Var2+'$=$'+format(self.Var2value,'.0f')
            ax1.set_title(title, weight = 'bold')
            
        for i,line in enumerate(lines):
            
            label = line[:3]+' '+line[3:]
            if magFac[i] != 1:
                label+=r'$\times$'+str(magFac[i])
            xaxis = ((H2_lines['distance']*u.cm).to(u.pc)).value
            if Normalize == True:
                ax1.plot(xaxis, H2_lines[line]/(max(H2_lines[line])), label = label, linewidth = 1)
            else:
                ax1.plot(xaxis, H2_lines[line]*magFac[i], label = label, linewidth = 1)
        ax1.plot([],[],color='grey', linestyle = '--', linewidth = 1,label = 'Temprature')
        ax1.legend(framealpha = 0)
        ax1.set_xlabel('z [pc]') 
        mhd_phys = self.mhd_phys()
        ax2 = ax1.twinx()
        ax2.plot((mhd_phys['distance']*u.cm).to(u.pc), mhd_phys['Tn'],color='grey', linestyle = '--', linewidth = 1)
        ax2.set_ylabel(r'Temperature [K]')
        if xmin:
            ax1.set_xlim(xmin,xmax)
            if logX == False:
                ax1.set_xlim(min(xaxis),max(xaxis))

            if logX == True:
                ax1.set_xlim(1e-12,max(xmax))
        if logX == True:        
            ax1.set_xscale('log')
            
        if ymin:
            ax1.set_ylim(ymin,ymax)

        if logY == True:
            ax1.set_yscale('log')
        if Normalize == True:
            ax1.set_ylabel(r'Normalized Emissivity') 
        else:
            ax1.set_ylabel(r'Emissivity [erg s$^{-1}$ cm$^{-3}$]') 
            
        
                         
        if save_PDF == True: 
            plt.savefig('H2_emissivity.pdf',bbox_inches = 'tight')
            
    def plot_population_diagram(self, v_levels = [0,1,2], J_levels = [], log10 = True, Fit = False):
        if Fit:
            def f_linear(x, a, b):
                    return a * x + b
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        H2_excit = self.H2_population()
        for i, level in enumerate(v_levels): 
            if J_levels == []:
                E = H2_excit['Energy(K)'][H2_excit['V']==level]
                Pop = H2_excit['log(N/g)'][H2_excit['V']==level]
            else:
                E = []
                Pop = []
                for Jlevel in J_levels[i]:
                    EJ = H2_excit['Energy(K)'][(H2_excit['V']==level)&(H2_excit['J']==Jlevel)][0]
                    PopJ = H2_excit['log(N/g)'][(H2_excit['V']==level)&(H2_excit['J']==Jlevel)][0]
                    E.append(EJ)
                    Pop.append(PopJ)
                E = np.array(E)
                Pop = np.array(Pop)
            
            if Fit:
                print(E,Pop)
                line_fir = curve_fit(f_linear, E, Pop)
                a = line_fir[0][0]
                b = line_fir[0][1]
                T = -1/a
                N = np.log10(np.exp(b + np.log(PartitionFunction(T))))
                
            
            if log10:
                Pop = Pop/np.log(10)
                if Fit:
                    a = a/np.log(10)
                    b = b/np.log(10)
                
            ax1.scatter(E, Pop, s=3, label='v='+str(int(level)))
            
            if Fit:
                ax1.plot(E, a*E+b, linewidth = 1, linestyle = '--', color ='grey', label = format(T, '.0f')+' K, '+format(N, '.2f')+' cm$^{-2}$')
                
        if log10:
            ax1.set_ylabel('log(N/g)')
        else:
            ax1.set_ylabel('ln(N/g)')
        plt.xlabel('Energy [K]')
        plt.legend(framealpha = 0) 
        


