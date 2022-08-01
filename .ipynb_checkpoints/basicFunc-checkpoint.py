#define the basic functions

import numpy as np
import os

#Merge two dictionary
def Merge(dict1, dict2): 
    return(dict2.update(dict1)) 

#create a dictionary to save the parameters
def createDefaultParams():
    inputParams = {}

    #input files
    inputParams['modele'] = ''                                  #output files radix
    inputParams['specfile'] = 'species_2022_05.in_depl'         #species file (spe/cies list / enthalpies / initial abundances)
    inputParams['chemfile'] = 'chemistry_2022_05_all.in_noadso' #chemistry file (reaction list / rates)
    inputParams['h2exfile'] = 'none'                            #h2* file (if none, population initialized depending on op_H2_in)
    inputParams['gridfile'] = 'none'                            #file containing the grid of position - radiation field

    #shock parameters
    inputParams['shock_type'] = 'S1'                            #'C' or 'J', Isoprotonic or isobaric steady state : 'S1' or 'S2', Isoprotonic or isobaric PDR: 'P1' or 'P2', Constant velocity or self-consistent Wind : 'W1' or 'W2'
    inputParams['Nfluids'] = '1'                                #1, 2 ou 3
    inputParams['Bbeta'] = '1.00E-01'                           #Bfield = Bbeta * sqrt(nH) (micro Gauss)
    inputParams['Vs_km'] = '1.00E+01'                           #shock speed (km/s)
    inputParams['DeltaVmin'] = '1.00E+03'                       #initial Vn - Vi (cm s-1)   
    inputParams['nH_init'] = '1.00E+04'                         #initial nH = n(H) + 2.0 n(H2) + n(H+) (cm-3)
    inputParams['Tn'] = '1.00E+01'                              #initial gas temperature (n,i,e) (K) (used only if shock_type = S1 or S2)
    inputParams['op_H2_in'] = '3'                               #initial H2 ortho/para ratio (used only if shock_type = S1 or S2)

    #environment
    inputParams['Zeta'] = '5.00E-17'                            #cosmic ray ionization rate (s-1)
    inputParams['F_ISRF'] = '1'                                 #radiation field spectrum - 1 = Mathis, 2 = Draine
    inputParams['RAD'] = '0.00E+00'                             #radiation field intensity (Habing units)
    inputParams['Av0'] = '1.00E-01'                             #initial extinction (magnitudes)
    inputParams['F_COUP_RAD'] = '0'                             #perform a full coupling with radiation field transfer (to compute dissociation rates, desorption, ...) - only if RAD ≠ 0
    inputParams['F_AV'] = '0'                                   #integrate Av or not
    inputParams['F_invAv'] = '0'                                #use grain coefficients to compute AV/NH (0) or scale grain coefficient to reproduce inv_Av_fac (1)
    inputParams['inv_Av_fac'] = '5.34D-22'                      #AV/NH (Galaxy : 5.34D-22) - only if F_invAv == 1
    inputParams['N_H2_0'] = '1.00E+06'                          #column density of H2 buffer (cm-2)
    inputParams['N_CO_0'] = '1.00E+06'                          #column density of CO buffer (cm-2)
    inputParams['vturb'] = '3.50E+00'                           #turbulent velocity (km s-1, used for Doppler broadening in FGK)

    #grain properties
    inputParams['F_TGR'] = '0'                                  #compute grain temperature (1) or keep it constant (0)
    inputParams['Tgrains'] = '15'                               #initial grain temperature (K)
    inputParams['amin_mrn'] = '1.00E-06'                        #grain MRN minimum radius (in cm)
    inputParams['amax_mrn'] = '3.00E-05'                        #grain MRN maximum radius (in cm)
    inputParams['alpha_mrn'] = '3.50E+00'                       #grain MRN index
    inputParams['rho_grc'] = '2.00E+00'                         #grain core volumic mass (g/cm3)
    inputParams['rho_grm'] = '1.00E+00'                         #grain mantle volumic mass (g/cm3)

    #excitation & cooling
    inputParams['ieqth'] = '1'                                  #thermal Balance (1 : solved, 0 : fixed T) - only for 'S' or 'P'
    inputParams['Cool_KN'] = '0'                                #Kaufman & Neufeld cooling (1) or analytical formula (0)
    inputParams['NH2_lev'] = '150'                              #number of H2 levels included
    inputParams['NH2_lines_out'] = '200'                        #maximum number of H2 lines in output file
    inputParams['H_H2_flag'] = 'BOTH'                           #H-H2 collisions : DRF, MM or BOTH (see documentation)
    inputParams['iforH2'] = '1'                                 #formation on grain model (option 1, 2, 3, or 4, see below)
    inputParams['ikinH2'] = '2'                                 #kinetic energy of H2 newly formed (option 1 or 2, see below)
    inputParams['pumpH2'] = '0'                                 #H2 pumping by UV photons taken into account (1) or not (0)
    inputParams['NCO_lev'] = '50'                               #number of CO levels included for FGK transfer (obsolete)
    
    #numerical parameters
    inputParams['integ_type'] = '0'                             #linear (0) or logarithmic (1)
    inputParams['Nstep_max'] = '30000'                          #maximum number of integration steps
    inputParams['timeJ'] = '9.99E+99'                           #shock age (years)
    inputParams['duration_max'] = '1.00E+09'                    #maximum shock duration (years)
    inputParams['Eps_V'] = '1.00E-07'                           #precision of computation
    inputParams['XLL'] = '1.00E+14'                             #characteristic viscous length (cm)
    
    
    #output specifications
    inputParams['F_W_HDF5_STD'] = '0'                           #write HDF5 standard output files (1) or not (0)
    inputParams['F_W_HDF5_CHE'] = '0'                           #write HDF5 chemical output files (1) or not (0)
    inputParams['F_W_ASCII'] = '1'                              #write ASCII final output files (1), save intermediary trajectories (2), or do none (0)
    inputParams['Npthdf5'] = '10000'                            #maximal number of points in HDF5 files
    inputParams['Nstep_w'] = '5'                                #number of steps between 2 outputs (for ascii and chemical HDF5 files)
    inputParams['speci_out'] = 'FD'                             #data in mhd_speci.out file - 'AD' (cm-3), 'CD' (cm-2) or 'FD' (n(x)/nH)
    inputParams['H2_out'] = 'AD'                                #data in H2_lev.out    file - 'AD' (cm-3), 'CD' (cm-2) or 'ln(N/g)'
    inputParams['line_out'] = 'local'                           #data in H2_line.out   file - 'local' (erg/s/cm3) or 'integrated' (erg/s/cm2/sr)
    inputParams['flag_analysis'] = 'N'                          #output chemical analysis (dominant reactions) (Y/N) (obsolete)

    #developer options
    inputParams['F_SORT'] = '1'                                 #sort reactions in increasing order before computing derivatives
    inputParams['F_CONS'] = '2'                                 #test for conservation laws and apply correction. (0) none, (1) only H2 and charge, (2) H2, charge, PAHs and all elements
    inputParams['F_CH'] = '0'                                   #compute CH velocity (1) or adopt neutral velocity (0)
    inputParams['F_S'] = '0'                                    #compute S  velocity (1) or adopt neutral velocity (0)
    inputParams['F_SH'] = '0'                                   #compute SH velocity (1) or adopt neutral velocity (0)
    inputParams['z0'] = '0.00E+00'                        #distance to the origin of the mass flow (used only for shock_type='W')

    #additional parameter description
    # iforH2 = 1                               ! Flag : H2 formation on grains
    #                                          !  -1: formation in the v,J = 0,0 and 0,1 levels only
    #                                          !   0: Proportional to Boltzmann distribution so that 1/3 of 4.4781 eV is in internal energy
    #                                          !   1: Proportional to Boltzmann distribution at 17249 K
    #                                          !   2: Dissociation limit : v = 14, J = 0,1 (4.4781 eV)
    #                                          !   3: v = 6, J = 0,1
    #                                          !   4: fraction = relative populations at t, initialised as H2_lev%density and changed during integration

    # ikinH2 = 2                               ! Flag : H2 formation energy released as kinetic energy
    #                                          !   1: 0.5 * (4.4781 eV - internal energy)
    #                                          !   2: MIN(1.4927 eV, 4.4781 eV - internal energy)
    
    return inputParams



#write parameter files
def writeInputFile(filename,Params):
    nStr = 41
    file = open(filename,'w')
    for key in list(Params.keys()):
        keyValue = Params[key]
        keyValueLength = len(keyValue)
        line = keyValue + (nStr - keyValueLength)*' ' + '! ' + key +'\n'
        if key == 'modele':
            file.write('!---- input files --------------------------------------------------------------\n')
        elif key == 'shock_type':
            file.write('!---- shock parameters ---------------------------------------------------------\n')
        elif key == 'Zeta':
            file.write('!---- environment --------------------------------------------------------------\n')
        elif key == 'F_TGR':
            file.write('!---- grain properties ---------------------------------------------------------\n')
        elif key == 'ieqth':
            file.write('!---- excitation & cooling -----------------------------------------------------\n')
        elif key == 'integ_type':
            file.write('!---- numerical parameters -----------------------------------------------------\n')
        elif key == 'F_W_HDF5_STD':
            file.write('!---- output specifications ----------------------------------------------------\n')
        elif key == 'F_SORT':
            file.write('!---- developer options --------------------------------------------------------\n')
        file.write(line)
    file.close()
    return True



#run a single model
def executeRun(masterDir):
    os.chdir(masterDir)
    os.system('./mhd_vode')
    return True



#Setting shock model, if UV_field_on = True, UV params should be set in specifiedParams and specifiedValues
def setModel(modelName, Params, UV_field_on = False):
    #default non-UV field, if UV_field_on = True, UV params should be set in specifiedParams and specifiedValues
    
    #set static model
    staticParams = createDefaultParams()
    staticParams
    for key in list(Params.keys()):
        if not isinstance(Params.get(key),str):
            staticParams[key] = format(Params.get(key),'.2E')
        else:
            staticParams[key] = Params.get(key)
    staticSetParams = ['modele', 'specfile', 'chemfile', 'h2exfile', 'shock_type', 'Nfluids', 'timeJ', 'duration_max']
    staticSetValues = [modelName+'-S', 'species_2022_05.in_depl', 'chemistry_2022_05_all.in_noadso', 'none', 'S1', '1', '9.99E+99', '1.00E+09']
    if UV_field_on == True:
        staticSetValues = [modelName+'-S', 'species_2022_05.in', 'chemistry_2022_05_all.in', 'none', 'S1', '1', '9.99E+99', '1.00E+09']
 
    for i in enumerate(staticSetParams):
        staticParams[i[1]] = staticSetValues[i[0]]
    
    
    #set product run model
    runParams = createDefaultParams()
    
    try:
        Nfluids = str(Params['Nfluids'])
    except:
        if Params['shock_type'] == 'C':
            Nfluids = '3'
        else:
            Nfluids = '1'
        
    for key in list(Params.keys()):
        if not isinstance(Params.get(key),str):
            runParams[key] = format(Params.get(key),'.2E')
        else:
            runParams[key] = Params.get(key)

    runSetParams = ['modele', 'specfile', 'chemfile', 'h2exfile', 'Nfluids']
    runSetValues = [modelName, '../output/'+modelName+'-S/species.out', 'chemistry_2022_05_all.in', '../output/'+modelName+'-S/h2levels.out', Nfluids]    
    for i in enumerate(runSetParams):
        runParams[i[1]] = runSetValues[i[0]]
    
    return staticParams, runParams



#run a single model
def runModel(staticParams, runParams, masterDir, Workdir):
    #static
    writeInputFile('input_mhd.in', staticParams)
    os.system('mv input_mhd.in '+ masterDir +'input/')
    executeRun(masterDir)
    
    #product 
    writeInputFile('input_mhd.in', runParams)
    os.system('mv input_mhd.in '+ masterDir +'input/')
    executeRun(masterDir)
    return staticParams['modele'], runParams['modele']