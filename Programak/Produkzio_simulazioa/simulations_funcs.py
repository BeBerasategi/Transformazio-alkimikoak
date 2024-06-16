# Molecular dynamics
import openmm as mm
from openmm import Vec3
from openmm import app, unit
from openmm.unit.quantity import Quantity
from openmmforcefields.generators import GAFFTemplateGenerator
from openmmtools.constants import ONE_4PI_EPS0
from openff.toolkit.topology import Molecule
from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion, AlchemicalState

import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Redefine print function to have flush=True
import functools
print = functools.partial(print, flush=True)

def load_molecule(results_path, molecule_path, molecule_name):
    # Open file where results will be saved:
    results_file = open(results_path+"/results_"+molecule_name+".txt", "a")
    print('Configuration of the system started at:', time.ctime())
    print('---------------------------------------------------------', file=results_file)
    print('Current time: ', time.ctime(), file=results_file)

    # Load the selected molecule
    print("Simulatutako molekula:", molecule_name, '\n')
    print("Simulatutako molekula:", molecule_name, '\n', file=results_file)
    molecule = Molecule(molecule_path+'/Molecule.sdf')
    molecule_pdb = app.PDBFile(molecule_path+'/Molecule.pdb')
    return molecule, molecule_pdb, results_file

def configure_system(molecule, molecule_pdb, system:str):
    forcefield = app.ForceField('amber14/tip3pfb.xml')
    gaff = GAFFTemplateGenerator(molecule, forcefield='gaff-2.11')
    forcefield.registerTemplateGenerator(gaff.generator)

    # Solvated system
    solvated_model = app.Modeller(molecule_pdb.topology, molecule_pdb.positions)
    solvated_model.addSolvent(forcefield, model='tip3p', padding=1.5*unit.nanometer)

    solvated_system = forcefield.createSystem(
        topology = solvated_model.topology,
        nonbondedMethod = app.PME,
        nonbondedCutoff = 0.9 * unit.nanometer,
        removeCMMotion = True,
        constraints = app.HBonds,
        rigidWater = True
    )

    # Vacuum system
    vacuum_model = app.Modeller(molecule_pdb.topology, molecule_pdb.positions)
    vacuum_model.topology.setPeriodicBoxVectors(Quantity(value=(
        Vec3(x=3.0, y=0.0, z=0.0),
        Vec3(x=0.0, y=3.0, z=0.0),
        Vec3(x=0.0, y=0.0, z=3.0)),
        unit = unit.nanometer
    ))

# Periodic Boundary Conditions for vaccum:
#    vacuum_system = forcefield.createSystem(
#        topology = vacuum_model.topology,
#        nonbondedMethod = app.PME,
#        nonbondedCutoff = 0.9 * unit.nanometer,
#        removeCMMotion = True,
#        constraints = app.HBonds,
#        rigidWater = True
#    )
    
    # Non-periodic boundary conditions and no Cutoff for vacuum:
    vacuum_system = forcefield.createSystem(
        topology = vacuum_model.topology,
        nonbondedMethod = app.NoCutoff, # Hau da aldaketa nagusia! Cutoffak eta PBCak kentzeko.
        # Remove: nonbondedCutoff = ..., 
        removeCMMotion = True,
        constraints = app.HBonds,
        rigidWater = True
    )

    # Record the indexes of the atoms that will undergo alchemical transformations
    ligand_index = set(range(vacuum_system.getNumParticles()))

    # This parameters are critical. Default values are being used.
    temperature = 300 * unit.kelvin
    friction_coefficient = 1 / unit.picosecond
    time_step = 4 * unit.femtosecond    

    if system == 'solvated':
        # Specify the alchemical region
        alchemical_atoms_solvated = AlchemicalRegion(ligand_index)
        # Create the alchemical version of the system
        factory_solvated = AbsoluteAlchemicalFactory()
        alchemical_system_solvated = factory_solvated.create_alchemical_system(solvated_system, alchemical_regions=alchemical_atoms_solvated)
        # Create the alchemical state
        alchemical_state_solvated = AlchemicalState.from_system(alchemical_system_solvated)
        # Select integrator
        integrator_solvated = mm.LangevinMiddleIntegrator(temperature, friction_coefficient,time_step)
        platform_solvated = mm.Platform.getPlatformByName('CUDA')
        simulation_solvated = app.Simulation(solvated_model.topology, alchemical_system_solvated, integrator_solvated, platform=platform_solvated)
        simulation_solvated.context.setPositions(solvated_model.positions)

        return simulation_solvated, solvated_model, alchemical_state_solvated

    elif system == 'vacuum':
        alchemical_atoms_vacuum = AlchemicalRegion(ligand_index)
        factory_vacuum = AbsoluteAlchemicalFactory()
        alchemical_system_vacuum = factory_vacuum.create_alchemical_system(vacuum_system, alchemical_regions=alchemical_atoms_vacuum)
        alchemical_state_vacuum = AlchemicalState.from_system(alchemical_system_vacuum)
        integrator_vacuum = mm.LangevinMiddleIntegrator(temperature, friction_coefficient,time_step)
        platform_vacuum = mm.Platform.getPlatformByName('CUDA')
        simulation_vacuum = app.Simulation(vacuum_model.topology, alchemical_system_vacuum, integrator_vacuum, platform=platform_vacuum)
        simulation_vacuum.context.setPositions(vacuum_model.positions)

        return simulation_vacuum, vacuum_model, alchemical_state_vacuum
    

def update_lambda(l_c, l_lj, alchemical_state):
    ''' 
    Change the values of the parameters of the alchemical state.
    Apply the changes to our system

    Note
    ----
    It returns the modifier state, not the system.
    '''
    alchemical_state.lambda_electrostatics = l_c     # Coulomb
    alchemical_state.lambda_sterics = l_lj           # Lennard-Jones
    return alchemical_state 


def plot_energy(file_name, data_file, report_interval):
    '''
    Function to plot the energy and volume fluctuations.
    It provides a way to see if the equilibration process has worked.
    '''
    state_data = pd.read_csv(data_file)
    state_data.columns = ['Potential Energy (kJ/mole)', 'Box Volume (nm^3)']
    time = np.arange(0, state_data.shape[0] * report_interval, report_interval)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))

    ax[0].plot(time, state_data['Potential Energy (kJ/mole)'])
    ax[0].set_xlabel('Step', fontweight='bold')
    ax[0].set_ylabel('Potential Energy (kJ/mole)', fontweight='bold')
    ax[1].plot(time, state_data['Box Volume (nm^3)'])
    ax[1].set_xlabel('Step', fontweight='bold')
    ax[1].set_ylabel('Box Volume (nm^3)', fontweight='bold')

    fig.tight_layout()
    fig.savefig(file_name, dpi=300)
    plt.close() # This is needed in order to manage memory efficiently.
    return state_data

def transformazio_alkimikoa_simulatu(simulation, model, alchemical_state, equilibration_path, NPT=False):
    # Parameters needed for the simulation in the NPT set:
    pressure = 1.0 * unit.atmosphere
    temperature = 300 * unit.kelvin
    mc_frequency = 25

    # Interpolazioa 16 baliotan egingo da:
    # 0-tik 1-ra ordenatuta, erroreak eman!
    lambda_coul = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    lambda_lj = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

    # Lambda kopurua:
    n_lambda = len(lambda_coul)
    # Egoera alkimiko bakoitzeko iterazio kopurua
    n_experiments = 5000 # 5ns/1ps =5000
    # Lagin bakoitzeko pausu kopurua:
    n_steps = 250 # 4fs * 250 = 1ps

    # Emaitzak gordetzeko, 3 heineko arraya
    U_kln = np.zeros([n_lambda,n_lambda,n_experiments], np.float64)

    # Balioak binaka aldatu:
    for k in range(0,n_lambda): 
        # This is only used in NPT simulations, to go back to NVT
        # It should ONLY BE USED for the SOLVATED SYSTEM.
        if NPT and k>0 :
            for i in range(simulation.system.getNumForces()):
                if simulation.system.getForce(i).getName() == 'MonteCarloBarostat':
                    simulation.system.removeForce(i)
                    break
            simulation.context.reinitialize(preserveState=False) # INPORTANTEA BOLUMENA BERRABIARAZTEKO

        # Posizioak berrezarri
        simulation.context.setPositions(model.positions) # Hobeto honekin edo hau gabe?

        # Lambdaren balio berriak jarri sisteman.
        update_lambda(lambda_coul[k], lambda_lj[k], alchemical_state).apply_to_context(simulation.context)
        
        # Energia minimizatu
        simulation.minimizeEnergy() 

        # Sistema orekatu    
        report_interval = 100   # 100 pausuero (4fs * 100 = 400fs) neurketa egin.
        equilibration_file = os.path.join(equilibration_path, f"NVT_Equilibration_Data_{k}.dat")
        data_reporter = app.StateDataReporter(file=equilibration_file, reportInterval=report_interval,potentialEnergy=True, volume=True)
        simulation.reporters.append(data_reporter)
        # 50ps-ko (4fs * 5000 = 20ps) simulazioa.
        simulation.step(5000)

        print(f"{k}. bikotea NVT-en orekatuta")
        simulation.reporters.pop(0)
        equilibration_figure = os.path.join(equilibration_path, f"NVT_Equilibration_Data_{k}.png")
        plot_energy(equilibration_figure, equilibration_file, report_interval)

        # This is only used in NPT simulations, to go NVT->NPT.
        # It should ONLY BE USED for the SOLVATED SYSTEM.
        if NPT:
            simulation.system.addForce(mm.MonteCarloBarostat(pressure, temperature, mc_frequency))
            simulation.context.reinitialize(preserveState=True)
            equilibration_file_NPT = os.path.join(equilibration_path, f"NPT_Equilibration_Data_{k}.dat")
            data_reporter_NPT = app.StateDataReporter(file=equilibration_file_NPT, reportInterval=report_interval, potentialEnergy=True, volume=True)
            simulation.reporters.append(data_reporter_NPT)
            # 50ps-ko (2* 4fs * 5000 = 2*20ps) simulazioa. Luzeagoa, ziurtatzeko.
            simulation.step(2*5000)  

            print(f"{k}. bikotea NPT-en orekatuta")
            simulation.reporters.pop(0)
            equilibration_figure = os.path.join(equilibration_path, f"NPT_Equilibration_Data_{k}.png")
            npt_data = plot_energy(equilibration_figure, equilibration_file_NPT, report_interval)
            npt_stats_file = os.path.join(equilibration_path, f"0_NPT_Equilibration_Data_Statistics.txt")
            with open(npt_stats_file, 'a') as stats_file:
                print("---------------------------------------------------------", file=stats_file)
                print(f"Bikotea: {k}. NPT-n orekatzea. Estatistikak: ", file=stats_file)
                print("[Batez-besteko balioak]", file=stats_file)
                print(npt_data.mean(axis=0), file=stats_file)
                print("\n[Desbideratze estandarrak]", file=stats_file)
                print(npt_data.std(axis=0), file=stats_file, end=2*'\n')      
                
        for iteration in range(n_experiments): 
            simulation.step(n_steps) # 1ps
            
            # Compute energies at 3 alchemical states (n-1,n,n+1)
            # Save them in the U_kln array in kJ/mol units.
            U_kln[k,k,iteration] = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

            if k>0:
                update_lambda(lambda_coul[k-1], lambda_lj[k-1], alchemical_state).apply_to_context(simulation.context)
                U_kln[k,k-1,iteration] = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            if k<(n_lambda-1):
                update_lambda(lambda_coul[k+1], lambda_lj[k+1], alchemical_state).apply_to_context(simulation.context)
                U_kln[k,k+1,iteration] = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            
            update_lambda(lambda_coul[k], lambda_lj[k], alchemical_state).apply_to_context(simulation.context)

    return U_kln

def emaitzen_tentsorea_gorde(molecule_name, model, U_kln, results_path):
    '''
    Save the tensor containing all the outputs of the simulations (in kJ/mol).

    Parameters
    ----------
    model : str
        'solvated' or 'vacuum'
    file_path : str
        Path of the 'result_files' folder
    '''

    # Sortutako tentsorea gorde:
    np.save(results_path+f'/measurements_{model}_{molecule_name}.npy', U_kln)
    # Lehenengo esperimentuko datuak, ikusi ahal izateko:
    df = pd.DataFrame(U_kln[:,:,0])
    df.to_csv(results_path+f'/U_kln_{model}_lehena_{molecule_name}.csv')
    # Azken esperimentuko datuak, ikusi ahal izateko:
    df = pd.DataFrame(U_kln[:,:,-1])
    df.to_csv(results_path+f'/U_kln_{model}_azkena_{molecule_name}.csv')

    # Tentsorea kargatzeko:
    # U_kln = np.load('measurements_solvated.npy')


def FEP_forward(U_kln, temperature=300* unit.kelvin):
    '''
    THIS FUNCTION NEEDS TO BE COMPARED WITH ITS NEWEST VERSION
    '''
    kT = unit.AVOGADRO_CONSTANT_NA*unit.BOLTZMANN_CONSTANT_kB*temperature
    kT = kT.value_in_unit(unit.kilojoule_per_mole) # Keep just the number!

    # 16. bikotearekin kasuan, kontuz, ezin dugu s+1 hartu.
    nlambda = U_kln.shape[0]
    niterations = U_kln.shape[2]
    u_diff = np.zeros((nlambda-1,niterations)) #-1, 16. bikotera ez iristeko.
    
    for k in range(nlambda-1):
        u_diff[k,:] = U_kln[k,k+1,:]-U_kln[k,k,:]

    exp = np.exp(-u_diff/kT)
    return np.sum((-kT)*np.log(np.mean(exp, axis=1))) * unit.kilojoules_per_mole