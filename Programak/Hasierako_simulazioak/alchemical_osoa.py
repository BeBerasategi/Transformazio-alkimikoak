# Beñat Berasategui Migueliz
# 2023-12-16

# Call it this way: python3 alchemical_osoa.py -m 'cyclopentanol'

import openmm as mm
from openmm import Vec3
from openmm import app, unit
from openmm.unit.quantity import Quantity
from openmmforcefields.generators import GAFFTemplateGenerator
from openmmtools.constants import ONE_4PI_EPS0
from openff.toolkit.topology import Molecule

# Datu fitxategiak irakurtzeko / irudikatzeko
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Openmm-ren bitartez potentzial alkimikoak definitzeko:
from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion, AlchemicalState

# Measure time:
import time
#import sys
start = time.time()
print("Simulation started", flush=True)

# Manage files and folders:
import os

# Manage command-line arguments.
import argparse

# Parser modulu hau ongi erabiltzen ikasi!!!
parser = argparse.ArgumentParser()

parser.add_argument('-m', '--molecule_name', default='cyclopentanol', type=str, help='Name of the molecule')
args = parser.parse_args()

PROJECT_ROOT_DIR = "."
OUTPUT_DIR = "outputs"
molecule_name = args.molecule_name

OUTPUT_PATH = os.path.join(PROJECT_ROOT_DIR, OUTPUT_DIR, molecule_name)
os.makedirs(OUTPUT_PATH, exist_ok=True)

RESULTS_PATH = os.path.join(OUTPUT_PATH, "result_files")
os.makedirs(RESULTS_PATH, exist_ok=True)

EQUILIBRATION_PATH = os.path.join(OUTPUT_PATH, "equilibration_files")
os.makedirs(EQUILIBRATION_PATH, exist_ok=True)
os.makedirs(EQUILIBRATION_PATH+'/solvated', exist_ok=True)
os.makedirs(EQUILIBRATION_PATH+'/vacuum', exist_ok=True)

print("Output path:" ,OUTPUT_PATH, flush=True)

# Open file where results will be saved:
results_file = open(RESULTS_PATH+"/results_"+molecule_name+".txt", "a")
print('Current Time:', time.ctime())
print('---------------------------------------------------------', file=results_file, flush=True)
print(time.ctime(), file=results_file, flush=True)


# Goian definitutako molekularekin egingo dugu lan:
print("Simulatutako molekula:", molecule_name, '\n', flush=True)
print("Simulatutako molekula:", molecule_name, '\n', file=results_file, flush=True)
molecule = Molecule('molekulak/'+molecule_name+'/Molecule.sdf')
molecule_pdb = app.PDBFile('molekulak/'+molecule_name+'/Molecule.pdb')

forcefield = app.ForceField('amber14/tip3pfb.xml')
gaff = GAFFTemplateGenerator(molecule, forcefield='gaff-2.11')
forcefield.registerTemplateGenerator(gaff.generator)

# Uretako sistema prestatu:
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

# Hutseko sistema prestatu:
vacuum_model = app.Modeller(molecule_pdb.topology, molecule_pdb.positions)
vacuum_model.topology.setPeriodicBoxVectors(Quantity(value=(
    Vec3(x=3.0, y=0.0, z=0.0),
    Vec3(x=0.0, y=3.0, z=0.0),
    Vec3(x=0.0, y=0.0, z=3.0)),
    unit = unit.nanometer
))

vacuum_system = forcefield.createSystem(
    topology = vacuum_model.topology,
    nonbondedMethod = app.PME,
    nonbondedCutoff = 0.9 * unit.nanometer,
    removeCMMotion = True,
    constraints = app.HBonds,
    rigidWater = True
)

# Atomoen indizeak gorde. Gero beharko ditugu transformazio alkimikoak zeini egin zehazteko.
ligand_index = set(range(vacuum_system.getNumParticles()))

###############################################################################
## AZKENEKO ALDAKETA: platfroms, factory eta alchemical_atoms bikoiztea...
###############################################################################

# Eskualde alkimikoa zehaztu:
alchemical_atoms_solvated = AlchemicalRegion(ligand_index)
alchemical_atoms_vacuum = AlchemicalRegion(ligand_index)

# Lehen sortutako sistemaren bertsio alkimikoa eratu.
# OpenMMToolsek beharrezko aldaketak egiten ditu potentzialean interpolazio alkimikoa egin ahal izateko.
factory_solvated = AbsoluteAlchemicalFactory()
factory_vacuum = AbsoluteAlchemicalFactory()
alchemical_system_solvated = factory_solvated.create_alchemical_system(solvated_system, alchemical_regions=alchemical_atoms_solvated)
alchemical_system_vacuum = factory_vacuum.create_alchemical_system(vacuum_system, alchemical_regions=alchemical_atoms_vacuum)

# Egoera alkimikoa sortu jadanik sortu dugun sistemaren parametroak kopiatuz.
alchemical_state_solvated = AlchemicalState.from_system(alchemical_system_solvated)
alchemical_state_vacuum = AlchemicalState.from_system(alchemical_system_vacuum)

# Egoera alkimikoaren parametroen balioak aldatu eta gure sistemari aldaketak ezarri.
def update_lambda(l_c, l_lj, alchemical_state):
    alchemical_state.lambda_electrostatics = l_c     # Coulomb
    alchemical_state.lambda_sterics = l_lj           # Lennard-Jones
    return alchemical_state # state-a aldatuta itzultzen du, ez system delakoa.

def plot_energy(file_name,data_file, report_interval):
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
    plt.close() # Esplizituki ixten ez badira, memoria asko kontsumitzen da.

# Simulazioa prestatu. Momentuz multzo kanonikoan egingo dut.
temperature = 300 * unit.kelvin
friction_coefficient = 1 / unit.picosecond
time_step = 4 * unit.femtosecond
integrator_solvated = mm.LangevinMiddleIntegrator(temperature, friction_coefficient,time_step)
integrator_vacuum = mm.LangevinMiddleIntegrator(temperature, friction_coefficient,time_step)

# Simulation objektuak sortu
platform_solvated = mm.Platform.getPlatformByName('CUDA')
platform_vacuum = mm.Platform.getPlatformByName('CUDA')
simulation_solvated = app.Simulation(solvated_model.topology, alchemical_system_solvated, integrator_solvated, platform=platform_solvated)
simulation_solvated.context.setPositions(solvated_model.positions)
simulation_vacuum = app.Simulation(vacuum_model.topology, alchemical_system_vacuum, integrator_vacuum, platform=platform_vacuum)
simulation_vacuum.context.setPositions(vacuum_model.positions)

def transformazio_alkimikoa_simulatu(simulation, model, alchemical_state, equilibration_path):
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
        # Posizioak berrezarri
        simulation.context.setPositions(model.positions) # Hobeto honekin edo hau gabe?
        # simulation.context.reinitialize(preserveState=False)

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

        print(f"{k}. bikotea orekatuta", flush=True)
        simulation.reporters.pop(0)
        equilibration_figure = os.path.join(equilibration_path, f"NVT_Equilibration_Data_{k}.png")
        plot_energy(equilibration_figure, equilibration_file, report_interval)

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

def emaitzen_tentsorea_gorde(model, U_kln, file_path):
    # Sortutako tentsorea gorde:
    np.save(file_path+f'/measurements_{model}_{molecule_name}.npy', U_kln)
    # Lehenengo esperimentuko datuak, ikusi ahal izateko:
    df = pd.DataFrame(U_kln[:,:,0])
    df.to_csv(file_path+f'/U_kln_{model}_lehena_{molecule_name}.csv')
    # Azken esperimentuko datuak, ikusi ahal izateko:
    df = pd.DataFrame(U_kln[:,:,-1])
    df.to_csv(file_path+f'/U_kln_{model}_azkena_{molecule_name}.csv')

    # Tentsorea kargatzeko:
    # U_kln = np.load('measurements_solvated.npy')

# FEP ekuazioaren inplementazioa:
def FEP_forward(U_kln):
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


U_kln = transformazio_alkimikoa_simulatu(simulation_solvated, solvated_model, alchemical_state_solvated, EQUILIBRATION_PATH+'/solvated')
emaitzen_tentsorea_gorde("solvated", U_kln, RESULTS_PATH)
print("Solvated system calculated and saved.", flush=True)
print("Time needed for calculating the solvated system: ", end='', file= results_file, flush=True)
end = time.time()
print((end - start)/60, "mins = ",(end - start)/3600, "h", file= results_file, flush=True)

# U_kln = np.load(RESULTS_PATH+f'/measurements_solvated_{molecule_name}.npy')
deltaG_solvated = FEP_forward(U_kln)
print(f"Solvated, free energy difference: {deltaG_solvated}", flush=True)
print(f"Solvated, free energy difference: {deltaG_solvated}", file= results_file, flush=True)

start = time.time()
U_kln = transformazio_alkimikoa_simulatu(simulation_vacuum, vacuum_model, alchemical_state_vacuum, EQUILIBRATION_PATH+'/vacuum')
emaitzen_tentsorea_gorde("vacuum", U_kln, RESULTS_PATH)
print("Vacuum system calculated and saved.", flush=True)
print("Time needed for calculating the vacuum system: ", end='', file= results_file, flush=True)
end = time.time()
print((end - start)/60, "mins = ",(end - start)/3600, "h", file= results_file, flush=True)

#U_kln = np.load(RESULTS_PATH+f'/measurements_vacuum_{molecule_name}.npy')
deltaG_vacuum = FEP_forward(U_kln)
print(f"Vacuum, free energy difference: {deltaG_vacuum}", flush=True)
print(f"Vacuum, free energy difference: {deltaG_vacuum}", file= results_file, flush=True)

deltaG = (deltaG_vacuum - deltaG_solvated) 
deltaG_cal = deltaG.in_units_of(unit.kilocalories_per_mole)

print(f"Solvatation free energy: {deltaG} = {deltaG_cal}", flush=True)
print(f"Solvatation free energy: {deltaG} = {deltaG_cal}", file= results_file, flush=True)

print('\n', file= results_file)

results_file.close()

