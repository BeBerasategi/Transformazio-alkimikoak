# Alchemical transformations: Code for production
# Beñat Berasategui. 2024-02-26.

# This code is intended to use in the Hyperion server of DIPC, to obtain more data and improve the statistical quantities 
# estimated from the simulations up to now. New molecules will be used.

# Call it this way: python3 alchemical_osoa.py -m '1,4-dioxane'

#############################################################################
##           IN THIS VERSION REQUIRED CODE FOR NPT IS ADDED                ##
#############################################################################

#---------------------------------------------------------------------------------------------------------------------------------------

import os
import sys
import argparse

# Silence command-line output temporarily
# sys.stdout, sys.stderr = os.devnull, os.devnull

# Unsilence command-line output
# sys.stdout, sys.stderr = sys.__stdout__, sys.__stderr__

# Read, plot and manipulate files and arrays
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Redefine print function to have flush=True
import functools
print = functools.partial(print, flush=True)

# Import the functions written for the simulation
from simulations_funcs import *

# Measure time:
import time

print("Configuring simulation------------------------------------------------")

# Use parser to get arguments from the command line
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--molecule', default='bromoform', type=str, help='Name of the molecule')
args = parser.parse_args()

# Configure file system
MOLEC_ROOT_DIR = "../../../dipc/bberasategui/Alchemical_transformations"
PROJECT_ROOT_DIR = "."
molecule_name = args.molecule

# Test values obtained in previous simulations for 1,4-dioxane
# Exp: -5.06 kcal/mol
# NTV_NVTnp: -5.2306 kcal/mol
# NPT_NVTnp: -3.7413 kcal/mol

OUTPUT_DIR = "outputs"
OUTPUT_SIMULATION = ["outputs_NVT", "outputs_NPT"]
RESULTS_PATH_DICT = dict.fromkeys(OUTPUT_SIMULATION, 0)
EQUILIBRATION_PATH_DICT = dict.fromkeys(OUTPUT_SIMULATION, 0)

for simulation_output in OUTPUT_SIMULATION:
    OUTPUT_PATH = os.path.join(PROJECT_ROOT_DIR, OUTPUT_DIR, simulation_output, molecule_name)
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    RESULTS_PATH = os.path.join(OUTPUT_PATH, "result_files")
    RESULTS_PATH_DICT[simulation_output] = RESULTS_PATH
    os.makedirs(RESULTS_PATH, exist_ok=True)

    EQUILIBRATION_PATH = os.path.join(OUTPUT_PATH, "equilibration_files")
    EQUILIBRATION_PATH_DICT[simulation_output] = EQUILIBRATION_PATH
    os.makedirs(EQUILIBRATION_PATH, exist_ok=True)
    os.makedirs(EQUILIBRATION_PATH+'/solvated', exist_ok=True)
    os.makedirs(EQUILIBRATION_PATH+'/vacuum', exist_ok=True)

    print(f"{simulation_output:14}: ", OUTPUT_PATH)

molecule_path = os.path.join(MOLEC_ROOT_DIR, 'molekulak', molecule_name)
molecule = Molecule(molecule_path + '/Molecule.sdf')
molecule_pdb = app.PDBFile(molecule_path + '/Molecule.pdb')

print(" ")

print("Starting simulation---------------------------------------------------")

#######################################
##            NPT-NVTp               ##
#######################################

# SOLVATED ---------------------------------------------------------------------------------------------
start = time.time()
molecule, molecule_pdb, results_file = load_molecule(RESULTS_PATH_DICT[OUTPUT_SIMULATION[1]], molecule_path, molecule_name)
simulation_solvated, solvated_model, alchemical_state_solvated = configure_system(molecule, molecule_pdb, system='solvated')

print("[Simulation in NPT-NVTnp]", file= results_file)

U_kln = transformazio_alkimikoa_simulatu(simulation_solvated, solvated_model, alchemical_state_solvated, EQUILIBRATION_PATH_DICT[OUTPUT_SIMULATION[1]]+'/solvated', NPT=True)
emaitzen_tentsorea_gorde(molecule_name, "solvated", U_kln, RESULTS_PATH_DICT[OUTPUT_SIMULATION[1]])
print("\nSolvated system calculated and saved.")
print("Time needed for calculating the solvated system: ", end='', file=results_file)
end = time.time()
print((end - start)/60, "mins = ",(end - start)/3600, "h", file= results_file)

# U_kln = np.load(RESULTS_PATH+f'/measurements_solvated_{molecule_name}.npy')
deltaG_solvated_NPT = FEP_forward(U_kln)
print(f"Solvated, free energy difference: {deltaG_solvated_NPT}\n")
print(f"Solvated, free energy difference: {deltaG_solvated_NPT}", file= results_file)

# VACUUM -----------------------------------------------------------------------------------------------
start = time.time()
simulation_vacuum, vacuum_model, alchemical_state_vacuum = configure_system(molecule, molecule_pdb, system='vacuum')

U_kln = transformazio_alkimikoa_simulatu(simulation_vacuum, vacuum_model, alchemical_state_vacuum, EQUILIBRATION_PATH_DICT[OUTPUT_SIMULATION[1]]+'/vacuum')
emaitzen_tentsorea_gorde(molecule_name, "vacuum", U_kln, RESULTS_PATH_DICT[OUTPUT_SIMULATION[1]])
emaitzen_tentsorea_gorde(molecule_name, "vacuum", U_kln, RESULTS_PATH_DICT[OUTPUT_SIMULATION[0]])
print("\nVacuum system calculated and saved.")
print("Time needed for calculating the vacuum system: ", end='', file= results_file)
end = time.time()
print((end - start)/60, "mins = ",(end - start)/3600, "h", file= results_file)

deltaG_vacuum = FEP_forward(U_kln)
print(f"Vacuum, free energy difference: {deltaG_vacuum}\n")
print(f"Vacuum, free energy difference: {deltaG_vacuum}", file= results_file)

deltaG = (deltaG_vacuum - deltaG_solvated_NPT) 
deltaG_cal = deltaG.in_units_of(unit.kilocalories_per_mole)

print("*******************************************************************************************")
print(f"Solvatation free energy NPT-NVTnp: {deltaG} = {deltaG_cal}")
print("*******************************************************************************************\n\n")
print(f"\nSolvatation free energy NPT-NVTnp: {deltaG} = {deltaG_cal}", file= results_file)

print('\n', file= results_file)
results_file.close()


#######################################
##            NVT-NVTp               ##
#######################################

# SOLVATED ---------------------------------------------------------------------------------------------
start = time.time()
molecule, molecule_pdb, results_file = load_molecule(RESULTS_PATH_DICT[OUTPUT_SIMULATION[0]], molecule_path, molecule_name)
simulation_solvated, solvated_model, alchemical_state_solvated = configure_system(molecule, molecule_pdb, system='solvated')

print("[Simulation in NVT-NVTnp]", file= results_file)

U_kln = transformazio_alkimikoa_simulatu(simulation_solvated, solvated_model, alchemical_state_solvated, EQUILIBRATION_PATH_DICT[OUTPUT_SIMULATION[0]]+'/solvated')
emaitzen_tentsorea_gorde(molecule_name, "solvated", U_kln, RESULTS_PATH_DICT[OUTPUT_SIMULATION[0]])
print("Solvated system calculated and saved.")
print("Time needed for calculating the solvated system: ", end='', file=results_file)
end = time.time()
print((end - start)/60, "mins = ",(end - start)/3600, "h", file= results_file)

# U_kln = np.load(RESULTS_PATH+f'/measurements_solvated_{molecule_name}.npy')
deltaG_solvated = FEP_forward(U_kln)
print(f"Solvated, free energy difference: {deltaG_solvated}")
print(f"Solvated, free energy difference: {deltaG_solvated}", file= results_file)

# VACUUM: SAME AS BEFORE! -----------------------------------------------------------------------------------------

# deltaG_vacuum SAME AS BEFORE 
print(f"Vacuum, free energy difference: {deltaG_vacuum}")
print(f"Vacuum, free energy difference: {deltaG_vacuum}", file= results_file)

deltaG = (deltaG_vacuum - deltaG_solvated) 
deltaG_cal = deltaG.in_units_of(unit.kilocalories_per_mole)

print("*******************************************************************************************")
print(f"Solvatation free energy NVT-NVTnp: {deltaG} = {deltaG_cal}")
print("*******************************************************************************************\n\n")
print(f"Solvatation free energy NVT-NVTnp: {deltaG} = {deltaG_cal}", file= results_file)


print('\n', file= results_file)
results_file.close()

