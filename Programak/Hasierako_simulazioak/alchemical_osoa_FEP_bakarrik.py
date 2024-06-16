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

temperature = 300 * unit.kelvin

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
    
    for k in range(0,nlambda-1):
        u_diff[k,:] = U_kln[k,k+1,:]-U_kln[k,k,:]

    exp = np.exp(-u_diff/kT)
    return np.sum((-kT)*np.log(np.mean(exp, axis=1))) * unit.kilojoules_per_mole


U_kln = np.load(RESULTS_PATH+f'/measurements_solvated_{molecule_name}.npy')
#U_kln = np.load(f'measurements_solvated_{molecule_name}.npy')
#U_kln = np.load(f'measurements_solvated.npy')
print("Solvated system calculated and saved.", flush=True)

deltaG_solvated = FEP_forward(U_kln)
print(f"Solvated, free energy difference: {deltaG_solvated}", flush=True)
print(f"Solvated, free energy difference: {deltaG_solvated}", file= results_file, flush=True)

U_kln = np.load(RESULTS_PATH+f'/measurements_vacuum_{molecule_name}.npy')
#U_kln = np.load(f'measurements_vacuum_{molecule_name}.npy')
#U_kln = np.load(f'measurements_vacuum.npy')
print("Vacuum system calculated and saved.", flush=True)

deltaG_vacuum = FEP_forward(U_kln)
print(f"Vacuum, free energy difference: {deltaG_vacuum}", flush=True)
print(f"Vacuum, free energy difference: {deltaG_vacuum}", file= results_file, flush=True)

deltaG = (deltaG_vacuum - deltaG_solvated) 
deltaG_cal = deltaG.in_units_of(unit.kilocalories_per_mole)

print(f"Solvatation free energy: {deltaG} = {deltaG_cal}", flush=True)
print(f"Solvatation free energy: {deltaG} = {deltaG_cal}", file= results_file, flush=True)

print('\n', file= results_file)

results_file.close()

