try:
  from pyPARSER import *
  from pyMol2 import *
  from pyWRITER import *
  from pyGrid import *
  from pyCOORD_MC import *
  from pyFindHB import *
  from pyEnergy2 import *
  from pyGaussian import *
  from pyConformer import *
  from pyRAND import *
  from pyMcEntropy import *
  from pySA import *
  from pyOptimizer import *
  from pyMC import *
  from pyFullSearch import *
  from pyDocker import *
  print('pyLiBELa is imported!')
except ImportError as err:
  raise err
import os

os.system("rm -rf dbfda.src.txt")

source_ligands_file = "http://files.docking.org/catalogs/source/dbfda.src.txt" #@param {type:"string"}

os.system(f'wget {source_ligands_file}')

os.system('mkdir output')

source_ligands = open(r'dbfda.src.txt','r')

ligand_smiles_list = [line.split('\t')[0] for line in source_ligands.readlines()]

source_ligands.close()


import timeit
import numpy as np

Input = PARSER()

# Some default options:
# 1. We are using two processors for Grid calculations;
#

Input.generate_conformers = True;
Input.dock_parallel = True;
Input.parallel_jobs = 2;
Input.write_grids = True
Input.use_grids = True
Input.write_mol2 = True
Input.atom_limit = 60 #@param {type:"number"}
    #atom limit not counting H

# Energy Calculation Parameters:
scoring_function = "0" #@param ["0", "1", "2", "3"]
Input.dielectric_model = "r" #@param ["r", "constant"]
Input.scoring_function = int(scoring_function)
grid_dimension = 30.0 #@param {type:"number"}
Input.grid_spacing = 0.5 #@param {type:"number"}
Input.solvation_alpha = 0.1 #@param {type:"number"}
Input.solvation_beta = -0.005 #@param {type:"number"}

# Optimization parameter:
Input.min_tol = 1E-10;
Input.min_delta = 1E-5;
Input.dock_min_tol = 1E-10;
search_box = 20.0 #@param {type:"number"}
Input.timeout = 30 #@param {type:"number"}
Input.min_timeout = 30 #@param {type:"number"}
Input.overlay_optimizer = "ln_auglag" #@param ["mma", "ln_auglag", "subplex", "none"]
Input.energy_optimizer = "direct" #@param ["direct", "isres", "crs", "esch", "stogo", "mma", "simplex", "none"]
if (Input.scoring_function < 3):
  delta = 2.5 #@param {type:"number"}
  Input.deltaij6 = (delta*delta*delta*delta*delta*delta)
  delta_es = 2.5 #@param {type:"number"}
  Input.deltaij_es6 = pow(delta_es, 6);
  Input.deltaij_es3 = (delta_es*delta_es*delta_es)

Input.search_box_x, Input.search_box_y, Input.search_box_z = search_box, search_box, search_box;
Input.x_dim, Input.y_dim, Input.z_dim = grid_dimension, grid_dimension, grid_dimension;

ligand_smiles = "NCC(O)=O" #@param {type:"string"}

os.system("rm -rf temp.smi")

with open("temp.smi", "w") as f:
    f.write(ligand_smiles)
#print(ligand_smiles)
#NCC(O)=O

os.system("rm -rf Lig2.mol2")
os.system("obabel -ismi temp.smi -omol2 -O Lig2.mol2")

Lig2 = Mol2(Input, "Lig2.mol2")


lig_src =  'test/lig.mol2.gz' #@param {type:"string"}
rec_src = 'test/rec.mol2.gz' #@param {type:"string"}

REC = Mol2(Input, rec_src)
RefLig = Mol2(Input, lig_src)

print('Receptor and reference ligand parsed successfully!')
print("Receptor has %4d atoms." % REC.N)
print('Reference Ligand has %4d atoms' % RefLig.N)
print('Search ligand has %4d atoms' % Lig2.N)

os.system("rm -rf McLiBELa_dock.mol2.gz")

Writer = WRITER(Input)
Coord = COORD_MC()
HB = FindHB()
Ene  = Energy2(Input)

for i in range(len(REC.residue_pointer)-1):
  HB.parse_residue(REC.residue_pointer[i]-1, REC.residue_pointer[i+1]-2, REC.resnames[i], REC, RefLig, 9.0);
HB.find_ligandHB(lig_src, RefLig);
print('The receptor has %5d / %5d HB donors/acceptors around the active site.' % (len(REC.HBdonors), len(REC.HBacceptors)));

Dock = Docker(Writer)
center = Coord.compute_com(RefLig)

print()
start_energy = Ene.compute_ene(REC, RefLig, RefLig.xyz);
print('Starting energy: %7.3f kcal/mol' % start_energy);
print()
print('Generating grids. This may take a while..')
Grids = Grid(Input, Writer, REC, center)
print('Grids computed!')
grid_energy = Ene.compute_ene(Grids, RefLig, RefLig.xyz);
print('Grid original energy: %7.3f kcal/mol' % grid_energy);
print('Grid error: %7.3f %%', abs(100.*(start_energy-grid_energy)/start_energy));
print()
print()
print('Starting docking calculation...')
Dock.run(REC, Lig2, RefLig, center, Input, Grids, 0)
print('Docking calculation finished!')

Writer.write_box(center, Grids.xbegin, Grids.ybegin, Grids.zbegin, Grids.xend, Grids.yend, Grids.zend)