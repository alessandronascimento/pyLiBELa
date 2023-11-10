from data.grid_loader import GridLoader
from data.ligand_loader import LigandLoader
import os

data_loader = GridLoader("/data/alex/workproj/grids", (7, 20, 20, 20))
arr, grid_names = data_loader.load_grids()

exclude_list = [file for file in os.listdir("/data/alex/down/targets_refined_set/") 
        if file not in grid_names]


ligand_loader = LigandLoader("/data/alex/down/targets_refined_set/")
ligs, names = ligand_loader.load_ligands(exclude_list=exclude_list)







