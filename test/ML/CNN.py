from data.grid_loader import GridLoader
from data.target_parser import Parser
from pathlib import Path
import numpy as np



def load_data():

    root_path : Path = 

    receptor_grids = GridLoader("/data/alex/workproj/rec_grids", (7, 20, 20, 20))
    rec_arr, rec_grid_names = receptor_grids.load_grids()

    ligand_grids = GridLoader("/data/alex/workproj/lig_grids", (7, 20, 20, 20))
    lig_arr, lig_grid_names = ligand_grids.load_grids()

    valid = [name for name in rec_grid_names if name in lig_grid_names]
    target_dict = Parser.parse_index_file("/data/alex/workproj/pyLiBELa/test/ML/INDEX_refined_set.2020")

    merged_grids = np.empty((len(valid), 2, 7, 20, 20, 20)) 
    targets = np.empty(len(valid))

    for i, name in enumerate(valid):
        merged_grids[i][0] = rec_arr[rec_grid_names.index(name)]
        merged_grids[i][1] = lig_arr[lig_grid_names.index(name)]

    return merged_grids

if __name__ == "__main__":
    grids = load_data() 
    print(grids.shape)
