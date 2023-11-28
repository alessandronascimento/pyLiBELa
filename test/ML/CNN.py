from data.grid_loader import GridLoader
from data.target_parser import Parser
from typing import Any
from pathlib import Path
import numpy as np



def load_data() -> tuple[np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]]]:

    root_path : Path = Path("/data/alex/workproj")

    receptor_grids = GridLoader(root_path / "rec_grids", (4, 20, 20, 20))
    rec_arr, rec_grid_names = receptor_grids.load_grids()

    ligand_grids = GridLoader( root_path / "lig_grids", (4, 20, 20, 20))
    lig_arr, lig_grid_names = ligand_grids.load_grids()

    valid = [name for name in rec_grid_names if name in lig_grid_names]
    target_dict = Parser.parse_index_file(root_path / "pyLiBELa" / "test" / "ML"/ "INDEX_refined_set.2020")

    merged_grids = np.empty((len(valid), 2, 4, 20, 20, 20)) 
    targets = np.empty(len(valid))

    for i, name in enumerate(valid):
        merged_grids[i][0] = rec_arr[rec_grid_names.index(name)]
        merged_grids[i][1] = lig_arr[lig_grid_names.index(name)]
        targets[i] = target_dict[name]

    return merged_grids, targets

if __name__ == "__main__":
    grids, targets = load_data() 
    print(grids.shape)
    print(targets.shape)
