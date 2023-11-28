from data.grid_loader import GridLoader
from data.target_parser import Parser
from typing import Any
from pathlib import Path
import numpy as np
import tensorflow as tf

GRID_CHANNELS = 4
GRID_DIM = 20

def load_data() -> tuple[np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]], int]:

    root_path : Path = Path("/data/alex/workproj")
    shape = (GRID_CHANNELS, GRID_DIM, GRID_DIM, GRID_DIM)

    receptor_grids = GridLoader(root_path / "rec_grids", shape)
    rec_arr, rec_grid_names = receptor_grids.load_grids()

    ligand_grids = GridLoader( root_path / "lig_grids", shape)
    lig_arr, lig_grid_names = ligand_grids.load_grids()

    valid = [name for name in rec_grid_names if name in lig_grid_names]
    total = len(valid)
    target_dict = Parser.parse_index_file(root_path / "pyLiBELa" / "test" / "ML"/ "INDEX_refined_set.2020")

    merged_grids = np.empty((total, 2, shape[0], shape[1], shape[2], shape[3])) 
    targets = np.empty(total)

    for i, name in enumerate(valid):
        merged_grids[i][0] = rec_arr[rec_grid_names.index(name)]
        merged_grids[i][1] = lig_arr[lig_grid_names.index(name)]
        targets[i] = target_dict[name]

    merged_grids = merged_grids.reshape((total, 2*shape[0], shape[1], shape[2], shape[3]))
    return merged_grids, targets, total

if __name__ == "__main__":
    grids, targets, total = load_data() 
    grids, targets = tf.convert_to_tensor(grids), tf.convert_to_tensor(targets)
    print(grids.shape)
    print(targets.shape)

    conv = tf.keras.layers.Conv3D(2, 3, activation='relu', data_format='channels_first', input_shape=grids.shape[1:], padding="same")
    fmap = conv(grids)
    print(fmap.shape)
