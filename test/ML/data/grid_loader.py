import numpy as np
import os
import struct
from typing import Generator, Any
from tqdm import tqdm
from pathlib import Path

class GridLoader():

    def __init__(self, filepath : str, grid_shape : tuple[int, ...]):

        self.filepath : Path = Path(filepath) 
        self.grid_shape : tuple[int, ...] = grid_shape

    def _iter_binary_file(self, grid_path : Path) -> Generator[tuple[np.float64], None, None]:

        with open(grid_path, "rb") as f:
            while ( byte_stream := f.read(8) ):
                val, = struct.unpack('d', byte_stream)
                yield val

    def _load_grid_to_ndarray(self, grid_path : Path) -> np.ndarray[Any, np.dtype[np.float64]]:

        arr = np.fromiter(self._iter_binary_file(grid_path), np.float64)
        arr = arr.reshape(*self.grid_shape)
        return arr

    
    def load_grids(self) -> tuple[np.ndarray , list[str]]:
        '''
        Loads and creates arrays from the grids 

        Returns:
        tuple[grids, names]
        'grids' is a ndarray of the loaded grid file as float64 
        'names' is a list of the respective names of each grid as string
        '''

        grid_list = [(grid, grid.name.replace(".grid", "")) for grid in self.filepath.rglob("*.grid")]
        grid_names = []
  
        num_of_files = len(grid_list)
        print(f"\n{num_of_files} .grid files found in {self.filepath}")

        grid_arr = np.empty((num_of_files, *self.grid_shape), dtype=np.float64)

        for i, (grid, name) in tqdm(enumerate(grid_list), 
                            desc="Loading Grids", total=num_of_files, colour="yellow"):
            grid_arr[i] = self._load_grid_to_ndarray(grid)
            grid_names.append(name)

        return (grid_arr, grid_names)




