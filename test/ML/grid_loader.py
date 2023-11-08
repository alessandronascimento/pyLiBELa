import numpy as np
import os
import struct
from typing import Generator, Tuple, Any, List
from tqdm import tqdm

class GridLoader():

    def __init__(self, filepath : str, grid_shape : Tuple[int, ...]):

        self.filepath : str = filepath
        self.grid_shape : Tuple[int, ...] = grid_shape

    def _iter_binary_file(self, grid_path : str) -> Generator[Tuple[np.float64], None, None]:

        with open(grid_path, "rb") as f:
            while ( byte_stream := f.read(8) ):
                val, = struct.unpack('d', byte_stream)
                yield val

    def _load_grid_to_ndarray(self, grid_path : str) -> np.ndarray[Any, np.dtype[np.float64]]:

        arr = np.fromiter(self._iter_binary_file(grid_path), np.float64)
        arr = arr.reshape(*self.grid_shape)
        return arr

    
    def create_grids_ndarray(self) -> Tuple[np.ndarray, List[str]]:
        '''
        Loads and creates arrays the Grids 

        Returns:
        Tuple[grids, names]
        'grids' is a ndarray of the loaded grids as float64 
        'names' is a list of the respective names of each grid as strings
        '''

        grid_list = [grid_file.replace('.grid', '') for grid_file in os.listdir(self.filepath) if ".grid" in grid_file]
        num_of_files = len(grid_list)
        print(f"{num_of_files} .grid files found in {self.filepath}")

        grid_arr = np.empty((num_of_files, *self.grid_shape), dtype=np.float64)

        for i, grid in tqdm(enumerate(grid_list), 
                            desc="Loading Grids", total=num_of_files):
            grid_arr[i] = self._load_grid_to_ndarray(f"{self.filepath}/{grid}.grid")

        return (grid_arr, grid_list)




