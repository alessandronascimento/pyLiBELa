import numpy as np
import pandas as pd
import os
import struct



# Define a function to read in the binary files and extract the data
def iter_binary_file(filepath : str) -> np.float64:
    with open(filepath, "rb") as f:
        while ( byte_stream := f.read(8) ):
            val, = struct.unpack('d', byte_stream)
            yield val

def load_grid_to_ndarray(filepath : str) -> np.ndarray:
    arr = np.fromiter(iter_binary_file(filepath), np.float64)
    arr = arr.reshape(7, 20, 20, 20)
    return arr


def main():

    file_path = ""

    _,_,file = next(os.walk(file_path))
    num_of_files = len(file)
    print(f"{num_of_files} files in {file_path}")

    grid_arr = np.empty((num_of_files, 7, 20, 20, 20))

    for i, grid in enumerate(os.listdir(file_path)):
        grid_arr[i] = load_grid_to_ndarray(f"{file_path}/{grid}")


if __name__ == "__main__":
    main()
