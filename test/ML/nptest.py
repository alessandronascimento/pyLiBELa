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

    file_path = "/data/alex/workproj/GPTtest"
    grid1 = load_grid_to_ndarray(f"{file_path}/1bgq.grid")
    
    columns = ["elec", "vdwA", "vdwB", "rec_solv_gauss", "solv_gauss", "hb_donor", "hb_acceptor"]
    print(grid1[0][19][19])

    print(np.mean(grid1))


if __name__ == "__main__":
    main()
