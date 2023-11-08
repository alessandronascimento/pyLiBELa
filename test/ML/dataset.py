import tensorflow as tf 
import os

def main():

    file_path = "/home/fac001/workingproj/grids"
    grids = [file.strip(".grids") for file in os.listdir(file_path)]

    generator = torch.Generator().manual_seed(2023)
    train, test, validation = torch.utils.data.random_split(grids, (0.8, 0.1, 0.1), generator=generator)
    print(f"{train}\n{test}\n{validation}")


if __name__ == "__main__":
    main()
