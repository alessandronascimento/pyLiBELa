import sys
import os

def main():
    
    # Set the source path to the directory containing all ligand/receptor directories.
    src = "/home/fac001/workingproj/targets_refined_set"
    # Set the output path for the grids.
    out = "/home/fac001/workingproj/grids"


    dirs = list(os.listdir(src))

    with open("tasks.txt", 'w+') as f:
        for dir in dirs:
            arg1 = f"{src}/{dir}/{dir}_rec.mol2.gz"
            arg2 = f"{src}/{dir}/{dir}_ligand.mol2.gz"
            arg3 = f"{out}/{dir}"
            f.write(f"{arg1} {arg2} {arg3}\n")

if __name__ == "__main__":
    main()

