import os

def main():
    
    # Set the source path to the directory containing all ligand/receptor directories.
    src = "/data/alex/down/targets_refined_set"
    # Set the output path for the grids.
    out = "/data/alex/workproj/lig_grids"


    dirs = list(os.listdir(src))

    with open("tasks.txt", 'w+') as f:
        for dir in dirs:
            arg1 = f"{src}/{dir}/{dir}_ligand.mol2.gz"     # receptor file
            arg2 = f"{src}/{dir}/{dir}_ligand.mol2.gz"  # ligand file
            arg3 = f"{out}/{dir}"                       # output directory
            f.write(f"{arg1} {arg2} {arg3}\n")

if __name__ == "__main__":
    main()

