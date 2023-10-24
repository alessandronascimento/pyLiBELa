import sys
import os

def main():
    
    src = "/home/alex/data/alex/down/targets_refined_set"
    out = "/home/alex/data/alex/workproj/grids" 

    dirs = list(os.listdir(src))

    with open("tasks.txt", 'w+') as f:
        for dir in dirs:
            arg1 = f"{src}/{dir}/{dir}_rec.mol2.gz"
            arg2 = f"{src}/{dir}/{dir}_ligand.mol2.gz"
            arg3 = f"{out}/{dir}"
            f.write(f"{arg1} {arg2} {arg3}\n")

if __name__ == "__main__":
    main()

