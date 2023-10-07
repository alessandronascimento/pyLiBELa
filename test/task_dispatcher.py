import sys
import os

def main():
    
    dir = sys.argv[1] 
    src = "/home/fac001/workingproj/targets_refined_set"
    out = "/home/fac001/workingproj/grids" 

    arg1 = f"{src}/{dir}/{dir}_rec.mol2.gz"
    arg2 = f"{src}/{dir}/{dir}_ligand.mol2.gz"
    arg3 = f"{out}/{dir}"

    os.system(f"./process_grids {arg1} {arg2} {arg3}")

if __name__ == "__main__":
    main()

