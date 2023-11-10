from rdkit import Chem, RDLogger
from pathlib import Path
from tqdm import tqdm


class LigandLoader():

    def __init__(self, filepath : str):

        self.filepath = Path(filepath) 

    def _load_ligand_as_rdkitmol(self, ligand_path : Path, suppress_warnings: bool = True) -> Chem.SDMolSupplier: 

        if suppress_warnings:
            RDLogger.DisableLog('rdApp.*')
        
        return Chem.SDMolSupplier(str(ligand_path))


    def load_ligands(self, exclude_list : list[str] = []) -> tuple[list[Chem.SDMolSupplier] , list[str]]:
        '''
        Loads and creates ligands lists from the defined filepath

        Returns:
        Tuple[ligands, names]
        'ligands' is a list of the loaded ligand files as a rdkit mol 3D representation 
        'names' is a list of the respective names of each ligand as string
        '''

        ligand_files = [ (lig, lig.name.replace("_ligand.sdf", "")) for lig in self.filepath.rglob("*ligand.sdf") 
            if lig.name.replace("_ligand.sdf", "") not in exclude_list]
        ligands_arr = []
        ligands_names = []
        size = len(ligand_files)
        print(f"{size} ligand.sdf files where found in {self.filepath}")

        for ligand_path, name in tqdm(ligand_files,
                           desc="Loading Ligands", total=size, colour='yellow'):
            file = self._load_ligand_as_rdkitmol(ligand_path)
            ligands_arr.append(next(file))
            ligands_names.append(name)

        return (ligands_arr, ligands_names)
