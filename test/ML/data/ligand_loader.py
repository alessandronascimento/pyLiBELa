import deepchem as dc
import os
from rdkit import Chem, RDLogger
from pathlib import Path
from tqdm import tqdm
from tempfile import NamedTemporaryFile

class InvalidLoadTypeException(Exception):
    pass

class LigandLoader():

    def __init__(self, filepath : str):

        self.filepath = Path(filepath) 

    def _load_ligand_as_MolGraphConv(self, ligand_path : Path):

        RDLogger.DisableLog("rdApp.*") 

        # The following is a hack to allow SDFLoader to read the file correctly
        tmp_file = NamedTemporaryFile(suffix=".sdf", delete=False)
        with open(tmp_file.name, 'w') as tmp:
            with open(ligand_path, "r") as sdf:
                for line in sdf.readlines():
                    tmp.write(line)
                    if "M  END" in line:
                        tmp.write("> <placeholder> (0)\n0\n\n")

        loader = dc.data.SDFLoader(['placeholder'], featurizer=dc.feat.MolGraphConvFeaturizer(use_edges=True), sanitize=True)
        disk_dataset = loader.create_dataset([f"{tmp_file.name}"])
        os.unlink(tmp_file.name)
        try: 
            data = disk_dataset.X
        except ValueError:
            return None
        else:
            return data[0] 

    def load_ligands(self, exclude_list : list[str] = []) -> tuple[list[Chem.SDMolSupplier | dc.feat.GraphData] , list[str]]:
        '''
        Loads and creates ligands lists from the defined filepath

        Returns:
        Tuple[ligands, names]
        'ligands' is a list of the loaded ligand files as a rdkit mol 3D representation 
        'names' is a list of the respective names of each ligand as string
        '''

        ligand_files = [ (lig, lig.name.replace("_ligand.sdf", "")) for lig in self.filepath.rglob("*ligand.sdf") 
            if lig.name.replace("_ligand.sdf", "") not in exclude_list]
        ligand_arr = []
        ligand_names = []
        size = len(ligand_files)
        print(f"\n{size} ligand.sdf files where found in {self.filepath}")

        for ligand_path, name in tqdm(ligand_files,
                           desc="Loading Ligands", total=size, colour='yellow'):

            file = self._load_ligand_as_MolGraphConv(ligand_path)
            ligand_arr.append(file)
            ligand_names.append(name)

        return (ligand_arr, ligand_names)
