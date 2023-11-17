import deepchem as dc
import tempfile
from deepchem.feat import ConvMolFeaturizer
from rdkit import Chem
import os

print("\n\n\n")

#file1 = Chem.SDMolSupplier("/home/alex/data/alex/down/targets_refined_set/1bgq/1bgq_ligand.sdf")
#file2 = Chem.SDMolSupplier("/home/alex/data/alex/down/targets_refined_set/6g98/6g98_ligand.sdf")
#mols = [next(file1), next(file2)]


paths = ["/home/alex/data/alex/down/targets_refined_set/1bgq/1bgq_ligand.sdf", "/home/alex/data/alex/down/targets_refined_set/6g98/6g98_ligand.sdf"]
tmp_file = tempfile.NamedTemporaryFile('w+', suffix=".sdf", prefix="HERE", delete=False)
for i, ligand in enumerate(paths):


    with open(tmp_file.name, 'w') as tmp:
        with open(ligand, "r") as sdf:
            for line in sdf.readlines():
                tmp.write(line)
                if "M  END" in line:
                    tmp.write("> <test> (1)\n0\n\n")
    
    

featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
loader = dc.data.SDFLoader(['test'], featurizer=featurizer, sanitize=True)

dataset = loader.create_dataset([f"{tmp_file.name}"])
os.unlink(tmp_file.name)

print(dataset.X)
'''
featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
out = featurizer.featurize(mols)

loader = dc.data.SDFLoader(tasks=["task"], featurizer=ConvMolFeaturizer())
loader.create_dataset("/home/alex/data/alex/down/targets_refined_set/1bgq/1bgq_ligand.sdf")

print(out)
'''

