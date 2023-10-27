import deepchem as dc
from rdkit import Chem

print("\n\n\n")

file = Chem.SDMolSupplier("/home/alex/data/alex/down/targets_refined_set/1bgq/1bgq_ligand.sdf")
mol = next(file)

featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
out = featurizer.featurize(mol)

print(out)
