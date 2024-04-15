import torch
from math import floor
from torch.utils.data import Dataset
import numpy as np
from pathlib import Path

path_to_pdbbind = "/data/home/alessandro/PDBbind_v2020"
LOG_DIR = "./log"

def read_pdbbind_data():
    #@title Getting list of targets with grids
    targets_ok = []
    target_ok_file = open(f'{path_to_pdbbind}/scripts/targets_ok.dat', 'r')
    for line in target_ok_file:
        line2 = line.split()
        targets_ok.append(line2[0])
    target_ok_file.close()

    #@title Reading PDBBind v2020 data
    binding_targets = []             # pdb id
    binding_years = []               # year of the structure
    binding_resolution = []          # resolution
    binding_score = []               # -logKd/Ki

    with open(f'{path_to_pdbbind}/index/INDEX_refined_data.2020', 'r') as binding_file:
        for line in binding_file:
            line2 = line.split()
            if line2[0][0] == '#': continue
            if line2[0] not in targets_ok: continue

            binding_targets.append(line2[0])
            if (line2[1] == 'NMR'):
                binding_resolution.append(0.00)
            else:
                binding_resolution.append(float(line2[1]))
                
            binding_years.append(int(line2[2]))
            binding_score.append(float(line2[3]))

    print('Binding data found for %d valid targets' % len(binding_targets))
    return binding_targets, binding_score

def process_file(file_path):
    as_path = Path(file_path)
    if not as_path.exists():
        raise IOError(f"File path '{as_path}' does not exist")

    with open(file_path, 'rb') as f:
        data = torch.frombuffer(f.read(), dtype=torch.float64).reshape(60, 60, 60, 3)
    return data


class GridDataset(Dataset):
    def __init__(self, targets_list, scores_list) -> None:
        self.targets = targets_list
        self.scores = scores_list 
        self._rec_idx = None
        # Sample weights
        self.weight_normal = 11/2
        self.weight_decoy = 11/20
    
    def __len__(self):
        # every targets has 10 decoys + 1 normal
        return len(self.targets) * 11

    def __getitem__(self, idx):
        grid_idx = floor(idx/11)
        inner_idx = idx % 11
        path = f"{path_to_pdbbind}/targets/{self.targets[grid_idx]}/grid_30_0.5_SF0"
        
        # Prevent from reloading rec grid every call
        if idx != self._rec_idx:
            self._rec = process_file(f"{path}/McGrid_rec.grid")
            self._rec_idx = idx

        if inner_idx == 0:
            self._file = f"{path}/McGrid_lig.grid"
            weight = self.weight_normal
        else:
            self._file = f"{path}/McGrid_dec_{inner_idx}.grid"
            weight = self.weight_decoy

        self._lig = process_file(self._file)
        return torch.cat((self._rec, self._lig), dim=-1), self.scores[grid_idx], weight


def test():
    #data = process_file(f'{path_to_pdbbind}/targets/8gpb/grid_30_0.5_SF0/McGrid_lig.grid')
    #print(data)
    targets, scores = read_pdbbind_data()
    dataset = GridDataset(targets[:23], scores[:23])
    for data in dataset:
        print(f"{dataset._file}:")
        print(data[0].dim(), data[1], data[2])

if __name__ == '__main__':
    test()