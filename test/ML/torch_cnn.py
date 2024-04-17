import torch
from math import floor
from torch.utils.data import Dataset
from torch import nn
from pathlib import Path
from sklearn.model_selection import train_test_split

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
        data = torch.frombuffer(f.read(), dtype=torch.float64).reshape(3, 60, 60, 60)
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
        return torch.cat((self._rec, self._lig), dim=0), self.scores[grid_idx], weight

def create_datasets(max_targets=None):
    
    binding_targets, binding_score = read_pdbbind_data()

    train_targets, test_targets, train_scores, test_scores = train_test_split(binding_targets[:max_targets],
                                                                          binding_score[:max_targets],
                                                                          train_size=0.8,
                                                                          shuffle=True)
    train_targets, valid_targets, train_scores, valid_scores = train_test_split(train_targets, train_scores, train_size=0.8)

    train_dataset = GridDataset(train_targets, train_scores) 
    test_dataset = GridDataset(test_targets, test_scores) 
    valid_dataset = GridDataset(valid_targets, valid_scores) 

    return train_dataset, test_dataset, valid_dataset

def get_device():
    device = ("cuda" if torch.cuda.is_available() 
        else "mps" if torch.backends.mps.is_available()
        else "cpu")

    return device

class CNNModel(nn.Module):

    def __init__(self) -> None:
        super().__init__()
        self.flatten = nn.Flatten()
        self._current_shape = (6,60,60,60) 
        self._old_out_channels = 0

        self.conv1 = nn.Conv3d(in_channels=6, out_channels=64, kernel_size=3)
        self.maxpool1 = nn.MaxPool3d(kernel_size=2)
        self.conv2 = nn.Conv3d(in_channels=64, out_channels=64, kernel_size=2)
        self.conv3 = nn.Conv3d(in_channels=64, out_channels=64, kernel_size=2)
        self.maxpool2 = nn.MaxPool3d(kernel_size=2)
        self.conv4 = nn.Conv3d(in_channels=64, out_channels=128, kernel_size=2)
        self.conv5 = nn.Conv3d(in_channels=128, out_channels=128, kernel_size=2)

    def make_conv(self, dims, out_channels, kernel_size, in_channels=None, padding=0, stride=1, dilation=1):
        if in_channels is None: in_channels = self._old_out_channels
        conv = nn.Conv3d(in_channels=in_channels, out_channels=out_channels, kernel_size=kernel_size)



        
    def _make_uniform(self, value):
        if isinstance(value, int):
            return (value, value, value)
        else: return value

    def _out_shape(self, dims, out_channels, kernel_size, padding=0, stride=1, dilation=1):
        """
        Helper function that returns a tuple with the shape
        information as outputted by a Conv3d or MaxPool3d 
        layer. Which can then be used as inputs for the next layer.
        """

        dims = self._make_uniform(dims)
        kernel_size = self._make_uniform(kernel_size)
        padding = self._make_uniform(padding)
        stride = self._make_uniform(stride)
        dilation = self._make_uniform(dilation) 
       
        # As specified in https://pytorch.org/docs/stable/generated/torch.nn.Conv3d.html
        f = lambda i : floor(1 + (dims[i] + 2*padding[i] - 
            dilation[i]*( kernel_size[i] - 1) - 1 ) / stride[i])

        return (out_channels, f(0), f(1), f(2))




    # model = torch.nn.Sequential(
    #     torch.nn.Conv3d(6, 32, 3, padding=1),
    #     torch.nn.MaxPool3d(2),
    #     torch.nn.Conv3d(32, 64, 3, padding=1),
    #     torch.nn.MaxPool3d(2),
    #     torch.nn.Conv3d(64, 128, 3, padding=1),
    #     torch.nn.MaxPool3d(2),
    #     torch.nn.Flatten(),
    #     torch.nn.Linear(128*3*3*3, 128),
    #     torch.nn.ReLU(),
    #     torch.nn.Linear(128, 1)

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
