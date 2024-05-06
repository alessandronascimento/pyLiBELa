import torch
import numpy as np
import wandb
from math import floor
from typing import OrderedDict, Literal, Optional
from torch.utils.data import Dataset, DataLoader
from torch.optim import Optimizer
from torch import nn
from pathlib import Path
from sklearn.model_selection import train_test_split

path_to_pdbbind = "/data/home/alessandro/PDBbind_v2020"
LOG_DIR = "./log"
USE_WANDB = True 

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
        data = torch.tensor(np.frombuffer(f.read(), dtype=np.float64).reshape((3,60,60,60)), dtype=torch.float32)
        
    return data


class GridDataset(Dataset):
    def __init__(self, targets_list, scores_list, device) -> None:
        self.targets = targets_list
        self.scores = scores_list 
        self.device = device
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

        out_tuple = (
            torch.cat((self._rec, self._lig), dim=0).to(self.device),
            torch.tensor(self.scores[grid_idx]).to(self.device),
            torch.tensor(weight).to(self.device)
        )

        # return torch.cat((self._rec, self._lig), dim=0).to(self.device), self.scores[grid_idx], weight
        return out_tuple

#TODO look into with num_workers and prefetch
def create_datasets(device : Literal['mps', 'cuda', 'cpu'],
                    batch_size : int,
                    max_targets : Optional[int] = None):
    '''
    Returns Dataloaders for train dataset, test dataset 
    and validation dataset respectively
    '''
    
    binding_targets, binding_score = read_pdbbind_data()

    train_targets, test_targets, train_scores, test_scores = train_test_split(binding_targets[:max_targets],
                                                                          binding_score[:max_targets],
                                                                          train_size=0.8,
                                                                          shuffle=True)
    train_targets, valid_targets, train_scores, valid_scores = train_test_split(train_targets, train_scores, train_size=0.8)

    train_dataset = DataLoader(GridDataset(train_targets, train_scores, device),
                               batch_size=batch_size)
    test_dataset = DataLoader(GridDataset(test_targets, test_scores, device),
                               batch_size=batch_size)
    valid_dataset = DataLoader(GridDataset(valid_targets, valid_scores, device),
                               batch_size=batch_size)

    return train_dataset, test_dataset, valid_dataset

def get_device():
    device = ("cuda" if torch.cuda.is_available() 
        else "mps" if torch.backends.mps.is_available()
        else "cpu")

    return device

class CNNModel(nn.Module):

    def __init__(self, input_shape : tuple[int, int, int, int]) -> None:
        super().__init__()
        self.conv_stack = nn.Sequential(OrderedDict([
            ('conv1',    nn.Conv3d(in_channels=input_shape[0], out_channels=32, kernel_size=3)),
            ('relu1',    nn.ReLU()),
            ('maxpool1', nn.MaxPool3d(kernel_size=2)),

            ('conv2',    nn.Conv3d(in_channels=32, out_channels=64, kernel_size=3)),
            ('relu2',    nn.ReLU()),

            ('conv3', nn. Conv3d(in_channels=64, out_channels=64, kernel_size=2)),
            ('relu3' , nn.ReLU()),
            ('maxpool2', nn.MaxPool3d(kernel_size=2)),

            ('conv4', nn.Conv3d(in_channels=64, out_channels=128, kernel_size=2)),
            ('relu4', nn.ReLU()),

            ('conv5', nn.Conv3d(in_channels=128, out_channels=128, kernel_size=2)),
            ('relu5', nn.ReLU()),

            ('flat', nn.Flatten()),
            ]))

        self.linear_size = self.conv_stack(torch.zeros([1]+list(input_shape))).shape[1]

        self.dense_stack = nn.Sequential(OrderedDict([
            ('linear1', nn.Linear(in_features=self.linear_size, out_features=128)),
            ('activ1', nn.ReLU()),

            ('batch_norm', nn.BatchNorm1d(128)),
            ('dropout', nn.Dropout(0.5)),

            ('linear2', nn.Linear(in_features=128, out_features=64)),
            ('activ2', nn.ReLU()),

            ('out', nn.Linear(in_features=64, out_features=1))
            ]))
    
    def forward(self, x):
        x = self.conv_stack(x)
        x = self.dense_stack(x)
        x = torch.flatten(x)
        return x


def train_model(model : nn.Module, 
                train_data: DataLoader, 
                valid_data: DataLoader,
                optimizer : Optimizer, 
                epochs : int,
                loss_fn):
    
    if USE_WANDB: wandb.init(
        # set the wandb project where this run will be logged
        project="Torch-grids",

        # track hyperparameters and run metadata
        config={
        "learning_rate": optimizer.state_dict()['param_groups'][0]['lr'],
        "architecture": "CNN",
        "dataset": "PDBbind",
        "epochs": epochs,
        }
    )
    for epoch in range(epochs):
        # TRAINING STEP -----------------------------
        running_loss = 0
        print(f"\nEPOCH {epoch+1}/{epochs}")
        model.train()
        pad = 15
        total_train = len(train_data)
        for i, (X, y, w) in enumerate(train_data):

            # progress bar
            pad_left = int(i*pad/total_train)
            print(f"\rRunning Epoch... |{'':{'='}<{pad_left}}>{'':{' '}>{pad-pad_left-1}}|", end='')

            pred = model(X)
            loss = loss_fn(pred, y, w)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            running_loss += loss.item()

        loss = running_loss/total_train
        if USE_WANDB: wandb.log({"train loss": loss})
        print('\n')
        print(f"Train Loss: {loss:.4f}")


        # VALIDATION STEP ---------------------------
        model.eval()
        with torch.no_grad():
            running_vloss = 0 
            for Xv, yv, wv in valid_data:
                vpred = model(Xv)
                vloss = loss_fn(vpred, yv, wv)
                running_vloss += vloss.item()

        vloss = running_loss/len(valid_data)
        if USE_WANDB: wandb.log({"valid loss": vloss})
        print(f"Validation Loss: {vloss:.4f}")

    if USE_WANDB: wandb.finish()

def WeightedMSELoss(x : torch.Tensor, y : torch.Tensor, w : torch.Tensor):
    l = (x - y)*(x - y) * w
    return l.sum()/l.size(dim=0)


def test():

    #data = process_file(f'{path_to_pdbbind}/targets/8gpb/grid_30_0.5_SF0/McGrid_lig.grid')
    #print(data)
    # dataset = GridDataset(targets[:23], scores[:23])


    # a = torch.tensor([1, 2, 3, 4], dtype=torch.float32)
    # b = torch.tensor([4, 3, 2, 1], dtype=torch.float32)
    # w = torch.tensor([1, 1, 1 ,1], dtype=torch.float32)
    # 
    # should give the same results with w vector set as ones
    # print(nn.MSELoss()(a, b))
    # print(WeightedMSELoss(a, b, w))

    device = get_device()
    cnn = CNNModel((6,60,60,60)).to(device)
    learning_rate = 1e-3
    batch_size = 20
    epochs = 20
    loss_fn = WeightedMSELoss 
    optimizer = torch.optim.Adam(cnn.parameters(), lr=learning_rate)
    train, _, valid = create_datasets(device, batch_size, max_targets=None)
    train_model(cnn, train, valid, optimizer, epochs, loss_fn)


if __name__ == '__main__':
    test()
