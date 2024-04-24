import torch
import numpy as np
from math import floor
from typing import OrderedDict, Literal, Optional
from torch.utils.data import Dataset, DataLoader
from torch.optim import Optimizer
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
        self.flatten = nn.Flatten()
        self._current_shape = input_shape 
        self._old_out_channels = 0
        self.model_stack = nn.Sequential(OrderedDict([
            ('conv1',    self.make_conv3d(out_channels=32,
                                        kernel_size=3,
                                        in_channels=input_shape[0])),
            ('maxpool1', self.make_maxpool3d(kernel_size=2)),
            ('conv2',    self.make_conv3d(out_channels=64, kernel_size=3)),
            ('conv3',    self.make_conv3d(out_channels=64, kernel_size=2)),
            ('maxpool2', self.make_maxpool3d(kernel_size=2)),
            ('conv4',    self.make_conv3d(out_channels=128, kernel_size=2)),
            ('conv5',    self.make_conv3d(out_channels=128, kernel_size=2)),
            ('flat',     self.make_flatten()),
            ('linear1',  self.make_linear(units=8)),
            ('activ1', nn.ReLU()),
            ('linear2',  self.make_linear(units=4)),
            ('activ2', nn.ReLU()),
            ('out',      self.make_linear(units=1))
            ]))

    def forward(self, x):
        x = self.model_stack(x)
        x = torch.flatten(x)
        return x

    def make_conv3d(self,
                  out_channels,
                  kernel_size,
                  in_channels=None,
                  padding=0,
                  stride=1,
                  dilation=1) -> nn.Conv3d:
        """
        Returns a Conv3d layer handling input shape automatically 
        """

        if in_channels is None: in_channels = self._old_out_channels
        conv = nn.Conv3d(in_channels=in_channels, 
                         out_channels=out_channels, 
                         kernel_size=kernel_size,
                         padding=padding,
                         stride=stride,
                         dilation=dilation)

        if not isinstance(self._current_shape, tuple): raise ValueError("Current shape must be multidimensional")

        dims = [self._current_shape[1], self._current_shape[2], self._current_shape[3]]
        out = self._out_dims(dims, 
                              kernel_size=kernel_size,
                              padding=padding,
                              stride=stride,
                              dilation=dilation)

        self._old_out_channels = out_channels
        self._current_shape = (out_channels, out[0], out[1], out[2])

        return conv

    def make_maxpool3d(self,
                       kernel_size,
                       stride=1,
                       padding=0,
                       dilation=1) -> nn.MaxPool3d:
        """
        Returns a MaxPool3d layer handling input shape automatically 
        """

        if not isinstance(self._current_shape, tuple): raise ValueError("Current shape must be multidimensional")

        maxpool = nn.MaxPool3d(kernel_size=kernel_size,
                               stride=stride,
                               padding=padding,
                               dilation=dilation)

        dims = [self._current_shape[1], self._current_shape[2], self._current_shape[3]]
        out = self._out_dims(dims, 
                              kernel_size=kernel_size,
                              padding=padding,
                              stride=stride,
                              dilation=dilation)
        self._current_shape = (self._old_out_channels, out[0], out[1], out[2])
        return maxpool 

    def make_flatten(self) -> nn.Flatten:
        """
        Returns a Flatten layer handling shape automatically 
        """
        if not isinstance(self._current_shape, tuple): raise ValueError("Current shape must be multidimensional")

        self._current_shape = self._current_shape[0] * self._current_shape[1] * self._current_shape[2] * self._current_shape[3]
        return nn.Flatten()

    def make_linear(self, units) -> nn.Linear:
        """
        Returns a Linear layer handling input shape automatically 
        """

        if not isinstance(self._current_shape, int): raise ValueError("Shape must be flat")

        lin = nn.Linear(in_features=self._current_shape, out_features=units)

        self._current_shape = units
        return lin

    def _make_uniform(self, value : int | list[int]) -> tuple[list[int], list[int]]:
        if isinstance(value, int):
            return [value,], [0, 0, 0]
        else: return value, [0, 1, 2]


    def _out_dims(self, dims : list[int],
                   kernel_size : int | list[int],
                   padding : int | list[int] = 0,
                   stride : int | list[int] = 1,
                   dilation : int | list[int] = 1) -> tuple[int, int, int]:
        """
        Helper function that returns a tuple with the shape 
        information (excluding channels) as outputted by a Conv3d 
        or MaxPool3d layer to be used by next layer.
        """

        kernel_size, kernel_shape = self._make_uniform(kernel_size)
        padding, padding_shape = self._make_uniform(padding)
        stride, stride_shape = self._make_uniform(stride)
        dilation, dilation_shape = self._make_uniform(dilation) 
       
        def reshape(i : int) -> int:
            return floor(1 + (dims[i] + 2*padding[padding_shape[i]] - 
            dilation[dilation_shape[i]]*( kernel_size[kernel_shape[i]] - 1) - 1 ) / stride[stride_shape[i]])

        return (reshape(0), reshape(1), reshape(2))

def train_model(model : nn.Module, 
                train_data: DataLoader, 
                valid_data: DataLoader,
                optimizer : Optimizer, 
                epochs : int,
                loss_fn : nn.MSELoss):

    batch_size = train_data.batch_size 
    if (batch_size := train_data.batch_size) is None:
        raise ValueError("Training DataLoader must have a batch size")

    for epoch in range(epochs):
        # TRAINING STEP -----------------------------
        avg_loss = 0
        running_loss = 0
        print(f"EPOCH {epoch+1}")
        model.train()
        for batch, (X, y, w) in enumerate(train_data):
            pred = model(X)
            loss = loss_fn(pred, y)
            # loss = loss * w.mean() # sample weight
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            running_loss += loss.item()

            if batch % 100 == 0:
                avg_loss = running_loss / batch_size # loss per batch
                print(f"batch #{batch} loss: {avg_loss:.4f}\r")
                running_loss = 0

        # VALIDATION STEP ---------------------------
        running_vloss = 0 
        model.eval()
        with torch.no_grad():
            for Xv, yv, wv in valid_data:
                vpred = model(Xv)
                vloss = loss_fn(vpred, yv)
                # vloss *= wv

                running_vloss += vloss.item()

        size = len(valid_data)
        avg_vloss = running_vloss / (size + 1)
        print(f"LOSS train {avg_loss} valid {avg_vloss}")


def test():

    #data = process_file(f'{path_to_pdbbind}/targets/8gpb/grid_30_0.5_SF0/McGrid_lig.grid')
    #print(data)
    # dataset = GridDataset(targets[:23], scores[:23])

    device = get_device()
    cnn = CNNModel((6,60,60,60)).to(device)
    learning_rate = 1e-3
    batch_size = 5
    epochs = 10
    loss_fn = nn.MSELoss()
    optimizer = torch.optim.Adam(cnn.parameters(), lr=learning_rate)
    train, _, valid = create_datasets(device, 1, max_targets=100)
    train_model(cnn, train, valid, optimizer, epochs, loss_fn)


if __name__ == '__main__':
    test()
