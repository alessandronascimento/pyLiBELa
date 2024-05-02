import torch
import torch.nn as nn
from math import floor
from typing import OrderedDict, Literal, Optional
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
import numpy as np

path_to_pdbbind = "/data/home/alessandro/PDBbind_v2020"
LOG_DIR = "logs"

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

    with open(f'{path_to_pdbbind}/index/INDEX_general_PL_data.2020', 'r') as binding_file:
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
  """
  Reads a binary file, decodes it as float64, reshapes, and normalizes the data.

  Args:
      file_path (str): Path to the binary file.

  Returns:
      torch.Tensor: The processed and normalized data.
  """
  with open(file_path, 'rb') as f:
    # Read the entire file content as a byte array
    data_bytes = f.read()
  # Convert the byte array to a PyTorch tensor of type float64
  data = torch.from_numpy(np.frombuffer(data_bytes, dtype=np.float64))

  # Reshape the data based on the expected format (channels, width, height, depth)
  data = data.reshape(3, 60, 60, 60)
  # Normalize the data (assuming you want values between 0 and 1)
  data_norm = (data - data.min()) / (data.max() - data.min())
  return data_norm

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
            torch.cat((self._rec, self._lig), dim=0).to(torch.float32).to(self.device),
            torch.tensor(self.scores[grid_idx]).float().to(self.device),
            torch.tensor(weight).float().to(self.device)
        )

        # return torch.cat((self._rec, self._lig), dim=0).to(self.device), self.scores[grid_idx], weight
        return out_tuple


def data_generator(name_list, score_list):
  weight_normal = 11 / 2
  weight_decoy = 11 / 20

  for i, name in enumerate(name_list):
    target = name #.decode()
    file_path = f"{path_to_pdbbind}/targets/{target}/grid_30_0.5_SF0"

    try:
      grid1 = process_file(f"{file_path}/McGrid_rec.grid")
      grid2 = process_file(f"{file_path}/McGrid_lig.grid")
    except IOError:
      # Handle file not found exceptions
      continue

    observable = score_list[i]
    combined_grid = torch.cat([grid1, grid2], dim=0)
    yield combined_grid, observable, weight_normal

    # Yield data for 10 related decoys
    for j in range(1, 11):
      grid2 = process_file(f"{file_path}/McGrid_dec_{j}.grid")
      observable = 0.00
      combined_grid = torch.cat([grid1, grid2], dim=0)
      yield combined_grid, observable, weight_decoy

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


class CNNModel(nn.Module):
  def __init__(self, in_channels=6, num_classes=1):  # Specify input and output channels
    super(CNNModel, self).__init__()
    self.model= nn.Sequential(
	#Conv layer 1    
        nn.Conv3d(in_channels, 32, kernel_size=(5, 5, 5), padding=1),
        nn.ReLU(),
        nn.MaxPool3d((2, 2, 2)), 
	# output shape = (32, 28, 28, 28)

	#Conv layer 2  
        nn.Conv3d(32, 64, kernel_size=(3, 3, 3), padding=1),
        nn.ReLU(),
#        nn.MaxPool3d((2, 2, 2)),
	# output shape = (64, 26, 26, 26)

	#Conv layer 3  
        nn.Conv3d(64, 128, kernel_size=(3, 3, 3), padding=1),
        nn.ReLU(),
        nn.MaxPool3d((2, 2, 2)),
	# output shape = (128, 12, 12, 12)
	
	#Conv layer 4  
        nn.Conv3d(128, 256, kernel_size=(3, 3, 3), padding=1),
        nn.ReLU(),
#        nn.MaxPool3d((2, 2, 2)),
	# output shape = (256, 10, 10, 10)


	#Conv layer 5  
        nn.Conv3d(256, 256, kernel_size=(3, 3, 3), padding=1),
        nn.ReLU(),
        nn.MaxPool3d((2, 2, 2)),
	#output shape =( 256, 4, 4, 4)
	
	## Flatten
	nn.Flatten(),

	# Dense Layer 1
        nn.Linear(87808, 128), 

	#Relu
        nn.ReLU(),
        #BatchNorm1d
        nn.BatchNorm1d(128),
        #Dropout
        nn.Dropout(p=0.15),
        #Linear 2
        nn.Linear(128, 1)
    )
  def forward(self, x):
        # Set 1
        out = self.model(x)
        return out


def train_model(model, train_loader, val_loader, epochs=30, learning_rate=1e-4):
  """
  Trains the PyTorch model.

  Args:
      model: The PyTorch model to train.
      train_loader: DataLoader for the training data.
      val_loader: DataLoader for the validation data.
      epochs: Number of training epochs (default: 30).
      learning_rate: Learning rate for the optimizer (default: 1e-4).

  Returns:
      A tuple containing the training history (fit_out) and the evaluation results (eval_out).
  """

  # Define loss function and optimizer
  criterion = nn.MSELoss()  # Mean squared error loss
  optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

  device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

  # Training loop
  fit_out = []
  for epoch in range(epochs):
    print(f"\nEpoch {epoch+1}/{epochs}")

    # Training phase
    model.train()
    train_loss = 0.0
    for data, target, weight in train_loader:
      optimizer.zero_grad()
      output = model(data)
      loss = criterion(output, target)
      loss.backward()
      optimizer.step()
      train_loss += loss.item() * weight.sum().item()

    # Validation phase
    model.eval()
    with torch.no_grad():
      val_loss = 0.0
      for data, target, _ in val_loader:  # Ignore weight during validation
        output = model(data)
        val_loss += criterion(output, target).item()

    # Print training and validation loss
    avg_train_loss = train_loss / len(train_loader)
    avg_val_loss = val_loss / len(val_loader)
    print(f"[Train]: Loss: {avg_train_loss:.4f}")
    print(f"[Val]: Loss: {avg_val_loss:.4f}")

    fit_out.append((avg_train_loss, avg_val_loss))

  # Evaluation on test seit (optional)
  # ... (similar logic using test_loader)

  return fit_out #, eval_out  # Replace eval_out with actual evaluation on test seti

def main():
  # Set default device (if CUDA is available)
  device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
  batch_size=5
  print(f"Using device: {device}")

  # Create and split datasets
  print('Creating and splitting the datasets...')
  train_loader, test_loader, val_loader = create_datasets(device, batch_size)  # Assuming you have a test loader now


  # Define and instantiate the model
  print('Creating the model')
  model = CNNModel().to(device)  # Move the model to the chosen device

  # Train the model
  print('Starting training...')
  fit_out = train_model(model, train_loader, val_loader)

  # ... (optional: evaluate on test set using test_loader)

  print("Training completed!")

if __name__ == "__main__":
  main()
