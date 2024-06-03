import sys
import tensorflow as tf
import numpy as np
import torch
from sklearn.model_selection import train_test_split
from typing import Literal, Optional
from diffusers import DDPMScheduler, UNet3DConditionModel

sys.path.append('/data/home/alessandro/workspace/pyLiBELa')
from pyPARSER import *
from pyMol2 import *

path_to_pdbbind = "/data/home/alessandro/PDBbind_v2020"
LOG_DIR = "logs"

# Reading PDBBind data
def read_pdbbind_data():
    targets_ok = []
    with open(f'{path_to_pdbbind}/scripts/targets_ok.dat', 'r') as target_ok_file:
        for line in target_ok_file:
            targets_ok.append(line.split()[0])

    binding_targets = []
    binding_years = []
    binding_resolution = []
    binding_score = []

    with open(f'{path_to_pdbbind}/index/INDEX_refined_data.2020', 'r') as binding_file:
        for line in binding_file:
            line2 = line.split()
            if line2[0][0] == '#' or line2[0] not in targets_ok:
                continue
            binding_targets.append(line2[0])
            binding_resolution.append(0.00 if line2[1] == 'NMR' else float(line2[1]))
            binding_years.append(int(line2[2]))
            binding_score.append(float(line2[3]))

    print('Binding data found for %d valid targets' % len(binding_targets))
    return binding_targets, binding_score

def process_grid_file(file_path):
    with open(file_path, 'rb') as f:
        data_bytes = f.read()
    data = np.frombuffer(data_bytes, dtype=np.float64).reshape(3, 60, 60, 60)
    data_norm = (data - data.min()) / (data.max() - data.min())
    return tf.convert_to_tensor(data_norm, dtype=tf.float32)

def get_ligand_xyz_charges(file_path):
    Input = PARSER()
    lig = Mol2(Input, file_path)
    xyz = np.array(lig.xyz)
    return tf.convert_to_tensor(xyz, dtype=tf.float32)

def dynamic_collate_fn(batch):
    rec_grids, ligand_coordinates, scores = zip(*batch)
    rec_grids = tf.stack(rec_grids)

    total_atoms = sum(coords.shape[0] for coords in ligand_coordinates)
    padded_coords = tf.zeros((len(batch), total_atoms, 3), dtype=tf.float32)

    start_idx = 0
    for i, coords in enumerate(ligand_coordinates):
        end_idx = start_idx + coords.shape[0]
        padded_coords = tf.tensor_scatter_nd_update(padded_coords, [[i, start_idx]], [coords])
        start_idx = end_idx

    scores = tf.stack(scores)
    return rec_grids, padded_coords, scores

class GridDataset(tf.data.Dataset):
    def __new__(cls, targets_list, scores_list):
        def generator():
            for target, score in zip(targets_list, scores_list):
                path = f"{path_to_pdbbind}/targets/{target}/grid_30_0.5_SF0"
                path_mol2 = f"{path_to_pdbbind}/targets/{target}/mol2"
                rec_grid = process_grid_file(f"{path}/McGrid_rec.grid")
                ligand_coordinates = get_ligand_xyz_charges(f"{path_mol2}/{target}_lig.mol2.gz")
                yield rec_grid, ligand_coordinates, tf.convert_to_tensor(score, dtype=tf.float32)
        return tf.data.Dataset.from_generator(generator, output_signature=(
            tf.TensorSpec(shape=(3, 60, 60, 60), dtype=tf.float32),
            tf.TensorSpec(shape=(None, 3), dtype=tf.float32),
            tf.TensorSpec(shape=(), dtype=tf.float32)
        ))

def create_datasets(device: Literal['gpu', 'cpu'], batch_size: int, max_targets: Optional[int] = None):
    binding_targets, binding_score = read_pdbbind_data()

    train_targets, test_targets, train_scores, test_scores = train_test_split(binding_targets[:max_targets],
                                                                              binding_score[:max_targets],
                                                                              train_size=0.8, shuffle=True)
    train_targets, valid_targets, train_scores, valid_scores = train_test_split(train_targets, train_scores, train_size=0.8)

    train_dataset = GridDataset(train_targets, train_scores).batch(batch_size).prefetch(tf.data.AUTOTUNE)
    test_dataset = GridDataset(test_targets, test_scores).batch(batch_size).prefetch(tf.data.AUTOTUNE)
    valid_dataset = GridDataset(valid_targets, valid_scores).batch(batch_size).prefetch(tf.data.AUTOTUNE)

    return train_dataset, test_dataset, valid_dataset

class CustomUNet3DConditionModel(tf.keras.Model):
    def __init__(self, in_channels, out_channels):
        super(CustomUNet3DConditionModel, self).__init__()

        self.ligand_layers = tf.keras.Sequential([
            tf.keras.layers.Dense(32, activation='relu'),
            tf.keras.layers.Dense(64, activation='relu'),
            tf.keras.layers.Dense(128, activation='relu')
        ])

        self.receptor_layers = tf.keras.Sequential([
            tf.keras.layers.Conv3D(16, kernel_size=3, padding='same', activation='relu'),
            tf.keras.layers.Conv3D(32, kernel_size=3, padding='same', activation='relu'),
            tf.keras.layers.Conv3D(64, kernel_size=3, padding='same', activation='relu'),
            tf.keras.layers.Conv3D(128, kernel_size=3, padding='same', activation='relu')
        ])

        self.pool = tf.keras.layers.GlobalAveragePooling3D()
        self.reduce_dim = tf.keras.layers.Dense(128)
        self.reduce_dim_320 = tf.keras.layers.Dense(320)

        self.unet = UNet3DConditionModel(in_channels=256, out_channels=out_channels)


    def call(self, ligand_coords, rec_grids, timestep):
        batch_size = tf.shape(ligand_coords)[0]
        num_atoms = tf.shape(ligand_coords)[1]

        ligand_coords_flat = tf.reshape(ligand_coords, [batch_size * num_atoms, -1])
        ligand_features = self.ligand_layers(ligand_coords_flat)
        ligand_features = tf.reduce_mean(tf.reshape(ligand_features, [batch_size, num_atoms, -1]), axis=1)

        receptor_features = self.receptor_layers(rec_grids)
        pooled_receptor_features = self.pool(receptor_features)

        receptor_features_reduced = self.reduce_dim(pooled_receptor_features)
#        print('receptor_features_reduced shape:', tf.shape(receptor_features_reduced))
        combined_features = tf.concat([ligand_features, receptor_features_reduced], axis=1)

        combined_features = tf.reshape(combined_features, [batch_size, 256, 1, 1, 1])
        combined_features = tf.tile(combined_features, [1, 1, 60, 60, 60])

        combined_features_torch = torch.tensor(combined_features.numpy())

        # Ensure rec_grids_torch matches the expected shape
        rec_grids_torch = self.reduce_dim_320(pooled_receptor_features)
        rec_grids_torch = tf.reshape(rec_grids_torch, [batch_size, 320, 1, 1, 1])
        rec_grids_torch = tf.tile(rec_grids_torch, [1, 1, 60, 60, 60])
        rec_grids_torch = torch.tensor(rec_grids_torch.numpy())
        
        timestep_torch = torch.tensor(timestep.numpy())


        return self.unet(combined_features_torch, timestep=timestep_torch, encoder_hidden_states=rec_grids_torch)



def add_noise(x, t, noise_scheduler):
    noise = tf.random.normal(tf.shape(x))
    alpha_t = tf.reshape(tf.gather(noise_scheduler.alphas_cumprod, t), [-1, 1, 1, 1])
    noisy_x = tf.sqrt(alpha_t) * x + tf.sqrt(1 - alpha_t) * noise
    return noisy_x, noise

def train_denoising_model(train_loader, model, noise_scheduler, optimizer, epochs=30, device="gpu", accumulation_steps=4):
    loss_fn = tf.keras.losses.MeanSquaredError()
    for epoch in range(epochs):
        for i, batch in enumerate(train_loader):
            rec_grids, ligand_coords, scores = batch

            with tf.GradientTape() as tape:
                t = tf.random.uniform([tf.shape(ligand_coords)[0]], minval=0, maxval=noise_scheduler.num_train_timesteps, dtype=tf.int32)
                noisy_coords, noise = add_noise(ligand_coords, t, noise_scheduler)

                # Debugging
                print("t (timestep) type:", type(t), "dtype:", t.dtype)
                print("noisy_coords type:", type(noisy_coords), "shape:", noisy_coords.shape)
                
                predicted_diffusion = model(noisy_coords, rec_grids, timestep=t)
                loss = loss_fn(predicted_diffusion, ligand_coords) / accumulation_steps

            gradients = tape.gradient(loss, model.trainable_variables)
            optimizer.apply_gradients(zip(gradients, model.trainable_variables))

            if (i + 1) % accumulation_steps == 0:
                print(f"Epoch {epoch + 1}, Step {i + 1}, Loss: {loss.numpy()}")

def main():
    device = "gpu" if tf.config.list_physical_devices('GPU') else "cpu"
    batch_size = 1
    print(f"Using device: {device}")

    train_loader, test_loader, val_loader = create_datasets(device, batch_size)
    print('Creating and splitting the datasets...')

    model = CustomUNet3DConditionModel(in_channels=3, out_channels=3)
    noise_scheduler = DDPMScheduler(num_train_timesteps=1000)
    optimizer = tf.keras.optimizers.Adam(learning_rate=1e-4)
    print('Creating the model...')

    print('Starting training...')
    train_denoising_model(train_loader, model, noise_scheduler, optimizer, epochs=10, device=device)
    print("Training completed!")

if __name__ == "__main__":
    main()
