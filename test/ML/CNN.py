from data.grid_loader import GridLoader
from data.target_parser import Parser
from typing import Any
from pathlib import Path
import numpy as np
import tensorflow as tf

GRID_CHANNELS = 4
GRID_DIM = 20
SEED = 2023
BATCH = 16

def load_data() -> tuple[np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]], int]:

    # directory containing ligand and receptor grids and target subdirectories
    root_path : Path = Path("/home/alexsouza/workingproj")

    shape = (GRID_CHANNELS, GRID_DIM, GRID_DIM, GRID_DIM)

    receptor_grids = GridLoader(root_path / "rec_grids", shape)
    rec_arr, rec_grid_names = receptor_grids.load_grids()

    ligand_grids = GridLoader( root_path / "lig_grids", shape)
    lig_arr, lig_grid_names = ligand_grids.load_grids()

    valid = [name for name in rec_grid_names if name in lig_grid_names]
    total = len(valid)
    target_dict = Parser.parse_index_file(root_path / "pyLiBELa" / "test" / "ML"/ "INDEX_refined_set.2020")

    merged_grids = np.empty((total, 2, shape[0], shape[1], shape[2], shape[3])) 
    targets = np.empty(total)

    for i, name in enumerate(valid):
        merged_grids[i][0] = rec_arr[rec_grid_names.index(name)]
        merged_grids[i][1] = lig_arr[lig_grid_names.index(name)]
        targets[i] = -np.log(target_dict[name])

    merged_grids = merged_grids.reshape((total, 2*shape[0], shape[1], shape[2], shape[3]))
    return merged_grids, targets, total

if __name__ == "__main__":
    grids, targets, total = load_data() 
    grids, targets = tf.convert_to_tensor(grids), tf.convert_to_tensor(targets)

    grids_norm_layer = tf.keras.layers.Normalization(axis=(1))
    grids_norm_layer.adapt(grids)
    grids = grids_norm_layer(grids)

    targets_norm_layer = tf.keras.layers.Normalization(axis=None)
    targets_revert_layer = tf.keras.layers.Normalization(axis=None, invert=True)
    targets_norm_layer.adapt(targets)
    targets_scaled = targets_norm_layer(targets)
    targets_revert_layer.adapt(targets)

    dataset = tf.data.Dataset.from_tensor_slices((grids, targets_scaled))
    
    train_val, test = tf.keras.utils.split_dataset(dataset, left_size=0.8, right_size=0.2,
                                                   shuffle=True, seed=SEED)
    train, valid = tf.keras.utils.split_dataset(train_val, left_size=0.8, right_size=0.2)

    cnn_model = tf.keras.Sequential([
        tf.keras.layers.Conv3D(64, 5, activation='relu', data_format='channels_first',
                               input_shape=grids.shape[1:], padding="same"),
        tf.keras.layers.MaxPool3D(data_format='channels_first'),
        tf.keras.layers.Conv3D(128, 3, activation='relu', data_format='channels_first',
                               padding="same"),
        tf.keras.layers.Conv3D(128, 3, activation='relu', data_format='channels_first',
                               padding="same"),
        tf.keras.layers.MaxPool3D(data_format='channels_first'),
        tf.keras.layers.Conv3D(256, 3, activation='relu', data_format='channels_first',
                               padding="same"),
        tf.keras.layers.Conv3D(256, 3, activation='relu', data_format='channels_first',
                               padding="same"),
        tf.keras.layers.Flatten(),
        tf.keras.layers.Dense(units=128, activation="relu"),
        tf.keras.layers.Dense(units=64, activation="relu"),
        tf.keras.layers.Dense(units=1),
    ])

    cnn_model.compile(loss="mse", 
                      optimizer=tf.keras.optimizers.Adam(learning_rate=1e-4), 
                      metrics=["RootMeanSquaredError"])

    out_info = cnn_model.fit(x=train.batch(BATCH), epochs=10, validation_data=valid.batch(BATCH))
    mse_test, rmse_test = cnn_model.evaluate(test.batch(BATCH))

    print(f"mse: {mse_test}")
    print(f"rmse_test: {rmse_test}")
