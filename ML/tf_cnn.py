import tensorflow as tf
# import keras_tuner as kt
from pathlib import Path
from sklearn.model_selection import train_test_split
from time import strftime

import keras.api._v2.keras as keras

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

path_to_pdbbind = "/data/home/alessandro/PDBbind_v2020"
LOG_DIR = "logs"

import logging

# get TF logger
log = logging.getLogger('tensorflow')
log.setLevel(logging.DEBUG)

# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# create file handler which logs even debug messages
fh = logging.FileHandler(LOG_DIR + '/tf_cnn.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
log.addHandler(fh)

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

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

#@title Auxiliary functions
def process_file(file_path):

  as_path = Path(file_path)
  if not as_path.exists():
    raise IOError(f"File path '{as_path}' does not exist")

  data = tf.io.read_file(file_path)
  data = tf.io.decode_raw(data, tf.float64)
  data = tf.reshape(data, (3, 60, 60, 60))
  data_norm = tf.keras.utils.normalize(data)
  return data_norm

def _log_file_not_found(file_path):
  with open("log_files_not_found", 'a') as f:
    f.write(f'{file_path}\n')


def data_generator(name_list, score_list):
  #ntargets=20
  weight_normal = 11/2
  weight_decoy = 11/20

  for i, name in enumerate(name_list):

    target = name.decode()
    file_path = f'{path_to_pdbbind}/targets/{target}/grid_30_0.5_SF0'

    try:
      grid1 = process_file(file_path + '/McGrid_rec.grid')
      grid2 = process_file(file_path + '/McGrid_lig.grid')

    except IOError:
      _log_file_not_found(file_path)
      continue

    observable = score_list[i]
    combined_grid = tf.concat([grid1, grid2], axis=0)

    # Yield the data as a tuple
    yield combined_grid, observable, weight_normal

    # Yields the 10 related decoys of the current grid
    for j in range(1,11):

      grid2 = process_file(file_path + '/McGrid_dec_{}.grid'.format(j))

      observable = 0.00
      combined_grid = tf.concat([grid1, grid2], axis=0)
      yield combined_grid, observable, weight_decoy

def create_datasets(max_targets=None):

    binding_targets, binding_score = read_pdbbind_data()

    train_names, test_names, train_scores, test_scores = train_test_split(binding_targets[:max_targets],
                                                                          binding_score[:max_targets],
                                                                          train_size=0.8,
                                                                          shuffle=True)
    train_names, valid_names, train_scores, valid_scores = train_test_split(train_names, train_scores, train_size=0.8)

    output_signature = (tf.TensorSpec(shape=(6, 60, 60, 60), dtype=tf.float64, name='combgrid'),
                        tf.TensorSpec(shape=(), dtype=tf.float32, name='ligscore'),
                        tf.TensorSpec(shape=(), dtype=tf.float32, name='weight'))
    batch_size = 5
    prefetch_size = 1

# Create the dataset from the generator function,
# with batch and prefetch sizes already determined
    train_dataset = tf.data.Dataset.from_generator(
        data_generator,
        output_signature=output_signature,
        args=(tf.convert_to_tensor(train_names, dtype=tf.string), 
              tf.convert_to_tensor(train_scores, dtype=tf.float32)),
        name="train_dataset_gen"
    ).batch(batch_size).prefetch(prefetch_size)

    test_dataset = tf.data.Dataset.from_generator(
        data_generator,
        output_signature=output_signature,
        args=(tf.convert_to_tensor(test_names, dtype=tf.string), 
              tf.convert_to_tensor(test_scores, dtype=tf.float32)),
        name="test_dataset_gen"
    ).batch(batch_size).prefetch(prefetch_size)

    valid_dataset = tf.data.Dataset.from_generator(
        data_generator,
        output_signature=output_signature,
        args=(tf.convert_to_tensor(valid_names, dtype=tf.string),
              tf.convert_to_tensor(valid_scores, dtype=tf.float32)),
        name="valid_dataset_gen"
    ).batch(batch_size).prefetch(prefetch_size)

    return train_dataset, test_dataset, valid_dataset

#@title Alex's model
def create_model():

    cnn_model = keras.Sequential([
            keras.layers.Conv3D(64, 5, activation='relu', data_format='channels_first', input_shape=(6,60,60,60), padding="same", name='conv1'),
            keras.layers.MaxPool3D(data_format='channels_first', name='maxpool1'),
            keras.layers.Conv3D(128, 3, activation='relu', data_format='channels_first', padding="same", name='conv2'),
            keras.layers.Conv3D(128, 3, activation='relu', data_format='channels_first', padding="same", name='conv3'),
            keras.layers.MaxPool3D(data_format='channels_first', name="maxpool2"),
            keras.layers.Conv3D(256, 3, activation='relu', data_format='channels_first', padding="same", name='conv4'),
            keras.layers.Conv3D(256, 3, activation='relu', data_format='channels_first', padding="same", name='conv5'),
            keras.layers.Flatten(name='flatten'),
            keras.layers.Dense(units=128, activation="relu", name='dense1'),
            keras.layers.Dense(units=128, activation="relu", name='dense2'),
            keras.layers.Dense(units=1, name='out_dense'),
        ])


    cnn_model.compile(loss="mse", 
                      optimizer=keras.optimizers.Adam(learning_rate=1e-4), 
                      metrics=[keras.metrics.MeanSquaredError()], 
                      weighted_metrics=[])

    return cnn_model


def tune_model():

    trials=2
    tuner = kt.BayesianOptimization(create_model,
                                    objective='mean_squared_error',
                                    max_trials=trials,
                                    directory=LOG_DIR,
                                    project_name='grids_cnn')

    early_stop_cb = keras.callbacks.EarlyStopping(restore_best_weights=True,
                                                  patience=2)

    
    train_dataset, _, valid_dataset = create_datasets(max_targets=100)
    tuner.search(train_dataset, epochs=2, validation_data=valid_dataset, callbacks=[early_stop_cb]) 

    hp_names = ['filter1', 'filter2', 'filter3', 'kernel1', 'kernel2', 'kernel3', 'units1', 'units2']

    with open(f"{LOG_DIR}/hp_tuning_results.txt", "a") as f:
        best_hps = tuner.get_best_hyperparameters(num_trials=trials)[0]

        for name in hp_names:
            f.write(f"{name} -> {best_hps.get(name)}\n")

    return tuner, best_hps 


def run_model(model: tf.keras.Model, hps=None):
    """
    Runs fitting and evaluation
    Args:
        model: May be a model created manually or a keras Tuner object. 
        hps: A keras Hyperprameters object to set the hyperparameters to 
            be used in the model. Defaults to None and will assume the 
            model already has pre-baked hyperparameters.
    Returns: The outputs of fitting and evaluation, respectively.

    """

    checkpoint_cb = keras.callbacks.ModelCheckpoint(filepath= LOG_DIR + '/{epoch:02d}-weights.keras',
                                                    save_weights_only=True,
                                                    verbose=1)

    early_stop_cb = keras.callbacks.EarlyStopping(restore_best_weights=True,
                                                  patience=1)

    def get_dir() -> str :
        return f"{LOG_DIR}/fit/{strftime('run_%m_%d_%H_%M')}"

    train_dataset, test_dataset, valid_dataset = create_datasets()

    print("\nNow starting model fitting\n")

    fit_out = model.fit(train_dataset, epochs=30, validation_data=valid_dataset, callbacks=[checkpoint_cb, early_stop_cb])
    eval_out = model.evaluate(test_dataset)

    return fit_out, eval_out

if __name__ == '__main__':

#    tuner, best_hps = tune_model()
#    run_model(tuner, best_hps)

    cnn_model = create_model()
    run_model(cnn_model)
