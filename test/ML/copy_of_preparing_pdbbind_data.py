import tensorflow as tf
import tensorflow_decision_forests as tfdf
import keras.api._v2.keras as keras
import keras_tuner as kt
from pathlib import Path
from sklearn.model_selection import train_test_split
from time import strftime


from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

path_to_pdbbind = "/data/home/alessandro/PDBbind_v2020"
LOG_DIR = "./log"

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)
physical_devices = tf.config.experimental.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(physical_devices[0], True)

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

#@title Auxiliary functions
def process_file(file_path):

  as_path = Path(file_path)
  if not as_path.exists():
    raise IOError(f"File path '{as_path}' does not exist")

  data = tf.io.read_file(file_path)
  data = tf.io.decode_raw(data, tf.float64)
  data = tf.reshape(data, (60, 60, 60, 3))
  return data

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
    combined_grid = tf.concat([grid1, grid2], axis=-1)

    # Yield the data as a tuple
    yield combined_grid, observable

    # Yields the 10 related decoys of the current grid
    for j in range(1,11):

      grid2 = process_file(file_path + '/McGrid_dec_{}.grid'.format(j))

      observable = 0.00
      combined_grid = tf.concat([grid1, grid2], axis=-1)
      yield combined_grid, observable

def slice_data(dataset):
  x = tf.concat([x for x, y in dataset], axis=-1)
  y = tf.concat([y for x, y in dataset], axis=-1)
  return x,y

def create_datasets():

    binding_targets, binding_score = read_pdbbind_data()

    train_names, test_names, train_scores, test_scores = train_test_split(binding_targets, binding_score, train_size=0.8, shuffle=True)
    train_names, valid_names, train_scores, valid_scores = train_test_split(train_names, train_scores, train_size=0.8)

    output_signature = (tf.TensorSpec(shape=(60, 60, 60, 6), dtype=tf.float64, name='combgrid'),
                        tf.TensorSpec(shape=(), dtype=tf.float32, name='ligscore'))
    batch_size = 5
    prefetch_size = 1

# Create the dataset from the generator function,
# with batch and prefetch sizes already determined
    train_dataset = tf.data.Dataset.from_generator(
        data_generator,
        output_signature=output_signature,
        args=(tf.convert_to_tensor(train_names, dtype=tf.string), tf.convert_to_tensor(train_scores, dtype=tf.float32)),
        name="train_dataset_gen"
    ).batch(batch_size).prefetch(prefetch_size)

    test_dataset = tf.data.Dataset.from_generator(
        data_generator,
        output_signature=output_signature,
        args=(tf.convert_to_tensor(test_names, dtype=tf.string), tf.convert_to_tensor(test_scores, dtype=tf.float32)),
        name="test_dataset_gen"
    ).batch(batch_size).prefetch(prefetch_size)

    valid_dataset = tf.data.Dataset.from_generator(
        data_generator,
        output_signature=output_signature,
        args=(tf.convert_to_tensor(valid_names, dtype=tf.string), tf.convert_to_tensor(valid_scores, dtype=tf.float32)),
        name="valid_dataset_gen"
    ).batch(batch_size).prefetch(prefetch_size)

    return train_dataset, test_dataset, valid_dataset

#@title Alex's model
def create_model(hp):

    # FIXME bad coding:
    # IF THIS CHANGES, THE NAMES LIST IN tune_model() MUST CHANGE 
    filters1 = hp.Int('filter1', min_value=40, max_value=80, step=2)
    kernel1 = hp.Int('kernel1', min_value=3, max_value=9, step=1)
    filters2 = hp.Int('filter2', min_value=90, max_value=200, step=5)
    kernel2 = hp.Int('kernel2', min_value=2, max_value=9, step=1)
    filters3 = hp.Int('filter3', min_value=90, max_value=200, step=5)
    kernel3 = hp.Int('kernel3', min_value=2, max_value=9, step=1)
    units1 = hp.Int('units1', min_value=50, max_value=200, step=10)
    units2 = hp.Int('units2', min_values=30, max_values=80, step=5)


    cnn_model = keras.Sequential([
            keras.layers.Conv3D(filters1, kernel1, activation='relu', data_format='channels_last', input_shape=(60,60,60,6), padding="same"),
            keras.layers.MaxPool3D(data_format='channels_last'),
            keras.layers.Conv3D(filters2, kernel2, activation='relu', data_format='channels_last', padding="same"),
            keras.layers.Conv3D(filters2, kernel2, activation='relu', data_format='channels_last', padding="same"),
            keras.layers.MaxPool3D(data_format='channels_last'),
            keras.layers.Conv3D(filters3, kernel3, activation='relu', data_format='channels_last', padding="same"),
            keras.layers.Conv3D(filters3, kernel3, activation='relu', data_format='channels_last', padding="same"),
            keras.layers.Flatten(),
            keras.layers.Dense(units=units1, activation="relu"),
            keras.layers.Dense(units=units2, activation="relu"),
            keras.layers.Dense(units=1),
        ])

    learning_rate = hp.Choice('learning_rate', values=[1e-2, 1e-3, 1e-4, 1e-5])

    cnn_model.compile(loss="mse", optimizer=keras.optimizers.Adam(learning_rate=learning_rate), metrics=["RootMeanSquaredError"])

    return cnn_model

def tune_model():
    tuner = kt.BayesianOptimization(create_model,
                                    objective='root_mean_squared_error',
                                    max_trials=10,
                                    directory=LOG_DIR,
                                    project_name='grids_cnn')

    early_stop_cb = keras.callbacks.EarlyStopping(restore_best_weights=True,
                                                  patience=2)

    
    train_dataset, test_dataset, valid_dataset = create_datasets()
    tuner.search(train_dataset, epochs=10, validation_data=valid_dataset) # validation created from train dataset so that the hp search doesn't overfit over the true validation

    hp_names = ['filter1', 'filter2', 'filter3', 'kernel1', 'kernel2', 'kernel3', 'units1', 'units2']

    with open(f"{LOG_DIR}/hp_tuning_results.txt", "a") as f:
        best_hps = tuner.get_best_hyperparameters(num_trials=1)[0]

        for name in hp_names:
            f.write(f"{name} -> {best_hps.get(name)}\n")


def run_model(model):

    checkpoint_cb = keras.callbacks.ModelCheckpoint(filepath= LOG_DIR + '/{epoch:02d}-{val_loss:.2f}.keras',
                                                    save_weights_only=True,
                                                    verbose=1)

    early_stop_cb = keras.callbacks.EarlyStopping(restore_best_weights=True,
                                                  patience=1)

    def get_dir() -> str :
        return f"{LOG_DIR}/{strftime("run_%m_%d_%H_%M")}"

    tensorboard_cb = keras.callbacks.TensorBoard(logdir=get_dir(), profile_batch=(100,200))
    

    model = create_model()
    train_dataset, test_dataset, valid_dataset = create_datasets()

    out_info = model.fit(train_dataset, epochs=10, validation_data=valid_dataset, callbacks=[checkpoint_cb, early_stop_cb, tensorboard_cb])
    mse_test = model.evaluate(test_dataset)

if __name__ == '__main__':

    train_dataset, test_dataset, valid_dataset = create_datasets()
    cnn_model = create_model()
    run_model(cnn_model)
