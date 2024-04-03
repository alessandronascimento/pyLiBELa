import tensorflow as tf
from pathlib import Path
from sklearn.model_selection import train_test_split

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession


config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)
path_to_pdbbind = "/data/home/alessandro/PDBbind_v2020"
physical_devices = tf.config.experimental.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(physical_devices[0], True)

# Commented out IPython magic to ensure Python compatibility.


#@title Conecting to Google Drive

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

binding_file = open(f'{path_to_pdbbind}/index/INDEX_refined_data.2020', 'r')
for line in binding_file:
  line2 = line.split()
  if line2[0][0] != '#':
    if line2[0] not in targets_ok: continue

    binding_targets.append(line2[0])
    if (line2[1] == 'NMR'):
      binding_resolution.append(0.00)
    else:
      binding_resolution.append(float(line2[1]))
    binding_years.append(int(line2[2]))
    binding_score.append(float(line2[3]))
binding_file.close()

print('Binding data found for %d valid targets' % len(binding_targets))

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
#@title Reading dataset


train_names, test_names, train_scores, test_scores = train_test_split(binding_targets, binding_score, train_size=0.8, shuffle=True)
train_names, valid_names, train_scores, valid_scores = train_test_split(train_names, train_scores, train_size=0.8)

def slice_data(dataset):
  x = tf.concat([x for x, y in dataset], axis=-1)
  y = tf.concat([y for x, y in dataset], axis=-1)
  return x,y

# Some dataset parameters

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

# train_val, test = tf.keras.utils.split_dataset(dataset, left_size=0.8, right_size=0.2, shuffle=True, seed=23)
# train, valid = tf.keras.utils.split_dataset(train_val, left_size=0.8, right_size=0.2)

# x_train, y_train = slice_data(train)
# x_test, y_test = slice_data(test)
# x_valid, y_valid = slice_data(valid)

# Repeat the dataset for multiple epochs (optional)
#dataset = dataset.repeat()

#@title Alex's model

cnn_model = tf.keras.Sequential([
        tf.keras.layers.Conv3D(64, 5, activation='relu', data_format='channels_last', input_shape=(60,60,60,6), padding="same"),
        tf.keras.layers.MaxPool3D(data_format='channels_last'),
        tf.keras.layers.Conv3D(128, 3, activation='relu', data_format='channels_last', padding="same"),
        tf.keras.layers.Conv3D(128, 3, activation='relu', data_format='channels_last', padding="same"),
        tf.keras.layers.MaxPool3D(data_format='channels_last'),
        tf.keras.layers.Conv3D(256, 3, activation='relu', data_format='channels_last', padding="same"),
        tf.keras.layers.Conv3D(256, 3, activation='relu', data_format='channels_last', padding="same"),
        tf.keras.layers.Flatten(),
        tf.keras.layers.Dense(units=128, activation="relu"),
        tf.keras.layers.Dense(units=64, activation="relu"),
        tf.keras.layers.Dense(units=1),
    ])

cnn_model.compile(loss="mse", optimizer=tf.keras.optimizers.Adam(learning_rate=1e-4), metrics=["RootMeanSquaredError"])

cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath='{epoch:02d}-{val_loss:.2f}.keras',
                                                     save_weights_only=True,
                                                          verbose=1)


out_info = cnn_model.fit(train_dataset, epochs=10, validation_data=valid_dataset, callbacks=[cp_callback])
mse_test = cnn_model.evaluate(test_dataset)

"""#Questions:

1. What if we used a batch normalization in the model? Something like:

```
tf.keras.layers.BatchNormalization()
```

in the begining of the model and after each hidden layer?

"""
