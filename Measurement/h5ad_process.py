import scanpy as sc
import pandas as pd

"""
READ H5AD FILE TO PANDAS DATAFRAME
"""
# read h5ad
anndata = sc.read_h5ad("D:\\SRTP2021\\12\\liver\\Liver_droplet.h5ad", backed='r')
# write to csv
anndata.write_csvs("D:\\SRTP2021\\12\\spleen\\spleen_csv.csv")
# to matrix
matrix = anndata.raw.X.to_backed()
array = matrix.toarray()
df = pd.DataFrame(array)
# get rows
ages = anndata.obs['age']
ages_int = [eval(age[:-1]) for age in ages]
# get cols
var = pd.read_csv("D:\\SRTP2021\\12\\liver\\liver_csv\\var.csv")
genes = var.iloc[:,0]

# set rows & cols
df.columns = genes
df['target'] = ages_int

# select gene subset
subset = df[['target','Cd52','Plaur','B2m','Gpnmb','Gbp2','Ifitm2','Cd27','Cd69','Ptprj','Cyba','Ifitm3','Itm2b','Ly6e','Icam2','Fxyd5','Rpsa']]
subset = df[['Ms4a1', 'Tnfrsf17', 'Tnfsf13b']]

"""
TEST THIS df VAR IN GSE100906 BEST 100 FEATURES
"""

import tensorflow as tf
from tensorflow import keras
import mysql.connector
import pandas as pd

best_100_features = pd.read_csv("D:\\SRTP2021\\7\\data_selected_final_100.csv")
best_genes = best_100_features.columns

# how many are not match?
def test_match(target, ref):
    mismatches = []
    for col in ref[1:]:
        if col not in target:
            mismatches.append(col)
    return mismatches

target_match = [gene for gene in genes]
mismatches = test_match(target_match, best_genes)

# seems like 'Bin2', 'Crebrf', 'mt-Co1', 'mt-Cytb', 'mt-Nd1', 'mt-Rnr1', 'mt-Rnr2' are missing
# we have to drop them from the original model
best_features_compatible = best_100_features.drop(mismatches, axis=1)

# build new model based on compatible features

# create feature_compatible dataset
sample_name = best_features_compatible.pop('sample_name')

# normalize
def normalize_dataframe(dataframe):
    arr = dataframe.to_numpy()
    nor_arr = tf.keras.utils.normalize(arr, axis=-1, order=2)
    dataframe_nor = pd.DataFrame(nor_arr,
                                 columns=dataframe.columns,
                                 index=dataframe.index)
    return dataframe_nor
best_features_compatible_nor = normalize_dataframe(best_features_compatible)
    # add aged column
aged_col = []
for rowname in sample_name:
    if rowname[0] == "Y":
        aged_col.append(0)
    else:
        aged_col.append(1)
best_features_compatible_nor["aged"] = aged_col

def create_tensorflow_datasets(data, train_n=120, test_n=120):
    target = data.pop("aged")
    dataset = tf.data.Dataset.from_tensor_slices((data.values, target.values))
    # shuffle
    dataset.shuffle(50000)
    # train_test split
    if train_n != -1:
        train_set = dataset.take(train_n)
        train_dataset = train_set.shuffle(len(train_set)).batch(1)
        test_set = dataset.skip(test_n)
        test_dataset = test_set.shuffle(len(test_set)).batch(1)
        return train_dataset, test_dataset, target
    else:
        return dataset.take(len(dataset)).shuffle(50000).batch(1)

train_X, test_X, target = create_tensorflow_datasets(best_features_compatible_nor)

def get_compiled_model():
  model = tf.keras.Sequential([
    tf.keras.layers.Dense(100, activation='relu',
        kernel_initializer='ones',
        #kernel_regularizer=tf.keras.regularizers.L1(1),
        #activity_regularizer=tf.keras.regularizers.L2(1)
                          ),
    #tf.keras.layers.Dropout(0.5),
    tf.keras.layers.Dense(10, activation='relu',
        kernel_initializer='ones',
        #kernel_regularizer=tf.keras.regularizers.L1(1),
        #activity_regularizer=tf.keras.regularizers.L2(1)
                          ),
    #tf.keras.layers.Dropout(0.5),
    tf.keras.layers.Dense(1, activation='sigmoid')
  ])
  model.compile(optimizer='adam',
                loss='binary_crossentropy',
                metrics=['accuracy'])
  return model
model = get_compiled_model()
model.fit(train_X, epochs=10)
model.evaluate(test_X)

# create df from h5ad dataset
df_compatible = df[best_features_compatible.columns]

# normalize
df_compatible_nor = normalize_dataframe(df_compatible)

    # add age column
df_compatible_nor_rownames = df_compatible_nor.index
aged_col_df = []
for row in df_compatible_nor_rownames:
    if row in [1,3]:
        aged_col_df.append(0)
    else:
        aged_col_df.append(1)
df_compatible_nor['aged'] = aged_col_df

dataset_nor = create_tensorflow_datasets(df_compatible_nor, -1)

model.evaluate(dataset_nor)

# normalization increased accuracy from 26% to 57.6%, but that's not good enough




"""
MODEL WITH MARROW DATA 41941*19862
"""
df = pd.read_csv("E:\\Jetbrains\\PyCharmPro\\projects\\SRTP21\\df.csv")
df.pop('Unnamed: 0')

# select best features from df

# 1st: variance threshold
from sklearn.feature_selection import VarianceThreshold
def variance_threshold_selector(data, threshold=0.5):
    selector = VarianceThreshold(threshold)
    selector.fit(data)
    return data[data.columns[selector.get_support(indices=True)]]

# plot variance vs feature count
feature_count = []
for i in range(10):
    df_selected_tmp = variance_threshold_selector(df, 0.1 * i)
    feature_count.append(len(df_selected_tmp.count()))

import matplotlib.pyplot as plt

plt.plot([0.1 * x for x in range(10)], feature_count)
plt.xlabel("Variance Threshold")
plt.ylabel("Number of Features")
plt.title("Marrow Feature Counts vs Variance Threshold")
plt.show()

df_selected = variance_threshold_selector(df, 0.2)

# 2nd: k best selector
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

def k_best_selector(data, target, k=100):
    selector = SelectKBest(chi2, k=k)
    selector.fit(data, target)
    return data[data.columns[selector.get_support(indices=True)]]

# order features
def write_order(name, df_selected, target):
    selector = SelectKBest(chi2, k=500)
    selector.fit(df_selected, target)
    scores = selector.scores_
    features = df_selected.columns.values
    feature_score = {}
    for i in range(len(features)):
        feature_score[features[i]] = scores[i]
        # sort by score
    sorted_f_s = {k:v for k, v in sorted(feature_score.items(), key=lambda item:item[1])}
        # find 500th score
    cutoff = list(sorted_f_s.values())[-500]
    best_500 = [k for k in sorted_f_s if sorted_f_s[k] >= cutoff]
        # check
    # for gene in best_500:
        # if gene not in df_selected_final.columns.values:
            # print(gene)
        # write
    with open('D:\\SRTP2021\\12\\'+name+'.csv', 'w') as out:
        for k in best_500:
            out.write(k+','+str(sorted_f_s[k])+'\n')
        out.close()

# pop target
df_selected['target'] = target
target = df_selected.pop('target')
df_selected_final = k_best_selector(df_selected, target, 500) # Plaur is the 1045th best gene
df_selected_final['target'] = target

# create dataset*
def create_tensorflow_datasets(data, target=[], train_n=120, test_n=120):
    if len(target) == 0:
        target = data.pop("target")
    dataset = tf.data.Dataset.from_tensor_slices((data.values, target.values))
    # shuffle
    dataset.shuffle(50000)
    # train_test split
    if train_n != -1:
        train_set = dataset.take(train_n)
        train_dataset = train_set.shuffle(len(train_set)).batch(1)
        test_set = dataset.skip(test_n)
        test_dataset = test_set.shuffle(len(test_set)).batch(1)
        return train_dataset, test_dataset
    else:
        return dataset.take(len(dataset)).shuffle(50000).batch(1)
# train_df_final_X, test_df_final_X = create_tensorflow_datasets(df_selected_final, target, 34000, 34000)

# create dataset
def create_dataset(df_selected_final):
    train_dataset = df_selected_final.sample(frac=0.7, random_state=42)
    test_dataset = df_selected_final.drop(train_dataset.index)

    train_labels = train_dataset.pop('target')
    test_labels = test_dataset.pop('target')
    return train_dataset, test_dataset, train_labels, test_labels

train_dataset, test_dataset, train_labels, test_labels = create_dataset(df_selected_final)

    # normalize
def normalize_dataframe(dataframe):
    arr = dataframe.to_numpy()
    nor_arr = tf.keras.utils.normalize(arr, axis=-1, order=2)
    dataframe_nor = pd.DataFrame(nor_arr,
                                 columns=dataframe.columns,
                                 index=dataframe.index)
    return dataframe_nor
nor_train_data = normalize_dataframe(train_dataset)
nor_test_data = normalize_dataframe(test_dataset)

# build model
def get_compiled_model():
  model = tf.keras.Sequential([
    tf.keras.layers.Dense(300, activation='relu',
                          input_shape=[len(train_dataset.keys())]),
    tf.keras.layers.Dense(100, activation='relu'),
    tf.keras.layers.Dense(50, activation='relu'),
    tf.keras.layers.Dense(10, activation='relu'),
    tf.keras.layers.Dense(1)
  ])
  optimizer = tf.keras.optimizers.RMSprop(0.001)
  model.compile(optimizer=optimizer,
                loss='mse',
                metrics=['mae', 'mse'])
  return model
model = get_compiled_model()

    # model more robust to outliers
def get_huber_compiled_model():
    model = tf.keras.Sequential([
        tf.keras.layers.Dense(300, activation='relu',
                              input_shape=[len(train_dataset.keys())]),
        tf.keras.layers.Dense(100, activation='relu'),
        tf.keras.layers.Dense(50, activation='relu'),
        tf.keras.layers.Dense(10, activation='relu'),
        tf.keras.layers.Dense(1)
    ])
    optimizer = tf.keras.optimizers.RMSprop(0.001)
    model.compile(optimizer=optimizer,
                  loss=tf.keras.losses.Huber(),
                  metrics=['mae', 'mse'])
    return model
    model = get_huber_compiled_model()
model = get_huber_compiled_model()

# start training
    # early stop
EPOCHS = 100
early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)

    # history
history = model.fit(nor_train_data, train_labels, epochs=EPOCHS,
                    validation_split=0.2, callbacks=[early_stop])
hist = pd.DataFrame(history.history)
hist['epoch'] = history.epoch

# plot history
import matplotlib.pyplot as plt
def plot_history(history, tissue):
    hist = pd.DataFrame(history.history)
    hist['epoch'] = history.epoch
    plt.figure()
    plt.xlabel('Epoch')
    plt.ylabel('Mean abs Error')
    plt.plot(hist['epoch'], hist['mae'],
             label='Train Error')
    plt.plot(hist['epoch'], hist['val_mae'],
             label='Val Error')
    plt.ylim([0,8])
    plt.title(tissue + ' Mean abs Error vs Epoch')
    plt.legend()
    plt.show()
plot_history(history, "Marrow")

test_predictions = model.predict(nor_test_data)

# plot labels vs predictions
plt.scatter(test_labels, test_predictions)
plt.xlabel('True Values')
plt.ylabel('Predictions')
plt.axis('equal')
plt.axis('square')

"""
SAME PROCESS WITH LUNG & LIVER DATA
"""
# LUNG
df_lung = pd.read_csv("E:\\Jetbrains\\PyCharmPro\\projects\\SRTP21\\df_lung.csv")
df_lung.pop('Unnamed: 0')

# select best features

# 1st: variance threshold
from sklearn.feature_selection import VarianceThreshold
def variance_threshold_selector(data, threshold=0.5):
    selector = VarianceThreshold(threshold)
    selector.fit(data)
    return data[data.columns[selector.get_support(indices=True)]]

# plot variance vs feature count
feature_count = []
for i in range(10):
    df_lung_selected_tmp = variance_threshold_selector(df_lung, 0.1 * i)
    feature_count.append(len(df_lung_selected_tmp.count()))

import matplotlib.pyplot as plt
plt.plot([0.1 * x for x in range(10)], feature_count)
plt.xlabel("Variance Threshold")
plt.ylabel("Number of Features")
plt.title("Lung Feature Counts vs Variance Threshold")
plt.show()

df_lung_selected = variance_threshold_selector(df_lung, 0.4)

# 2nd: k best selector
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
def k_best_selector(data, target, k=100):
    selector = SelectKBest(chi2, k=k)
    selector.fit(data, target)
    return data[data.columns[selector.get_support(indices=True)]]

# order
lung_target = df_lung_selected.pop('target')
write_order('df_lung_best_500_score', df_lung_selected, lung_target)

# pop target
lung_target = df_lung_selected.pop('target')
df_lung_selected_final = k_best_selector(df_lung_selected, lung_target, 500) # Plaur 1155th best gene
df_lung_selected_final['target'] = lung_target

# create dataset
train_lung_dataset, test_lung_dataset, train_lung_labels, test_lung_labels = create_dataset(df_lung_selected_final)

    # normalize
nor_lung_train_data = normalize_dataframe(train_lung_dataset)
nor_lung_test_data = normalize_dataframe(test_lung_dataset)

# build model
lung_model = get_compiled_model()

    # early stop
EPOCHS = 100
early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)
    # history
lung_history = lung_model.fit(nor_lung_train_data, train_lung_labels, epochs=EPOCHS,
                    validation_split=0.2, callbacks=[early_stop])
lung_hist = pd.DataFrame(lung_history.history)
lung_hist['epoch'] = lung_history.epoch

# plot history
import matplotlab.pyplot as plt
plot_history(lung_history, 'Lung')

# predict
test_lung_predictions = lung_model.predict(nor_lung_test_data)





# LIVER
df_liver = pd.read_csv("E:\\Jetbrains\\PyCharmPro\\projects\\SRTP21\\df_liver.csv")
df_liver.pop('Unnamed: 0')

# select best features

# 1st: variance threshold
from sklearn.feature_selection import VarianceThreshold
def variance_threshold_selector(data, threshold=0.5):
    selector = VarianceThreshold(threshold)
    selector.fit(data)
    return data[data.columns[selector.get_support(indices=True)]]

# plot variance vs feature count
feature_count = []
for i in range(10):
    df_liver_selected_tmp = variance_threshold_selector(df_liver, 0.1 * i)
    feature_count.append(len(df_liver_selected_tmp.count()))

import matplotlib.pyplot as plt
plt.plot([0.1 * x for x in range(10)], feature_count)
plt.xlabel("Variance Threshold")
plt.ylabel("Number of Features")
plt.title("Liver Feature Counts vs Variance Threshold")
plt.show()

df_liver_selected = variance_threshold_selector(df_liver, 0.4)

# 2nd: k best selector
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
def k_best_selector(data, target, k=100):
    selector = SelectKBest(chi2, k=k)
    selector.fit(data, target)
    return data[data.columns[selector.get_support(indices=True)]]

# order
liver_target = df_liver_selected.pop('target')
write_order('df_liver_best_500_score', df_liver_selected, liver_target)

# pop target
liver_target = df_liver_selected.pop('target')
df_liver_selected_final = k_best_selector(df_liver_selected, liver_target, 500) # Plaur 1926th best gene
df_liver_selected_final['target'] = liver_target

# create dataset
train_liver_dataset, test_liver_dataset, train_liver_labels, test_liver_labels = create_dataset(df_liver_selected_final)

    # normalize
nor_liver_train_data = normalize_dataframe(train_liver_dataset)
nor_liver_test_data = normalize_dataframe(test_liver_dataset)

# build model
liver_model = get_compiled_model()

    # early stop
EPOCHS = 100
early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)
    # history
liver_history = liver_model.fit(nor_liver_train_data, train_liver_labels, epochs=EPOCHS,
                    validation_split=0.2, callbacks=[early_stop])
liver_hist = pd.DataFrame(liver_history.history)
liver_hist['epoch'] = liver_history.epoch

# plot history
import matplotlab.pyplot as plt
plot_history(liver_history, "Liver")

# predict
test_liver_predictions = liver_model.predict(nor_liver_test_data)

"""
CODES FOR PERSONAL SURFACE
"""
# read df_selected
df_selected_final = pd.read_csv("D:\\SRTP2021\\12\\model_marrow\\df_selected_final_11.30.csv")
df_selected_final.pop('Unnamed: 0')

# create dataset
train_dataset, test_dataset, train_labels, test_labels = create_dataset(df_selected_final)
# normalize dataset
nor_train_data = normalize_dataframe(train_dataset)
nor_test_data = normalize_dataframe(test_dataset)

# build model
def get_compiled_model():
  model = tf.keras.Sequential([
    tf.keras.layers.Dense(300, activation='relu',
                          input_shape=[len(train_dataset.keys())]),
    tf.keras.layers.Dense(100, activation='relu'),
    tf.keras.layers.Dense(50, activation='relu'),
    tf.keras.layers.Dense(10, activation='relu'),
    tf.keras.layers.Dense(1)
  ])
  optimizer = tf.keras.optimizers.RMSprop(0.001)
  model.compile(optimizer=optimizer,
                loss='mse',
                metrics=['mae', 'mse'])
  return model
model = get_compiled_model()

    # model more robust to outliers
def get_huber_compiled_model():
    model = tf.keras.Sequential([
        tf.keras.layers.Dense(300, activation='relu',
                              input_shape=[len(train_dataset.keys())]),
        tf.keras.layers.Dense(100, activation='relu'),
        tf.keras.layers.Dense(50, activation='relu'),
        tf.keras.layers.Dense(10, activation='relu'),
        tf.keras.layers.Dense(1)
    ])
    optimizer = tf.keras.optimizers.RMSprop(0.001)
    model.compile(optimizer=optimizer,
                  loss=tf.keras.losses.Huber(),
                  metrics=['mae', 'mse'])
    return model
    model = get_huber_compiled_model()
model = get_huber_compiled_model()

# start training
    # early stop
EPOCHS = 50
early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=100)

    # history
history = model.fit(nor_train_data, train_labels, epochs=EPOCHS,
                    validation_split=0.2, callbacks=[early_stop])
hist = pd.DataFrame(history.history)
hist['epoch'] = history.epoch

# plot history
import matplotlib.pyplot as plt
def plot_history(history, tissue):
    hist = pd.DataFrame(history.history)
    hist['epoch'] = history.epoch
    plt.figure()
    plt.xlabel('Epoch')
    plt.ylabel('Mean abs Error')
    plt.plot(hist['epoch'], hist['mae'],
             label='Train Error')
    plt.plot(hist['epoch'], hist['val_mae'],
             label='Val Error')
    plt.ylim([0,8])
    plt.title(tissue + ' Mean abs Error vs Epoch')
    plt.legend()
    plt.show()
plot_history(history, "Marrow")




### repetitive learning
import numpy as np

# try doing in python
    # plot first: pretty similar to test set result
plt.scatter(train_predictions, train_labels, alpha=0.05)
plt.show()
    # concatenate
label_predictions = np.concatenate((train_predictions, train_labels), axis=1)
label_predictions_df = pd.DataFrame(label_predictions, columns=['prediction', 'actual'])
# never mind

# read from R
new_train_data = pd.read_csv('D:\\SRTP2021\\12\\repetitive_learning\\new_train_data.csv')
new_train_data.pop('Unnamed: 0')
new_train_data_arr = new_train_data.to_numpy()

nor_train_data['new_1'] = new_train_data_arr
nor_train_data['target'] = train_labels

nor_train_data_new1 = nor_train_data[nor_train_data['new_1'] == True]

nor_train_data_new1.pop('new_1')
train_labels_new1 = nor_train_data_new1.pop('target')


# now train on new data
model_new1 = get_compiled_model()

# learn once
nor_train_data.pop('new_1')
nor_train_data.pop('target')

EPOCHS = 50
early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=100)
    # history
history = model_new1.fit(nor_train_data, train_labels, epochs=EPOCHS,
                    validation_split=0.2, callbacks=[early_stop])
hist = pd.DataFrame(history.history)
hist['epoch'] = history.epoch

# learn twice
EPOCHS = 50
early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=100)
    # history
history = model_new1.fit(nor_train_data_new1, train_labels_new1, epochs=EPOCHS,
                    validation_split=0.2, callbacks=[early_stop])
hist = pd.DataFrame(history.history)
hist['epoch'] = history.epoch

# now predict
test_predictions_new1 = model.predict(nor_test_data)



# test prediction
test_predictions = model.predict(nor_test_data)

# plot labels vs predictions
plt.scatter(test_labels, test_predictions)
plt.xlabel('True Values')
plt.ylabel('Predictions')
plt.axis('equal')
plt.axis('square')

"""
Select markers & controls
"""

# Marrow
with open('D:\\SRTP2021\\12\\genes_to_select.txt','r') as infile:
    lines = infile.readlines()
    genes = [line[:-1] for line in lines]
    infile.close()

genes.remove('H2-eb1')
genes.remove('Cd20')
genes.remove('Bcma')

df_genes = df[genes]

# Lung
df_lung_genes = df_lung[genes]

# Liver
df_liver_genes = df_liver[genes]

# F
K84_T24 = pd.read_csv('D:\\SRTP2021\\15_human_rnaseq\\K_84_T_24.csv')

"""
Order & find membrane protein
"""
from uniprot_web_scraping import *
# Liver
liver_mem = main_function(inpath='D:\\SRTP2021\\12\\scores\\df_liver_best_500_gene.csv',
                         outpath='D:\\SRTP2021\\12\\scores\\df_liver_best_500_mem.csv')

# lung
lung_mem = main_function(inpath='D:\\SRTP2021\\12\\scores\\df_lung_best_500_gene.csv',
                         outpath='D:\\SRTP2021\\12\\scores\\df_lung_best_500_mem.csv')
for i in lung_mem:
    print(i)

# marrow
marrow_mem = main_function(inpath='D:\\SRTP2021\\12\\scores\\df_best_500_gene.csv',
                           outpath='D:\\SRTP2021\\12\\scores\\df_best_500_mem.csv')
for i in marrow_mem:
    print(i)

# add scores
    # marrow
with open('D:\\SRTP2021\\12\\scores\\df_liver_best_500_score.csv') as dfscore:
    lines = dfscore.readlines()
    gene_score = [line.split(',') for line in lines]
    dfscore.close()

with open('D:\\SRTP2021\\12\\scores\\df_liver_best_500_mem.csv') as dfmem:
    lines = dfmem.readlines()
    df_mem = [line[:-1] for line in lines]
    dfmem.close()

with open('D:\\SRTP2021\\12\\scores\\df_liver_best_500_memscore.csv','w') as dfout:
    for arr in gene_score:
        if arr[0] in df_mem:
            dfout.write(arr[0]+','+arr[1])
    dfout.close()
