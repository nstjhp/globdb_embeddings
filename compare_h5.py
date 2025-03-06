import h5py
import pandas as pd
#with h5py.File("Ecoli/protein_embeddings.h5") as hf:
#    print(len(hf.keys()))
#with h5py.File("Ecoli/protein_embeddings2.h5") as hf:
#    print(len(hf.keys()))

# see https://stackoverflow.com/a/74127100
dictionary = {}
with h5py.File("Ecoli/protein_embeddings.h5", "r") as f:
    for key in f.keys():
        print(key)

        ds_arr = f[key][()]   # returns as a numpy array
        dictionary[key] = ds_arr # appends the array in the dict under the key
df = pd.DataFrame.from_dict(dictionary)
df
df.mean(axis=0) # for each protein, 1 would be across each of 1024 positions of the embeddings

dictionary2 = {}
with h5py.File("Ecoli/protein_embeddings2.h5", "r") as f:
    for key in f.keys():
        print(key)

        ds_arr = f[key][()]   # returns as a numpy array
        dictionary2[key] = ds_arr # appends the array in the dict under the key
df2 = pd.DataFrame.from_dict(dictionary2)
df2.mean(axis=0)

