import h5py
import pandas as pd
import numpy as np

#with h5py.File("Ecoli/protein_embeddings.h5") as hf:
#    print(len(hf.keys()))
#with h5py.File("Ecoli/protein_embeddings2.h5") as hf:
#    print(len(hf.keys()))

# see https://stackoverflow.com/a/74127100
dictionary = {}
with h5py.File("4448170_node-a06_CPU.h5", "r") as f:
    for key in f.keys():
        print(key)

        ds_arr = f[key][()]   # returns as a numpy array
        dictionary[key] = ds_arr
dfcpu = pd.DataFrame.from_dict(dictionary)
dfcpu.mean(axis=0) # for each protein, 1 would be across each of 1024 positions of the embeddings

dictionary2 = {}
with h5py.File("4448208_node-d01_GPU.h5", "r") as f:
    for key in f.keys():
        print(key)

        ds_arr = f[key][()]   # returns as a numpy array
        dictionary2[key] = ds_arr
dfgpu = pd.DataFrame.from_dict(dictionary2)

mean_cpu = dfcpu.mean(axis=0)
mean_gpu = dfgpu.mean(axis=0)

# Calculate the absolute difference in means for each protein
abs_diff_means = (mean_cpu - mean_gpu).abs()
np.sort(abs_diff_means)

# Calculate the relative difference in means for each protein
rel_diff_means = abs_diff_means / mean_gpu.abs()
np.sort(rel_diff_means)

def cosine_similarity(a, b):
    """Calculate cosine similarity between two vectors."""
    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

# Calculate cosine similarity for each protein (for each column)
cos_sim = {}
for col in dfcpu.columns:
    vec_cpu = dfcpu[col].values
    vec_gpu = dfgpu[col].values
    cos_sim[col] = cosine_similarity(vec_cpu, vec_gpu)
cos_sim = pd.Series(cos_sim)
np.sort(cos_sim)

# Pearson correlation coefficient and Mean Squared Error (MSE) for each protein
pearson_corr = {}
mse = {}
for col in dfcpu.columns:
    vec_cpu = dfcpu[col].values
    vec_gpu = dfgpu[col].values
    # Pearson correlation: np.corrcoef returns a 2x2 matrix
    pearson_corr[col] = np.corrcoef(vec_cpu, vec_gpu)[0, 1]
    mse[col] = np.mean((vec_cpu - vec_gpu) ** 2)
pearson_corr = pd.Series(pearson_corr)
mse = pd.Series(mse)

np.sort(pearson_corr)
np.sort(mse)


def compare_embedding(protein, dfcpu, dfgpu, index=None, new_value=None):
    """
    Compare embeddings for a given protein from CPU and GPU dataframes.

    Parameters:
      protein (str): The protein identifier (column name in the dataframes).
      dfcpu (pd.DataFrame): DataFrame containing CPU embeddings.
      dfgpu (pd.DataFrame): DataFrame containing GPU embeddings.
      index (int, optional): Index in the CPU embedding to modify.
      new_value (float, optional): New value to set at the specified index in the CPU embedding.

    The function prints:
      - The full CPU and GPU embeddings for the protein.
      - The cosine similarity between them.
      - If index and new_value are provided, the updated CPU embedding and the new cosine similarity.
    """
    # Check if the protein exists in both DataFrames
    if protein not in dfcpu.columns or protein not in dfgpu.columns:
        print(f"Protein {protein} not found in both DataFrames.")
        return

    # Extract embeddings as numpy arrays (copy to allow modifications)
    vec_cpu = dfcpu[protein].values.copy()
    vec_gpu = dfgpu[protein].values.copy()

    # Print original embeddings
    print(f"Original CPU embedding for protein '{protein}':\n{vec_cpu}\n")
    print(f"GPU embedding for protein '{protein}':\n{vec_gpu}\n")

    # Calculate and print cosine similarity
    orig_cos_sim = cosine_similarity(vec_cpu, vec_gpu)
    print(f"Original Cosine Similarity: {orig_cos_sim:.4f}\n")

    # If index and new_value are provided, update the CPU embedding
    if index is not None and new_value is not None:
        # Check index bounds
        if index < 0 or index >= len(vec_cpu):
            print(f"Index {index} is out of bounds for the embedding (length {len(vec_cpu)}).")
            return

        old_val = vec_cpu[index]
        vec_cpu[index] = new_value
        print(f"Updated CPU embedding for protein '{protein}':")
        print(f" - Changed position {index} from {old_val} to {new_value}")
        print(vec_cpu, "\n")

        # Recalculate and print new cosine similarity
        new_cos_sim = cosine_similarity(vec_cpu, vec_gpu)
        print(f"New Cosine Similarity after update: {new_cos_sim:.4f}")

# Example usage:
# compare_embedding("P12345", dfcpu, dfgpu)
# To update the 1st position (index 0) from its current value to 10:
# compare_embedding("P12345", dfcpu, dfgpu, index=0, new_value=10)

