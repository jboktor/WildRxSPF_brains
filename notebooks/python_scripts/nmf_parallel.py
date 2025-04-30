import polars as pl
import pandas as pd
from sklearn.impute import SimpleImputer
# from sklearn.decomposition import NMF as sklearn_nmf
import numpy as np
import os
import sys
sys.path.append('/central/groups/mthomson/jboktor/spatial_genomics/osNMF/osNMF')
from knn_impute import kNN_imputation
from morans_index import MoranIforPPs
from other_functions import show_graph_with_labels, filter_genes, weighted_correlation
from sklearn_nmf import sklearn_nmf
from stability import findcorrelation, HungarianError, amariMaxError, instability

# Function to run NMF on bootstrapped data
def run_nmf_replicate(K, b, filtered_data, path):
    print(f"Working on replicate {b} for K={K}...")
    np.random.seed(b)
    X = np.maximum(filtered_data, 0)
    n_samples = X.shape[0]
    random_indices = np.random.choice(n_samples, n_samples, replace=True)
    bootstrap_X = X.iloc[random_indices]

    nmf = sklearn_nmf(n_components=K, l1_ratio=1, alpha=0, max_iter=5000, random_state=1)
    nmf.fit(bootstrap_X)
    PPs_tmp = nmf.components_
    np.savez(os.path.join(path, f'DG_PP_{b}'), PPs_tmp=PPs_tmp)

# Read and preprocess the data
def load_and_preprocess_data():
    try:
        filtered_data_polars = pl.read_csv(
            '/central/groups/mthomson/jboktor/spatial_genomics/jess_2024-01-23/data/interim/seurat_data_2024-06-15.csv',
            has_header=True,
            ignore_errors=True
        )
        
        first_col_name = filtered_data_polars.columns[0]
        filtered_data_polars = filtered_data_polars.rename({first_col_name: "rownames"})
        filtered_data = filtered_data_polars.to_pandas()
        filtered_data.set_index('rownames', inplace=True)
        
        if filtered_data.isna().any().any():
            imputer = SimpleImputer(strategy='mean')
            filtered_data[:] = imputer.fit_transform(filtered_data)
        
        return filtered_data

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

if __name__ == "__main__":
    K = int(sys.argv[1])
    b = int(sys.argv[2])
    folder_path = sys.argv[3]

    filtered_data = load_and_preprocess_data()
    if filtered_data is None:
        exit(1)

    path = os.path.join(folder_path, f'K={K}')
    os.makedirs(path, exist_ok=True)

    run_nmf_replicate(K, b, filtered_data, path)
