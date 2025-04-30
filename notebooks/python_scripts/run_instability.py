import sys
import pickle
sys.path.append('/central/groups/mthomson/jboktor/spatial_genomics/osNMF/osNMF')
from knn_impute import kNN_imputation
from morans_index import MoranIforPPs
from other_functions import show_graph_with_labels, filter_genes, weighted_correlation
from sklearn_nmf import sklearn_nmf
from stability import findcorrelation, HungarianError, amariMaxError, instability


def run_instability(folder_path, K1, K2, B, shape_flat, name_DG):
    instability_DG, instability_std_DG, distMat_DG = instability(_folder_path=folder_path,
                                                                 _k1=K1,
                                                                 _k2=K2,
                                                                 _numReplicates=B,
                                                                 _n_features=shape_flat,
                                                                 _name=name_DG)
    print("\n DG instability values:", instability_DG)
    results = {
        "instability_DG": instability_DG,
        "instability_std_DG": instability_std_DG,
        "distMat_DG": distMat_DG
    }
    return results

if __name__ == "__main__":
    folder_path = sys.argv[1]
    K1 = int(sys.argv[2])
    K2 = int(sys.argv[3])
    B = int(sys.argv[4])
    shape_flat = int(sys.argv[5])
    name_DG = sys.argv[6]
    output_file = sys.argv[7]
    
    results = run_instability(folder_path, K1, K2, B, shape_flat, name_DG)
    
    # Save results to a pickle file
    with open(output_file, 'wb') as f:
        pickle.dump(results, f)
