import os
import sys
import popalign as PA
import numpy as np
import pandas as pd
import pickle
from datetime import datetime

def popalign_workflow(input_folder, output_folder):
    print(f"Processing input folder: {input_folder}")
    print(f"Output folder: {output_folder}")

    try:
        mymatrix = os.path.join(input_folder, 'matrix.mtx')
        mybarcodes = os.path.join(input_folder, 'barcodes.tsv')
        mygenes = os.path.join(input_folder, 'genes.tsv')
        mymetadata = os.path.join(input_folder, 'metadata.csv')

        pop = PA.load_multiplexed(matrix=mymatrix,
                                  barcodes=mybarcodes,
                                  genes=mygenes,
                                  metafile=mymetadata,
                                  controlstring='SPF',
                                  outputfolder=output_folder,
                                  only=[], # list of sample names to only load the specified samples
                                  col=None, # either None or a column name from the meta data
                                  value=None, # if col != None, specify value in column to filter samples
                                  existing_obj=None)

        PA.print_ncells(pop)
        PA.normalize(pop, scaling_factor=1000)
        PA.plot_gene_filter(pop, offset=0.9)
        PA.filter(pop, remove_ribsomal=False, remove_mitochondrial=False)
        print(f"Number of filtered genes: {len(pop['filtered_genes'])}")

        # osNMF analysis
        PA.onmf(pop, ncells=20000, nfeats=list(range(3,20)), nreps=2, niter=500)

        # Save files for one featureset
        PA.choose_featureset(pop, alpha=3, multiplier=3)

        PA.build_gmms_by_celltypes(
            pop, 
            ks=(3,15), 
            niters=3, 
            training=0.7, 
            nreplicates=10, 
            reg_covar='auto', 
            types=None, 
            criteria='aic', 
            only=None
        )

        # Visualize GMMs
        PA.render_models(pop, figsizegrouped=(30,30), samples=pop['order'], figsizesingle=(6,5), mode='grouped')

        # Calculate LLR ranking
        PA.rank(pop,
                ref='SPF run4 roi4', # label of the reference sample
                k=20000, # number of cells per bootstrapping sample
                niter=500, # number of iterations
                method='LLR', # LLR for log-likelihood ratio or LL for log-likelihood
                mincells=50, # sample's minimum number of cells to be included in ranking
                figsize=(10,5)) # plot figure size

        PA.align(pop, ref='SPF run4 roi4',
                 method='conservative', # one of: test2ref, ref2test, conservative
                 figsizedeltas=(10,10),
                 figsizeentropy=(10,10))

        PA.plot_deltas(pop, figsize=(10,10), sortby='mu', pthresh=2)

        # Collect delta stats
        df_list = []
        for k in pop['deltas'].keys():
            df = pop['deltas'][k]['combined'].transpose()
            df['cell_type'] = k
            df_list.append(df)

        # Save the combined dataframe as a CSV file
        combined_df = pd.concat(df_list, ignore_index=True)
        delta_stats_df = os.path.join(output_folder, 'delta_stats.csv')
        combined_df.to_csv(delta_stats_df, index=False)
        print(f"Combined dataframe saved to {delta_stats_df}")

        PA.samples_grid(pop, method='tsne', figsize=(20,20), size_background=.1, size_samples=.3, samplecolor='magenta', showplot=True)

        # Save pop object as a pickle file
        current_date = datetime.now().strftime('%Y-%m-%d')
        filename = os.path.join(output_folder, f'popalign_{current_date}.pkl')

        with open(filename, 'wb') as file:
            pickle.dump(pop, file)
        print(f"Popalign object saved to {filename}")

    except Exception as e:
        print(f"Error processing folder {input_folder}: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python popalign_workflow.py <input_folder> <output_folder>")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    
    popalign_workflow(input_folder, output_folder)
