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

        pop['ncores'] = 32
        PA.print_ncells(pop)
        PA.normalize(pop, scaling_factor=1000)
        PA.plot_gene_filter(pop, offset=0.9)
        PA.filter(pop, remove_ribsomal=False, remove_mitochondrial=False)
        print(f"Number of filtered genes: {len(pop['filtered_genes'])}")

        # osNMF analysis
        PA.onmf(pop, ncells=20000, nfeats=list(range(3,20)), nreps=3, niter=500)

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

        # Save pop object as a pickle file
        current_date = datetime.now().strftime('%Y-%m-%d')
        filename = os.path.join(output_folder, f'popalign_{current_date}.pkl')
        with open(filename, 'wb') as file:
            pickle.dump(pop, file)
        print(f"Popalign object saved to {filename}")

        # analysis_root = os.path.dirname(output_folder)
        # dirs = [d for d in os.listdir(analysis_root) if os.path.isdir(os.path.join(analysis_root, d))]
        
        all_dfs = []
        quality_slices = pop['order']
        
        # Loop through slices and do all-vs-all comparisons for det
        for qs in quality_slices:
            try:
                print(qs + os.path.basename(output_folder))

                output_dir_sub = f'{output_folder}/multicomp/{qs}'
                PA.mkdir(output_dir_sub) 
                pop['output'] = output_dir_sub

                PA.align(pop, ref=qs,
                          method='celltype',
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
                delta_stats_df = os.path.join(output_dir_sub, 'delta_stats.csv')
                combined_df.to_csv(delta_stats_df, index=False)
                print(f"Combined dataframe saved to {delta_stats_df}")

                # # Calculate LLR ranking
                # PA.rank(pop,
                #         ref=qs, # label of the reference sample
                #         k=10000, # number of cells per bootstrapping sample
                #         niter=200, # number of iterations
                #         method='LLR', # LLR for log-likelihood ratio or LL for log-likelihood
                #         mincells=50, # sample's minimum number of cells to be included in ranking
                #         figsize=(10,5)) # plot figure size
                
                # # saving LLR data
                # ranking_stats_df = os.path.join(output_dir_sub, 'ranking_stats.csv')
                # pop['rankings'].to_csv(ranking_stats_df, index=True)

            except Exception as e:
                print(f"Error processing quality slice '{qs}' in directory '{output_folder}': {e}")

        # return output dir 
        pop['output'] = output_folder

        # Extracting delta mu vectors for each comparison
        for ref in quality_slices:
            for sample in quality_slices:
                # if sample == ref: 
                #     continue
            
                print(ref + ' - ' + sample)

                PA.align(pop, ref=ref,
                    method='celltype', # one of: test2ref, ref2test, conservative
                    figsizedeltas=(10,10),
                    figsizeentropy=(10,10))

                samplist = pop['order']
                celltypes = pop['samples'][ref]['gmm_types']

                # Initialize lists to store data for the DataFrame
                ref_cell_types = []
                rep_id = []
                aligned_cell_types = []
                delta_mu_vectors = []

                for i, currtype in enumerate(celltypes):  # for each reference subpopulation
                    mu_ref = PA.get_gmm_means(pop, ref, rep = 0)[i]
                    itest = PA.getalignedcompnum(pop, i, sample, rep = 0)
                    if len(itest) > 0:
                        curr_del_mu = PA.compute_delta_mu_prenorm(pop, mu_ref, itest, sample, rep = 0)
                        ref_cell_types.append(currtype)
                        aligned_cell_types.append([celltypes[j] for j in itest])
                        delta_mu_vectors.append(curr_del_mu)

                # Create the DataFrame
                go_labels = [pop['top_feat_labels']] * len(aligned_cell_types)
                df = pd.DataFrame({
                    'Reference': ref,
                    'Test': sample,
                    'Reference_Cell_Type': ref_cell_types,
                    'Aligned_Cell_Types': aligned_cell_types,
                    'Delta_Mu_Vectors': delta_mu_vectors,
                    'top_go_labels': go_labels
                })

                # Append the dataframe to the list
                all_dfs.append(df)
                
        # Concatenate all dataframes into one
        final_df = pd.concat(all_dfs, ignore_index=True)

        # Save the final dataframe to a CSV file
        csv_path =  f'{output_folder}/delta_mu_vectors_all.csv'
        final_df.to_csv(csv_path, index=False)

        # Visualize GMMs
        PA.render_models(pop, figsizegrouped=(30,30), samples=pop['order'], figsizesingle=(6,5), mode='grouped')

        # PA.samples_grid(pop, method='tsne', figsize=(20,20), size_background=.1, size_samples=.3, samplecolor='magenta', showplot=True)

    except Exception as e:
        print(f"Error processing folder {input_folder}: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python popalign_workflow.py <input_folder> <output_folder>")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    
    popalign_workflow(input_folder, output_folder)
