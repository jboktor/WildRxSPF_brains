#!/bin/bash
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=24
#SBATCH --mem=300G
#SBATCH --time=12:00:00   # time limit
#SBATCH --partition=expansion
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL   # Notify on failure.
#SBATCH --output=/resnick/groups/MazmanianLab/jboktor/snRNAseq_analysis/WILDRxSPF_brains/.cluster_runs/mapmycells/spipe_all_%j.log
#SBATCH --error=/resnick/groups/MazmanianLab/jboktor/snRNAseq_analysis/WILDRxSPF_brains/.cluster_runs/mapmycells/spipe_all_%j.err

source /home/${USER}/.bashrc
source activate mapmycells
module load cuda/12.2.1-gcc-11.3.1-sdqrj2e

ref_path="/resnick/groups/MazmanianLab/jboktor/snRNAseq_analysis/WILDRxSPF_brains/data/input/ref_genomes/MapMyCellsRefs"
output_path="/resnick/groups/MazmanianLab/jboktor/snRNAseq_analysis/WILDRxSPF_brains/data/interim/MapMyCells_trial"
anndata_path="/resnick/groups/MazmanianLab/jboktor/snRNAseq_analysis/WILDRxSPF_brains/data/input/parse_bio_anndata/anndata_ensembl.h5ad"

mkdir -p $output_path

python -m cell_type_mapper.cli.from_specified_markers \
--query_path $anndata_path \
--query_markers.serialized_lookup $ref_path/mouse_markers_230821.json \
--precomputed_stats.path $ref_path/precomputed_stats_ABC_revision_230821.h5 \
--extended_result_path $output_path/output.json \
--csv_result_path $output_path/output.csv \
--log_path $output_path/log.txt \
--drop_level CCN20230722_SUPT \
--cloud_safe False \
--type_assignment.normalization raw \
--type_assignment.n_processors 24 \
--type_assignment.chunk_size 25000 \
--type_assignment.rng_seed 42 \
--tmp_dir /resnick/scratch/jbok/tmp
