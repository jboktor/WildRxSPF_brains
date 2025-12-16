#!/bin/bash
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=1-00:00:00   # time limit
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL   # Notify on failure.
#SBATCH --output=/central/groups/MazmanianLab/joeB/snRNAseq_analysis/WILDRxSPF_brains/.cluster_runs/spipe_all_concat/spipe_all_%j.log
#SBATCH --error=/central/groups/MazmanianLab/joeB/snRNAseq_analysis/WILDRxSPF_brains/.cluster_runs/spipe_all_concat/spipe_all_%j.err

source /home/${USER}/.bashrc
source activate spipe

split-pipe --mode comb \
--nthreads 32 \
--output_dir /resnick/groups/mthomson/jboktor/split-pipe_workflow/combined_parallel \
--sublibraries /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool1_S1 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool10_S9 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool11_S10 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool12_S11 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool13_S12 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool14_S13 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool15_S14 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool16_S15 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool2_S2 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool3_S3 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool4_S4 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool5_S5 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool6_S6 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool8_S7 /resnick/groups/mthomson/jboktor/split-pipe_workflow/mode-all/Mazman_subpool9_S8
