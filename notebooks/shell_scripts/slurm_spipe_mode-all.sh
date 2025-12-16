#!/bin/bash
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=36
#SBATCH --mem=300G
#SBATCH --time=2-00:00:00   # time limit
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL   # Notify on failure.
#SBATCH --output=/central/groups/MazmanianLab/joeB/snRNAseq_analysis/WILDRxSPF_brains/.cluster_runs/spipe_all_concat/spipe_all_%j.log
#SBATCH --error=/central/groups/MazmanianLab/joeB/snRNAseq_analysis/WILDRxSPF_brains/.cluster_runs/spipe_all_concat/spipe_all_%j.err

# Parse command line arguments
while getopts "1:2:o:" opt; do
  case $opt in
    1) FASTQ1="$OPTARG";;
    2) FASTQ2="$OPTARG";;
    o) OUTPUT_DIR="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1;;
  esac
done

# Check if required arguments are provided
if [ -z "$FASTQ1" ] || [ -z "$FASTQ2" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 -1 <fastq1> -2 <fastq2> -o <output_dir>"
    exit 1
fi

source /home/${USER}/.bashrc
source activate spipe

# Run pipeline
split-pipe \
--mode all \
--chemistry v3 \
--kit WT_mega \
--nthreads 36 \
--fq1 "$FASTQ1" \
--fq2 "$FASTQ2" \
--output_dir "$OUTPUT_DIR" \
--genome_dir /resnick/groups/mthomson/jboktor/split-pipe_workflow/GRCm39_v113 \
--sample "WT11" A1-A4 \
--sample "ASO11" A5-A7 \
--sample "WT12" A8-A10 \
--sample "ASO13" A11-B1 \
--sample "WT13" B2-B4 \
--sample "ASO15" B5-B7 \
--sample "WT14" B8-B10 \
--sample "ASO16" B11-C1 \
--sample "WT15" C2-C4 \
--sample "ASO17" C5-C7 \
--sample "WT17" C8-C10 \
--sample "ASO20" C11-D1 \
--sample "WildR 11 AMY" D2-D4 \
--sample "SPF 9 AMY" D5-D7 \
--sample "WildR 10 AMY" D8-D10 \
--sample "SPF 11 AMY" D11-E1 \
--sample "WildR 11 HYP" E2-E4 \
--sample "SPF 10 HYP" E5-E7 \
--sample "WildR 7 AMY" E8-E10 \
--sample "SPF 7 AMY" E11-F1 \
--sample "WildR 7 HYP" F2-F4 \
--sample "SPF 8 HYP" F5-F7 \
--sample "SPF 8 AMY" F8-F10 \
--sample "WildR 6 AMY" F11-F12 \
--sample "SPF 6 HYP" G1-G2 \
--sample "WildR 6 HYP" G3-G4 \
--sample "WildR 9 AMY" G5-G6 \
--sample "SPF 6 AMY" G7-G8 \
--sample "WildR 5 AMY" G9-G10 \
--sample "WildR 10 HYP" G11-G12 \
--sample "SPF 3D4 HYP" H1-H2 \
--sample "SPF 9 HYP" H3-H4 \
--sample "WildR 9 HYP" H5-H6 \
--sample "WildR 3D4 HYP" H7-H8 \
--sample "SPF 1D4 HYP" H9-H10 \
--sample "SPF 3D4 AMY" H11-H12
