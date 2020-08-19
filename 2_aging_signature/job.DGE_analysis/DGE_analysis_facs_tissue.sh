#!/bin/bash
#SBATCH -c 1                               # Request 4 core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-6:00                          # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --array=1-23                    # Partition to run in
#SBATCH --mem=32000                          # Memory total in MB (for all cores)
#SBATCH -o //home/jz286/code/tabula-muris-senis/2_aging_signature/job_info_debug/hostname_%j.out  
#SBATCH -e /home/jz286/code/tabula-muris-senis/2_aging_signature/job_info_debug/hostname_%j.err                 
#SBATCH --mail-type=None                    # Type of email notification- BEGIN,END,FAIL,ALL

ROOT_FOLDER="/n/groups/price/martin/tms_gene_data"

# i_tissue=4
i_tissue=$SLURM_ARRAY_TASK_ID
tissue=$(head -n ${i_tissue} "${ROOT_FOLDER}/DGE_result/facs.tissue_list" | tail -1 )
echo $tissue

OUTPUT_FOLDER="${ROOT_FOLDER}/DGE_result/DGE_facs_tissue.tpm"
n_gene="F"
DATA_NAME="facs"

Rscript DGE_analysis.R $OUTPUT_FOLDER $DATA_NAME $n_gene $tissue