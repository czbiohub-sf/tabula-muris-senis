#!/bin/bash
#SBATCH -c 4                               # Request 4 core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                          # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --array=3,4,6,8,9,11                  # Partition to run in
#SBATCH --mem=96000                          # Memory total in MB (for all cores)
#SBATCH -o //home/jz286/code/tabula-muris-senis/2_aging_signature/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/code/tabula-muris-senis/2_aging_signature/job_info/hostname_%j.err                 
#SBATCH --mail-type=None                    # Type of email notification- BEGIN,END,FAIL,ALL

ROOT_FOLDER="/n/groups/price/martin/tms_gene_data"

# i_tissue=3
i_tissue=$SLURM_ARRAY_TASK_ID
tissue=$(head -n ${i_tissue} "${ROOT_FOLDER}/DGE_result/droplet.tissue_list" | tail -1 )
echo $tissue
# tissue='toy'

OUTPUT_FOLDER="${ROOT_FOLDER}/DGE_result/DGE_droplet_tissue.tpm"
n_gene="F"
DATA_NAME="droplet"

Rscript DGE_analysis.R $OUTPUT_FOLDER $DATA_NAME $n_gene $tissue