#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=3-00:00:00                # Max time for job to run
#SBATCH --job-name=cpg                  # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=qcai17@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=20                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=20G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-90]                 # array job

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /hb/groups/kelley_lab/tina/mytilus/04_methylation_extractor/final_minscore0.6/sample.txt)
sample=$(echo ${LINE} | awk '{ print $1; }')

find "${sample}"_R1_merged_bismark_bt2_pe.bismark.cov \
| xargs basename -s _R1_merged_bismark_bt2_pe.bismark.cov \
| xargs -I{} coverage2cytosine \
--genome_folder /hb/groups/kelley_lab/tina/mytilus/ref_genome/GCF_021869535.1/ \
-o {}.CpG_report.merged_CpG_evidence.cov \
--merge_CpG \
--zero_based \
"${sample}"_R1_merged_bismark_bt2_pe.bismark.cov
