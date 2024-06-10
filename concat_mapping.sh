#!/bin/bash

#SBATCH --partition=128x24                # Partition/queue to run on
#SBATCH --time=0-12:00:00                # Max time for job to run
#SBATCH --job-name=bismark                  # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=qcai17@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=20                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=40G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-200]                 # array job

module load parallel miniconda3.9
conda activate bismark

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/sample_final.txt)
sample=$(echo ${LINE} | awk '{ print $2; }')
tissue=$(echo ${LINE} | awk '{ print $3; }')
replicate=$(echo ${LINE} | awk '{ print $4; }')
strand=$(echo ${LINE} | awk '{ print $5; }')

echo "running Bismark for sample: ${sample} (${tissue}, ${replicate})"

bowtie2_dir="/hb/groups/kelley_lab/tina/mytilus/bowtie2/bowtie2-2.4.2-sra-linux-x86_64"
genome_folder="/hb/groups/kelley_lab/tina/mytilus/ref_genome/GCF_021869535.1/"
trimmed_dir="/hb/groups/kelley_lab/tina/mytilus/02_trim/final_trim0603/concat"

mkdir -p final_map/${tissue}/${sample}
cd final_map/${tissue}/${sample}

bismark --p 4 ${genome_folder} \
--gzip -score_min L,0,-0.6 \
-1 ${trimmed_dir}/${sample}_R1_merged.fq.gz \
-2 ${trimmed_dir}/${sample}_R2_merged.fq.gz


