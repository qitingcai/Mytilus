#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=3-00:00:00                # Max time for job to run
#SBATCH --job-name=concat                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=qcai17@ucsc.edu      # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=20               # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Amount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # Don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-700]                     # Array job, adjust the range as needed

# Load the sample name based on the SLURM array task ID
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/sample.txt)
sample=$(echo ${LINE} | awk '{ print $2; }')

# Check if sample is empty
if [ -z "$sample" ]; then
  echo "No sample name found for task ID ${SLURM_ARRAY_TASK_ID}."
  exit 1
fi

echo "Processing sample: $sample"

# Define the input and output file paths
input_path="/hb/groups/kelley_lab/tina/mytilus/02_trim/final_trim0603"
output_path="/hb/groups/kelley_lab/tina/mytilus/02_trim/final_trim0603/concat"

# Ensure the output directory exists
mkdir -p $output_path

# Check if all expected files exist
for suffix in L001 L002 L003 L004; do
  for read in R1 R2; do
    val_suffix=${read}_001.65bp_5prime.fq.gz
    if [ "$read" == "R2" ]; then
      val_suffix=${read}_001.65bp_5prime.fq.gz
    fi
    if [ ! -f ${input_path}/${sample}_${suffix}_${val_suffix} ]; then
      echo "File ${input_path}/${sample}_${suffix}_${val_suffix} does not exist."
      exit 1
    fi
  done
done

# Concatenate R1 files
cat ${input_path}/${sample}_L001_R1_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L002_R1_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L003_R1_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L004_R1_001.65bp_5prime.fq.gz > ${output_path}/${sample}_R1_merged.fq.gz

# Check if concatenation of R1 files was successful
if [ $? -ne 0 ]; then
  echo "Concatenation of R1 files for sample ${sample} failed."
  exit 1
fi

# Concatenate R2 files
cat ${input_path}/${sample}_L001_R2_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L002_R2_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L003_R2_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L004_R2_001.65bp_5prime.fq.gz > ${output_path}/${sample}_R2_merged.fq.gz

# Check if concatenation of R2 files was successful
if [ $? -ne 0 ]; then
  echo "Concatenation of R2 files for sample ${sample} failed."
  exit 1
fi

echo "Finished processing sample: $sample"
