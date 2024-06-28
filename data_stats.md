# To understand the sequenced reads more, below are scripts used for extracting library size of the data before and after trimming, as well as cytosine coverage from bam and coverage files

## getting library size (number of bases) from raw data and reads after 65bp trimmed:

```
https://bioinf.shenwei.me/seqkit/download/
stats tool

install seqkit: 

https://bioinf.shenwei.me/seqkit/download/

seqkit stats *.gz -T > output.tsv

head output.tsv (sum_len   number of bases or residues, with gaps or spaces counted)
file	format	type	num_seqs	sum_len	min_len	avg_len	max_len
12O-F_S66_L001_R1_001.fastq.gz	FASTQ	DNA	178506	26954406	151	151.0	151
12O-F_S66_L001_R2_001.fastq.gz	FASTQ	DNA	178506	26954406	151	151.0	151
12O-F_S66_L002_R1_001.fastq.gz	FASTQ	DNA	176036	26581436	151	151.0	151
12O-F_S66_L002_R2_001.fastq.gz	FASTQ	DNA	176036	26581436	151	151.0	151
12O-F_S66_L003_R1_001.fastq.gz	FASTQ	DNA	183650	27731150	151	151.0	151
12O-F_S66_L003_R2_001.fastq.gz	FASTQ	DNA	183650	27731150	151	151.0	151
12O-F_S66_L004_R1_001.fastq.gz	FASTQ	DNA	175265	26465015	151	151.0	151
12O-F_S66_L004_R2_001.fastq.gz	FASTQ	DNA	175265	26465015	151	151.0	151
12O-G_S22_L001_R1_001.fastq.gz	FASTQ	DNA	179405	27090155	151	151.0	151

###sum the four lanes

```

# getting the positions of CGs in reference genome:

```
FIRST Getting the number of CGs in the reference genome using bismark bam2nuc
bam2nuc --genome_folder "$genome_folder" "$bam_file"

genome_folder="/hb/groups/kelley_lab/tina/mytilus/ref_genome/GCF_021869535.1/"

```



```

from Bio import SeqIO

# Define the window size
window_size = 2

# Open the reference genome FASTA file
with open("GCF_021869535.1_xbMytCali1.0.p_genomic.fa", "r") as fasta_file:
    # Iterate through each sequence record in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).lower()  # Convert sequence to lowercase
        chromosome = record.id

        # Iterate over the sequence to find CpG sites and create windows around them
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i+window_size]
            if "cg" in window:  # Check for lowercase "cg"
                # Annotate the window as a CpG site and print to BED format
                start = i + 1
                end = i + window_size
                print(f"{chromosome}\t{start}\t{end}\tCpG_site")

#Running the python script

#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=3-00:00:00                # Max time for job to run
#SBATCH --job-name=extract_cpg                  # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=qcai17@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=20                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=30G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

python extract_cpg.py  > annotated_new_CpG_sites.bed


###sort the bed file:
sort -k1,1 -k2,2n -k3,3n annotated_new_CpG_sites.bed -o sorted_annotated_new_CpG_sites.bed

head sorted_annotated_new_CpG_sites.bed 
NW_026262581.1	171	172	CpG_site
NW_026262581.1	270	271	CpG_site
NW_026262581.1	346	347	CpG_site
NW_026262581.1	353	354	CpG_site
NW_026262581.1	356	357	CpG_site
NW_026262581.1	402	403	CpG_site
NW_026262581.1	413	414	CpG_site
NW_026262581.1	425	426	CpG_site
NW_026262581.1	432	433	CpG_site
NW_026262581.1	436	437	CpG_site


#Using bedtools and samtools to get intersection with bed file and the read depth at each position
 nano cpg.sh 
#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=3-00:00:00                # Max time for job to run
#SBATCH --job-name=extract_cpg                  # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=qcai17@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=20                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=30G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

# Define the base directory where all folders are located
BASE_DIR="/hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603/final_map_min0.6/F"
ANNOTATION_FILE="$BASE_DIR/sorted_annotated_new_CpG_sites.bed"

# Find all BAM files with the specific naming pattern within the subdirectories of the base directory
find "$BASE_DIR" -type f -name "*_merged_bismark_bt2_pe.bam" | while read -r BAM_FILE; do
  # Get the directory of the current BAM file
  DIR=$(dirname "$BAM_FILE")

  # Define the output file names based on the current directory
  SORTED_BAM_FILE="$DIR/$(basename "$BAM_FILE" .bam)_sorted.bam"
  OVERLAP_OUTPUT="$DIR/cpg_reads_overlap_sorted.bed"
  COVERAGE_OUTPUT="$DIR/cpg_coverage.txt"

  # Sort the BAM file if not already sorted
  if [[ ! -f "$SORTED_BAM_FILE" ]]; then
    samtools sort "$BAM_FILE" -o "$SORTED_BAM_FILE"
    samtools index "$SORTED_BAM_FILE"
  fi

  # Run bedtools intersect
  bedtools intersect -a "$ANNOTATION_FILE" -b "$SORTED_BAM_FILE" -wa -wb | sort -k1,1 -k2,2n > "$OVERLAP_OUTPUT"

  # Run samtools depth
  samtools depth -a -b "$OVERLAP_OUTPUT" "$SORTED_BAM_FILE" > "$COVERAGE_OUTPUT"

  echo "Processed directory: $DIR"
done

echo "All directories processed."

###output is cpg_coverage.txt for each sample:
head cpg_coverage.txt 
NW_026262581.1	97827	2
NW_026262581.1	97834	2
NW_026262581.1	97838	2
NW_026262581.1	97845	4
NW_026262581.1	97858	4
NW_026262581.1	97861	4
NW_026262581.1	97867	4
NW_026262581.1	97890	4
NW_026262581.1	97903	2
NW_026262581.1	278953	1

##Counting the number of CG sites in each sample--> counting the number of lines in the txt file
nano extract_cp.sh
#!/bin/bash

# Define the base directory
BASE_DIR="/hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603/final_map_min0.6/F"

# Define the output CSV file
OUTPUT_CSV="cpg_coverage_stats_new.tsv"

# Find all cpg_coverage.txt files and count lines
echo "file_path,line_count" > $OUTPUT_CSV
find $BASE_DIR -type f -name "cpg_coverage.txt" | while read -r file; do
subdirectory=$(dirname "$file" | sed "s|$BASE_DIR/||")
line_count=$(wc -l < "$file")
    echo -e "$subdirectory\t$file\t$line_count" >> $OUTPUT_CSV
done

####head cpg_coverage.txt 
NW_026262581.1	97827	2
NW_026262581.1	97834	2
NW_026262581.1	97838	2
NW_026262581.1	97845	4
NW_026262581.1	97858	4
NW_026262581.1	97861	4
NW_026262581.1	97867	4
NW_026262581.1	97890	4
NW_026262581.1	97903	2
NW_026262581.1	278953	1


#counting the number of reads aligned to each CG sites --> getting the sum of the third column in the txt file
nano count_reads.sh
#!/bin/bash

# Define the directory where cpg_coverage.txt files are located
directory="/hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603/final_map_min0.6/F"

# Output file name for summaries
output_file="coverage_summaries.tsv"

# Find all cpg_coverage.txt files and process each
find "$directory" -name "cpg_coverage.txt" | while read -r file; do
    # Extract sample name from the directory name
    sample=$(basename "$(dirname "$file")")
    sum=$(awk '{ sum += $3 } END { print sum }' "$file")
    echo -e "$sample\t$sum" >> "$output_file"
done



###counting number of CG sites in the coverage files --> counting the number of lines in each of the coverage file
#### example of coverage file after coverager2cysosine
head PLB7-G_S55.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov
NW_026262581.1	97825	97827	0.000000	0	3
NW_026262581.1	97832	97834	0.000000	0	3
NW_026262581.1	97836	97838	0.000000	0	3
NW_026262581.1	97843	97845	0.000000	0	3
NW_026262581.1	97856	97858	0.000000	0	3
NW_026262581.1	97859	97861	0.000000	0	3
NW_026262581.1	97865	97867	0.000000	0	3
NW_026262581.1	97888	97890	0.000000	0	3
NW_026262581.1	97901	97903	0.000000	0	3
NW_026262581.1	278951	278953	0.000000	0	4

#### extract number of lines and saving as tsv
find . -type f -name "*.merged_CpG_evidence.cov" -exec wc -l {} + | awk '{print $1 "\t" $2}' > wc_results.tsv

### counting number of CG sites in the coverage files --> summing the coverage (methylated + unmethylated occurances)
R code:

file_list <- list.files(input_dir, pattern = ".CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov", full.names = TRUE)

compute_coverage_statistics <- function(file) {
df<- read_tsv(file, col_names = FALSE)
colnames(df) <- c("chromosome", "start position", "end position", "methylation percentage", "count methylated", "count unmethylated", "coverage")
# Calculate coverage as the sum of methylated and unmethylated counts
compute_coverage_statistics <- function(file) {
df <- df %>%
mutate(coverage = `count methylated` + `count unmethylated`)
coverage_stats <- list( sum_reads =sum(df$coverage))
  return(coverage_stats)
}

results <- lapply(file_list, compute_coverage_statistics)
# Combine results into a single data frame
results_df <- bind_rows(results)

# Output CSV file path
output_csv <- "/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/cg_coverage_statistics.csv"

# Write results to CSV
write_csv(results_df, output_csv)


```

## all stats output: https://docs.google.com/spreadsheets/d/1Dk6OWD8RvfzZCB_5r3vdp2GrWqkdKjiXGEAURNU_Eu0/edit?usp=sharing

