# Downloading fastq screen
wget https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.3.tar.gz
tar -xzf v0.15.3.tar.gz

### downloading Bowtie2 and adding bowtie to path
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip/download
unzip download
`export PATH=$HOME/tools/bowtie2/bowtie2-2.4.2-sra-linux-x86_64:$PATH`

#### adding bismark to path, for bisulfite data analysis only
`export PATH=/hb/groups/kelley_lab/tina/mytilus/Bismark-master/bismark:$PATH`

### running fastq screen
#### one example /hb/groups/kelley_lab/tina/mytilus/02_trim/final_trim/concat_lanes/12O-F_S66_R1_merged.fq.gz

### need a configuration file (indicating path of genomes) to run Fastq screen
cd /hb/groups/kelley_lab/tina/mytilus/fastqc_screen/FastQ-Screen-0.15.3
cat fastq_screen.conf

`FastQ-Screen-0.15.3/fastq_screen --bisulfite  /hb/groups/kelley_lab/tina/mytilus/02_trim/final_trim/concat_lanes/12O-F_S66_R1_merged.fq.gz`
