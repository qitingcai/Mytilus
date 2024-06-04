
### downloading bowtie2 from here:

wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip/download

### decompress
unzip download

### add location to system PATH
export PATH=$HOME/tools/bowtie2/bowtie2-2.4.2-sra-linux-x86_64:$PATH

### downloading bismark
wget https://github.com/FelixKrueger/Bismark/archive/master.zip
unzip master.zip

 ./bismark --version
Bismark Version: v0.24.2
