# Running multiqc on fastqc files

## create conda environment according to https://bioconda.github.io/recipes/multiqc/README.html#package-multiqc

> module load miniconda3.9

### can omit
> conda config --add channels defaults
> conda config --add channels bioconda
> conda config --add channels conda-forge
> conda config --set channel_priority strict

# start here
> conda create -n multiqc 
> conda activate multiqc

### whole workflow installing newest version of multiqc
> mamba create --name myenvname multiqc
> mamba install multiqc
> mamba update multiqc

### getting the help document to see running options
> multiqc --help 

# running multiqc for samples after trimming in the directory
> multiqc -v -o multiqc/ fastqc/ 