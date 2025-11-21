#!/bin/bash
# install minimap2, seqkit, and blast
sudo apt-get install -y ncbi-blast+ seqkit minimap2
Rscript -e "install.packages(c('R.utils','data.table','dplyr','reshape2'), repos='https://cloud.r-project.org/')"

# install cutadapt
pip install cutadapt

# Download and extract samtools
mkdir tools
cd tools
curl -L -O https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2
tar -jxf samtools-1.22.tar.bz2
cd samtools-1.22
make
