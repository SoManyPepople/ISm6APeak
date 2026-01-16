#create conda enviroment and install required packages
conda create -n MyPackage -c conda-forge -c bioconda -c defaults macs2 bedtools rseqc samtools \
bioconda::bioconductor-exomepeak conda-forge::r-devtools
r-base=4.2 conda-forge::r-devtools bioconda::bioconductor-rtracklayer bioconda::bioconductor-genomicfeatures
conda activate MyPackage
conda install bioconda::bioconductor-exomepeak
