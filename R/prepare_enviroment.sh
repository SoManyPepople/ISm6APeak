
conda create -n TestMyPackage -c conda-forge -c bioconda -c defaults macs2 bedtools rseqc samtools \
r-base r-devtools bioconda::bioconductor-exomepeak2 conda-forge::r-tidyverse conda-forge::r-vcfr
conda activate TestMyPackage
# conda install bioconda::bioconductor-bsgenome.mmusculus.ucsc.mm10 bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38 \
# bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg19 bioconda::bioconductor-bsgenome.mmusculus.ucsc.mm9
R

options(repos = c( CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_ver <- BiocManager::version()       # e.g. "3.17"
bioc_ver_short <- sub("^Bioc", "", bioc_ver) # if BiocManager::version() returns "3.17" already, use it directly

repos_bioc_base <- sprintf("https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/%s", bioc_ver)
options(repos = c(
  CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
  BioCsoft = paste0(repos_bioc_base, "/bioc"),
  BioCann  = paste0(repos_bioc_base, "/data/annotation"),
  BioCexp  = paste0(repos_bioc_base, "/data/experiment")
))
getOption("repos")

remotes::install_github("SoManyPepople/mysterypackage")
library(mysterypackage)
# install.packages("~/software/BSgenome.Hsapiens.UCSC.hg19_1.4.3.tar.gz", repos=NULL,type="source")
# install.packages("~/software/BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz", repos=NULL,type="source")
# install.packages("~/software/BSgenome.Mmusculus.UCSC.mm39_1.4.3.tar.gz", repos=NULL,type="source")
install_mystery_packages(pkgs = c("TRESS","MeTPeak","exomePeak2","doParallel","foreach","reshape","qvalue","Guitar","ggsci","tidyverse","data.table",
"BSgenome.Hsapiens.UCSC.hg38","BSgenome.Mmusculus.UCSC.mm39"))
install_mystery_packages(pkgs=c("exomePeak","MeRIPtools"),github_map=c(exomePeak="ZW-xjtlu/exomePeak", MeRIPtools="scottzijiezhang/MeRIPtools"))

#test require packages
for(p in c("exomePeak", "exomePeak2", "TRESS", "MeTPeak","MeRIPtools","foreach","tidyverse","data.table")){
  library(p,character.only=T)
}

