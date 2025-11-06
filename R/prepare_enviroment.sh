source("~/software/mysterypackage/mysterypackage/R/install_required_r_packages.R")
install_mystery_packages(pkgs = c("exomePeak", "exomePeak2", "TRESS", "MeTPeak","MeRIPtools","foreach","tidyverse","data.table")[1])
install_mystery_packages(pkgs = c("exomePeak", "exomePeak2", "TRESS", "MeTPeak","MeRIPtools","foreach","tidyverse","data.table")[2])
install_mystery_packages(pkgs = c("exomePeak", "exomePeak2", "TRESS", "MeTPeak","MeRIPtools","foreach","tidyverse","data.table")[3])
install_mystery_packages(pkgs = c("exomePeak", "exomePeak2", "TRESS", "MeTPeak","MeRIPtools","foreach","tidyverse","data.table")[4])
install_mystery_packages(pkgs = c("exomePeak", "exomePeak2", "TRESS", "MeTPeak","MeRIPtools","foreach","tidyverse","data.table")[5])
install_mystery_packages(pkgs = c("exomePeak", "exomePeak2", "TRESS", "MeTPeak","MeRIPtools","foreach","tidyverse","data.table")[6])
install_mystery_packages(pkgs = c("exomePeak", "exomePeak2", "TRESS", "MeTPeak","MeRIPtools","foreach","tidyverse","data.table")[7])
install_mystery_packages(pkgs = c("exomePeak", "exomePeak2", "TRESS", "MeTPeak","MeRIPtools","foreach","tidyverse","data.table")[8])




conda create -n TestMyPackage -c conda-forge -c bioconda -c defaults macs2 bedtools rseqc samtools \
r-base bioconda::bioconductor-exomepeak2 conda-forge::r-tidyverse conda-forge::r-vcfr

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

install_mystery_packages(pkgs = c("doParallel","foreach","reshape","qvalue","Guitar","ggsci","tidyverse","data.table","TRESS"))
install.packages("~/software/exomePeak/exomePeak-master/",repos = NULL,type="source")
install.packages("~/software/MeTPeak/MeTPeak-master/",repos = NULL,type="source")
install_mystery_packages(pkgs="MeRIPtools",github_map=c(MeRIPtools="scottzijiezhang/MeRIPtools"))

install_mystery_packages(pkgs = c("exomePeak", "exomePeak2", "TRESS", "MeTPeak","MeRIPtools","foreach","tidyverse","data.table"))


