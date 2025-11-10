options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = F)
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)

dt.parameter.combo <- foreach(pvalue=c(NA,1,0.05,0.01,0.001),.combine='rbind')%do%{#default is NA
  foreach(qvalue=c(1,0.05,0.01,0.001), .combine='rbind')%do%{#default is 5e-2
    foreach(nomodel=c(TRUE,FALSE), .combine='rbind')%do%{#default is FALSE
      foreach(nolambda=c(TRUE,FALSE),.combine='rbind')%do%{#default is FALSE
        foreach(mfold=c("-m 5 50", "-m 2 50"), .combine='rbind')%do%{#default is -m 5 50
          data.table(pvalue_cutoff=pvalue, qvalue_cutoff=qvalue, nomodel=nomodel, nolambda=nolambda, mfold=mfold)
        }
      }
    }
  }
} %>% dplyr::filter(!(!is.na(pvalue_cutoff) & qvalue_cutoff != 1)) %>% mutate(ComboName=paste0("Combo",1:nrow(.)))


bin_path <- args[1]
InputBAM <- args[2]
RIPBAM <- args[3]
outdir <- args[4]
prefix <- args[5]
SelectedCombo <- args[6]
Organism <- args[7]
Annot.gtf <- args[8]
Annot.sqlite <- args[9]
Annot.genome <- args[10]

t1<- Sys.time()

bin_path
InputBAM
RIPBAM
outdir
prefix
SelectedCombo
Organism
Annot.gtf
Annot.sqlite
Annot.genome

setwd(dir = outdir)

if(SelectedCombo %in% dt.parameter.combo$ComboName){
  RunMACS2(MACS2_path = bin_path,
           org=ifelse(Organism=="Human","hs","mm"),
           InputBAM = InputBAM,
           RIPBAM = RIPBAM,
           outdir = paste0(outdir,"/MACS2_", SelectedCombo),
           prefix=prefix,
           MACS2.options=" --keep-dup all ",
           nomodel=dt.parameter.combo[ComboName==SelectedCombo,nomodel],#deafault is FALSE, i.e, build shifting model
           mfold = dt.parameter.combo[ComboName==SelectedCombo,mfold],#Select the regions within MFOLD range of high-confidence enrichment ratio against background to build model
           #extsize=c(200,150),#default is 200
           inputformat=c("BAM","BAMPE")[1],
           pvalue_cutoff=dt.parameter.combo[ComboName==SelectedCombo,pvalue_cutoff],#default not set,mutually exclusive with qvalue,if set pvalue, qvalue will not be calculated
           qvalue_cutoff=dt.parameter.combo[ComboName==SelectedCombo,qvalue_cutoff],#default set as 0.05,ignored when pvalue is setted
           nolambda=dt.parameter.combo[ComboName==SelectedCombo,nolambda],#default is FALSE
           fe_cutoff=1,
           verbose = TRUE
           )
}
Sys.time() - t1

rm(list=ls())
q(save="no")
