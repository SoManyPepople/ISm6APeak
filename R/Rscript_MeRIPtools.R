options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = F)
setwd("/data/m6A_calling_strategy")
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)
source("/data/m6A_calling_strategy/Script/Wrapper_function_for_m6A_calling_using_MACS2SMePeak.R")

dt.parameter.combo <-  foreach(binSize=c(50),.combine = 'rbind')%do%{#default is 50
  foreach(test.method=c("Fisher","Binomial"),.combine='rbind')%do%{#default is unknown
    foreach(min_counts=c(5,15),.combine='rbind')%do%{#default is 15
      foreach(fdr_cutoff=c(0.1,0.05,1e-2),.combine='rbind')%do%{#default is 0.05
        foreach(oddratio_cutoff=c(0,1),.combine='rbind')%do%{#default is 1
          data.table(binSize=binSize, test.method=test.method, min_counts=min_counts, fdr_cutoff=fdr_cutoff, oddratio_cutoff=oddratio_cutoff)
        }
      }
    }
  }
} %>% mutate(ComboName=paste0("Combo",1:nrow(.))) %>% mutate(IsDefaultCombo = binSize==50 & min_counts==15 & fdr_cutoff==0.05 & oddratio_cutoff==1) %>%
  dplyr::arrange(desc(IsDefaultCombo))


bin_path <- args[1]#"/home/huangp/anaconda3/envs/m6A_seq/bin/"
InputBAM <- args[2]#"/data/m6A_calling_strategy/7_MACS2SPeak/BAM//HEK_Abcam_mRNA1_Input.sortedByCoord.UN.bam"
RIPBAM <- args[3]#"/data/m6A_calling_strategy/7_MACS2SPeak/BAM//HEK_Abcam_mRNA1_RIP.sortedByCoord.UN.bam"
outdir <- args[4]#"/data/m6A_calling_strategy/RIPPeakS/MeRIPtools_Abcam/"
prefix <- args[5]#"HEK_Abcam_mRNA1"
SelectedCombo <- args[6]#"Combo13"
Organism <- args[7]
Annot.gtf <- args[8]
Annot.sqlite <- args[9]
Annot.genome <- args[10]

Strandness <- args[11]
nThread <- as.integer(args[12])

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
Strandness
nThread

t1 <- Sys.time()
if(SelectedCombo %in% dt.parameter.combo$ComboName){
  RunMeRIPtools(
    Samples=prefix,
    bamFolder=InputBAM %>% strsplit(split=prefix,fixed=T) %>% sapply("[",1),
    stranded.bamFolder=paste0(outdir %>% strsplit(split="/MeRIPtools",fixed=T) %>% sapply("[",1),"/tmp"),#stranded bam directory
    OutputDir=paste0(outdir,"/MeRIPtools_", SelectedCombo,"/", prefix),
    gtf.file=Annot.gtf,
    samtools_path = bin_path,
    bedtools_path = bin_path,
    fragmentLength=150,
    strandToKeep=ifelse(Strandness=="s","same","opposite"),
    binSize=dt.parameter.combo[ComboName==SelectedCombo,binSize],
    threads = nThread,
    test.method=dt.parameter.combo[ComboName==SelectedCombo,test.method],
    min_counts=dt.parameter.combo[ComboName==SelectedCombo,min_counts],#
    fdr_cutoff=dt.parameter.combo[ComboName==SelectedCombo,fdr_cutoff],
    oddratio_cutoff=dt.parameter.combo[ComboName==SelectedCombo,oddratio_cutoff],#
    genome_file=Annot.genome
  )           
}
Sys.time() - t1

rm(list=ls())
q(save="no")
