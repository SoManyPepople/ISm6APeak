options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = F)
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)


dt.parameter.combo <-  foreach(pvalue=c(NA,1,1e-2,1e-5,1e-8),.combine = 'rbind')%do%{#default is NA
  foreach(fdr=c(NA,1,5e-2,1e-3,1e-5),.combine='rbind')%do%{#default is 5e-2
    foreach(window.size=c(50,100),.combine='rbind')%do%{#default is 50
      foreach(sliding.step=c(30,50),.combine='rbind')%do%{#default is 50
        data.table(pvalue_cutoff=pvalue, fdr_cutoff=fdr, window.size=window.size,sliding.step=sliding.step)
      }
    }
  }
} %>% filter((is.na(pvalue_cutoff) & !is.na(fdr_cutoff)) | (is.na(fdr_cutoff) & !is.na(pvalue_cutoff))) %>%
  #filter(window.size==50,sliding.step==50) %>%
  mutate(ComboName=paste0("Combo",1:nrow(.)))


bin_path <- args[1]#"/home/huangp/anaconda3/envs/m6A_seq/bin/"
InputBAM <- args[2]#"/data/m6A_calling_strategy/7_MACS2SPeak/BAM//HEK_Abcam_mRNA1_Input.sortedByCoord.UN.bam"
RIPBAM <- args[3]#"/data/m6A_calling_strategy/7_MACS2SPeak/BAM//HEK_Abcam_mRNA1_RIP.sortedByCoord.UN.bam"
outdir <- args[4]#"/data/m6A_calling_strategy/RIPPeakS/MeTPeak_Abcam/"
prefix <- args[5]#"HEK_Abcam_mRNA1"
SelectedCombo <- args[6]#"Combo19"
Organism <- args[7]
Annot.gtf <- args[8]
Annot.sqlite <- args[9]
Annot.genome <- args[10]

t1 <- Sys.time()

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
  RunMeTPeak(InputBAM=InputBAM,
             RIPBAM=RIPBAM,
             outdir = paste0(outdir,"/MeTPeak_", SelectedCombo),
             prefix=prefix,
             gtf=Annot.gtf,
             window.width=dt.parameter.combo[ComboName==SelectedCombo,window.size],
             sliding.step=dt.parameter.combo[ComboName==SelectedCombo,sliding.step],
             #mininal_peak_length=50,
             peak_cutoff_pvalue=dt.parameter.combo[ComboName==SelectedCombo,pvalue_cutoff],#default 1e-5
             peak_cutoff_fdr =dt.parameter.combo[ComboName==SelectedCombo,fdr_cutoff] ,#If it is specified, then use "fdr" instead of "p" in peak calling
             fold_enrichment=1)
}
Sys.time() - t1

rm(list=ls())
q(save="no")
