options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = F)
library(mysterypackage)
getNamespaceExports("mysterypackage")
suppressMessages(library(GenomicFeatures))
detach(package:GenomicFeatures,force=T)
suppressMessages(library(txdbmaker))
suppressMessages(library(GenomicFeatures))


require(exomePeak)
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)

# RunexomePeak <- function(InputBAM="~/m6A_calling_strategy/6_alignment_star/R2_bam/mESC_NEB_WT_rep2_Input_STAR.Aligned.out.SortedByCoord.unique.deduplicate.R2.bam",
#                          RIPBAM="~/m6A_calling_strategy/6_alignment_star/R2_bam/mESC_NEB_WT_rep2_RIP_STAR.Aligned.out.SortedByCoord.unique.deduplicate.R2.bam",
#                          outdir="~/m6A_calling_strategy/exomePeak",
#                          prefix="mESC_NEB_WT_rep2_UD_test",
#                          gtf="~/genome_db/mm39/gencode.GRCm39.vM30.annotation.gtf",
#                          window.width=50,
#                          sliding.step=50,
#                          #mininal_peak_length=50,
#                          peak_cutoff_pvalue=c(1,0.05,0.01,1e-5,1e-6),#default 1e-5
#                          peak_cutoff_fdr = c(NA,1,0.05,0.01,0.001),#If it is specified, then use "fdr" instead of "p" in peak calling
#                          fold_enrichment=1
# ){
#   if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
#   if(!file.exists(paste0(outdir,"/",prefix,"/peak.xls"))){
#
#     suppressMessages(library(GenomicFeatures))
#     detach(package:GenomicFeatures,force=T)
#     suppressMessages(library(txdbmaker))
#     suppressMessages(library(GenomicFeatures))
#
#     suppressMessages(library(AnnotationDbi))
#     detach(package:AnnotationDbi,force=T)
#     suppressMessages(library(AnnotationDbi))
#     suppressMessages(library(exomePeak))
#     res <- exomePeak::exomepeak(GENE_ANNO_GTF=gtf,
#                                 IP_BAM = RIPBAM,
#                                 INPUT_BAM = InputBAM,
#                                 EXPERIMENT_NAME=prefix,
#                                 OUTPUT_DIR = outdir,
#                                 WINDOW_WIDTH=window.width,
#                                 SLIDING_STEP=sliding.step,
#                                 #MINIMAL_PEAK_LENGTH=mininal_peak_length,
#                                 PEAK_CUTOFF_PVALUE=peak_cutoff_pvalue,
#                                 PEAK_CUTOFF_FDR=peak_cutoff_fdr,
#                                 FOLD_ENRICHMENT=fold_enrichment
#     )
#   }else{message(paste0("Detected existed peak results, hence skip ", prefix))}
#
# }

dt.parameter.combo <- foreach(pvalue=c(NA,1,1e-2,1e-5,1e-8),.combine = 'rbind')%do%{#default is NA
  foreach(fdr=c(NA,1,5e-2,1e-3,1e-5),.combine='rbind')%do%{#default is 5e-2
    foreach(window.size=c(50,100,200),.combine='rbind')%do%{#default is 200
      foreach(sliding.step=c(30,50),.combine='rbind')%do%{#default is 30
        data.table(pvalue_cutoff=pvalue, fdr_cutoff=fdr, window.size=window.size,sliding.step=sliding.step)
      }
    }
  }
} %>% filter((is.na(pvalue_cutoff) & !is.na(fdr_cutoff)) | (is.na(fdr_cutoff) & !is.na(pvalue_cutoff))) %>%
  #filter(window.size==50,sliding.step==50) %>%
  mutate(ComboName=paste0("Combo",1:nrow(.)))

t1 <- Sys.time()
bin_path <- args[1]#"/home/huangp/anaconda3/envs/m6A_seq/bin/"
InputBAM <- args[2]#"/data/m6A_calling_strategy/7_MACS2SPeak/BAM//HEK_Abcam_mRNA1_Input.sortedByCoord.UN.bam"
RIPBAM <- args[3]#"/data/m6A_calling_strategy/7_MACS2SPeak/BAM//HEK_Abcam_mRNA1_RIP.sortedByCoord.UN.bam"
outdir <- args[4]#"/data/m6A_calling_strategy/RIPPeakS/exomePeak_Abcam/"
prefix <- args[5]#"HEK_Abcam_mRNA1"
SelectedCombo <- args[6]#"Combo29"
Organism <- args[7]
Annot.gtf <- args[8]
Annot.sqlite <- args[9]
Annot.genome <- args[10]


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

if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
setwd(dir = outdir)
if(SelectedCombo %in% dt.parameter.combo$ComboName){
  mysterypackage::RunexomePeak(InputBAM=InputBAM,
               RIPBAM=RIPBAM,
               outdir = paste0(outdir,"/exomePeak_", SelectedCombo),
               prefix=prefix,
               gtf=Annot.gtf,
               window.width=dt.parameter.combo[ComboName==SelectedCombo,window.size],
               sliding.step=dt.parameter.combo[ComboName==SelectedCombo,sliding.step],
               #mininal_peak_length=50,
               peak_cutoff_pvalue=dt.parameter.combo[ComboName==SelectedCombo,pvalue_cutoff],#default 1e-5
               peak_cutoff_fdr = dt.parameter.combo[ComboName==SelectedCombo,fdr_cutoff],#If it is specified, then use "fdr" instead of "p" in peak calling
               fold_enrichment=1
  )
}
Sys.time() - t1
rm(list=ls())
q(save="no")

