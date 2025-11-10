options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = F)
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)

#function to run exomePeak2
RunexomePeak2 <- function(InputBAM="/data2/HCC_m6A/7_RIPPeakS/BAM//HCC1NT_Input.sortedByCoord.UN.bam",
                          RIPBAM="/data2/HCC_m6A/7_RIPPeakS/BAM//HCC1NT_RIP.sortedByCoord.UN.bam",
                          outdir="/data2/HCC_m6A/RIPPeakS//exomePeak2_SYSY",
                          prefix="HCC1NT",
                          strandness=c("unstrand","1st_strand","2nd_strand")[3],
                          gtf="~/genome_db/gencode.vM33.annotation.gtf",
                          mode=c("exon","full_transcript", "whole_genome")[1],
                          genome=c("mm39","hg38")[1],
                          bin.size=100,
                          step.size=25,
                          test.method=c("Poisson","DESeq2")[1],#only Poisson could work
                          p.cutoff=c(1,5e-2,1e-3,1e-5,1e-10,1e-12)[3],#default 1e-10
                          parallel=1){
  if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
  if(!file.exists(paste0(outdir,"/",paste0(prefix,"__",mode, "__",test.method),"/peaks.bed"))){
    suppressMessages(library(AnnotationDbi))
    detach(package:AnnotationDbi,force=T)
    suppressMessages(library(AnnotationDbi))
    suppressMessages(library(exomePeak2))
    res <- exomePeak2(gff=gtf,
                      bam_ip = RIPBAM,
                      bam_input = InputBAM,
                      genome=genome,
                      mode = mode,
                      strandness = strandness,#determine the strandness option using all three parameters
                      experiment_name=paste0(prefix,"__",mode, "__",test.method),
                      save_dir = outdir,
                      bin_size = bin.size,
                      step_size = step.size,
                      test_method = test.method,
                      p_cutoff = p.cutoff
    )
  }else{message(paste0("Detected existed peak results, hence skip ", paste0(prefix,"__",mode, "__",test.method)))}
}

dt.parameter.combo <- foreach(mode=c("exon","full_transcript"),.combine='rbind')%do%{#default is exon
  foreach(bin.size=c(25,50,100),.combine='rbind')%do%{#default is 25
    foreach(step.size=c(25,50),.combine='rbind')%do%{#default is 25
      foreach(pvalue=c(1,1e-3,1e-10,1e-12),.combine='rbind')%do%{#default is 1e-10
        data.table(mode=mode,bin.size,step.size,pvalue_cutoff=pvalue)
      }
    }
  }
} %>% dplyr::filter(bin.size>=step.size) %>%
  #filter(mode=="exon" & bin.size==100 & step.size==30) %>%
  dplyr::arrange(mode,bin.size,step.size,pvalue_cutoff) %>%
  mutate(ComboName=paste0("Combo",1:nrow(.)))
dt.parameter.combo <- dt.parameter.combo %>% dplyr::filter(mode != "whole_genome")#too slow to run

#
bin_path <- as.character(args[1])#"/home/huangp/anaconda3/envs/m6A_seq/bin/"
InputBAM <- as.character(args[2])#"/data/m6A_calling_strategy/7_MACS2SPeak/BAM//HEK_Abcam_mRNA1_Input.sortedByCoord.UN.bam"
RIPBAM <- as.character(args[3])#"/data/m6A_calling_strategy/7_MACS2SPeak/BAM//HEK_Abcam_mRNA1_RIP.sortedByCoord.UN.bam"
outdir <- args[4]#"/data/m6A_calling_strategy/RIPPeakS/exomePeak2_Abcam/"
prefix <- args[5]#"HEK_Abcam_mRNA1"
SelectedCombo <- args[6]#"Combo35"
Organism <- args[7]
Annot.gtf <- args[8]
Annot.sqlite <- args[9]
Annot.genome <- args[10]
Strandness <- as.character(args[11])
nThread <- as.integer(args[12])
# bin_path <- "/home/huangp/anaconda3/envs/m6A_seq/bin/"
# InputBAM <- "/data2/HCC_m6A/7_RIPPeakS/BAM//HCC1NT_Input.sortedByCoord.UN.bam"
# RIPBAM <- "/data2/HCC_m6A/7_RIPPeakS/BAM//HCC1NT_RIP.sortedByCoord.UN.bam"
# outdir <- "/data2/HCC_m6A/RIPPeakS//exomePeak2_SYSY"
# prefix <- "HCC1NT"
# SelectedCombo <- "Combo35"
# Organism <- "Human"
# Annot.gtf <- "/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf"
# Annot.sqlite <- "/data/m6A_calling_strategy/RIPPeakS_resource/UCSC.hg38.knownGene.sqlite"
# Annot.genome <- "/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt"
# Strandness <- "S"

stranded.bamfolder <- paste0(outdir %>% strsplit(split="exomePeak2",fixed=T) %>% sapply("[",1),"/tmp")


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
stranded.bamfolder
nThread

if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
setwd(dir = outdir)
t1 <- Sys.time()
if(SelectedCombo %in% dt.parameter.combo$ComboName){
  RunexomePeak2(
    InputBAM=list.files(stranded.bamfolder,pattern="(.bam)$",full.names = T) %>% grep(pattern="_Input",value=T) %>% grep(pattern=prefix,value=T),#stranded bam reads same to mRNA
    RIPBAM=list.files(stranded.bamfolder,pattern="(.bam)$",full.names = T) %>% grep(pattern="_RIP",value=T) %>% grep(pattern=prefix,value=T),
    outdir=paste0(outdir,"/exomePeak2_",SelectedCombo),
    prefix=prefix,
    strandness=ifelse(Strandness=="s", "2nd_strand", "1st_strand"),
    gtf=Annot.gtf,
    mode=dt.parameter.combo[ComboName==SelectedCombo,mode],
    genome=ifelse(Organism=="Human","hg38","mm39"),
    bin.size=dt.parameter.combo[ComboName==SelectedCombo,bin.size],
    step.size=dt.parameter.combo[ComboName==SelectedCombo,step.size],
    test.method=c("Poisson","DESeq2")[1],
    p.cutoff=dt.parameter.combo[ComboName==SelectedCombo,pvalue_cutoff] ,
    parallel=nThread
  )
}
Sys.time() - t1

rm(list=ls())
q(save="no")
