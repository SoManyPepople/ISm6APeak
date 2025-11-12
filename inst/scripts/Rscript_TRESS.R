options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = F)
suppressMessages(require(TRESS))

require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)

#function to run TRESS
RunTRESS <- function(InputBAM="/data/m6A_calling_strategy/7_MACS2SPeak/BAM/HEK_NEB_IVT2_Input.sortedByCoord.UN.bam",
                     RIPBAM="/data/m6A_calling_strategy/7_MACS2SPeak/BAM/HEK_NEB_IVT2_RIP.sortedByCoord.UN.bam",
                     AnnotFile="/data/m6A_calling_strategy/TRESS/UCSC.hg38.knownGene.sqlite",
                     lowcount=0,
                     binsize=50,
                     WhichThreshold = c("pval","fdr", "lfc", "pval_lfc", "fdr_lfc")[3],
                     pval.cutoff = c(NA,1,0.05,1e-3,1e-5,1e-7,1e-9),#default 1e-5
                     fdr.cutoff = c(NA,1,0.05,1e-3,1e-5),#default 0.05
                     lfc.cutoff = 0,
                     OutputDir = "/data/m6A_calling_strategy/TRESS",
                     experiment_name="HEK_NEB_IVT2",
                     IncludeIntron = FALSE){
  suppressMessages(require(TRESS))
  if(!dir.exists(OutputDir)){dir.create(OutputDir,recursive = T)}
  if(!file.exists(paste0(OutputDir,"/",experiment_name,"_peaks.xls"))){
    Input.file <- InputBAM %>% strsplit(split="/",fixed=T) %>% sapply(tail, 1)
    IP.file <- RIPBAM %>% strsplit(split="/",fixed=T) %>% sapply(tail,1)
    InputDir <- InputBAM %>% str_replace(pattern=Input.file, replacement = "")
    TRESS_peak(IP.file = IP.file,
               Input.file = Input.file,
               Path_To_AnnoSqlite = AnnotFile,
               InputDir = InputDir,
               OutputDir = OutputDir, # specify a directory for output
               experiment_name = experiment_name,
               binsize = binsize,
               lowcount=lowcount,
               WhichThreshold = WhichThreshold,
               pval.cutoff0 = pval.cutoff,
               fdr.cutoff0 = fdr.cutoff,
               lfc.cutoff0 = lfc.cutoff,
               IncludeIntron = IncludeIntron
    )
  }else{
    message(paste0("Detected existed peak file, skip TRESS peak calling for ", prefix))
  }
}

t1 <- Sys.time()
dt.parameter.combo <- foreach(IncludeIntron=c(FALSE,TRUE),.combine ='rbind')%do%{#default is FALSE
  foreach(binsize=c(50,100),.combine='rbind')%do%{#default is 50
    foreach(lowcount=c(0,10,30),.combine = 'rbind')%do%{#default is 30
      foreach(pvalue=c(NA,1,1e-2,1e-5,1e-8),.combine='rbind')%do%{#default is 1e-5
        foreach(fdr=c(NA,1,5e-2,1e-3,1e-5),.combine='rbind')%do%{#default is 0.05
          foreach(Threshold=c("pval","fdr"),.combine='rbind')%do%{#default is fdr
            data.table(IncludeIntron=IncludeIntron,binsize=binsize,lowcount=lowcount,
                       Threshold=Threshold,pval.cutoff=pvalue,fdr.cutoff=fdr)
          }
        }
      }
    }
  }
} %>% dplyr::filter((Threshold=="pval" & is.na(fdr.cutoff) & !is.na(pval.cutoff)) | (Threshold=="fdr" & is.na(pval.cutoff) & !is.na(fdr.cutoff))) %>%
  #filter(IncludeIntron==FALSE & binsize==100 & lowcount==50) %>%
  dplyr::arrange(IncludeIntron,binsize,lowcount,Threshold, pval.cutoff, fdr.cutoff) %>%  mutate(ComboName=paste0("Combo",1:nrow(.)))


bin_path <- args[1]#"/home/huangp/anaconda3/envs/m6A_seq/bin/"
InputBAM <- args[2]#"/data/m6A_calling_strategy/7_MACS2SPeak/BAM//HEK_Abcam_mRNA1_Input.sortedByCoord.UN.bam"
RIPBAM <- args[3]#"/data/m6A_calling_strategy/7_MACS2SPeak/BAM//HEK_Abcam_mRNA1_RIP.sortedByCoord.UN.bam"
outdir <- args[4]#"/data/m6A_calling_strategy/RIPPeakS/TRESS_Abcam/"
prefix <- args[5]#"HEK_Abcam_mRNA1"
SelectedCombo <- args[6]#"Combo48"
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
  RunTRESS(InputBAM=InputBAM,
           RIPBAM=RIPBAM,
           AnnotFile=Annot.sqlite,
           lowcount=dt.parameter.combo[ComboName==SelectedCombo,lowcount],
           binsize=dt.parameter.combo[ComboName==SelectedCombo,binsize],
           WhichThreshold = dt.parameter.combo[ComboName==SelectedCombo,Threshold],
           pval.cutoff = dt.parameter.combo[ComboName==SelectedCombo,pval.cutoff],#default 1e-5
           fdr.cutoff = dt.parameter.combo[ComboName==SelectedCombo,fdr.cutoff],#default 0.05
           lfc.cutoff = 0,
           OutputDir = paste0(outdir,"/TRESS_",SelectedCombo),
           experiment_name=prefix,
           IncludeIntron = dt.parameter.combo[ComboName==SelectedCombo,IncludeIntron]
           )
}
Sys.time() - t1
rm(list=ls())
q(save="no")

