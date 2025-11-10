options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = F)
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)

#function to run MeRIPtools
RunMeRIPtools <- function(Samples="HCC1NT",#one sample
                          bamFolder="/data2/HCC_m6A/7_RIPPeakS/BAM/",#input BAM folder
                          stranded.bamFolder="/data2/HCC_m6A/RIPPeakS/tmp/",#folder that contained stranded bam that reverse to mRNA
                          OutputDir="/data2/HCC_m6A//RIPPeakS/MeRIPtools_SYSY/HCC1NT",
                          gtf.file="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",
                          samtools_path="~/anaconda3/envs/m6A_seq/bin/",
                          bedtools_path = "~/anaconda3/envs/m6A_seq/bin",
                          fragmentLength=150,
                          strandToKeep="opposite",#Illumina Truseq library kit or other similar protocol resulted in sequencing opposite strand of the RNA molecules in read1
                          binSize=c(50,100)[1],#Usually we do 50 bp bins for typical library of 20-30 Million reads. If you have deeper coverage, a smaller bin size could increase resolution of the analysis
                          threads = 60,
                          test.method=c("Fisher","Binomial")[2],
                          min_counts=c(5,15,30)[1],#minimal number of total reads (Input + IP) presented in a bin to call this bin peak
                          fdr_cutoff=c(1,0.1,0.05,1e-3)[2],
                          oddratio_cutoff=c(0,1,2)[1],#if method is Binomial, set the IP/Input ratio threshold to call a bin peak; else, represented by the (IP/gene-median-IP)/(Input/gene-median-Input),
                          genome_file="/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt"
){
  suppressMessages(require(MeRIPtools))
  suppressMessages(require(foreach))
  suppressMessages(require(tidyverse))
  suppressMessages(require(data.table))
  if(!dir.exists(OutputDir)){dir.create(OutputDir,recursive = T)}
  dir.tmp <- paste0(OutputDir,"/tmp")
  if(!dir.exists(dir.tmp)){dir.create(dir.tmp)}
  if(!file.exists(paste0(OutputDir,"/",unique(Samples),"_peak.tsv"))){
    # prepare bam folder
    Input.BAM.old <- list.files(stranded.bamFolder,full.names = T,pattern = "(.bam)$") %>% grep(pattern=unique(Samples),value=T) %>% grep(pattern="_Input",value=T)
    RIP.BAM.old <- list.files(stranded.bamFolder,full.names = T,pattern = "(.bam)$") %>% grep(pattern=unique(Samples),value=T) %>% grep(pattern="_RIP",value=T)
    InputBAM.new <- paste0(dir.tmp,"/",paste0(unique(Samples),"_REP1"),".input.bam")
    RIPBAM.new <- paste0(dir.tmp,"/",paste0(unique(Samples),"_REP1"),".m6A.bam")
    system(command = paste0("cp ", Input.BAM.old, " ", InputBAM.new),wait = T)
    system(command = paste0("cp ", RIP.BAM.old, " ", RIPBAM.new),wait = T)
    # system(command = paste0(samtools_path,"samtools view --threads 16 -s 0.5 -b  ",Input.BAM.old, " >", InputBAM.new),wait = T)
    # system(command = paste0(samtools_path,"samtools view --threads 16 -s 0.5 -b  ",RIP.BAM.old, " >", RIPBAM.new),wait = T)
    system(command = paste0("cp ", InputBAM.new, " ", str_replace(InputBAM.new,pattern="REP1",replacement="REP2")),wait = T)
    system(command = paste0("cp ", RIPBAM.new, " ", str_replace(RIPBAM.new,pattern="REP1",replacement="REP2")),wait = T)
    #count reads
    MeRIP <- countReads(samplenames = paste0(Samples,"_REP",1:2),
                        gtf = gtf.file,
                        bamFolder = dir.tmp,
                        outputDir = NA,
                        modification = "m6A",
                        strandToKeep = strandToKeep,
                        fragmentLength = fragmentLength,
                        binSize = binSize,
                        #paired = F,
                        threads = threads,
                        saveOutput = FALSE)
    #call peak
    if(test.method=="Fisher"){
      peak <- callPeakFisher(MeRIP = MeRIP, min_counts = min_counts,
                             peak_cutoff_fdr = fdr_cutoff,
                             peak_cutoff_oddRatio = oddratio_cutoff,
                             threads = threads)
    }else{
      peak <- callPeakBinomial(MeRIP = MeRIP, min_counts = min_counts,
                               peak_cutoff_fdr = fdr_cutoff,
                               peak_cutoff_oddRatio = oddratio_cutoff,
                               threads = threads)
    }
    #define joint peak
    peak <- reportJointPeak(MeRIPdata = peak, joint_threshold = 2, threads = threads)
    dt.peak <- jointPeak(peak) %>% as.data.table() %>% mutate(name=paste(chr,start,end,strand,sep="_"))
    fwrite(dt.peak[,1:6],file=paste0(dir.tmp,"/", unique(Samples),".peak.bed"),row.names = F,col.names = F,sep="\t")

    #calculate peak intensity
    dt.BAM.Depth <- data.table(Sample=paste0(unique(Samples),c("_Input","_RIP")),Depth=NA,Mapped=NA)
    dt.BAM.Depth$Depth <- c(as.integer(system(command = paste0(samtools_path,"samtools view -F 256 -c ", InputBAM.new), wait = T, intern = T)),
                            as.integer(system(command = paste0(samtools_path,"samtools view -F 256 -c ", RIPBAM.new), wait = T, intern = T)))
    dt.BAM.Depth <- dt.BAM.Depth %>% mutate(Mapped=Depth/1000000)
    dt.peak.intensity <- PeakIntensity.parallel(n.cores=threads,
                                                peakbeds=paste0(dir.tmp,"/", unique(Samples),".peak.bed"),
                                                names.peakbeds=unique(Samples),
                                                InputBAM=list.files(path=stranded.bamFolder,pattern = "(.bam)$",full.names = T) %>% grep(pattern=unique(Samples),value=T) %>%
                                                  grep(pattern="_Input",value=T),
                                                names.Input=unique(Samples),
                                                RIPBAM=list.files(path=stranded.bamFolder,pattern = "(.bam)$",full.names = T) %>% grep(pattern=unique(Samples),value=T) %>%
                                                  grep(pattern="_RIP",value=T),
                                                names.RIP=unique(Samples),
                                                Depth=dt.BAM.Depth$Mapped,
                                                names.Depth=dt.BAM.Depth$Sample,
                                                out.tmp=dir.tmp,
                                                bedtools_path = bedtools_path,
                                                bt_coverage.options=" -split -mean ",
                                                genome_file=genome_file,
                                                strandness = ifelse(strandToKeep=="same","s", "S"),
                                                pseudo_count=0.5,
                                                foldchange_cutoff=0,
                                                RIP_coverage_cutoff=0
    )
    dt.peak.intensity <- dt.peak.intensity %>% mutate(score=intensity)
    #save peak result
    fwrite(dt.peak.intensity,file=paste0(OutputDir,"/",unique(Samples),"_peak.tsv"),row.names=F,col.names=T,sep="\t")
    #remove intermediate files/dirs
    system(command = paste0("rm -rf ", dir.tmp))
  }else{message(paste0("Detect peak file, skip ", unique(Samples)))}
}

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

if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
setwd(dir = outdir)
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
