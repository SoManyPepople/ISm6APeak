options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = F)
suppressMessages(require(data.table))
suppressMessages(require(tidyverse))
require(foreach)
require(doParallel)
require(IRanges)
require(bedtoolsr)

#bedtools intensity for one sample
BedtoolsIntensity <- function(input.m6a.sample = "/data/m6A_calling_strategy/M6APeakS_Abcam/tmp/MACS2_Combo5_bedtools_intensity_input_m6A_HEK_Abcam_mRNA1.tsv",
                              prefix="HEK_SYSY_mRNA1",
                              RIPLibraryScaleFactor=paste0(dt.stranded.BAM.Depth.scale.factor[Replicate=="HEK_SYSY_mRNA1",round(RIP,2)],"/",
                                                           dt.stranded.BAM.Depth.scale.factor[Replicate=="HEK_SYSY_mRNA1",round(Input,2)]),
                              stranded.bam.dir="/data/m6A_calling_strategy/RIPPeakS/tmp/",
                              bedtools_path = "/home/huangp/anaconda3/envs/m6A_seq/bin/",
                              dir.tmp=paste0("/data/m6A_calling_strategy/RIPPeakS/tmp/","HEK_SYSY_mRNA1_SYSY_exomePeak2_intensity_tmp"),
                              #threads = 10,
                              genome_file="/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt",#org specific parameter
                              intron.bed="/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.hg38.bed",#org specific parameter
                              intronremoved = c(TRUE,FALSE)[1],
                              bedtools_option= c("counts ", "mean ")[1],
                              binsize=500,
                              Pseudocount=1,
                              SummaryMethod="max"
){
  if(!dir.exists(dir.tmp)){dir.create(dir.tmp)}
  #collect both strand m6a
  dt.input.m6a <- fread(input.m6a.sample,header=T)
  # dt.input.m6a <- rbind(dt.input.m6a %>% mutate(strand = "+"), dt.input.m6a %>% mutate(strand = "-")) %>% dplyr::arrange(seqnames,start,end,strand)
  # dt.m6a.sample.sorted <- bedtoolsr::bt.sort(dt.input.m6a) %>% as.data.table()
  # colnames(dt.m6a.sample.sorted) <- colnames(dt.input.m6a)
  # dt.m6a.sample.merge <- rbind(
  #   bedtoolsr::bt.merge(dt.m6a.sample.sorted %>% dplyr::filter(strand=="+"),s = T,d = 5, c=5, o="mean") %>% as.data.table() %>%
  #     mutate(seqnames=V1,start=V2,end=V3,strand="+",score=V4) %>% mutate(name=paste(seqnames,start,end,strand,sep="_")) %>%
  #     dplyr::select(seqnames,start,end,name,score,strand),
  #   bedtoolsr::bt.merge(dt.m6a.sample.sorted %>% dplyr::filter(strand=="-"),s = T,d = 5, c=5, o="mean") %>% as.data.table() %>%
  #     mutate(seqnames=V1,start=V2,end=V3,strand="-",score=V4) %>% mutate(name=paste(seqnames,start,end,strand,sep="_")) %>%
  #     dplyr::select(seqnames,start,end,name,score,strand)
  # )
  dt.m6a.both.strand <- rbind(dt.input.m6a %>% mutate(strand = "+"), dt.input.m6a %>% mutate(strand = "-")) %>% dplyr::arrange(seqnames,start,end,strand)
  #intron removed if yes
  if(intronremoved){
    message(paste0("intron remove for ", prefix))
    dt.m6a.sample.intron.removed <- LongPeakRemoveIntron(dt.peak = dt.m6a.both.strand,
                                                         genome_file = genome_file,
                                                         dt.intron = fread(intron.bed,header = F))
    dt.m6a.sample.intron.removed <- dt.m6a.sample.intron.removed %>% dplyr::arrange(seqnames,start,end,strand)
    #divide peaks into small bins
    if(binsize != 0){
      gr.peaks <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.sample.intron.removed,keep.extra.columns = T)
      gr.tiles <- IRanges::tile(gr.peaks,width = binsize)
      names(gr.tiles) <- gr.peaks$name
      gr.tiles <- unlist(gr.tiles,use.names = T)
      dt.tiles <- as.data.table(gr.tiles)
      dt.tiles$name <- names(gr.tiles)
      dt.tiles <- dt.tiles %>% left_join(x=.,y=dt.m6a.sample.intron.removed %>% dplyr::distinct(name,score))
      dt.tiles <- dt.tiles %>% dplyr::select(seqnames,start,end,name,score,strand)
    }else{
      dt.tiles <- dt.m6a.sample.intron.removed %>% dplyr::distinct(seqnames,start,end,name,score,strand)
    }
  }else{
    #intron-kept
    dt.m6a.sample.intron.kept <- dt.m6a.both.strand
    #divide peaks into small bins
    if(binsize != 0){
      gr.peaks <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.sample.intron.kept,keep.extra.columns = T)
      gr.tiles <- IRanges::tile(gr.peaks,width = binsize)
      names(gr.tiles) <- gr.peaks$name
      gr.tiles <- unlist(gr.tiles,use.names = T)
      dt.tiles <- as.data.table(gr.tiles)
      dt.tiles$name <- names(gr.tiles)
      dt.tiles <- dt.tiles %>% left_join(x=.,y=dt.m6a.sample.intron.kept %>% dplyr::distinct(name,score))
      dt.tiles <- dt.tiles %>% dplyr::select(seqnames,start,end,name,score,strand) %>% distinct()
    }else{
      dt.tiles <- dt.m6a.sample.intron.kept %>% dplyr::select(seqnames,start,end,name,score,strand) %>% distinct()
    }
  }
  #write bed into dir tmp
  fwrite(dt.tiles[,1:6],file=paste0(dir.tmp,"/", prefix,".peak.both.strand.bed"),row.names = F, col.names=F,sep="\t")
  #intensity calculation
  message(paste0("Peak intensity estimation  of ", prefix))
  dt.parameter.bedtools.intensity <- data.table(peakbeds=rep(paste0(dir.tmp,"/", prefix,".peak.both.strand.bed"),2),
                                                BAMs=c(list.files(stranded.bam.dir,pattern = "(Input).+(.bam)$",full.names = T) %>% grep(pattern=prefix,value=T),
                                                       list.files(stranded.bam.dir,pattern = "(RIP).+(.bam)$",full.names = T) %>% grep(pattern=prefix,value=T)),
                                                coverage.txt=paste0(dir.tmp,"/",prefix,c("_Input","_RIP"), ".both.strand.peak.coverage.txt"),
                                                log.file=paste0(dir.tmp,"/",prefix,c("_Input","_RIP"), ".both.strand.peak.coverage.log")
  )
  fwrite(dt.parameter.bedtools.intensity, file = paste0(dir.tmp,"/parameter.bedtools.intensity.parallel.txt"), sep="\t",row.names=F,col.names=F)
  parallel.bedtools.intensity.cmd <- paste0("parallel -j 2  --will-cite -a ", paste0(dir.tmp,"/parameter.bedtools.intensity.parallel.txt") ," --colsep '\t' '", paste0(bedtools_path,"bedtools coverage "),
                                            bedtools_option, " -a {1} -b {2} -g ", genome_file," 1>{3} 2>{4}'")
  system(command = parallel.bedtools.intensity.cmd, wait = T)
  #calculate the pos and negative strand (intron-removed intensity)
  message(paste0("new intensity "))
  dt.peak.intensity <- list(Input=fread(paste0(dir.tmp,"/",prefix,c("_RIP"), ".both.strand.peak.coverage.txt")) %>% mutate(Library="RIP") %>% distinct(),
                            RIP=fread(paste0(dir.tmp,"/",prefix,c("_Input"), ".both.strand.peak.coverage.txt")) %>% mutate(Library="Input") %>% distinct()) %>% Reduce(x=.,rbind)
  colnames(dt.peak.intensity) <- c("seqnames","start","end","name","score","strand","depth","Library")
  dt.peak.intensity <- dt.peak.intensity %>% dplyr::arrange(name,start,end,strand) %>% pivot_wider(values_from = c("depth"),names_from = c("Library")) %>%
    as.data.table() %>% mutate(RIP=RIP/as.numeric(RIPLibraryScaleFactor %>% strsplit(split="/",fixed=T) %>% sapply("[",1)),
                               Input=Input/as.numeric(RIPLibraryScaleFactor %>% strsplit(split="/",fixed=T) %>% sapply("[",2))) %>%
    mutate(Intensity=(RIP+Pseudocount)/(Input+Pseudocount)) %>%
    group_by(name,strand) %>% mutate(Intensity.max=max(Intensity)) %>% as.data.table() %>% dplyr::arrange(name,strand)
  dt.intensity  <- dt.peak.intensity %>% dplyr::rename(score.old=score,score=Intensity.max) %>% dplyr::filter(Intensity==score) %>%
    group_by(name,strand) %>% slice_head(n=1) %>% as.data.table() %>%
    group_by(name) %>% mutate(strand.new=strand[which.max(score)], score.new=score[which.max(score)]) %>% as.data.table()
  dt.intensity <- dt.intensity %>% mutate(strand = dplyr::case_when(strand == "+" ~ "pos", strand == "-" ~ "neg")) %>%
    pivot_wider(id_cols = c("name","score.old","score.new","strand.new"), names_from = "strand", values_from = c("score","Input","RIP")) %>% as.data.table() %>%
    dplyr::rename(pos=score_pos,neg=score_neg)
  system(command = paste0("rm -rf ", dir.tmp), wait = T)
  return(dt.intensity)
}

input.m6a.sample.file <- args[1]
prefix <- args[2]
RIPLibraryScaleFactor <- args[3]
stranded.bam.dir <- args[4]
bedtools_path <- args[5]
dir.tmp <- args[6]
genome_file <- args[7]
intron.bed.file <- args[8]
intronremoved <- as.logical(args[9])
bedtools_mode <- args[10]
binsize <- as.numeric(args[11])
Pseudocount <- as.numeric(args[12])
SummaryMethod <- args[13]
ComboName <- args[14]
Strandness <- args[15]

options(bedtools.path = bedtools_path)
bedtools_option <- ifelse(bedtools_mode=="counts", paste0("-",Strandness, " -split -counts "), paste0("-",Strandness, " -split -mean "))


input.m6a.sample.file
prefix
RIPLibraryScaleFactor
stranded.bam.dir
bedtools_path
dir.tmp
genome_file
intron.bed.file
intronremoved
bedtools_mode
binsize
Pseudocount
SummaryMethod
ComboName
Strandness
bedtools_option

if(!dir.exists(dir.tmp)){dir.create(dir.tmp,recursive = T)}
setwd(dir = dir.tmp)
t1 <- Sys.time()
if(!file.exists(paste0(stranded.bam.dir, "/", prefix, "__", ComboName, "_bedtools_intensity.tsv"))){
  dt.intensity <- BedtoolsIntensity(input.m6a.sample = input.m6a.sample.file,
                                    prefix=prefix,
                                    RIPLibraryScaleFactor=RIPLibraryScaleFactor,
                                    stranded.bam.dir=stranded.bam.dir,#stranded bam same as mRNA
                                    bedtools_path = bedtools_path,
                                    dir.tmp=dir.tmp,
                                    genome_file=genome_file,#org specific parameter
                                    intron.bed=intron.bed.file,#org specific parameter
                                    intronremoved = intronremoved,
                                    bedtools_option= bedtools_option,
                                    binsize=binsize,
                                    Pseudocount=Pseudocount,
                                    SummaryMethod=SummaryMethod)
  fwrite(dt.intensity, file=paste0(stranded.bam.dir, "/", prefix, "__", ComboName, "_bedtools_intensity.tsv"),sep="\t",col.names = T,row.names=F)
}

Sys.time() - t1
rm(list=ls())
q(save="no")
