

#function to load m6As from each method
LoadPeakMethod <- function(Method=c("MACS2","MeRIPtools", "MeTPeak", "TRESS", "exomePeak2", "exomePeak")[5],
                           Method.peak.dir=paste0("/data/m6A_calling_strategy/RIPPeakS/",c("MACS2","MeRIPtools", "MeTPeak", "TRESS", "exomePeak2", "exomePeak")[5],"_Abcam"))
{
  if(Method=="MACS2"){
    message(paste0("start load  m6a of ", Method))
    Samples <- list.files(path=Method.peak.dir, pattern="peaks.xls",full.names = F,recursive = T) %>%
      str_replace(pattern="_peaks.xls",replacement = "")
    m6As <- list.files(path=Method.peak.dir, pattern="peaks.xls",full.names = T,recursive = T)
    names(m6As) <- Samples
    foreach(i = 1:length(m6As),.combine='rbind')%do%{
      fread(file = m6As[i]) %>% dplyr::select(seqnames=chr,start,end,name,score=fold_enrichment) %>% mutate(strand="*") %>%
        mutate(name=paste(names(m6As)[i],"MACS2", seqnames, start, end, strand,sep="_")) %>%
        dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(Sample=names(m6As)[i])
    } %>%  mutate(Method=Method)
  }else if(Method=="exomePeak"){
    message(paste0("start load  m6a of ", Method))
    Samples <- list.files(path=Method.peak.dir, pattern="peak.xls",full.names = T,recursive = T) %>% grep(pattern="con_peak.xls",invert = T,value=T) %>%
      strsplit(split="/") %>% sapply(tail,2,simplify = F) %>% sapply(head,1)
    m6As <- list.files(path=Method.peak.dir, pattern="peak.xls",full.names = T,recursive = T) %>% grep(pattern="con_peak.xls",invert = T,value=T)
    names(m6As) <- Samples
    foreach(i = 1:length(m6As),.combine='rbind')%do%{
      fread(file = m6As[i],header = T) %>% dplyr::mutate(seqnames=chr, start=chromStart, end=chromEnd, score=fold_enrchment) %>%
        mutate(name=paste(names(m6As)[i],"exomePeak", seqnames, start, end, strand,sep="_")) %>%
        dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(Sample=names(m6As)[i])
    } %>% mutate(Method=Method)
  }else if(Method=="exomePeak2"){
    message(paste0("start load  m6a of ", Method))
    Samples <- list.files(path=Method.peak.dir, pattern="peaks.csv",full.names = T,recursive = T) %>%
      strsplit(split="/") %>% sapply(tail,2,simplify = F) %>% sapply(head,1) %>% strsplit(split="__",fixed=T) %>% sapply("[",1)
    m6As <- list.files(path=Method.peak.dir, pattern="peaks.csv",full.names = T,recursive = T)
    names(m6As) <- Samples
    foreach(i = 1:length(m6As),.combine='rbind')%do%{
      fread(file = m6As[i],header = T) %>% dplyr::mutate(seqnames=chr, start=chromStart, end=chromEnd, score=(RPM.IP+0.5)/(RPM.input+0.5)) %>%#add pseudocount 0.1
        dplyr::filter(RPM.IP>RPM.input) %>%
        mutate(name=paste(names(m6As)[i],"exomePeak2", seqnames, start, end, strand,sep="_")) %>%
        dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(Sample=names(m6As)[i])
    } %>% mutate(Method=Method)
  }else if(Method=="MeTPeak"){
    message(paste0("start load  m6a of ", Method))
    Samples <- list.files(path=Method.peak.dir, pattern="peak.xls",full.names = T,recursive = T) %>%
      strsplit(split="/") %>% sapply(tail,2,simplify = F) %>% sapply(head,1)
    m6As <- list.files(path=Method.peak.dir, pattern="peak.xls",full.names = T,recursive = T)
    names(m6As) <- Samples
    foreach(i = 1:length(m6As),.combine='rbind')%do%{
      fread(file = m6As[i],header = T) %>% dplyr::mutate(seqnames=chr, start=chromStart, end=chromEnd, score=fold_enrchment) %>%
        mutate(name=paste(names(m6As)[i],"MeTPeak", seqnames, start, end, strand,sep="_")) %>%
        dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(Sample=names(m6As)[i])
    } %>% mutate(Method=Method)
  }else if(Method=="TRESS"){
    message(paste0("start load  m6a of ", Method))
    Samples <- list.files(path=Method.peak.dir, pattern="peaks.xls",full.names = T,recursive = T) %>%
      strsplit(split="/") %>% sapply(tail,1) %>% str_replace(pattern="_peaks.xls",replacement="")
    m6As <- list.files(path=Method.peak.dir, pattern="peaks.xls",full.names = T,recursive = T)
    names(m6As) <- Samples
    foreach(i = 1:length(m6As),.combine='rbind')%do%{
      fread(file = m6As[i],header = T) %>% dplyr::mutate(seqnames=chr, score=2^(as.numeric(lg.fc))) %>%
        mutate(name=paste(names(m6As)[i],"TRESS", seqnames, start, end, strand,sep="_")) %>%
        dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(Sample=names(m6As)[i])
    } %>% mutate(Method=Method)
  }else if(Method=="MeRIPtools"){
    message(paste0("start load  m6a of ", Method))
    Samples <- list.files(path=Method.peak.dir, pattern="peak.tsv",full.names = T,recursive = T) %>%
      strsplit(split="/") %>% sapply(tail,1) %>% str_replace(pattern="_peak.tsv",replacement="")
    m6As <- list.files(path=Method.peak.dir, pattern="peak.tsv",full.names = T,recursive = T)
    names(m6As) <- Samples
    foreach(i = 1:length(m6As),.combine='rbind')%do%{
      fread(file = m6As[i],header = T) %>% dplyr::mutate(seqnames=chr) %>%
        mutate(name=paste(names(m6As)[i],"MeRIPtools", seqnames, start, end, strand,sep="_")) %>%
        dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(Sample=names(m6As)[i])
    } %>% mutate(Method=Method)
  }
}
#function to calculate peak intensity parallel
PeakIntensity.parallel <- function(n.cores=2,
                                   peakbeds=list.files(path="~/Fetus_m6A/MACS2SMePeak/MACS2SMePeak",pattern = '^(FP).+(bed)$',full.names = T),
                                   names.peakbeds=strsplit(list.files(path="~/Fetus_m6A/MACS2SMePeak/MACS2SMePeak",pattern = '^(FP).+(bed)$',full.names = T), split="/", fixed=T) %>%
                                     sapply(tail,1) %>% strsplit(split="_MACS2SMePeak",fixed=T) %>% sapply("[",1),
                                   InputBAM=list.files(path="/data/m6A_calling_strategy/7_MACS2SPeak/BAM", pattern = '^(HEK).+(Input).+(UN).+(bam)$',full.names = T),
                                   names.Input=list.files(path="/data/m6A_calling_strategy/7_MACS2SPeak/BAM", pattern = '^(HEK).+(Input).+(UN).+(bam)$',full.names = F) %>%
                                     strsplit(split=".",fixed=T) %>% sapply("[",1) %>% gsub(pattern="_Input",replacement=""),
                                   RIPBAM=list.files(path="/data/m6A_calling_strategy/7_MACS2SPeak/BAM", pattern = '^(HEK).+(RIP).+(UN).+(bam)$',full.names = T),
                                   names.RIP=list.files(path="/data/m6A_calling_strategy/7_MACS2SPeak/BAM", pattern = '^(HEK).+(RIP).+(UN).+(bam)$',full.names = F) %>%
                                     strsplit(split=".",fixed=T) %>% sapply("[",1) %>% gsub(pattern="_RIP",replacement=""),
                                   Depth=dt.BAM.Depth$Mapped,
                                   names.Depth=dt.BAM.Depth$Sample,
                                   out.tmp="/data/m6A_calling_strategy/SlidingWindow/tmp",
                                   bedtools_path = "~/anaconda3/envs/m6A_seq/bin",
                                   bt_coverage.options=" -split -mean ",
                                   genome_file="~/genome_db/STAR_index/STAR_hg38/chrNameLength.txt",
                                   strandness = c("s", "S", "unstranded"),
                                   pseudo_count=0.1,
                                   foldchange_cutoff=1,
                                   RIP_coverage_cutoff=2){
  #write commandline
  require(foreach)
  require(data.table)
  require(dplyr)
  if(!dir.exists(out.tmp)){dir.create(out.tmp,recursive = T)}
  #obtain the intersect sample name
  Samples <- Reduce(intersect, list(names.peakbeds, names.Input, names.RIP,
                                    gsub(names.Depth,pattern="_Input",replacement=""), gsub(names.Depth,pattern="_RIP",replacement="")))
  n.cores <- min(c(n.cores, length(Samples)*2))
  peakbeds <- peakbeds[match(Samples,names.peakbeds)]
  InputBAM <- InputBAM[match(Samples,names.Input)]
  RIPBAM <- RIPBAM[match(Samples,names.RIP)]
  Input.Depth <- Depth[match(paste0(Samples,"_Input"), names.Depth)]
  RIP.Depth <- Depth[match(paste0(Samples,"_RIP"), names.Depth)]
  #determine bt_coverage.options
  if(strandness != "unstranded"){ bt_coverage.options <- paste0(bt_coverage.options, "-", strandness, " ")}
  #sort input bed files

  if(RIP_coverage_cutoff>0){
    #bash command to perform bedtools coverage in RIP samples first (keep bins with at lest RIP_coverage_cutoff=2)
    #build parameter table for GNU parallel
    dt.parameter.RIP <- foreach(i = 1:length(Samples), .combine='rbind')%do%{
      data.table(a=peakbeds[i], b=RIPBAM[i], o=paste0(out.tmp,"/",Samples[i], "_RIP_coverage.txt"))
    }
    fwrite(dt.parameter.RIP, file=paste0(out.tmp,"/RIP_parallel_bedtools_parameter_table.txt"),sep="\t",row.names = F,col.names=F)
    parallel_cmd.RIP <- paste0("parallel -j ", n.cores," --will-cite -a ", out.tmp, "/RIP_parallel_bedtools_parameter_table.txt --colsep '\t' '",
                               bedtools_path,"/bedtools coverage  -g ", genome_file, bt_coverage.options, " -a {1} -b {2} 1>{3} 2>/dev/null '")
    system(command = parallel_cmd.RIP, wait = T)
    #determine the accomplishment of the bedtools coverage
    no.RIP <- system(command = paste0("ls  ",out.tmp, "/*_RIP_coverage.txt | wc -l"),wait = T,intern = T) %>% as.integer()
    if(no.RIP == length(Samples)){
      message(paste0("Finished RIP coverage calculation!\nfilter bins according to RIP coverage...\n",  Sys.time()))
      for(k in 1:length(Samples)){
        RIP_cov <- fread(paste0(out.tmp,"/", Samples[k],"_RIP_coverage.txt"),sep="\t",header=F,nThread = floor(72/n.cores))
        colnames(RIP_cov) <- c("chr","start","end","name","score","strand","Cov")
        filtered_bins <- RIP_cov %>% dplyr::filter(Cov>=RIP_coverage_cutoff)
        if(nrow(filtered_bins)>0){
          fwrite(filtered_bins[,1:6],file=paste0(out.tmp, "/", Samples[k], "_filtered_bins.bed"), row.names = F, col.names=F, sep="\t",nThread = floor(72/n.cores))}else{
            Samples <- Samples[Samples != Samples[k]]
          }
      }
    }
    rm(k)
    #estimate intensity for filtered bins
    #build parameter table
    dt.parameter.InputRIP <- foreach(k1=1:length(Samples),.combine ='rbind')%do%{
      if(file.exists(paste0(out.tmp, "/", Samples[k1], "_filtered_bins.bed"))){
        data.table(a=rep(paste0(out.tmp, "/", Samples[k1], "_filtered_bins.bed"),2),
                   b=c(InputBAM[k1],RIPBAM[k1]), c=paste0(out.tmp, "/", Samples[k1], "_", c("Input","RIP"),"_coverage.txt"))
      }
    }
  }else{
    #estimate intensity for filtered bins
    #build parameter table
    dt.parameter.InputRIP <- foreach(k1=1:length(Samples),.combine ='rbind')%do%{
      data.table(a=rep(peakbeds[k1],2), b=c(InputBAM[k1],RIPBAM[k1]), c=paste0(out.tmp, "/", Samples[k1], "_", c("Input","RIP"),"_coverage.txt"))
    }
  }
  fwrite(dt.parameter.InputRIP, file=paste0(out.tmp,"/InputRIP_parallel_bedtools_parameter_table.txt"),sep="\t",row.names = F,col.names=F)
  parallel_cmd.InputRIP <- paste0("parallel -j ", n.cores," --will-cite -a ", out.tmp, "/InputRIP_parallel_bedtools_parameter_table.txt --colsep '\t' '",
                                  bedtools_path,"/bedtools coverage -g ", genome_file, bt_coverage.options, " -a {1} -b {2} 1>{3} 2>/dev/null '")
  system(command = parallel_cmd.InputRIP, wait = T)

  #determine the accomplishment of the bedtools coverage
  no.Input <- system(command = paste0("ls  ",out.tmp, "/*_Input_coverage.txt | wc -l"),wait = T,intern = T) %>% as.integer()
  no.RIP <- system(command = paste0("ls  ",out.tmp, "/*_RIP_coverage.txt | wc -l"),wait = T,intern = T) %>% as.integer()
  if(min(no.Input,no.RIP) == length(Samples)){
    message(paste0("Finished coverage calculation!\nStart peak intensity calculation...\n",  Sys.time()))
    peak_level_all <- foreach(j = 1:length(Samples),.combine='rbind')%do%{
      Input_cov <- fread(paste0(out.tmp,"/", Samples[j],"_Input_coverage.txt"),sep="\t",header=F,nThread = floor(72/n.cores)) %>% mutate(Cov.norm=V7/Input.Depth[j])
      RIP_cov <- fread(paste0(out.tmp,"/", Samples[j],"_RIP_coverage.txt"),sep="\t",header=F,nThread = floor(72/n.cores)) %>% mutate(Cov.norm=V7/RIP.Depth[j])
      colnames(Input_cov) <- colnames(RIP_cov) <- c("chr","start","end","name","score","strand","Cov","Cov.norm")
      peak_level <- inner_join(x=Input_cov,y=RIP_cov,by=c("chr","start","end","name","score","strand"),suffix=c(".Input",".RIP")) %>%
        as.data.table() %>% mutate(intensity=(Cov.norm.RIP+pseudo_count)/(Cov.norm.Input+pseudo_count))
      message(paste0("Finished the ",j, "/", length(Samples)))
      peak_level %>% mutate(Sample = Samples[j]) %>% dplyr::filter(intensity > foldchange_cutoff)
    }
    message(paste0("Finished intensity calculation!\n",  Sys.time()))
  }
  return(peak_level_all)
}
#function to remove intron for long peak
#function to remove intron for long peak
LongPeakRemoveIntron <- function(dt.peak = dt.MeTPeak.m6As %>% dplyr::filter(ComboName=="Combo2" & Sample=="HEK_NEB_mRNA1" & end-start>=5000),
                                 genome_file="~/genome_db/STAR_index/STAR_hg38/chrNameLength.txt",
                                 dt.intron=readRDS("/data/m6A_calling_strategy/Analysis/dt.intron.hg38.RDS")
){
  #sort input and intron
  dt.input <- bedtoolsr::bt.sort(dt.peak,g = genome_file) %>% as.data.table()
  dt.intron <- bedtoolsr::bt.sort(dt.intron,g = genome_file) %>% as.data.table()
  dt.peak.rm.intron <- bedtoolsr::bt.subtract(a=dt.input,b=dt.intron,s=T) %>% as.data.table()
  dt.peak.rm.intron <- dt.peak.rm.intron %>% group_by(V4) %>% mutate(SumWidth=sum(V3-V2)) %>% as.data.table()
  attr(dt.peak.rm.intron,"groups") <- NULL
  colnames(dt.peak.rm.intron) <- c(colnames(dt.peak),"SumWidth")
  return(dt.peak.rm.intron)
}

#main function to run M6APeakS
runM6APeakS  <- function(
    InputBAMs=list.files(path="/data/m6A_calling_strategy/7_MACS2SPeak/BAM/",pattern = "HEK",full.names = T) %>%
      grep(pattern="Input",value=T) %>% grep(pattern="Abcam", value=T, invert = F) %>% grep(pattern="test",value=T,invert = F) %>%
      grep(pattern="bai",value=T,invert=T),
    RIPBAMs=list.files(path="/data/m6A_calling_strategy/7_MACS2SPeak/BAM/",pattern = "HEK",full.names = T) %>%
      grep(pattern="RIP",value=T) %>% grep(pattern="Abcam", value=T, invert = F) %>% grep(pattern="test",value=T,invert = F) %>%
      grep(pattern="bai",value=T,invert=T),
    # tmp.dir="/data/m6A_calling_strategy/M6APeakS_Abcam/tmp/",
    bin.dir="/home/huangp/anaconda3/envs/m6A_seq/bin/",#absolute path
    # log.dir="/data/m6A_calling_strategy/M6APeakS_Abcam/log",
    out.dir="/data/m6A_calling_strategy/M6APeakS_Abcam/",
    script.dir="/data/m6A_calling_strategy/Script/",
    Organism="Human",
    org.gtf_file="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",
    org.sqlite_file="/data/m6A_calling_strategy/RIPPeakS_resource/UCSC.hg38.knownGene.sqlite",
    org.genome_file="/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt",
    org.intron.bed="/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.hg38.bed",
    org.genebed="/data/m6A_calling_strategy/RIPPeakS_resource/hg38_gencode.gene.bed",
    ModelCutoffTable="/data/m6A_calling_strategy/RIPPeakS_resource/dt.cutoff.at.all.FPR.M6APeakS.tsv",
    selectedModel="NEB",
    Model.FPR.cutoff="Moderate",#could be NULL (default, which is Moderate mode), character values refer to recommanded three mode:Strict, Sensitive, Moderate, or numeric value between 0.02-1, except Millipore, 0.04-1
    org.benchmark.bed="/data/m6A_calling_strategy/RIPPeakS_resource/dt.benchmark.GLORI_HEK293.m6As.bed",
    n.cores=40,
    rm.inter=TRUE
){
  options(scipen = 9)
  options(bedtools.path = bin.dir)
  if(!dir.exists(out.dir)){dir.create(out.dir,recursive = T)}
  tmp.dir <- paste0(out.dir,"/tmp")
  log.dir <- paste0(out.dir,"/log")
  if(!dir.exists(tmp.dir)){dir.create(tmp.dir,recursive = T)}
  if(!dir.exists(log.dir)){dir.create(log.dir,recursive = T)}
  t1 <- Sys.time()
  #step 0. obtain the samples, infer and calculate depth
  message(paste0("[",Sys.time(),"] ","step 0.1 keep paired samples"))
  #0.1 paired samples with input and rip
  Sample.Input <- InputBAMs[file.exists(InputBAMs)] %>% strsplit(split="/",fixed=T) %>% sapply(tail,1) %>% strsplit(split="_Input") %>% sapply("[",1)
  names(InputBAMs) <- Sample.Input
  Sample.RIP <- RIPBAMs[file.exists(RIPBAMs)] %>% strsplit(split="/",fixed=T) %>% sapply(tail,1) %>% strsplit(split="_RIP") %>% sapply("[",1)
  names(RIPBAMs) <- Sample.RIP
  Samples <- intersect(Sample.Input, Sample.RIP)
  #obtain samples finished peak calling
  peakfile.pattern <- paste0("M6APeakS_peak_",selectedModel,"_", ifelse(is.null(Model.FPR.cutoff),"default",Model.FPR.cutoff))
  finished.samples <- list.files(out.dir, pattern=peakfile.pattern,full.names = F) %>% strsplit(split="_M6APeakS_peak") %>% sapply(head,1)
  Samples <- Samples[!Samples %in% finished.samples]
  if(length(Samples)>0){
    #keep new input rip bams
    InputBAMs <- InputBAMs[Samples]
    RIPBAMs <- RIPBAMs[Samples]
    #0.2 infer strandness for input and rip
    message(paste0("[",Sys.time(),"] ","step 0.2 infer strandness"))
    #build parameter table for GNU parallel
    dt.parameter.infer.strandness.parallel <- foreach(i = 1:length(Samples), .combine='rbind')%do%{
      data.table(arg1=rep(2000000,2),#nReads to use
                 arg2=rep(org.genebed,2),#gene bed to use
                 arg3=c(InputBAMs[i],RIPBAMs[i]),
                 arg4=paste0(tmp.dir,"/",Samples[i],c("_Input","_RIP"),"_infer.strandness.tsv")#output
      )
    }
    fwrite(dt.parameter.infer.strandness.parallel, file=paste0(tmp.dir,"/parameter.infer.strandness.parallel.txt"),sep="\t",row.names = F,col.names=F)
    parallel.infer.strandness.cmd <- paste0("parallel -j ", floor(n.cores/4)," --will-cite -a ", paste0(tmp.dir,"/parameter.infer.strandness.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"infer_experiment.py"),
                                            " -s {1} -r {2} -i {3} 1>{4} 2>&1 '")
    system(command = parallel.infer.strandness.cmd, wait = T)
    #0.3 bam index as input for exomePeak and MeTPeak (run bam index to wait forinfer strandness finished )
    message(paste0("[",Sys.time(),"] ","step 0.3 bam index"))
    #build parameter table for parallel
    dt.parameter.bai.parallel <- foreach(i=1:length(Samples),.combine='rbind')%do%{
      data.table(arg1=c(InputBAMs[i],RIPBAMs[i]))
    }
    fwrite(dt.parameter.bai.parallel, file=paste0(tmp.dir,"/parameter.bam.index.parallel.txt"),sep="\t",row.names = F,col.names=F)
    parallel.bai.cmd <- paste0("parallel -j ", floor(n.cores/4)," --will-cite -a ", paste0(tmp.dir,"/parameter.bam.index.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"samtools"),
                               " index {1} 2>&1 '")
    system(command = parallel.bai.cmd, wait = T)
    #0.4 pull out stranded reads which is reverse to the mRNA
    #0.4 pull out stranded reads which is reverse to the mRNA
    message(paste0("[",Sys.time(),"] ","step 0.4 pull out stranded reads which is reverse to the mRNA"))
    dt.infer.strandness.res <- foreach(i=1:length(Samples),.combine='rbind')%do%{
      foreach(j=c("Input","RIP"),.combine='rbind')%do%{
        dt.res <- fread(paste0(tmp.dir,"/",Samples[i],"_", j, "_infer.strandness.tsv"))
        if(sum(str_detect(dt.res$V6,pattern=fixed("2++,2--")))>0){
          Library <- "PE"
          if(dt.res[str_detect(V6,fixed("2++,2--")),V7]>=0.5){Strandness <- "R2"}else{
            if(dt.res[str_detect(V6,fixed("1++,1--")),V7]>=0.5){Strandness <- "R1"}else{
              Strandness <- "unstranded"
            }
          }
          data.table(Sample=Samples[i], Condition=j, Library=Library, Strandness=Strandness)
        }else{
          Strandness <- "unstranded"
          message(paste0("Error! BAM files is not PE"))
        }
      }
    }

    if(length(unique(dt.infer.strandness.res$Strandness))==1 & unique(dt.infer.strandness.res$Strandness)!="unstranded"){
      if(unique(dt.infer.strandness.res$Strandness == "R1")){
        dt.parameter.stranded.bam.parallel <- foreach(i=1:length(Samples),.combine='rbind')%do%{
          data.table(arg1=rep(4,2),#thread for samtools view
                     arg2=rep("-F",2),
                     arg3=rep("64", 2),
                     arg4=c(InputBAMs[i], RIPBAMs[i]),
                     arg5=paste0(tmp.dir,"/",Samples[i], c("_Input","_RIP"), ".R2.bam"),
                     arg6=paste0(log.dir,"/",Samples[i], c("_Input","_RIP"), ".R2.log")
          )
        }
      }else if(unique(dt.infer.strandness.res$Strandness == "R2")){
        dt.parameter.stranded.bam.parallel <- foreach(i=1:length(Samples),.combine='rbind')%do%{
          data.table(arg1=rep(4,2),#thread for samtools view
                     arg2=rep("-f",2),
                     arg3=rep("64", 2),
                     arg4=c(InputBAMs[i], RIPBAMs[i]),
                     arg5=paste0(tmp.dir,"/",Samples[i], c("_Input","_RIP"), ".R1.bam"),
                     arg6=paste0(log.dir,"/",Samples[i], c("_Input","_RIP"), ".R1.log")
          )
        }
      }else{
        message(paste0("Error! failed to pull out stranded reads from PE BAMs"))
      }
      fwrite(dt.parameter.stranded.bam.parallel, file=paste0(tmp.dir,"/parameter.stranded.bam.parallel.txt"),sep="\t",row.names = F,col.names=F)
      parallel.stranded.bam.cmd <- paste0("parallel -j ", floor(n.cores/5)," --will-cite -a ", paste0(tmp.dir,"/parameter.stranded.bam.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"samtools"),
                                          " view --threads {1} -b -h {2} {3} {4} -o {5} >{6} 2>&1 '")
      system(command = parallel.stranded.bam.cmd, wait = T)

      #calculate stranded bam depth
      registerDoParallel(cl=min(floor(n.cores/4), length(Samples)))
      dt.stranded.BAM.Depth <- foreach(s=Samples,.combine='rbind')%dopar%{
        sample.input.s.bam <- list.files(tmp.dir,pattern="(bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_Input",value=T)
        sample.ip.s.bam <- list.files(tmp.dir,pattern="(bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_RIP",value=T)
        Depth <- c(as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.input.s.bam), wait = T, intern = T)),
                   as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.ip.s.bam), wait = T, intern = T)))
        data.table(Sample=paste0(s,c("_Input","_RIP")), Depth=Depth) %>% mutate(Replicate=s)
      }
      stopImplicitCluster()
      dt.stranded.BAM.Depth  <- dt.stranded.BAM.Depth %>% mutate(Mapped=Depth/1000000)
      dt.stranded.BAM.Depth <- dt.stranded.BAM.Depth %>% mutate(Library=ifelse(str_detect(Sample,"_Input"),"Input","RIP")) %>%
        tidyr::pivot_wider(id_cols = c("Replicate"),names_from = c("Library"),values_from = c("Mapped")) %>% as.data.table() %>%
        mutate(RIPLibraryScaleFactor=paste0(round(RIP,2),"/",round(Input,2)))
      print(dt.stranded.BAM.Depth)

      #step 1. Run each method to call peaks
      #according to selected Model and FPR cutoff, determine the selected method combo and option
      if(!file.exists(ModelCutoffTable)){message(paste0("Error! Can't find model cutoff table file at ", ModelCutoffTable))}
      dt.method.combo.model <- fread(ModelCutoffTable,header = T) %>% dplyr::filter(M6APeakModel==selectedModel)
      if(is.null(Model.FPR.cutoff)){
        dt.method.combo.selected <- dt.method.combo.model %>% dplyr::filter(Mode=="Moderate")}else{
          if(is.numeric(Model.FPR.cutoff)){
            dt.method.combo.selected <- dt.method.combo.model %>% dplyr::filter(FPR==Model.FPR.cutoff)
          }else{
            dt.method.combo.selected <- dt.method.combo.model %>% dplyr::filter(Mode==Model.FPR.cutoff)
          }
        }
      if(nrow(dt.method.combo.selected)<=0){message(paste0("Error! Please choose appropriate model and FPR cutoff! "))}
      print(dt.method.combo.selected)
      #1.1 run non-exomePeak2 method combo
      dt.method.combo.selected.nonexomepeak2 <- dt.method.combo.selected %>% dplyr::filter(Method != "exomePeak2")
      if(nrow(dt.method.combo.selected.nonexomepeak2)>0){
        dt.parameter.nonexomePeak2.parallel <- foreach(M=unique(dt.method.combo.selected.nonexomepeak2$Method),.combine = 'rbind')%do%{
          foreach(combo = unique(dt.method.combo.selected.nonexomepeak2[Method==M,ComboName]), .combine='rbind')%do%{
            data.table(Method=M,
                       bin_path=rep(bin.dir,length(Samples)),
                       InputBAM=InputBAMs,
                       RIPBAM=RIPBAMs,
                       outdir=rep(paste0(out.dir,"/", M),length(Samples)),
                       prefix=Samples,
                       SelectedCombo=rep(combo,length(Samples)),
                       Organism=rep(Organism, length(Samples)),
                       Annot.gtf=rep(org.gtf_file, length(Samples)),
                       Annot.sqlite=rep(org.sqlite_file, length(Samples)),
                       Annot.genome=rep(org.genome_file, length(Samples)),
                       Strandness=rep("S", length(Samples)),
                       nthread=2,#thread for MeRIPtools and exomePeak2 method and ignore for other method
                       Logfile=paste0(log.dir,"/", M, "_",combo,"_", Samples,".log")
            )
          }
        }
        message(paste0("[",Sys.time(),"] ","step 1.0 call peaks use non-exomePeak2 methods: MACS2, MeTPeak, MeRIPtools, exomePeak, exomePeak2, and TRESS"))
        fwrite(dt.parameter.nonexomePeak2.parallel, file=paste0(tmp.dir,"/parameter.Method.except.exomePeak2.parallel.txt"),sep="\t",row.names = F,col.names=F)
        parallel.Method.except.exomePeak2.cmd <- paste0("parallel -j ", n.cores," --will-cite -a ", paste0(tmp.dir,"/parameter.Method.except.exomePeak2.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"/Rscript "),
                                                        paste0(script.dir,"/Rscript_"), "{1}.R {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} 1>{14}  2>&1 '")
        system(command = parallel.Method.except.exomePeak2.cmd,wait = T)
      }

      dt.method.combo.selected.exomepeak2 <- dt.method.combo.selected %>% dplyr::filter(Method == "exomePeak2")
      if(nrow(dt.method.combo.selected.exomepeak2)>0){
        dt.parameter.exomePeak2.parallel <- foreach(M=unique(dt.method.combo.selected.exomepeak2$Method),.combine = 'rbind')%do%{
          foreach(combo = unique(dt.method.combo.selected.exomepeak2[Method==M,ComboName]), .combine='rbind')%do%{
            data.table(Method=M,
                       bin_path=rep(bin.dir,length(Samples)),
                       InputBAM=InputBAMs,
                       RIPBAM=RIPBAMs,
                       outdir=rep(paste0(out.dir,"/", M),length(Samples)),
                       prefix=Samples,
                       SelectedCombo=rep(combo,length(Samples)),
                       Organism=rep(Organism, length(Samples)),
                       Annot.gtf=rep(org.gtf_file, length(Samples)),
                       Annot.sqlite=rep(org.sqlite_file, length(Samples)),
                       Annot.genome=rep(org.genome_file, length(Samples)),
                       Strandness=rep("S", length(Samples)),
                       nthread=2,#thread for MeRIPtools and exomePeak2 method and ignore for other method
                       Logfile=paste0(log.dir,"/", M, "_",combo,"_", Samples,".log")
            )

          }
        }
        fwrite(dt.parameter.exomePeak2.parallel, file=paste0(tmp.dir,"/parameter.Method.exomePeak2.parallel.txt"),sep="\t",row.names = F,col.names=F)
        parallel.Method.exomePeak2.cmd <- paste0("parallel -j ", min(max(floor(n.cores/(quantile(c(dt.stranded.BAM.Depth$Input,dt.stranded.BAM.Depth$RIP),0.75)/30)/8),1),10)," --will-cite -a ", paste0(tmp.dir,"/parameter.Method.exomePeak2.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"/Rscript "),
                                                 paste0(script.dir,"/Rscript_"), "{1}.R {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} 1>{14}  2>&1 '")
        message(paste0("step 1.1 call peaks use exomePeak2 method"))
        system(command = parallel.Method.exomePeak2.cmd,wait = T)
      }

      #load and save m6As
      dt.Method.m6a <- foreach(M = unique(dt.method.combo.selected$Method),.combine='rbind')%do%{
        foreach(combo = unique(dt.method.combo.selected[Method==M,ComboName]), .combine='rbind')%do%{
          LoadPeakMethod(Method=M, Method.peak.dir = paste0(out.dir,"/", M, "/",M,"_", combo)) %>% mutate(Method=M, ComboName=combo) %>%
            dplyr::filter(Sample %in% Samples)
        }
      }
      dt.Method.m6a %>% distinct(Method,ComboName,Sample)
      #step 2. calculate bedtools-based intensity for each method's peaks
      message(paste0("[",Sys.time(),"] ","step 2 calculate bedtools-based intensity for peaks from selected method combo"))
      #bedtools intensity for exomePeak2|Combo35_depth_Ex_bin50_P2
      dt.parameter.bedtools.parallel <- foreach(M=unique(dt.Method.m6a$Method), .combine='rbind')%do%{
        foreach(combo=unique(dt.Method.m6a[Method==M, ComboName]), .combine='rbind')%do%{
          foreach(s=Samples, .combine = 'rbind')%do%{
            fwrite(dt.Method.m6a %>% dplyr::filter(ComboName==combo & Method==M & Sample==s), file=paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),row.names=F,col.names=T,sep="\t")
            selected.option <- dt.method.combo.selected[Method==M & ComboName==combo,optionID]
            data.table(input.m6a = paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),
                       prefix=s,
                       RIPLibraryScaleFactor=dt.stranded.BAM.Depth[match(s,Replicate),RIPLibraryScaleFactor],
                       stranded.bam.dir = tmp.dir,
                       bedtools_path = bin.dir,
                       dir.tmp = paste0(tmp.dir,"/", M, "_", combo, "_bedtools_intensity_tmp_",s),
                       genome_file = org.genome_file,
                       intron.bed = org.intron.bed,
                       intronremoved = str_detect(selected.option, "Ex"),
                       bedtools_mode = ifelse(str_detect(selected.option, "mean"),"mean","counts"),
                       binsize = as.integer(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",3) %>% str_replace(pattern="bin",replacement = "")),
                       Pseudocount = as.numeric(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",4) %>% str_replace(pattern="P",replacement = "")),
                       SummaryMethod = "max",
                       ComboName = paste0(M, "__", combo, "_",selected.option),
                       Strandness = "S",
                       Logfile=paste0(log.dir,"/",M, "_", combo,"_bedtools_intensity_",s,".log")
            )
          }
        }
      }
      fwrite(dt.parameter.bedtools.parallel, file=paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt"),sep="\t",row.names = F,col.names=F)
      intensity.nthread <- max(floor(n.cores/(quantile(c(dt.stranded.BAM.Depth$Input,dt.stranded.BAM.Depth$RIP),0.75)/30)/8),1)
      parallel.Method.bedtools.cmd <- paste0("parallel -j ", min(intensity.nthread,10)," --will-cite -a ", paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"/Rscript "),
                                             paste0(script.dir,"/Rscript_bedtools_intensity.R")," {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} 1>{16}  2>&1 '")
      system(command = parallel.Method.bedtools.cmd, wait = T)

      #step 3. filter and merged at cutoff
      #obtain m6a at FPR cutoff
      dt.m6a.FPR <- foreach(M=unique(dt.method.combo.selected$Method),.combine='rbind')%do%{
        foreach(combo = dt.method.combo.selected[Method==M,unique(ComboName)], .combine='rbind')%do%{
          message(paste0("Obtain peak at FPR cutoff at ",ifelse(is.null(Model.FPR.cutoff), "default", Model.FPR.cutoff)," for ", M, " ", combo, " ", dt.method.combo.selected[Method==M & ComboName==combo,optionID]))
          intensity.files <- list.files(tmp.dir,pattern = paste0(M,"__",combo),full.names = T) %>% grep(pattern="(_bedtools_intensity.tsv)$",value=T)
          names(intensity.files) <- list.files(tmp.dir,pattern = paste0(M,"__",combo),full.names = F) %>% grep(pattern="(_bedtools_intensity.tsv)$",value=T) %>%
            strsplit(split="__",fixed=T) %>% sapply("[",1)
          intensity.files <- intensity.files[Samples]
          dt.intensity <- foreach(s=names(intensity.files),.combine='rbind')%do%{fread(file = intensity.files[s]) %>% mutate(Sample=s)}#load intensity for all sample
          dt.intensity <- dt.intensity %>% dplyr::select(name,Sample,pos,neg)
          #pull out m6a with intensity higher than cutoff at FPR20
          #filtered and un overlapped  m6a at pos strand
          dt.m6a.filtered.pos <- dt.intensity %>% dplyr::filter(pos >= dt.method.combo.selected[Method==M & ComboName==combo,cutoff])
          if(nrow(dt.m6a.filtered.pos)>0){
            dt.m6a.filtered.pos <- dt.m6a.filtered.pos %>% dplyr::mutate(name = name %>% strsplit(split=paste0(M,"_"),fixed=T) %>% sapply(tail,1)) %>%
              dplyr::mutate(seqnames=name %>% strsplit(split="_",fixed=T) %>% sapply("[",1),
                            start=name %>% strsplit(split="_",fixed=T) %>% sapply("[",2),
                            end=name %>% strsplit(split="_",fixed=T) %>% sapply("[",3),
                            strand="+") %>%
              dplyr::select(seqnames,start,end,name,score=pos,strand,Sample)
            dt.m6a.filtered.pos.gr <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.filtered.pos,keep.extra.columns = T)
            dt.m6a.filtered.pos.merged <- foreach(sam = unique(dt.m6a.filtered.pos$Sample), .combine='rbind')%do%{
              GenomicRanges::reduce(dt.m6a.filtered.pos.gr[dt.m6a.filtered.pos.gr$Sample==sam,]) %>% as.data.table() %>%
                mutate(score=1000,name=paste(seqnames,start,end,strand,sep="_"),Sample=sam) %>%
                dplyr::select(seqnames,start,end,name,score,strand,Sample)
            }
          }else{
            dt.m6a.filtered.pos.merged <- data.table(seqnames=NULL,start=NULL,end=NULL,name=NULL,score=NULL,strand=NULL,Sample=NULL)
          }
          #filtered  and un overlapped m6a at neg strand
          dt.m6a.filtered.neg <- dt.intensity %>% dplyr::filter(neg >= dt.method.combo.selected[Method==M & ComboName==combo,cutoff])
          if(nrow(dt.m6a.filtered.neg)>0){
            dt.m6a.filtered.neg <- dt.m6a.filtered.neg %>% dplyr::mutate(name = name %>% strsplit(split=paste0(M,"_"),fixed=T) %>% sapply(tail,1)) %>%
              dplyr::mutate(seqnames=name %>% strsplit(split="_",fixed=T) %>% sapply("[",1),
                            start=name %>% strsplit(split="_",fixed=T) %>% sapply("[",2),
                            end=name %>% strsplit(split="_",fixed=T) %>% sapply("[",3),
                            strand="-") %>%
              dplyr::select(seqnames,start,end,name,score=neg,strand,Sample)
            dt.m6a.filtered.neg.gr <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.filtered.neg,keep.extra.columns = T)
            dt.m6a.filtered.neg.merged <- foreach(sam = unique(dt.m6a.filtered.neg$Sample), .combine='rbind')%do%{
              GenomicRanges::reduce(dt.m6a.filtered.neg.gr[dt.m6a.filtered.neg.gr$Sample==sam,]) %>% as.data.table() %>%
                mutate(score=1000,name=paste(seqnames,start,end,strand,sep="_"),Sample=sam) %>%
                dplyr::select(seqnames,start,end,name,score,strand,Sample)
            }
          }else{
            dt.m6a.filtered.neg.merged <- data.table(seqnames=NULL,start=NULL,end=NULL,name=NULL,score=NULL,strand=NULL,Sample=NULL)
          }
          dt.m6a.filtered <- rbind(dt.m6a.filtered.pos.merged, dt.m6a.filtered.neg.merged)
          if(nrow(dt.m6a.filtered)>0){
            #calculate exonic sum width to filter extreme long peak (<=6kb)
            dt.m6a.filtered.sumwidth <- LongPeakRemoveIntron(dt.peak = dt.m6a.filtered %>% distinct(seqnames,start,end,name,score,strand),#unique m6a peak
                                                             dt.intron = fread(org.intron.bed),
                                                             genome_file = org.genome_file)
            dt.m6a.filtered <- dt.m6a.filtered %>% left_join(x=.,y=dt.m6a.filtered.sumwidth  %>% dplyr::distinct(name,SumWidth), by=c("name")) %>%
              dplyr::mutate(SumWidth=case_when(is.na(SumWidth) ~ 0, .default = SumWidth))
            #the overlap gene
            dt.m6a.filtered.overlapped.gene <- bedtoolsr::bt.intersect(a=dt.m6a.filtered, b=fread(org.genebed),
                                                            s=T, f=0.5, F=0.8, c=T,e = T) %>% as.data.table()
            dt.m6a.filtered <- dt.m6a.filtered %>% left_join(x=.,y=dt.m6a.filtered.overlapped.gene %>% dplyr::select(name=V4,Sample=V7,OverlappedNoGene=V9) %>%
                                                               dplyr::distinct(Sample,name,OverlappedNoGene), by=c("Sample","name"))
            dt.m6a.filtered %>% mutate(Method=M, ComboName=combo, optionID=dt.method.combo.selected[Method==M & ComboName==combo,optionID])
          }
        }
      }
      #merge m6a
      message(paste0("[",Sys.time(),"] ","step 3. Merge peak at FPR at cutoff ",ifelse(is.null(Model.FPR.cutoff), "default ", Model.FPR.cutoff)))
      if(is.null(dt.m6a.FPR)){
        dt.unmerged.m6a <- data.table(seqnames=NULL,start=NULL,end=NULL,name=NULL,score=NULL,strand=NULL,Sample=NULL)
      }else{dt.unmerged.m6a <- dt.m6a.FPR %>% dplyr::filter(!(SumWidth>3000 & OverlappedNoGene>2))}
      if(nrow(dt.unmerged.m6a)>0){
        #dt.topn.m6a[SumWidth<=6000,.N,by=.(Sample,Method,ComboName,optionID)]
        #merge together
        dt.m6a.filtered.pos <- dt.unmerged.m6a %>% dplyr::filter(strand=="+") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% distinct()
        if(nrow(dt.m6a.filtered.pos)>0){
          dt.m6a.filtered.pos.gr <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.filtered.pos,keep.extra.columns = T)
          dt.m6a.filtered.pos.merged <- foreach(sam = unique(dt.m6a.filtered.pos$Sample), .combine='rbind')%do%{
            GenomicRanges::reduce(dt.m6a.filtered.pos.gr[dt.m6a.filtered.pos.gr$Sample==sam,]) %>% as.data.table() %>%
              mutate(score=1000,name=paste(seqnames,start,end,strand,sep="_"),Sample=sam) %>%
              dplyr::select(seqnames,start,end,name,score,strand,Sample)
          }
        }else{
          dt.m6a.filtered.pos.merged <- data.table(seqnames=NULL,start=NULL,end=NULL,name=NULL,score=NULL,strand=NULL,Sample=NULL)
        }
        #filtered  and un overlapped m6a at neg strand
        dt.m6a.filtered.neg <- dt.unmerged.m6a %>% dplyr::filter(strand=="-") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% distinct()
        if(nrow(dt.m6a.filtered.neg)>0){
          dt.m6a.filtered.neg.gr <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.filtered.neg,keep.extra.columns = T)
          dt.m6a.filtered.neg.merged <- foreach(sam = unique(dt.m6a.filtered.neg$Sample), .combine='rbind')%do%{
            GenomicRanges::reduce(dt.m6a.filtered.neg.gr[dt.m6a.filtered.neg.gr$Sample==sam,]) %>% as.data.table() %>%
              mutate(score=1000,name=paste(seqnames,start,end,strand,sep="_"),Sample=sam) %>%
              dplyr::select(seqnames,start,end,name,score,strand,Sample)
          }
        }else{
          dt.m6a.filtered.neg.merged <- data.table(seqnames=NULL,start=NULL,end=NULL,name=NULL,score=NULL,strand=NULL,Sample=NULL)
        }
        dt.m6a.filtered <- rbind(dt.m6a.filtered.pos.merged, dt.m6a.filtered.neg.merged)
        #calculate exonic sum width to filter extreme long peak (<=6kb)
        dt.m6a.filtered.sumwidth <- LongPeakRemoveIntron(dt.peak = dt.m6a.filtered %>% distinct(seqnames,start,end,name,score,strand),
                                                         dt.intron = fread(org.intron.bed),
                                                         genome_file = org.genome_file)
        dt.m6a.filtered <- dt.m6a.filtered %>% left_join(x=.,y=dt.m6a.filtered.sumwidth  %>% dplyr::distinct(name,SumWidth), by=c("name")) %>%
          dplyr::mutate(SumWidth=case_when(is.na(SumWidth) ~ 0, .default = SumWidth))
        #check the overlapped gene count of merged m6a
        # gene.bed.file <- "/data/m6A_calling_strategy/RIPPeakS_resource/mm39_gencode.gene.bed"
        dt.m6a.filtered.overlapped.gene <- bedtoolsr::bt.intersect(a=dt.m6a.filtered, b=fread(org.genebed),s=T, f=0.5, F=0.8, c=T,e = T) %>% as.data.table()
        #dt.m6a.filtered.overlapped.gene[,.N,by=.(Sample=V7,nGene=V9)] %>% dplyr::arrange(Sample,nGene)
        #filter out merged peak with more than two genes and exonic width over 3kb
        dt.m6a.notkept <- dt.m6a.filtered.overlapped.gene[V9>2 & V8>3000,]
        #dt.m6a.notkept <- dt.m6a.filtered.overlapped.gene[V9>1,]
        dt.m6a.filtered <- dt.m6a.filtered %>% dplyr::filter(!name %in% dt.m6a.notkept$V4) %>% dplyr::filter(str_detect(seqnames,"chr")) #key filtering
        dt.m6a.overlapgene <- bedtoolsr::bt.intersect(a=dt.m6a.filtered, b=fread(org.genebed),s=T, f=0.5, F=0.8, wo=T,e = T) %>% as.data.table()
        #obtain the intersect of GLORI_mESC m6a site
        if(!is.null(org.benchmark.bed)){
          if(file.exists(org.benchmark.bed)){
            dt.m6a.overlap <- bedtoolsr::bt.intersect(a=dt.m6a.filtered, b=fread(org.benchmark.bed)[,1:6], s=T, wo=T) %>% as.data.table() %>% distinct()
            dt.benchmark.sites <- fread(org.benchmark.bed)[,c(1:6,9)]
            dt.benchmark.gene <- fread(org.benchmark.bed)[,c(1:6,9)] %>% dplyr::rename(OverlappedGenes=V9) %>%
              tidyr::separate_longer_delim(OverlappedGenes,delim = "|") %>% as.data.table()
            benchmark.gene <- unique(dt.benchmark.gene[OverlappedGenes!="",OverlappedGenes])
            #calculate TPR and FPR
            dt.TPR.FPR <- foreach(width.cut = c(Inf),.combine='rbind')%do%{
              dt.m6a.overlap %>% dplyr::filter(V8<=width.cut) %>%  group_by(Sample=V7) %>% mutate(nTP=n_distinct(V4), nSite=n_distinct(V12), TPR=n_distinct(V12)/nrow(dt.benchmark.sites)) %>%
                as.data.table() %>% distinct(Sample,nTP,nSite,TPR) %>% left_join(x=.,y=dt.m6a.filtered[SumWidth<=width.cut,list(nPeak=.N),by=.(Sample)],by=c("Sample")) %>% mutate(FPR=1-nTP/nPeak) %>%
                mutate(Strategy="m6APeakS", InvolvedMethod=dt.method.combo.selected$Method %>% unique() %>% paste(collapse = ","),
                       FPR=round(FPR,2), TPR=round(TPR,2)) %>% mutate(LongPeak=paste0("<=",round(width.cut/1000),"kb"))
            } %>% mutate(MergedFPR=Model.FPR.cutoff)
            message("The average FPR of merged m6a at FPR<=",Model.FPR.cutoff, " is ", mean(dt.TPR.FPR$FPR))
            dt.genelevel.TPR.FPR <- foreach(width.cut = c(Inf),.combine='rbind')%do%{
              dt.m6a.overlapgene %>% dplyr::filter(V8<=width.cut) %>%  group_by(Sample=V7) %>% mutate(nTP=n_distinct(V12[V12 %in% benchmark.gene]), nGene=n_distinct(V12)) %>%
                mutate(FPR=round(1-nTP/nGene,2), TPR=round(nTP/length(benchmark.gene),2)) %>%
                as.data.table() %>% distinct(Sample,nTP,nGene,TPR,FPR) %>%
                mutate(Strategy=paste0("M6APeakS_",selectedModel), InvolvedMethod=dt.method.combo.selected$Method %>% unique() %>% paste(collapse = ","),
                       FPR=round(FPR,2), TPR=round(TPR,2)) %>% mutate(LongPeak=paste0("<=",round(width.cut/1000),"kb"))
            } %>% mutate(MergedFPR=Model.FPR.cutoff)
            message("The average TPR of merged m6a at FPR<=",Model.FPR.cutoff, " is ", mean(dt.TPR.FPR$TPR))
            dt.TPR.FPR <- dt.TPR.FPR %>% left_join(x=.,y=dt.genelevel.TPR.FPR %>% dplyr::select(Sample,nTP_gene=nTP,nGene,TPR_gene=TPR,FPR_gene=FPR),by=c("Sample"))
            print(dt.TPR.FPR)
          }else{message("Could not find benchmark m6A sites!")}
        }
        dt.m6a.filtered <- dt.m6a.filtered %>% left_join(x=.,y=dt.m6a.overlapgene %>% dplyr::select(name=V4,Sample=V7,OverlappedGene=V12),by=c("name","Sample")) %>% group_by(name,Sample) %>%
          mutate(nOverlappedGene=n_distinct(OverlappedGene[!is.na(OverlappedGene)]), OverlappedGenes=paste(OverlappedGene,collapse = "|")) %>% as.data.table() %>%
          dplyr::distinct(seqnames,start,end,name,score,strand,Sample,SumWidth,nOverlappedGene,OverlappedGenes)
      }

      #step4 calculate the bedtools intensity of final m6a peak (depth_Ex_bin100_P2)
      #save m6a peaks
      for(s in Samples){fwrite(dt.m6a.filtered %>% dplyr::filter(Sample==s),file=paste0(out.dir,"/",s,"_M6APeakS_peak_", selectedModel, "_", ifelse(is.null(Model.FPR.cutoff),"default",Model.FPR.cutoff),".tsv"))}
    }
  }
  Sys.time()-t1
  #remove intermediate files
  if(rm.inter==TRUE){system(command = paste0("rm -rf ",tmp.dir), wait = T)}
}
