#curation of demo BAM files for library
require(data.table)
require(foreach)
require(tidyverse)
#load M6APeakS peaks form /data2/mESC_m6A/M6APeakS/NEB_GSE*M6APeakS_peak_NEB_default.tsv
NEB_peak_files <- list.files("/data2/mESC_m6A/M6APeakS/",pattern = "^(NEB_GSE).+(M6APeakS_peak_NEB_default.tsv)$",full.names = T)
names(NEB_peak_files) <- list.files("/data2/mESC_m6A/M6APeakS/",pattern = "^(NEB_GSE).+(M6APeakS_peak_NEB_default.tsv)$",full.names = F) %>%
  strsplit(split="_M6APeakS") %>% sapply("[",1)
NEB_bam_files <- list.files("/data2/mESC_m6A/7_RIPPeakS/BAM/",pattern = "^(NEB_GSE).+(.sortedByCoord.UN.bam)$",full.names = T)
names(NEB_bam_files) <- list.files("/data2/mESC_m6A/7_RIPPeakS/BAM/",pattern = "^(NEB_GSE).+(.sortedByCoord.UN.bam)$",full.names = F) %>%
  strsplit(split=".sortedByCoord.UN.bam") %>% sapply("[",1)

#calculate mapped read in each chromosome for each BAM files
dt.readratio.bam <- foreach(i = 1:length(NEB_peak_files),.combine = 'rbind')%do%{
  idxstats.input <- read.table(pipe(paste0("~/anaconda3/envs/TestMyPackage/bin/samtools idxstats ", NEB_bam_files[paste0(names(NEB_peak_files)[i],"_Input")])), sep="\t",
                         col.names = c("chrom", "length", "mapped", "unmapped"))
  idxstats.input$ratio <- idxstats.input$mapped / sum(idxstats.input$mapped)*100
  idxstats.input <- idxstats.input[str_detect(idxstats.input$chrom,"chr"), c("chrom", "ratio")] %>% mutate(condition="Input")
  idxstats.rip <- read.table(pipe(paste0("~/anaconda3/envs/TestMyPackage/bin/samtools idxstats ", NEB_bam_files[paste0(names(NEB_peak_files)[i],"_RIP")])), sep="\t",
                               col.names = c("chrom", "length", "mapped", "unmapped"))
  idxstats.rip$ratio <- idxstats.rip$mapped / sum(idxstats.rip$mapped)*100
  idxstats.rip <- idxstats.rip[str_detect(idxstats.rip$chrom,"chr"), c("chrom", "ratio")] %>% mutate(condition="RIP")
  rbind(idxstats.input,idxstats.rip) %>% mutate(ratio=round(ratio,2)) %>% pivot_wider(id_cols = "chrom",values_from = "ratio",names_from = "condition") %>%
    as.data.table() %>% mutate(sample=names(NEB_peak_files)[i])
}
dt.readratio.bam %>% dplyr::filter(chrom=="chr1")
dt.readratio.bam %>% dplyr::filter(chrom=="chr13")

dt.peakratio <- foreach( i = 1:length(NEB_peak_files),.combine = 'rbind')%do%{
  # dt.peak <- fread(NEB_peak_files[i]) %>% slice_sample(n = 500) %>% dplyr::select(seqnames,start,end)
  # dt.peak <- fread(NEB_peak_files[i]) %>% filter(seqnames=="chr21")
  dt.peak <- fread(NEB_peak_files[i])
  dt.peak[,.(nPeak=.N),by=.(seqnames)]
  # fwrite(dt.peak, file=paste0("~/software/mysterypackage/mysterypackage/inst/extdata/",names(NEB_peak_files)[i],".subset.regions.bed"),row.names = F,col.names = F,sep="\t")
}
dt.peakratio %>% dplyr::filter(seqnames=="chr1")
dt.peakratio %>% dplyr::filter(seqnames=="chr13")


# samtools view -b -L regions.bed input.bam > output.selected.bam
#use chr13
selected_chr="chr13"
selected_length=(read.table(pipe(paste0("~/anaconda3/envs/TestMyPackage/bin/samtools idxstats ", NEB_bam_files[paste0(names(NEB_peak_files)[i],"_Input")])), sep="\t",
           col.names = c("chrom", "length", "mapped", "unmapped")) %>% filter(chrom==selected_chr))$length
fwrite(data.table(seqnames=selected_chr,start=0,end=selected_length), file=paste0("~/software/mysterypackage/mysterypackage/inst/extdata/mm39_",selected_chr,".bed"),row.names=F,col.names = F,sep="\t")
dt.bam.subset.cmd <- foreach(i = 1:length(NEB_peak_files), .combine = 'c')%do%{
  cmd.input <- paste0("samtools view --threads 6 -b -L ~/software/mysterypackage/mysterypackage/inst/extdata/mm39_", selected_chr, ".bed ",
                      NEB_bam_files[paste0(names(NEB_peak_files)[i],"_Input")], " >~/software/mysterypackage/mysterypackage/inst/extdata/",names(NEB_peak_files)[i],selected_chr, "_Input.bam &")
  cmd.rip <- paste0("samtools view --threads 6 -b -L ~/software/mysterypackage/mysterypackage/inst/extdata/mm39_", selected_chr,".bed ",
                      NEB_bam_files[paste0(names(NEB_peak_files)[i],"_RIP")], " >~/software/mysterypackage/mysterypackage/inst/extdata/",names(NEB_peak_files)[i],selected_chr,"_RIP.bam &")
  c(cmd.input,cmd.rip)
}
cat(dt.bam.subset.cmd,sep="\n")

#rename the demo bam
# mv ~/software/mysterypackage/mysterypackage/inst/extdata/NEB_GSE146465_mESC_M3WT_rep1chr13_Input.bam ~/software/mysterypackage/mysterypackage/inst/extdata/NEB_GSE4_mESC1_Input.bam &
# mv ~/software/mysterypackage/mysterypackage/inst/extdata/NEB_GSE146465_mESC_M3WT_rep1chr13_RIP.bam ~/software/mysterypackage/mysterypackage/inst/extdata/NEB_GSE4_mESC1_RIP.bam &
# mv ~/software/mysterypackage/mysterypackage/inst/extdata/NEB_GSE183399_mESC_M3Ctr_rep1chr13_Input.bam ~/software/mysterypackage/mysterypackage/inst/extdata/NEB_GSE8_mESC1_Input.bam
# mv ~/software/mysterypackage/mysterypackage/inst/extdata/NEB_GSE183399_mESC_M3Ctr_rep1chr13_RIP.bam ~/software/mysterypackage/mysterypackage/inst/extdata/NEB_GSE8_mESC1_RIP.bam


