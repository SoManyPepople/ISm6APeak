#' RunexomePeak: Process Exome Peak Data
#'
#' @description Brief description.
#' @param ... Arguments.
#' @return Processed data.
#' @export
#function to run exomePeak
RunexomePeak <- function(InputBAM="~/m6A_calling_strategy/6_alignment_star/R2_bam/mESC_NEB_WT_rep2_Input_STAR.Aligned.out.SortedByCoord.unique.deduplicate.R2.bam",
                         RIPBAM="~/m6A_calling_strategy/6_alignment_star/R2_bam/mESC_NEB_WT_rep2_RIP_STAR.Aligned.out.SortedByCoord.unique.deduplicate.R2.bam",
                         outdir="~/m6A_calling_strategy/exomePeak",
                         prefix="mESC_NEB_WT_rep2_UD_test",
                         gtf="~/genome_db/mm39/gencode.GRCm39.vM30.annotation.gtf",
                         window.width=50,
                         sliding.step=50,
                         #mininal_peak_length=50,
                         peak_cutoff_pvalue=c(1,0.05,0.01,1e-5,1e-6),#default 1e-5
                         peak_cutoff_fdr = c(NA,1,0.05,0.01,0.001),#If it is specified, then use "fdr" instead of "p" in peak calling
                         fold_enrichment=1
){
  if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
  if(!file.exists(paste0(outdir,"/",prefix,"/peak.xls"))){

    suppressMessages(library(GenomicFeatures))
    detach(package:GenomicFeatures,force=T)
    suppressMessages(library(txdbmaker))
    suppressMessages(library(GenomicFeatures))

    suppressMessages(library(AnnotationDbi))
    detach(package:AnnotationDbi,force=T)
    suppressMessages(library(AnnotationDbi))
    suppressMessages(library(exomePeak))
    res <- exomePeak::exomepeak(GENE_ANNO_GTF=gtf,
                                IP_BAM = RIPBAM,
                                INPUT_BAM = InputBAM,
                                EXPERIMENT_NAME=prefix,
                                OUTPUT_DIR = outdir,
                                WINDOW_WIDTH=window.width,
                                SLIDING_STEP=sliding.step,
                                #MINIMAL_PEAK_LENGTH=mininal_peak_length,
                                PEAK_CUTOFF_PVALUE=peak_cutoff_pvalue,
                                PEAK_CUTOFF_FDR=peak_cutoff_fdr,
                                FOLD_ENRICHMENT=fold_enrichment
    )
  }else{message(paste0("Detected existed peak results, hence skip ", prefix))}

}
