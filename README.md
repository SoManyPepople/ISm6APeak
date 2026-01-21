## Manual for ISm6APeak

## Prepare require software via conda

### create conda enviroment and install required packages

```         
conda create -n ISm6APeak_env -c conda-forge -c bioconda -c defaults macs2 bedtools rseqc samtools bioconda::bioconductor-exomepeak conda-forge::r-devtools r-base=4.2 conda-forge::r-devtools bioconda::bioconductor-rtracklayer bioconda::bioconductor-genomicfeatures
```

```         
conda activate ISm6APeak_env
```

```         
conda install bioconda::bioconductor-exomepeak
```

### Install ISm6APeak from github

```         
devtools::install_github("SoManyPepople/ISm6APeak")
```

### Collected resource

ISm6APeak require five types of genomic annotation files as input:\
1.GTF file required by the integrated tools like exomePeak, MeTPeak, and exomePeak2 \
2. knownGene sqlite file required by the integrated tool: TRESS \
3. genomic chrNameLength file required by the integrated tools and ISm6APeak \
4. genomic position of all introns(bed file) required by ISm6APeak \
5. genomic position of all genes (bed file) required by ISm6APeak

There were some downloadable annotation files for hg38 and mm39 in [here](https://doi.org/10.6084/m9.figshare.31077001).

### QC and alignment to prepare BAM files

We recommond STAR for alignment and bamtools to filtered for unique mapped reads as input BAM files for ISm6APeak \
To test the installtion of ISm6APeak, demo BAM files could be downloaded from [here](https://doi.org/10.6084/m9.figshare.31077001)

### Perform m6A peak detection

```         
runM6APeakS(
  InputBAMs="full_path_to_input_BAM_files",
  RIPBAMs="full_path_to_RIP_BAM_files",
  bin.dir="/path/to/conda/envs/bin/",
  out.dir="/path/to/outputdir",
  script.dir=system.file("scripts",package ="ISm6APeak"),
  Organism=c("Mouse","Human")[1],
  org.gtf_file="/path_to_annotation_dir/gencode.vM33.annotation.gtf",
  org.sqlite_file="/path_to_annotation_dir/UCSC.mm39.knownGene.sqlite",
  org.genome_file="/path_to_annotation_dir/chrNameLength_mm39.txt",
  org.intron.bed="/path_to_annotation_dir/dt.intron.mm39.bed",
  org.genebed="/path_to_annotation_dir/mm39_gencode.gene.bed",
  ModelCutoffTable=system.file("extdata","dt.cutoff.at.all.FPR.M6APeakS.tsv",package ="ISm6APeak"),
  selectedModel=c("NEB","SYSY","Millipore","Abcam")[1],
  Model.FPR.cutoff=NULL,#could be numeric value such as 0.01-1, or one of character variables including "Strict",  "Moderate", "Sensitive" 
  org.benchmark.bed=NULL,
  n.cores=75,
  rm.inter=TRUE
)
```
