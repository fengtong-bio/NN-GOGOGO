######################################################################
# install DiffSegR 20241105
mamba create -n DiffSegR -c conda-forge mamba r-base=4.4.1

-- New R --
install.packages("remotes")
#1st install
remotes::install_github("aLiehrmann/DiffSegR")

install.packages("BiocManager")
BiocManager::install("GenomicAlignments")
BiocManager::install("Rsamtools")
BiocManager::install("DESeq2")
BiocManager::install("Rsubread")
BiocManager::install("rtracklayer")

install.packages("pak")
pak::pkg_install("sanssouci-org/sanssouci")

#2nd install
remotes::install_github("aLiehrmann/DiffSegR")

######################################################################
# use DiffSegR 20241105
library(DiffSegR)
cd "working_directory"

sample_info=read.table("DiffSegR_F_sample_info.txt", header = T, sep="\t", check.names = F)
knitr::kable(sample_info, row.names = FALSE)
nb_threads_tot = 40  #threads
nb_threads_locus = 40  #threads
working_directory <- "/pub3/Liu-group/fengtong/project_ruan/fish_40_MF/diff_contig/DiffSegR_F"  #working_directory
data <- newExperiment(
  sampleInfo   = sample_info,
  loci         = data.frame(
    seqid      = "CHR_FF",   #chr name
    chromStart = 1, 
    chromEnd   = 743556,  #chr length
    locusID    = "FF"  #ID
  ),
  referenceCondition = "Group_F",
  otherCondition     = "Group_M",
  nbThreads          = nb_threads_tot,
  nbThreadsByLocus   = nb_threads_locus,
  coverage           = working_directory
)
print(data)
coverage(
  data         = data,
  coverageType = "average"
)
features <- segmentationLFC(
  data  = data, 
  alpha = 2
)
SExp <- counting(
  data     = data, 
  features = features
)
SummarizedExperiment::strand(SExp)
segment_coordinates <- as.data.frame(SummarizedExperiment::rowRanges(SExp))
knitr::kable(segment_coordinates[1:5, colnames(segment_coordinates)%in%c(
  "seqnames",
  "start",
  "end",
  "width",
  "strand")
]) 
counts <- SummarizedExperiment::assay(SExp)
knitr::kable(counts[1:5,])
dds <- dea(
  SExp              = SExp,
  design            = ~ condition,
  sizeFactors       = rep(1,40),  #all samples 40
  significanceLevel = 1e-1  #P 1e-2 = 0.01;1e-1 = 0.1
)
DERs <- dds[SummarizedExperiment::mcols(dds)$DER,]
length(DERs)
DERs <- as.data.frame(SummarizedExperiment::rowRanges(DERs))
knitr::kable(DERs[1:5, colnames(DERs)%in%c(
  "seqnames",
  "start",
  "end",
  "width",
  "strand",
  "log2FoldChange",
  "padj")
])
