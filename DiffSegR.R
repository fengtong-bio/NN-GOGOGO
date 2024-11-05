######################################################################
# install DiffSegR
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
# use DiffSegR
