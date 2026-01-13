mamba create -n metaphlan4 -c conda-forge -c bioconda metaphlan mamba

cd PATH_TO_metaphlan4_database/mpa_vJan25_CHOCOPhlAnSGB_202503/
wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_latest ./
cat mpa_latest

wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vJan25_CHOCOPhlAnSGB_202503_bt2.tar ./
wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan25_CHOCOPhlAnSGB_202503.md5 ./
wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan25_CHOCOPhlAnSGB_202503.nwk ./
wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan25_CHOCOPhlAnSGB_202503.tar ./
wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan25_CHOCOPhlAnSGB_202503_marker_info.txt.bz2 ./
wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan25_CHOCOPhlAnSGB_202503_species.txt.bz2 ./

metaphlan L1EIE3101239-HX1.R1.raw.fastq.gz,L1EIE3101239-HX1.R2.raw.fastq.gz --input_type fastq -o L1EIE3101239-HX1.txt --nproc 10 --mapout L1EIE3101239-HX1.btout.bz2 --db_dir PATH_TO_metaphlan4_database/mpa_vJan25_CHOCOPhlAnSGB_202503
