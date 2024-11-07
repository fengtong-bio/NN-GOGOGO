######################################################################
# install sylph 20241107
mamba create -n sylph -c conda-forge -c bioconda mamba sylph

# https://www.nature.com/articles/s41587-024-02412-y
# https://github.com/bluenote-1577/sylph
######################################################################
# Sketching a large database

export PATH=PATH/envs/sylph/bin:$PATH
sylph sketch -l rhi_1377_genome_list.txt -o rhi_1377_sylph -t 40

time for line in $(cat DJ_file_name.txt)
do

sylph profile rhi_1377_sylph.syldb -1 $line\_meta_clean_R1.fq.gz -2 $line\_meta_clean_R2.fq.gz -u -t 80 -o $line\_.results.tsv

done
