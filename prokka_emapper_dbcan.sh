#rm -rf MAG_prokka MAG_eggnog MAG_cazy
#mkdir MAG_prokka
#mkdir MAG_eggnog
#mkdir MAG_cazy

for i in $(cat fa_file.list)
do

export PATH=/PATH/mamba/mambaforge/bin:$PATH
export PATH=/PATH/mamba/mambaforge/envs/prokka/bin:$PATH
prokka /PATH/$i --outdir /PATH/MAG_prokka/$i --prefix $i --metagenome --force --cpus 2 --centre X --compliant

export PATH=/PATH/mamba/mambaforge/envs/eggnog/bin:$PATH
mkdir /PATH/MAG_eggnog/$i
emapper.py -i /PATH/MAG_prokka/$i/$i.faa --output_dir /PATH/MAG_eggnog/$i --cpu 2 -o $i -d bact
emapper.py  -m diamond -i /PATH/MAG_prokka/$i/$i.faa --output_dir /PATH/MAG_eggnog/$i -o $i --cpu 2 --data_dir /PATH/software/SK_conda/eggNOG_database --pfam_realign

export PATH=/PATH/mamba/mambaforge/envs/cazy/bin:$PATH
mkdir /PATH/MAG_cazy/$i
run_dbcan --out_dir /PATH/MAG_cazy/$i --db_dir /PATH/software/cazy_db --hmm_cpu 2 /PATH/MAG_prokka/$i/$i.faa protein -t hmmer

export PATH=/PATH/mamba/mambaforge/envs/dbcan4/bin:$PATH
mkdir /PATH/MAG_dbCAN/$i
run_dbcan --out_dir /PATH/MAG_dbCAN/$i --db_dir /PATH/MAG_prokka/$i/$i.faa protein --tools dbcansub --dbcan_thread 1

run_dbcan easy_substrate --input_raw_data /PATH/MAG_prokka/$i/$i.faa --mode protein --output_dir /PATH/MAG_dbCAN/$i --db_dir /PATH/software/SK_conda/CAZy_dbcan_database --input_gff /PATH/MAG_prokka/$i/$i.faa --gff_type prodigal --hmm_cpu 2 --dbcan_thread 2 --log-level INFO

done
