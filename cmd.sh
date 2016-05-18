# using pre-fix of fastq file 
python ATACseq.py  -P 3 -M 100 -O test_out -R -S liver -G mm9 -Q paired  -C Configure_ATACseq.yaml  -I test_data/liver-CD31_S19_L001_R1_001.fastq.gz -I2 test_data/liver-CD31_S19_L001_R2_001.fastq.gz 
