#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
# out and err files are in your home folder. Organize according to your preferences. 
#$ -e home/dave/lmamr44_megahit/humann2_err
#$ -o home/dave/lmamr44_megahit/humann2_out
#$ -q Single.q
unset SGE_ROOT
#echo $JOB_ID
#NSLOTS=4


cat /home/dave/lmamr44_analysis_ready/$1.fw.gz /home/dave/lmamr44_analysis_ready/$1.rw.gz /home/dave/lmamr44_analysis_ready/$1.se.gz > /home/dave/lmamr44_analysis_ready/$1.all.gz

/home/dave/.local/bin/humann2  --input /home/dave/lmamr44_analysis_ready/$1.all.gz --output /state/partition1/$1\_humann2 --metaphlan /home/dave/tools/metaphlan_dir/ --bowtie2 /share/apps/bin/ --diamond /home/dave/tools/ --protein-database /home/dave/tools/humann2-0.11.1/ --threads 4

#or with the uniref50db
# /home/dave/.local/bin/humann2  --input $path/concat_fastqs/$1.all.fastq.gz --output /state/partition1/$1\_humann2 --metaphlan /home/dave/tools/metaphlan_dir/ --bowtie2 /share/apps/bin/ --diamond /home/dave/tools/ --protein-database /home/dave/tools/humann2-0.11.1/uniref/ --threads 4

rm /home/dave/lmamr44_analysis_ready/$1.all.gz

scp -C -l 500000 /state/partition1/$1\_humann2/* 10.1.1.1:/home/dave/lmamr44_megahit/$1\_humann2/

rm -r /state/partition1/$1\_humann2/
