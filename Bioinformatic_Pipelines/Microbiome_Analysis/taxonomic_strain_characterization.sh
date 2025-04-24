#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
# out and err files are in your home folder. Organize according to your preferences.
#$ -e /home/dave/lmamr44_megahit/metaphlan.err
#$ -o /home/dave/lmamr44_megahit/metaphlan.out
#$ -q Single.q
unset SGE_ROOT
#echo $JOB_ID
#NSLOTS=1

path=/home/dave/lmamr44_analysis_ready


/home/dave/MBIO5623/python2/Python-2.7.14/python /home/dave/tools/metaphlan_dir/metaphlan2.py  <(zcat $path/$1.fw.fastq.gz $path/$1.rw.fastq.gz $path/$1.se.fastq.gz)  --bowtie2_build /share/apps/bin/bowtie2-build --bowtie2_exe /share/apps/bin/bowtie2 --input_type fastq --mpa_pkl /home/dave/tools/metaphlan_dir/db_v20/mpa_v20_m200/mpa_v20_m200.pkl --bowtie2db /home/dave/tools/metaphlan_dir/db_v20 --bowtie2out /state/partition1/$1.bowtie2.bz2 -s /state/partition1/$1.sam.bz2 > /state/partition1/$1\_metaphlan_profile.txt

/home/dave/MBIO5623/python2/Python-2.7.14/python /home/dave/tools/metaphlan_dir/strainphlan_src/sample2markers.py --ifn_samples /state/partition1/$1.sam.bz2 --input_type sam --output_dir /home/dave/lmamr44_megahit/sample2markers/ --samtools_exe /home/dave/tools/metaphlan_dir/metaphlan_bin/samtools --bcftools_exe /home/dave/tools/metaphlan_dir/metaphlan_bin/bcftools

/home/dave/MBIO5623/python2/Python-2.7.14/python /home/dave/tools/metaphlan_dir/strainphlan.py --mpa_pkl /home/dave/tools/metaphlan_dir/db_v20/mpa_v20_m200/mpa_v20_m200.pkl --ifn_samples /home/dave/lmamr44_megahit/sample2markers/*.markers --output_dir /home/dave/lmamr44_megahit/sample2markers/  &> /home/dave/lmamr44_megahit/sample2markers/strain.log

scp -C -l 500000 /state/partition1/$1*.txt 10.1.1.1:/home/dave/lmamr44_megahit/

rm /state/partition1/$1.*

