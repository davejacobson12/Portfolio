#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
# out and err files are in your home folder. Organize according to your preferences.
#$ -e /home/dave/assembly.err
#$ -o /home/dave/assembly.out
#$ -q all.q
unset SGE_ROOT
#echo $JOB_ID
### provide these arguments again

path=/home/dave/burkina_analysis_ready
outpath=/home/dave/burkina_assembly

/home/dave/MBIO5623/megahit_v1.1.3_LINUX_CPUONLY_x86_64-bin/megahit -r $path/$1.se.fastq.gz -1 $path/$1.fw.fastq.gz -2 $path/$1.rw.fastq.gz -t 32 --k-list 21,31,41,61,81,101 -o /state/partition1/$1\_megahit

echo "Sample N50 TotalBP AlignmentRate"; cat /home/dave/burkina_samps.txt  | while read line; do  echo $line | tr "\n" " "; grep "N50" $path/$line\_megahit/log | grep "avg" | awk '{print$(NF-1)" "$(NF-13)}'| tr "\n" " ";grep "alignment rate" $path/$line\_bowtie2/$line.log | awk '{print$1}'; done > assembly_stats.txt


mkdir /home/dave/burkina_assembly/$1\_bowtie2
mkdir /home/dave/burkina_assembly/$1\_metabat
mkdir /home/dave/burkina_assembly/$1\_checkm

/share/apps/bin/bowtie2-build $outpath/$1\_megahit/final.contigs.fa $outpath/$1\_bowtie2/$1\_bowtie_ref

(/share/apps/bin/bowtie2 -p 5 -x $outpath/$1\_bowtie2/$1\_bowtie_ref -U $path/$1.se.fastq.gz -1 $path/$1.fw.fastq.gz -2 $path/$1.rw.fastq.gz --no-unal -S /state/partition1/$1.sam) 2>/state/partition1/$1.log

/share/apps/bin/samtools view -@ 5 -S -b -o /state/partition1/$1.bam /state/partition1/$1.sam

/share/apps/bin/samtools sort -@ 5 -o /state/partition1/$1.sorted.bam /state/partition1/$1.bam

/share/apps/bin/samtools rmdup -S /state/partition1/$1.sorted.bam /state/partition1/$1.rmdup.bam

scp -C -l 500000 /state/partition1/$1.sorted.bam 10.1.1.1:/home/dave/burkina_assembly/$1\_bowtie2/
scp -C -l 500000 /state/partition1/$1.log 10.1.1.1:/home/dave/burkina_assembly/$1\_bowtie2/

rm /state/partition1/$1.*

/home/dave/tools/metabat/jgi_summarize_bam_contig_depths --outputDepth /state/partition1/$1\_depth.txt --pairedContigs /state/partition1/$1\_paired.txt --minContigLength 750 --minContigDepth 2 $outpath/$1\_bowtie2/$1.sorted.bam

/home/dave/tools/metabat/metabat2 -a /state/partition1/$1\_depth.txt -i $outpath/$1\_megahit/final.contigs.fa -o /state/partition1/$1\_bins -m 1500

scp -C -l 500000 /state/partition1/$1* 10.1.1.1:/home/dave/burkina_assembly/$1\_metabat/

rm /state/partition1/$1*

/home/dave/tools/CheckM-1.0.12/bin/checkm lineage_wf $outpath/$1\_metabat/ $outpath/$1\_checkm/ -x fa -f $outpath/$1\_checkm/$1\_checkm.txt


