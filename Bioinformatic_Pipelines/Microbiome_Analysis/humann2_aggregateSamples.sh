#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
# out and err files are in your home folder. Organize according to your preferences.
#$ -e /home/dave/humann2_populations/diverse_pops_uniref50/humann2_downstream.err
#$ -o /home/dave/humann2_populations/diverse_pops_uniref50/humann2_downstream.out
#$ -q Single.q
unset SGE_ROOT
#echo $JOB_ID
#NSLOTS=1

path=/home/dave/humann2_populations/diverse_pops_uniref50

/home/dave/.local/bin/humann2_join_tables -i $path/humann2_tsv -o $path/humann2_downstream/diversePops_pathabund.tsv --file_name pathabund

#/home/dave/.local/bin/humann2_join_tables -i $path/humann2_tsv -o $path/humann2_downstream/diversePops_geneFam.tsv --file_name genefamalies

/home/dave/.local/bin/humann2_renorm_table -i $path/humann2_downstream/diversePops_pathabund.tsv -o $path/humann2_downstream/diversePops_pathabund_relab.tsv --units relab

#/home/dave/.local/bin/humann2_renorm_table -i $path/humann2_downstream/diversePops_geneFam.tsv -o $path/humann2_downstream/diversePops_geneFam_relab.tsv --units relab

cat $path/humann2_downstream/pathways.txt | while read line; do awk -F ":" -v var=$line 'NR==1 || $1==var' $path/humann2_downstream/diversePops_pathabund_relab.tsv > $path/humann2_downstream/$line\_humann2.txt ; done

