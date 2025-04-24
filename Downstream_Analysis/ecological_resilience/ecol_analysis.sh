#!/bin/sh
source activate qiime1

#analyze ecological resilience metrics using a suite of rscripts
cat /Users/dave/Desktop/reference_based_mapping/ddig.txt | while read pathway name; do

    #run path rich and unclass
    awk -F "\t" 'NR == 1 || $1 ~ /g__/ || $1 ~ /unclass/' $pathway\_humann2.txt > $pathway\_processed.txt
    sed 's/# Pathway/Pathway|Species/' $pathway\_processed.txt | tr "|" "\t" > $pathway\_forR.txt

    Rscript pathway_richness.r $pathway\_forR.txt $pathway\_a.txt $pathway $name
    Rscript unclassAbund.R $pathway\_forR.txt $pathway $name 


    #run faithpd
    sed 's/g__.*s__//g' $pathway\com_matrix.txt | awk -F "\t" '{print$1}' | while read line; do grep "$line" -m 1 eztaxon_id_taxonomy.txt | awk -v a="$line" '{print$1"\t"a}' ; done | while read number taxa ; do grep "$number" -A 1 eztaxon_qiime_full.fasta | sed "s/${number}/${taxa}/g" ;done  > $1\_taxa.fasta
    sed 's/g__.*s__//' $pathway\com_matrix.txt > $pathway\edited_comm.txt

    mafft $pathway\_taxa.fasta > $pathway\_align.fasta
    make_phylogeny.py -i $pathway\_align.fasta -o $pathway\_R.tre

    Rscript faithPD.R $pathway\_R.tre $pathway $pathway\edited_comm.txt $name

done



Rscript /Users/dave/Desktop/reference_based_mapping/ecol_plot.r 

rm /Users/dave/Desktop/reference_based_mapping/delete.txt
rm /Users/dave/Desktop/reference_based_mapping/tempRich.txt
rm /Users/dave/Desktop/reference_based_mapping/*R.tre
rm /Users/dave/Desktop/reference_based_mapping/*align.fasta
rm /Users/dave/Desktop/reference_based_mapping/*taxa.fasta

source deactivate qiime1

#



