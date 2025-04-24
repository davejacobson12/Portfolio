#Generalized script for analyzing 16S rRNA V4 data.

#If the run needs to be downloaded from basepsace and converted to fastq, you will need to do the following steps. 
#Download BaseSpace Run
	#Attain BaseSpace token
	#	Login to developer.basespace.illumina.com
	#	Select MyApps
	#	Select Create New Application and follow prompts
	#		Description box can just write “temp”
	#	On applications page, select your app
	#	Navigate to Credentials tab and make sure you copy access token

	#Download BaseSpce run downloader (download .zip file)
	#	unzip the file
	#	create file called run_downloader.sh in nano
	#		On first line, type the following	
	#			python BaseSpaceRunDownloader_v2.py -r -a
	#			after the -a paste in your access token
	#		save the nano file
	#	you need to give permission to run the file
	#		chmod +x /path/to/file/run_downloader.sh
	#On BaseSpace, find the run ID
	#	illumina.com/run/numericrunID/runName
	#open run_downloader.sh in nano
	#	copy the numericrunID in space after -r
	#	Exit and save nano file
	#Execute file
	#	in folder with .sh file, type ./run_downloader.sh
	#		A new folder will be created and the files will be downloaded from BaseSpace

#Convert bcl files to fastq
	#within the downloaded run directory, you may need to remove the SampleSheet.csv file

# -R is the whole folder that was just downloaded
# -o output folder that will have forward, reverse, and index fastqs
# --use-bases-mask format of how your forward and reverse reads are generated. Y*,I*,
bcl2fastq -R bcl_folder -o output_folder --use-bases-mask Y*,I*,I*,Y* --create-fastq-for-index-reads

#Fastq to fasta 
#gunzip -c collapsed.gz | paste ---- | cut -f 1,2 | sed 's/^/>/' | tr "\t" "\n"

### Use adapterremoval2

AdapterRemoval --file1 sample_R1.fastq.gz --file2 sample_R2.fastq.gz --trimns --trimqualities --minquality 30 --minalignmentlength 40  --maxns 0 --collapse --gzip --basename sample 

### Determine sequence header names for reads that merged

grep "^@M" lmamr33.merge.fastq | awk '{print$1}' | sed 's/^@M_//' > seqnames.final

### Use sequence names to pull out the associated barcodes with each read
# -f input file of barcodes from MiSeq
# -s file of sequence names to be retained
# -o output barcode file

filter_fasta.py -f index_file.fastq -s seqnames.final -o filtered_barcodes.fastq

### demultiplex merged reads. Requires a properly formatted qiime metadata file - tab seperated
###   #SampleID BarcodeSequence LinkePrimerSequence Description
### Sample names cannot have underscores, LinkerPrimerSequence column can be empty, Description column has to be last column - can be filled with NAs. 
### Can make columns with metadata categories between LinkerPrimerSequence and Description
# -i merged output from pear
# -m properly formatted metadata file with barcode sequences for each sample
# -b output of filter_fasta for barcode files
# -o output directory
# --rev_comp_mapping_barcodes must be used if the barcodes in the metadata file are not reverse complemented. Otherwise, do not use this flag

split_libraries_fastq.py -i lmamr33.merge.fastq -m metadata.txt -b filtered_barcodes.fastq -o demultiplexed_directory

### Perl script written by Krithi to get list of unique sequences and how many times they occur
### This will be used for creating a de novo clustering database
# seqs.fna is directory from output of split_libraries_fastq.py
# derep.fna is fasta file of dereplicated sequences

#Rename fasta
ls *R1.fastq.gz | awk -F "_" '{print$1}' | while read line; do awk -F ":" '{print$1}' $line.fa| sed  "s/^>@M.*/>${line}/g" | awk '/>{print$0"_"(++I)}!/>/' $line\_renamed.fa; done
#Another option for the awk part of the statement
awk '/^>/{print$0"_" ++i; next}{print}'

perl Make_unique.pl seqs.fna > derep.fna

### Eliminate sequences that occur less than x number of times. The x depends on the data
# -sortbysize output of make_unique.pl command
# -fastaout new output for filtered sequences
# -minsize get rid of sequences that occur less than this number of times

usearch -sortbysize derep.fna -fastaout filterabundance.fa -minsize 5

### Cluster the filtered sequences by sequence similarity to create OTUs at 97%. This creates the de novo database
# -cluster_otus sequences filtered by abundance, output from sortbysize command
# -otus output file that is the de novo database
# -relabel used to relabel the OTUs in the de novo database, we use Parse by convention

vsearch -cluster_otus filterabundance.fa -otus rep_db.fa -relabel Parse

### Cluster the demultiplexed sequences against the newly created database
# -usearch_global seqs.fna from demultiplexing output
# -db de novo database, output from cluster_otus command
# -strand direction of reads, should always be plus
# -id clustering percent identity, by convention use 97%
# -uc output file in usearch format

vsearch -usearch_global seqs.fna -db rep_db.fa -strand plus -id 0.97 -uc otu_table.uc

### convert output of clustering to a otu table in tab seperated format, use script Krithi wrote

perl Convert_uc_otu_map.pl otu_table.uc otu_table.txt

### assign taxonomy from green genes to the Parse sequences in de novo database. Two part command
### print qiime config to get path to green genes reference taxonomy

print_qiime_config.py 

### assign taxonomy command requires path to green genes db on computer, which is output from print qiime config
# -i de novo reference databse
# -r reference sequence file path from print_qiime_config, on the line that starts "assign_taxonomy_reference_seqs_fp"
# -t reference taxonomy file path from print_qiime_config, on the line that starts "assign_taxonomy_id_to_taxonomy_fp"
# -o output directory

assign_taxonomy.py -i rep_db.fa -r path_to_green_genes_reference_sequences -t path_to_green_genes_taxonomy -o assigned_taxonomy

assign_taxonomy.py -i otus.fa -r /home/dave/Desktop/eztaxon_.fasta -t /home/.txt -o assigned_taxa 

### make the otu table with taxonomies
# -i otu table in text 
# -t path to output from assign_taxonomy, use the text file
# -o output biom file

make_otu_table.py -i otu_table.txt -t taxonomy_assigned/rep_db_tax_assignments.txt -o otus.biom 

### Rarify your biom file, typically to 10,000 reads
### Do this in three sets, from 10 to 100 by 10s, 200 to 1000 by 100s, 2000 to 10000 by 1000s
# -i biom file
# -x max number of sequences per sample
# -m min number of sequences per sample
# -s steps between min and max
# -n number of iterations at each step, typically 10
# -o output folder
# -O number of threads

-x 1000 -m 100 -s 100 n -10
-x 10000 -m 1000 -s 1000 -n 10

parallel_multiple_rarefactions.py -i otus.biom -x 100 -m 10 -s 10 -n 10 -o rarefactions -O 16

### create tree that will be used for phylogenetic analysis. Two step command


mafft rep_db.fa > file_for_tree.fa

# -i input file, representative database
# -o output newick tree file in 

make_phylogeny.py -i file_for_tree.fa -o rep_tree.tre

### Generate alpha diversity metrics for rarefied sequences
# -i rarefactions folder
# -m alpha diversity metrics to generate, seperated by commas. Typically species richness (observed OTUs) and Faith's phylogenetic diversity (PD_whole_tree)
# -t generated tree from make phylogeny
# -o output directory
# -O number of threads

parallel_alpha_diversity.py -i rarefactions -m observed_otus,PD_whole_tree -t rep_tree.tre -o alpha_diversity -O 8

### Collate information from alpha diversity output into one file
# -i directory with alpha diversity files
# -o output directory

collate_alpha.py -i alpha_diversity -o collate_alpha

add_alpha_to_mapping_file.py -m norman+BFTM_samps.txt -i collate_alpha/observed_otus.txt,collate_alpha/PD_whole_tree.txt --depth 9000 --collated_input -o alpha+meta.txt


### Make rarefaction plots with alpha diversity metrics
# -i collate alpha folder
# -o output folder
# -m metadata file


make_rarefaction_plots.py -i collate_alpha -m metadata.txt -o rareplots

### Generate beta diversity metrics (weighted and unweighted UniFrac distances) that will be used to create PCoA plots. The depth of sequences used depends on the rarefaction plots.
# -i input rarefaction biom file of desired depth
# -t representative tree
# -o output folder

beta_diversity.py -i rarefactions/rarefaction_10000_0.biom -t rep_tree.tre -o beta_diversity

### Generate PCoA scores. Can use either weighted unifrac or unweighted unifrac measurements, or both
# -i weighted or unweighted unifrac text file
# -o output matrix of PCoA scores

principal_coordinates.py -i weighted_unifrac_rarefaction_10000_0.txt -o weighted.pc

### Generate 3D beta diversity plot with unifrac PCoA scores
# -i output of principal_coordinates command
# -m metadata file in qiime format
# -o output folder
# --ignore_missing_samples is required if there are samples in the metadata file that are not in biom file, or vice versa

make_emperor.py -i weighted.pc -m metadata.txt -o weighted_3D_plot --ignore_missing_samples

### Generate a summary of taxa abundance for each sample
# -i rarefaction biom file of desired depth
# -a return numeric values, rather than proportion
# -o output directory

summarize_taxa.py -i jaco3003/lmamr20_goldberg_analysis/rarefactions/rarefaction_10000_0.biom -a -L 1,2,3,5,6,7 -o temp_taxa_summarized




