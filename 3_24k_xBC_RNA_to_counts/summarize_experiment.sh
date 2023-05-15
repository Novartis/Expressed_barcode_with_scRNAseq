#!/usr/bin/env bash

# Run this in the project level output directory

# Get 10x TP barcodes from cell ranger assay. Code may vary depending on 10x run
#input_folder: Cell ranger assay. Eg: /dlab/NGS/ONC/AAC142/aggr_141_142_mapped/outs

#python $TILING_HOME/get_10x_barcodes_from_CR.py <input_folder> > CR_10x_barcodes.csv

# Append cDNA and lib to the sample name
#head -n1 CR_10x_barcodes.csv > CR_10x_barcodes_cDNA_lib.csv
#awk -F"," '(NR>1){print $1"_cDNA,"$2}' CR_10x_barcodes.csv >> CR_10x_barcodes_cDNA_lib.csv
#awk -F"," '(NR>1){print $1"_library,"$2}' CR_10x_barcodes.csv >> CR_10x_barcodes_cDNA_lib.csv

# Let us now assume that the CR_10x_barcodes_cDNA_lib.csv file exist

# Summary stats
(echo -e "sample\tcategory\tcount" && for file in `find -type f -name "*stats.txt" | grep -v half`; do sample=`echo $file | cut -f4 -d"/" | cut -f1 -d"." | sed 's/_stats//g'`; awk -v sample=$sample '($0 ~ /^#/){split($0,arr,":"); print sample"\t"substr(arr[1],2,100)"\t"arr[2]}' $file; done;) > stats_summary.txt

# tenx counts file:
(head -n1 `ls */*tenx_associations.txt | head -n1` && ls */*tenx_associations.txt | xargs -I file sed 1d file) > tenx_counts.txt

# cw file
for folder in `ls -d * | grep "-"`; do cd $folder; python /da/onc/krishvi7/bitbucket/oncp-tiling/compute_cw.py *tenx_associations.txt *_stats.txt > tenx_exp_cw.txt& cd ..; done;
(echo -e "sample\tpercent_barcodes_read\tunique_barcode_reqd"  && for folder in `ls | grep "-"`; do sample=`echo $folder | cut -f2 -d"/"`; awk -v sample=$sample '(NR > 1){print $0}' $folder/tenx_exp_cw.txt ; done;) > tenx_cw.txt
(echo -e "sample\tpercent_barcodes_read\tunique_barcode_reqd\tcellline\tday\tsource" && awk '(NR>1){split($1,arr,"_"); print $0"\t"arr[1]"\t"arr[2]"\t"arr[3]}' tenx_cw.txt) > tenx_cw_annt.txt

# tenx_recounts and co-occurence
# recounts merges clonal barcodes that look similar to one other and the 1st one is seen at a much higher frequency than the second
# co-occurence checks if a 10x barcode is associated with multiple clones. If so, it reports both a matrix per sample and tall skilly summary
python /da/onc/krishvi7/bitbucket/oncp-expressed/barcode_cleanup/recluster_barcodes.py tenx_counts.txt CR_10x_barcodes_cDNA_lib.csv > tenx_recounts.txt

# Recount recounts
# Uses TP file and expression distances between 10x barcodes to choose TP sets from tenx recounts. Write tenx_recounts2 file
python /da/onc/krishvi7/bitbucket/oncp-expressed/barcode_cleanup/recount_recounts.py tenx_recounts.txt  CR_10x_barcodes_cDNA_lib.csv

# Sorts recount 2 file
(head -n1 tenx_recounts2.txt && sed 1d tenx_recounts2.txt | sort -k1,2 -k4nr) > tenx_recounts2_sorted.txt

# Created summaries for xBC<->10x mapping
python /da/onc/krishvi7/bitbucket/oncp-expressed/summarize_tenx_xbc_counts.py tenx_counts.txt tenx_recounts2_sorted.txt ../CR_10x_barcodes_cDNA_lib.csv


# R commands to get matrix of recounts2_sorted and merging with DNA data
# library(data.table)
# library(dplyr)

# tenx_recounts<-fread("tenx_recounts2_sorted.txt")
# recounts_mat<-dcast(tenx_recounts, exp_bc ~ sample, value.var="cpm", fun.aggregate=sum, fill=0.0)
# write.table(recounts_mat, "recounts_mat.txt",row.names=F,quote=F,sep="\t")

# recounts_mat_subset <- recounts_mat %>% select(exp_bc, HCC827_0, HCC827_14)
# names(recounts_mat_subset) <- c("exp_bc", "Day0_RNA", "Treated_RNA")

# dna_counts<-fread("../../20190923_HCC827_drug_DNA/output/counts_annt_Day0_EGF_single.txt")
# dna_counts_mat<-dcast(dna_counts, bc~sample, value.var="cpm", fun.aggregate=mean)
# names(dna_counts_mat) <- c("exp_bc", "Treated_DNA", "Day0_DNA")

# rna_dna_HCC827<- merge(dna_counts_mat, recounts_mat_subset, by="exp_bc")
# write.table(rna_dna_HCC827,"dna_rna_comparison_v2.txt", sep="\t", quote=F, row.names=F)