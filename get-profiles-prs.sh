#!/bin/bash
####script to generate profiles in plink
################
###set parametres###
################
threshold=(0.00000005 0.00005 0.05 0.01 0.1 0.5) #p-value thresholds 
inputfile= gwasfile
tempfile= tempgwasfile
prunedfile=gwasfile-ld0.8-maf0.5
gwasstat=sumstatsfile
################
###prune and QC SNPs###
################
cov=clinical-data-file.txt
plink1.9 --bfile $inputfile --keep  $cov  --maf 0.01 --hwe 0.000001 --indep-pairwise 50 5 0.8 -make-bed --out $tempfile ##LD 0.8, cab ne chancges 
plink1.9 --bfile $inputfile --keep  $cov   --extract $tempfile.prune.in --make-bed --out $prunedfile
################
###get thresholds and plink profiles###
################
for i in ${!threshold[*]}
do
    awk '{ if ($10 <= '${threshold[$i]}') print $1,$5,$9}' $gwasstat.txt > $gwasstat-A2-thres${threshold[$i]} ###get all SNPs associated at selected level
    plink1.9 --bfile $prunedfile  --score $gwasstat-A2--thres${threshold[$i]} --out IRL-profileA2-$gwasstat-thres${threshold[$i]}
done
#####################
####clean a space####
#####################
rm *.*.nopred
rm $inputfile.bed
rm $inputfile.bim
rm $inputfile.fam
rm $tempfile.prune.out
rm $tempfile.prune.in
rm *.log

