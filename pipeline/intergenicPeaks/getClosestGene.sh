#!/bin/bash

module load bedtools/2.27.1-foss-2017a


##### sort the bed files 

## this is the bed file of all ensembl and refSeq transcripts. 

sort -k1,1 -k2,2n $QUANT_INTERGENIC/allExons_refSeq_ensembl.bed > $QUANT_INTERGENIC//allExons_refSeq_ensembl_sorted.bed

### this is the bed file of the most distal 3' position per gene. 

sort -k1,1 -k2,2n $QUANT_INTERGENIC/toExtend_longestEnsembl_refSeq_n100.bed > $QUANT_INTERGENIC/toExtend_longestEnsembl_refSeq_n100_sorted.bed


##### we want to calculate the distance between the most distal 3' end per gene and the next annotation (refSeq + ensembl), to prevent considering RNAseq singal coming from another annotation. 

bedtools closest -d -s -io -iu -D a -a $QUANT_INTERGENIC/toExtend_longestEnsembl_refSeq_n100_sorted.bed -b $QUANT_INTERGENIC/allExons_refSeq_ensembl_sorted.bed > $QUANT_INTERGENIC/toExtend_longestEnsembl_refSeq_n100_sorted_distances.bed



