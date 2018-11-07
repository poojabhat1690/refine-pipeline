#!/bin/bash
#$ -S /bin/bash
#$ -q public.q
#$ -cwd
# -pe smp 6-24
#$ -l hostname=compute-*
#$ -v GENOME
#$ -v QUANT_ALIGN
#$ -v PIPELINE
#$ -v LOG
#### after mapping 

PARAMETER="$QUANT_ALIGN/sampleInfo.txt"

INPUT=$QUANT_ALIGN/n_100_global_a0/
OUTDIR=$INPUT
DIR=$OUTDIR

touch $OUTDIR/log_afterMapping.txt

####################
# counting positions where reads end :
####################

module load bedtools 
module load R/3.2.2 #ml R/3.3.0 #??

echo "loading bedtools" >>$OUTDIRlog_afterMapping.txt

bedtools_version=$(bedtools --version)
echo the version of bedtools being used is "$bedtools_version" >>"$DIR"/log_afterMapping.txt

for BAM in $DIR/*_slamdunk_mapped_filtered.bam; do

    echo Processing $BAM >>"$DIR"/log_afterMapping.txt
    echo "the first step is to convert the bam to bed to find the end of reads this is done using bamtobed from bedtools" >>"$DIR"/log_afterMapping.txt

    bedtools bamtobed -i  $BAM > ${BAM}_bamTobed.bed 

    bamTobedlines=$(wc -l ${BAM}_bamTobed.bed )

    echo "the number of reads that are in the bed file- as a result of the reads mapping in the bam file are" "$bamTobedlines"  >>"$DIR"/log_afterMapping.txt

    awk -vFS="\t" '$6 == "-"' ${BAM}_bamTobed.bed > ${BAM}_bamTobed_minusStrand.bed
    awk -vFS="\t" '$6 == "+"' ${BAM}_bamTobed.bed > ${BAM}_bamTobed_plusStrand.bed
done

####################
# minus and plus strand . For the minus strand I will consider the start and chr as the id and for the plus strand I will consider the end and the chr as the id.  
####################

cat "$DIR"/*_bamTobed_minusStrand.bed | awk -vFS="\t" '{print $1,$2}'  | sort | uniq -c | perl -pe 's#^\s+##'  > "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed
cat "$DIR"/*_bamTobed_plusStrand.bed | awk  -vFS="\t" '{print $1,$3}' | sort | uniq -c | perl -pe 's#^\s+##' > "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed

onPositive=$(wc -l  "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed)
onNegative=$(wc -l  "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed)

echo the numeber of unique positions of the genome that are covered by reads on the plus strand is "$onPositive" and on the minus strand is "$onNegative" >>"$DIR"/log_afterMapping.txt

###############
# Filter above THRESHOLD
###############


###### this is the thresholding based on the number of mapped instances

#MINUSSUM=`awk '{ sum += $1 } END { print sum } '  "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed`
#PLUSSUM=`awk '{ sum += $1 } END { print sum } '  "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed`
#TOTALSUM=`expr $PLUSSUM + $MINUSSUM`
#THRESHOLD_1=`expr  $"$TOTALSUM" / 60000000`
#THRESHOLD=`expr $"$TOTALSUM" / 1000000`

######## this is the thresholding based on the number of polyA-size filtered reads 
#ml samtools
#FILE=`ls $DIR/*_slamdunk_mapped_filtered.bam`
#a=`samtools view -H "$FILE" | grep @RG`
#TOTALSUM=`echo "$a" | cut -f2 -d { | cut -f2 -d : | cut -f1 -d ,`
#THRESHOLD=`expr $"$TOTALSUM" / 1000000`
### the value of threshold above takes into account the number of reads that are mapped (including multimapping) 

######## creating a threshold based on the number of priming sites


#cat "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed > "$DIR"/totalPrimingSites.bed
#a=`cut -f 1 -d ' ' "$DIR"/totalPrimingSites.bed`
#b=`echo "$a" | sort -n | uniq -c`
#cumulative=`echo "$b" | awk '{total += $0; $0 = total}1' -`
#totalEntries=`echo "$a" | wc -l` 
#cumFrac=`echo "$cumulative" | awk '{$2=$1/a;print}' a="$totalEntries" - | cut -f2 -d ' '`
#lessThan=`echo "$cumFrac" | awk '$1<=0.85' -`
#nlines=`echo "$lessThan" | wc -l`
#THRESHOLD=`expr $nlines + 1`


###### creating a threshold based on TPM.... 


#awk '{print $1, $2, $3, "+"}' "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed > "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_includingStrand.bed
#awk '{print $1, $2, $3, "-"}' "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed > "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_includingStrand.bed

#cat  "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_includingStrand.bed "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_includingStrand.bed > "$DIR"/allSites_countsUnique.bed

#awk  '{ print $1, $2, $3, $4, $1*1000 }' "$DIR"/allSites_countsUnique.bed > "$DIR"/allSites_countsUnique_RPK.bed
#scaleFac=`awk '{sum+=$5} END {print sum/1000000}' "$DIR"/allSites_countsUnique_RPK.bed`

#awk '{print $1, $2, $3, $4, $5, $5/ a}' a="$scaleFac" "$DIR"/allSites_countsUnique_RPK.bed | awk '$6>=1' - > "$DIR"/allSites_TPMfiltered.bed

#THRESHOLD=1
#awk -vFS=" " '$4 == "-"' "$DIR"/allSites_TPMfiltered.bed | awk '{print $1, $2, $3}' -  > "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$THRESHOLD".bed
#awk -vFS=" " '$4 == "+"' "$DIR"/allSites_TPMfiltered.bed | awk '{print $1, $2, $3}' -  > "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$THRESHOLD".bed






THRESHOLD=2

awk -vT=$THRESHOLD '{ if ($1 >=T) print  }' "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed >"$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$THRESHOLD".bed
awk -vT=$THRESHOLD '{ if ($1 >=T) print  }' "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed >"$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$THRESHOLD".bed

echo "taking the positions that have greater than 10 reads ending at the position"  >>"$DIR"/log_afterMapping.txt
positive_10=$(wc -l "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$THRESHOLD".bed)
negative_10=$(wc -l "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$THRESHOLD".bed)

echo "the number of positions that occur that have greater than 10 reads ending at the particular position are" "$positive_10" on the plus strand and "$negative_10" on the minus strand >>"$DIR"/log_afterMapping.txt



### writing a bed file : 

awk -vOFS="\t" '{print $2, $3-1, $3, $1}' $DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$THRESHOLD".bed" | sort -k1,1 -k2,2n  > "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$THRESHOLD".bed_sorted.bed_changedCoordinates.bed

awk -vOFS="\t" '{print $2, $3, $3+1, $1}' $DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$THRESHOLD".bed" | sort -k1,1 -k2,2n  > "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$THRESHOLD".bed_sorted.bed_changedCoordinates.bed





bedtools merge -d 0 -i "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$THRESHOLD".bed_sorted.bed_changedCoordinates.bed -c 4 -o count,collapse >   "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$THRESHOLD".bed_sorted_merged.bed

bedtools merge  -d 0  -i "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$THRESHOLD".bed_sorted.bed_changedCoordinates.bed -c 4 -o count,collapse > "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$THRESHOLD".bed_sorted_merged.bed 
 
 
 

#### firstly getting a bed file of 120 nt around the provided utr annotation file :


ml  R/3.3.0

####### then getting the fasta sequence from these coordinates :


ml bedtools 



### not merging this with the polyA peaks 

## first getting the bed file +/- 60nt of the polyA peaks
Rscript  --vanilla $PIPELINE/primingSites/overlappingPolyApeaks.R "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$THRESHOLD".bed_sorted_merged.bed "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$THRESHOLD".bed_sorted_merged.bed "$DIR"/peaks_"$THRESHOLD"_120bps.bed 

bedtools getfasta -s -fi $GENOME -bed "$DIR"/peaks_"$THRESHOLD"_120bps.bed > "$DIR"/peaks_"$THRESHOLD"_120bps.fa
