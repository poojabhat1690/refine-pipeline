#!/bin/bash







#### this is a script that can be used to obtain polyA reads from multiple fastQ files. 
### a parameter file is required that has information about the file and one column with the index with the file 


# 8. Parse parameter file to get variables.
#module load bedtools
#module load  bam2fastq


### reading in the variables passed in the command line... 

while getopts 'a: i: o:' OPTION; do
	  case "$OPTION" in
		      a)
			      avalue="$OPTARG"
			  echo "the adapter trimmed has the sequence $OPTARG"
				          ;;

					      i)
						      ivalue="$OPTARG"
						        echo "the input directory containing QuantSeq data is $OPTARG"
							          ;;
										
							o) 
								ovalue="$OPTARG"
								echo "the output dorectory is $OPTARG"
									;;

					  			  ?)
										 echo "script usage: $(basename $0) [-a adapter] [-i input directory]" >&2
														        exit 1
															      ;;
															        esac
															done
															shift "$(($OPTIND -1))"



#### i value is the input directory ... this should consist of two directories... 

	### quantseq - containing all the zipped files of the quantSeq data
	### rnaseq - containing all the rnaseq data for the condition, mapeed and indexed. 

QUANT_ALIGN=$ovalue/polyAmapping_allTimepoints
mkdir -p "$QUANT_ALIGN"

ls "$ivalue"/quantseq/*.{gz,fastq,fq} | perl -pe "s#$ivalue/quantseq/##" > $QUANT_ALIGN/sampleInfo.txt


INPUT="$ivalue"/quantseq/

PARAMETER="$QUANT_ALIGN/sampleInfo.txt"

#set dirs

OUTDIR=$QUANT_ALIGN
LOG=$QUANT_ALIGN/logs/
mkdir -p $LOG

#get input file (index)

#COUNTER=0
#for index in `cat $PARAMETER`; do #parameter file = list of files
 #   COUNTER=$((COUNTER+1))
 #	  if [ $COUNTER = $SGE_TASK_ID ]; then
#	break
 #  fi
#done

#init files
touch "$OUTDIR"/logFile.txt ## initialise a log file
touch "$OUTDIR"/readStats.txt ## initialise a file to store your read stats at different processing steps
touch "$OUTDIR"/processingSteps.txt ### just another file to list the processing steps 

touch "$LOG"/stderr_"$index".txt
touch "$LOG"/stdo_"$index".txt


date >"$OUTDIR"/logFile.txt
 date >>"$LOG"/stderr_"$index".txt
 date >>"$LOG"/stdo_"$index".txt

#################
#ADAPTER trimming
#################


module load cutadapt/1.9.1-foss-2017a-python-2.7.13
module load fastx-toolkit/0.0.14-foss-2017a
module load fastqc/0.11.5-java-1.8.0_121
cutadapt_version=$(cutadapt --version)
echo running cutadapt version "$cutadapt_version" >> "$LOG"/stdo_"$index".txt

## the step below trims the adapter : you should replace this based on your library prep method

while read index; do
	 echo "processing $index"
  
	  ### trimming adapter and adding the relavant information. 

  	cutadapt -a $avalue -o "$OUTDIR"/"$index"_trimmed.fastq  --trim-n $INPUT/"$index" 2>"$LOG"/stderr_"$index".txt 1>>"$LOG"/stdo_"$index".txt
  	adapterTrimming=$(cat  "$OUTDIR"/"$index"_trimmed.fastq | echo $((`wc -l`/4)))
  	echo the number of reads after adapter trimming is "$adapterTrimming" >>"$LOG"/stdo_"$index".txt
	

	### QC  
	
	module load fastqc/0.11.5-java-1.8.0_121
	fastqc "$OUTDIR"/"$index"_trimmed.fastq 
	
	### trimming from the 5' end....

	echo fastx_trimmer to remove 12 nts from the 5 end of the reads >>"$LOG"/stdo_"$index".txt  
	echo version of fastx_trimmer used is Version: `which fastx_trimmer` >>"$LOG"/stdo_"$index".txt
	fastx_trimmer -Q33 -f 13 -m 1 -i "$OUTDIR"/"$index"_trimmed.fastq  > "$OUTDIR"/"$index"_5primetrimmed_trimmed.fastq 2>"$LOG"/stderr_"$index".txt 
	fivePrimeTrimmed=$(cat  "$OUTDIR"/"$index"_5primetrimmed_trimmed.fastq | echo $((`wc -l`/4)))                                                                                    echo number of reads after 5 trimmig "$fivePrimeTrimmed" >>"$LOG"/stdo_"$index".txt
	echo "retaining reads that have >=5 As the the 3 end" >>"$LOG"/stdo_"$index".txt


	#### removing A's from the 3' end

	#AAAAA$ could theoretically also happen in PHRED score
	egrep -A2 -B1 'AAAAA$' "$OUTDIR"/"$index"_5primetrimmed_trimmed.fastq  | sed '/^--$/d' > "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads.fastq  2>"$LOG"/stderr_"$index".txt
	polyAreads=$(cat "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads.fastq | echo $((`wc -l`/4)))
	echo readd withs polyA "$polyAreads"  >>"$LOG"/stdo_"$index".txt
	echo remove the polyAs at the end of polyA reads >>"$LOG"/stdo_"$index".txt
	
	#echo this is done using cutadapat. At this step we also perform a size filteirng of minimum 18 nucleotides >>"$LOG"/stdo_"$index".txt
	#cut super long polyA to avoid internal polyA cutting of 5 or more As
	cutadapt --no-indels -m 18 -e 0 -a "A{1000}"  -o "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq  "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads.fastq 1>>"$LOG"/stdo_"$index".txt 2>"$LOG"/stderr_"$index".txt
	finalReads=$(cat "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq | echo $((`wc -l`/4)))
	echo final reads after size filtering "$finalReads"  >>"$LOG"/stdo_"$index".txt

	################
	# write stats
	#################

	touch "$LOG"/preProcessingNumbers_"$index".txt

	echo initialFile:"$initialrReads" >>"$LOG"/preProcessingNumbers_"$index".txt
	echo adapterTrimmed:"$adapterTrimming" >>"$LOG"/preProcessingNumbers_"$index".txt
	echo fivePrimeTrimming:"$fivePrimeTrimmed" >>"$LOG"/preProcessingNumbers_"$index".txt
	echo polyAcontaining:"$polyAreads" >>"$LOG"/preProcessingNumbers_"$index".txt
	echo finalFile:"$finalReads" >>"$LOG"/preProcessingNumbers_"$index".txt
	echo the pre processing before mapping has been completed.  >>"$LOG"/stdo_"$index".txt


	echo initialFile:"$initialrReads"  >>"$LOG"/stdo_"$index".txt
	echo adapterTrimmed:"$adapterTrimming"  >>"$LOG"/stdo_"$index".txt
	echo fivePrimeTrimming:"$fivePrimeTrimmed"  >>"$LOG"/stdo_"$index".txt
	echo polyAcontaining:"$polyAreads"  >>"$LOG"/stdo_"$index".txt
	echo finalFile:"$finalReads"  >>"$LOG"/stdo_"$index".txt


	#################
	# getting the read length distribution 
	#################

	awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq >"$OUTDIR"/"$index"_lengthDistribution.txt


  done <"$QUANT_ALIGN/sampleInfo.txt"






































































