#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=100G






ml   python/2.7.13-foss-2017a
pip install --user slamdunk



while getopts 'a: i: o: g:' OPTION; do
	          case "$OPTION" in
			    a) avalue="$OPTARG"
				echo "the adapter trimmed has the sequence $OPTARG"
				;;

			 	i) ivalue="$OPTARG"
				   echo "the input directory containing QuantSeq data is $OPTARG"
				;;

				o)  ovalue="$OPTARG"
				    echo "the output dorectory is $OPTARG"
				;;		
				
				g) genome="$OPTARG"
				   echo " the genome: $OPTARG"
				;;

				?)	echo "script usage: $(basename $0) [-a adapter] [-i input directory] [-o output directory] [-g genome]" >&2
					exit 1
				;;
		 esac
	done
shift "$(($OPTIND -1))"







OUTDIR=$ovalue/n_100_global_a0/
mkdir -p $OUTDIR


QUANT_ALIGN=$ovalue/polyAmapping_allTimepoints
PARAMETER="$QUANT_ALIGN/sampleInfo.txt"



while read index; do
	         echo "processing $index"
		      slamdunk map -r $genome -o $OUTDIR/ -n 100 -5 0 -a 0 -e $QUANT_ALIGN/$index_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq
		      slamdunk filter -o $OUTDIR -mq 0 -mi 0.95  $OUTDIR/$index_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved_slamdunk_mapped.bam         
 done <"$QUANT_ALIGN/sampleInfo.txt"





