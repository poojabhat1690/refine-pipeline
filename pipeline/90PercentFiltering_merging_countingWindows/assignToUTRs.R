









## this script creates a file for the accepted refSeq and ensembl 3'UTR overlapping ends. 


library(checkmate)
library(GenomicRanges)
### the input should be info about the peaks that overlap with refSeq and ensembl 

refSeq_total = read.table(paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/refSeq_total.bed"),sep="\t",stringsAsFactors = F)

#### creating an addit

ensembl_total = read.table(paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/ensembl_total.bed"),sep="\t",stringsAsFactors = F)
ensembl_intronOverlapping =  read.table(paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/intron_total.bed"),sep="\t",stringsAsFactors = F)
ensembl_intronOverlapping$V14 = ensembl_intronOverlapping$V17 
ensembl_exonOverlapping =  read.table(paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/exon_total.bed"),sep="\t",stringsAsFactors = F)
ensembl_exonOverlapping$V14 = ensembl_exonOverlapping$V17

assertDataFrame(x = refSeq_total,ncols = 21)
assertDataFrame(x = ensembl_total,ncols = 21)
assertDataFrame(x = ensembl_intronOverlapping,ncols = 21)
assertDataFrame(x = ensembl_exonOverlapping,ncols = 21)

refSeq_total$overlap = "refSeq"
ensembl_total$overlap = "ensembl"
ensembl_intronOverlapping$overlap = "ensembl_intron"
ensembl_exonOverlapping$overlap = "ensembl_exon"

refSeq_pas = refSeq_total[-grep("noPAS",refSeq_total$V21),]
refSeq_noPas = refSeq_total[grep("noPAS",refSeq_total$V21),]


ensembl_pas = ensembl_total[-grep("noPAS",ensembl_total$V21),]
ensembl_noPas = ensembl_total[grep("noPAS",ensembl_total$V21),]

ensembl_intronOverlapping_pas = ensembl_intronOverlapping[-grep("noPAS",ensembl_intronOverlapping$V21),]
ensembl_intronOverlapping_noPas = ensembl_intronOverlapping[grep("noPAS",ensembl_intronOverlapping$V21),]

ensembl_exonOverlapping_pas = ensembl_exonOverlapping[-grep("noPAS",ensembl_exonOverlapping$V21),]
ensembl_exonOverlapping_noPas = ensembl_exonOverlapping[grep("noPAS",ensembl_exonOverlapping$V21),]


total_pas = rbind(refSeq_pas,ensembl_pas,ensembl_intronOverlapping_pas,ensembl_exonOverlapping_pas)

total_nopas = rbind(refSeq_noPas, ensembl_noPas,ensembl_intronOverlapping_noPas,ensembl_exonOverlapping_noPas)

## filtering based on the A threshold . here V10 is contains the A thresholds


total_pas_accepted =  total_pas[which(total_pas$V10<0.36),]
total_nopas_accepted = total_nopas[which(total_nopas$V10<0.24),]


total_pas_notAccetped = total_pas[which(total_pas$V10>=0.36),]
total_noPas_notAccepted =   total_nopas[which(total_nopas$V10>=0.24),]

allRejected = rbind(total_noPas_notAccepted, total_pas_notAccetped)
allRejected_granges = makeGRangesFromDataFrame(allRejected,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T )


allPeaks = rbind(total_pas_accepted, total_nopas_accepted)


peaks_accepted_nonAccepted = rbind(total_pas, total_nopas)

acceptedPeaks_sum = sum(as.numeric(unlist(strsplit(allPeaks$V5,split = ",",fixed = T))))
peaks_accepted_nonAccepted_sum = sum(as.numeric(unlist(strsplit(peaks_accepted_nonAccepted$V5,split = ",",fixed = T))))


cat(paste0("The fraction of counts accepted:",acceptedPeaks_sum/peaks_accepted_nonAccepted_sum))


## writing these tables 

allPeaks$V4 = unlist(lapply(strsplit(allPeaks$V14,",",T),function(x) x[[1]]))
write.table(allPeaks,paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/ends_all_10threshold_n100.txt"),sep="\t",quote = F,row.names = F)

### writng a part of this to a bed file : 

ensemblIdsOfMissing = allPeaks[which(allPeaks$V4 == ""),]$V17
allPeaks$V4[which(allPeaks$V4 == "")] <-ensemblIdsOfMissing

allPeaks_short = cbind.data.frame(allPeaks[,c(1:6)],allPeaks$overlap)
write.table(allPeaks_short,paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/ends_all_10threshold_n100.bed"),sep = "\t",quote = F,row.names = F,col.names = F)



