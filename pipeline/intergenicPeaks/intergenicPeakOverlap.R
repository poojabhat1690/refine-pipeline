

extendedFiles = list.files ("//clustertmp/bioinfo/pooja/SLAMannotation/dr_allData/output/parallel/extensions/")
extendedFilesPath = paste0("//clustertmp/bioinfo/pooja/SLAMannotation/dr_allData/output/parallel/extensions/",extendedFiles)
fileWithLines = c()
for (file in list.files ("//clustertmp/bioinfo/pooja/SLAMannotation/dr_allData/output/parallel/extensions/"))
  {
  file=paste0("//clustertmp/bioinfo/pooja/SLAMannotation/dr_allData/output/parallel/extensions/",file)
  if (file.size(file) == 0) next
  fileWithLines = c(fileWithLines,file)
}


interGenic_thresholds = lapply(fileWithLines,function(x) read.table(x,stringsAsFactors = F))

library(reshape)
numberPerBin = melt(lapply(interGenic_thresholds, function(x) { if (is.null(x)) { return(0) } else {return(nrow(x))}}))
library(ggplot2)

#### NEED TO CHANGE FOLDER PATH

pdf(paste0(BOut, "/coverage/numberOfBins_rnaseq_new.pdf"),height=5,width=10)
ggplot(numberPerBin,aes(y=value,x=c(1:nrow(numberPerBin)))) + geom_point() + theme_bw() + xlab("Bin number") + ylab("Number of peaks with RNAseq signal")  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#### NEED TO CHANGE FOLDER PATH


#pdf("/Users/pooja.bhat/Dropbox/UTRannotation/mESC/intergenicPeaks/plots/numberOfBins_rnaseq_new.pdf",height=5,width=10)
#ggplot(numberPerBin,aes(y=value,x=c(1:nrow(numberPerBin)))) + geom_point() + theme_bw() + xlab("Bin number") + ylab("Number of peaks with RNAseq signal")  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dev.off()

#### NEED TO CHANGE FOLDER PATH
write.table(numberPerBin, paste0(BOut, "/coverage/numberOfRNAseqContaingBins.txt"),sep="\t")

# for(i in 1:length(interGenic_thresholds)){
#   
#   interGenic_thresholds[[i]] = interGenic_thresholds[[i]][which(interGenic_thresholds[[i]]$V3 - interGenic_thresholds[[i]]$V2 == 200),]
# }


totalBins = do.call(rbind, interGenic_thresholds)
library(GenomicRanges)
### reading in the nonOverlapping peaaks we defined by hierarchical overlapping of ends with refSeq UTRs (protein coding), ensembl UTRs (protein coding), exons (refSeq protein coding), introns (refSeq protein coding)

nonOverlappning = read.table(paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/nonOverlapping_total.bed"),sep="\t",stringsAsFactors = F)

#START PLUS ONE??????
nonOverlappingPas = nonOverlappning[-grep("noPAS",nonOverlappning$V12),]
nonOverlappingNoPas = nonOverlappning[grep("noPAS",nonOverlappning$V12),]

### filtering based on the A threshold we defined to distinguish true priming event from internal priming events. 

nonOverlappingNoPas = nonOverlappingNoPas[which(nonOverlappingNoPas$V10<0.24),]
nonOverlappingPas = nonOverlappingPas[which(nonOverlappingPas$V10<0.36),]

nonOverlapping = rbind(nonOverlappingPas, nonOverlappingNoPas)

nonOverlappingPeaks_plus = nonOverlapping[which(nonOverlapping$V6 == "+" ),]
nonOverlappingPeaks_minus = nonOverlapping[which(nonOverlapping$V6 == "-" ),]


creatingGranges  = function(counts_offSet100_10PercentThreshold,nonOverlappingPeaks_plus,nonOverlappingPeaks_minus){
  
  if(is.null(counts_offSet100_10PercentThreshold))
  {
    return(NULL)
  }
  
  # nonOverlappingPeaks_plus = nonOverlapping[which(nonOverlapping$V6 == "+" ),]
  # nonOverlappingPeaks_minus = nonOverlapping[which(nonOverlapping$V6 == "-" ),]
  nonOverlappingPeaks_plus_Granges = with(nonOverlappingPeaks_plus,GRanges(V1,IRanges(start=V2,end=V3),strand = V6))
  nonOverlappingPeaks_minus_Granges = with(nonOverlappingPeaks_minus,GRanges(V1,IRanges(start=V2,end=V3),strand = V6))
  
  start(nonOverlappingPeaks_plus_Granges) <- start(nonOverlappingPeaks_plus_Granges)+1
  start(nonOverlappingPeaks_minus_Granges) <- start(nonOverlappingPeaks_minus_Granges)+1
  
  counts_offSet100_10PercentThreshold_positive = counts_offSet100_10PercentThreshold[which(counts_offSet100_10PercentThreshold$V6=="+"),]
  counts_offSet100_10PercentThreshold_negative = counts_offSet100_10PercentThreshold[which(counts_offSet100_10PercentThreshold$V6=="-"),]
  
  counts_offSet100_10PercentThreshold_positive_granges  = with(counts_offSet100_10PercentThreshold_positive,GRanges(V1,IRanges(start=V2,end=V3),strand = V6))
  counts_offSet100_10PercentThreshold_negative_granges =  with(counts_offSet100_10PercentThreshold_negative,GRanges(V1,IRanges(start=V2,end=V3),strand = V6))
  
  start(counts_offSet100_10PercentThreshold_positive_granges) <- start(counts_offSet100_10PercentThreshold_positive_granges)+1
  start(counts_offSet100_10PercentThreshold_negative_granges) <- start(counts_offSet100_10PercentThreshold_negative_granges)+1
  
  
  overlaps_positive = findOverlaps(nonOverlappingPeaks_plus_Granges,counts_offSet100_10PercentThreshold_positive_granges)
  overlaps_negative = findOverlaps(nonOverlappingPeaks_minus_Granges,counts_offSet100_10PercentThreshold_negative_granges)
  
  
  querySubject_positive =cbind(  nonOverlappingPeaks_plus[queryHits(overlaps_positive),],counts_offSet100_10PercentThreshold_positive[subjectHits(overlaps_positive),])
  querySubject_negative =cbind(  nonOverlappingPeaks_minus[queryHits(overlaps_negative),],counts_offSet100_10PercentThreshold_negative[subjectHits(overlaps_negative),])
  
  subjectQuery = rbind(querySubject_positive,querySubject_negative)
  
  return(subjectQuery)
}

### overlap the offset counting windows with non-overlapping peaks. 

overlapPeaks = lapply(interGenic_thresholds,function(x) creatingGranges(counts_offSet100_10PercentThreshold = x,nonOverlappingPeaks_plus = nonOverlappingPeaks_plus,nonOverlappingPeaks_minus = nonOverlappingPeaks_minus))

numberOverlappingPeaks = melt(lapply(overlapPeaks, function(x) { if (is.null(x)) { return(0) } else {return(nrow(x))}}))
pdf(paste0(BOut, "/coverage/numberOfPeaks_polyAsupport_new.pdf"),height = 5,width = 10)
ggplot(numberOverlappingPeaks,aes(x=c(1:nrow(numberOverlappingPeaks)),y=value)) + geom_point() + theme_bw() + xlab("Bin number") + ylab("Number of regions with polyA support")  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

overlapPeaks = do.call(rbind,overlapPeaks)
overlapPeaks = overlapPeaks[!duplicated(overlapPeaks[,c(1:3)]),]

colnames(overlapPeaks) = c("chr",	"start","end"	,"peakName"	,"numberOfpolyAends"	,"strand",	"120ntSeq"	,"internalId(chr_start)",	"20nts (+1 to +20)",	"Acontent"	,"upStreamSequence(-40 to -5)",	"Acontent_motif",	"chr",	"RNAseqbinStart",	"RNAseqbinEnd",	"geneName","score",	"strand","refSeq/ensemblid","closestDistancetonextexon",	"mean(RNAseq)",	"nBIN","RNAseq(bin)/rnaSeq(inUTR)")										

######### WRINTING OUT THE TABLES####################


write.table(overlapPeaks,paste0(BOut, "/coverage/allIntergenicPeaks_n100_new.txt"),sep="\t",quote = F, row.names = F,col.names = F)

write.table(numberOverlappingPeaks,paste0(BOut, "/coverage/numberOfOverlappingPeaksPerBin.txt"),sep="\t",quote = F, row.names = F,col.names = F)




colnames(overlapPeaks) = paste("V",1:ncol(overlapPeaks),sep="")


overlapPeaks$V4 = overlapPeaks$V16

overlapPeaks = overlapPeaks[,c(1:6)]

overlapPeaks_plus = overlapPeaks[which(overlapPeaks$V6 == "+"),]
overlapPeaks_minus = overlapPeaks[which(overlapPeaks$V6 == "-"),]


write.table(overlapPeaks,paste0(BOut, "/coverage/overlappingPeaks_threshold_n100.bed"),sep="\t", quote = F, row.names = F,col.names = F)

write.table(overlapPeaks_plus,paste0(BOut, "/coverage/overlappingPeaks_threshold_n100_plus.bed"),sep="\t", quote = F, row.names = F,col.names = F)
write.table(overlapPeaks_minus,paste0(BOut, "/coverage/overlappingPeaks_threshold_n100_minus.bed"),sep="\t", quote = F, row.names = F,col.names = F)


