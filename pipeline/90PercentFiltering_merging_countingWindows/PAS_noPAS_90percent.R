









## this script creates a file for the accepted refSeq and ensembl 3'UTR overlapping ends. 


library(checkmate)
library(GenomicRanges)
### the input should be info about the peaks that overlap with refSeq and ensembl 

refSeq_total = read.table(paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/refSeq_total.bed"),sep="\t",stringsAsFactors = F)

ensembl_total = read.table(paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/ensembl_total.bed"),sep="\t",stringsAsFactors = F)

assertDataFrame(x = refSeq_total,ncols = 21)
assertDataFrame(x = ensembl_total,ncols = 21)

refSeq_total$overlap = "refSeq"
ensembl_total$overlap = "ensembl"

refSeq_pas = refSeq_total[-grep("noPAS",refSeq_total$V21),]
refSeq_noPas = refSeq_total[grep("noPAS",refSeq_total$V21),]


ensembl_pas = ensembl_total[-grep("noPAS",ensembl_total$V21),]
ensembl_noPas = ensembl_total[grep("noPAS",ensembl_total$V21),]

total_pas = rbind(refSeq_pas,ensembl_pas)

total_nopas = rbind(refSeq_noPas, ensembl_noPas)

## filtering based on the A threshold . here V10 is contains the A thresholds


total_pas_accepted =  total_pas[which(total_pas$V10<0.36),]
total_nopas_accepted = total_nopas[which(total_nopas$V10<0.24),]

#total_pas_accepted$type = "PAS"
#total_nopas_accepted$type = "noPAS"

total_pas_notAccetped = total_pas[which(total_pas$V10>=0.36),]
total_noPas_notAccepted =   total_nopas[which(total_nopas$V10>=0.24),]


allRejected = rbind(total_noPas_notAccepted, total_pas_notAccetped)
##allRejected_granges = makeGRangesFromDataFrame(allRejected,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T )


#### I want to now overlap this with ends from 3pSeq to identify possibly missed ends... 

##ulitskyData = read.table("//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/3primeEnds_convertedTobed_dr10.bed",sep="\t",stringsAsFactors = F)

##ulitskyData_ranges = makeGRangesFromDataFrame(df = ulitskyData,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)
# 
# #### adding the data from white et.al 
# whiteData = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/whiteEtal/elife-30860-supp18-v1.tsv",stringsAsFactors = F,sep="\t")
# whiteData = whiteData[whiteData$Chr %in% c(1:19),]



##allRejected_granges = makeGRangesFromDataFrame(allRejected,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T )

##allOverlaps_rejection = findOverlaps(allRejected_granges,ulitskyData_ranges)

##rejectedAndOverlapping = allRejected[queryHits(allOverlaps_rejection),]
##rejectedAndOverlapping = rejectedAndOverlapping[!duplicated(rejectedAndOverlapping),]


allPeaks = rbind(total_pas_accepted, total_nopas_accepted)

### now also including ends that are overlapping compared to the 3pseq data

##allPeaks = rbind(allPeaks,rejectedAndOverlapping)


peaks_accepted_nonAccepted = rbind(total_pas, total_nopas)

acceptedPeaks_sum = sum(as.numeric(unlist(strsplit(allPeaks$V5,split = ",",fixed = T))))
peaks_accepted_nonAccepted_sum = sum(as.numeric(unlist(strsplit(peaks_accepted_nonAccepted$V5,split = ",",fixed = T))))


library(checkmate)

### keep total_pas_accepted and total_nopas_accepted separately
#### write a function to get the 90 % ends


######### 


get90Percent = function(ends_all){


ends_all$V4 = ends_all$V14
ends_all_rearranged = cbind.data.frame(ends_all$V1, ends_all$V2, ends_all$V3, ends_all$V4, ends_all$V5, ends_all$V6, ends_all$overlap, "UTRoverlapping")
colnames(ends_all_rearranged) = paste("V",c(1:ncol(ends_all_rearranged)),sep="")
assertDataFrame(ends_all_rearranged,nrows = nrow(ends_all),ncols = 8)
querySubject = rbind(ends_all_rearranged)

newSums = c()
for(i in 1:nrow(querySubject)){
  takeThis = sum(as.numeric(strsplit(as.character(querySubject[i,]$V5),split = ",",fixed = T)[[1]]))
  
  newSums = c(newSums,takeThis)
}
querySubject$sumCounts = newSums

### splitting data based on the gene name to compare per gene polyA read ends

querySubject_split = split(x = querySubject,f = querySubject$V4,drop = T)

## ordering this based on the number of counts - decreasing number of counts

querySubject_split_order = lapply(querySubject_split,function(x) x[order(x$sumCounts,decreasing=T),])

counts_sum = lapply(querySubject_split_order,function(x) x$sumCounts/sum(x$sumCounts))
a  =Map(cbind,querySubject_split_order,counts_sum)
a = lapply(a,function(x) x[order(x[,10],decreasing=F),])
b = lapply(a,function(x) cumsum(x[,10]))
b = Map(cbind,a,b)
b = lapply(b, function(x) x[which(x[,11]>0.1),])
b  =lapply( b , setNames , nm = c("chr","startPeak","endPeak","external_gene_name","counts","strand","origin","peakKind","sumCounts","fractionCounts","cumSumCounts") )
totalPeaks = do.call(rbind,b)
totalPeaks_bed = totalPeaks[,c(1:6)]
return(totalPeaks)
}



percent90_noPAS = get90Percent(ends_all = total_nopas_accepted)
percent90_PAS = get90Percent(ends_all = total_pas_accepted)


###############ggetting counting windows from these 
library(GenomicRanges)

getCountingWindows = function(annotation_custom){
  colnames(annotation_custom) = paste0("V",c(1:ncol(annotation_custom)))
  annotation_custom_positive = annotation_custom[which(annotation_custom$V6 == "+"),]
  annotation_custom_negative = annotation_custom[which(annotation_custom$V6 == "-"),]
  annotation_custom_positive = rbind(annotation_custom_positive)
  annotation_custom_positive$V2 = annotation_custom_positive$V3 -250
  annotation_custom_positive$V2 = annotation_custom_positive$V2 + 1
  annotation_custom_positive_split = split(annotation_custom_positive,f = annotation_custom_positive$V4,drop = T )
  positive_ranges = lapply(annotation_custom_positive_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4)))
  allAnnotations_plus_ranges_reduced = lapply(positive_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )
  
  
  reducedToDf = function(reduced){
    reduced <- data.frame(seqnames=seqnames(reduced),
                          starts=start(reduced),
                          ends=end(reduced),
                          names=c(names(reduced)),
                          scores=0,strand = strand(reduced))
    return(reduced)
  }
  
  allAnnotations_plus_ranges_reduced_df = lapply(allAnnotations_plus_ranges_reduced,function(x) reducedToDf(x))
  
  
  
  ##### minus strand 
  
  
  annotation_custom_negative = rbind(annotation_custom_negative)
  
  annotation_custom_negative$V3 = annotation_custom_negative$V2 + 250
  
  ## changing to 1 based
  
  annotation_custom_negative$V2 = annotation_custom_negative$V2 + 1
  
  annotation_custom_negative_split = split(annotation_custom_negative,f = annotation_custom_negative$V4,drop = T )
  
  negative_ranges = lapply(annotation_custom_negative_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4)))
  
  
  allAnnotations_minus_ranges_reduced = lapply(negative_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )
  #allAnnotations_plus_ranges_reduced = unlist(reduce(split(a[[1]], elementMetadata(a[[1]])$names))) 
  
  
  
  allAnnotations_minus_ranges_reduced_df = lapply(allAnnotations_minus_ranges_reduced,function(x) reducedToDf(x))
  
  
  allAnnotations_plus_ranges_reduced_df = do.call(rbind,allAnnotations_plus_ranges_reduced_df)
  allAnnotations_minus_ranges_reduced_df = do.call(rbind,allAnnotations_minus_ranges_reduced_df)
  
  allAnnotations = rbind(allAnnotations_plus_ranges_reduced_df,allAnnotations_minus_ranges_reduced_df)
  
  allAnnotations$starts = allAnnotations$starts -1
  return(allAnnotations)
  
  }

countingWindows_noPAS = getCountingWindows(percent90_noPAS)

countingWindows_PAS = getCountingWindows(percent90_PAS)



write.table(countingWindows_noPAS,paste0(BOut,"/polyAmapping_allTimepoints/n_100_global_a0/countingWindows_noPAS.bed"),sep="\t",quote = F,row.names = F,col.names=F)

write.table(countingWindows_PAS,paste0(BOut,"/polyAmapping_allTimepoints/n_100_global_a0/countingWindows_PAS.bed"),sep="\t",quote = F,row.names = F,col.names = F)

write.table(percent90_noPAS,paste0(BOut,"/polyAmapping_allTimepoints/n_100_global_a0/highest90percent_noPAS.bed"),sep="\t",quote = F,row.names = F,col.names = F)

write.table(percent90_PAS,paste0(BOut,"/polyAmapping_allTimepoints/n_100_global_a0/highest90percent_PAS.bed"),sep="\t",quote = F,row.names = F,col.names = F)

countingWindows_PAS$type = "PAS"
countingWindows_noPAS$type = "noPAS"

allCountingWindows = rbind(countingWindows_noPAS,countingWindows_PAS)

write.table(allCountingWindows,paste0(BOut,"/polyAmapping_allTimepoints/n_100_global_a0/allCountingWindows.bed"),sep="\t",quote = F,row.names = F,col.names = F)