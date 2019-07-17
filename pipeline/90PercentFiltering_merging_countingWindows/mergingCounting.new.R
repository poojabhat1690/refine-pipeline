### creating counting windows for counting T>C conversion and for transcriptional output
library(checkmate)
library(dplyr)
#### counting windows for half life calculations
#BOut="//clustertmp/bioinfo/pooja/SLAMannotation/dr/output/"
annotation_custom = read.table(paste0(BOut, "/final90percent/ends_greater90percent_intergenic_n100.txt"),stringsAsFactors = F,sep="\t",header = T)

assertDataFrame(annotation_custom,ncols=11)

annotation_custom = annotation_custom[,c(1:8)]

colnames(annotation_custom) = paste("V",c(1:8),sep="")


### refSeq annotation



### adding refSeq anntoation here 

#refSeq = read.table("/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation/annotations_mm10/refSeq_mm10_GRCm38_06-march-2017/processed/refSeq_mrna_utrsPresent.bed",stringsAsFactors = F)
#ucscDir="///groups/ameres/Pooja/Projects/zebrafishAnnotation//dr11/refSeq_dr11_GRCz11_2019Sep20//"

refFlat <- read.table(paste0(ucscDir, "/refFlat.txt"), stringsAsFactors = F)
refSeqTranscr <- refFlat[grep("NM_", refFlat$V2),]


#add ncRNA
#ensemblDir="//groups/ameres/Pooja/Projects/zebrafishAnnotation//dr11/ensembl_dr11_Ensembl_Genes_93//"


allEnsembl = read.table(paste0(ensemblDir,"transcriptStartsAndEnds_all.txt"),sep="\t",header=T,stringsAsFactors = F)

### removing the rRNas and snoRNA
# allEnsembl = allEnsembl[-which(allEnsembl$transcript_biotype =="snRNA"),]
# allEnsembl = allEnsembl[-which(allEnsembl$transcript_biotype =="rRNA"),]

classesToInclude = c("antisense", "bidrectional_promoter_lncRNA", "lincRNA", "macro_lncRNA", "processed_transcript", "sense_intronic", "sense_overlapping","protein_coding")
allEnsembl = allEnsembl[allEnsembl$transcript_biotype %in% classesToInclude,]

### not including these anymore but just including long non-coding RNAs from ensembl

ncRNA <- read.table(paste0(ensemblDir, "/ncRNA_refSeq.txt"), stringsAsFactors=F, header=T,sep="\t")

write.table(ncRNA, paste0(BOut, "/final90percent/included.ncRNA.txt"),sep="\t",quote = F,row.names = F,col.names = F)

#refSeqTranscr <- rbind(refSeqTranscr, ncRNA)
#refSeqTranscr <- merge(refFlat, refSeq, by.x = "V2", by.y = "V7")
#nrow(refSeq)
#nrow(refSeqTranscr)
allEnsembl <- unique(allEnsembl[,c("chromosome_name", "transcript_start", "transcript_end", "external_gene_name", "strand", "ensembl_transcript_id")])
colnames(allEnsembl) <- c("chr", "start", "end", "gid", "strand", "tid")

allEnsembl <- cbind(allEnsembl[,1:4], 0,  allEnsembl[,5:6])

allEnsembl$V8 = "EnsemblOriginal"
colnames(allEnsembl) <- paste0("V", 1:8)



annotation_custom_positive = annotation_custom[which(annotation_custom$V6 == "+"),]
annotation_custom_negative = annotation_custom[which(annotation_custom$V6 == "-"),]

# noEnsemblEntries = allEnsembl$V7[which(allEnsembl$V4 == "")]
# allEnsembl$V4[which(allEnsembl$V4 == "")]<-noEnsemblEntries

allEnsembl_positive = allEnsembl[which(allEnsembl$V6 == "+"),]
allEnsembl_negative = allEnsembl[which(allEnsembl$V6 == "-"),]


#### FIX THIS....

# mtTranscripts = read.table(paste0(ensemblDir,"/chrMT.txt"),sep="\t",stringsAsFactors = F)
# mtTranscripts = mtTranscripts[,c("V1","V2","V3","V7","V4","V6")]
# colnames(mtTranscripts) <- c("chr", "start", "end", "gid", "strand", "tid")
# 
# mtTranscripts_positive = mtTranscripts[which(mtTranscripts$strand == "+"),]
# mtTranscripts_negative = mtTranscripts[which(mtTranscripts$strand == "-"),]

######################## positive strand ###########################
library(GenomicRanges)
annotation_custom_positive = rbind(annotation_custom_positive,allEnsembl_positive)

annotation_custom_positive$V2 = annotation_custom_positive$V3 -250
annotation_custom_positive$V2 = annotation_custom_positive$V2 + 1

annotation_custom_positive_split = split(annotation_custom_positive,f = annotation_custom_positive$V4,drop = T )

positive_ranges = lapply(annotation_custom_positive_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4)))

#positive_ranges <- with(annotation_custom_positive,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4))


allAnnotations_plus_ranges_reduced = lapply(positive_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )
#allAnnotations_plus_ranges_reduced = unlist(reduce(split(a[[1]], elementMetadata(a[[1]])$names))) 


reducedToDf = function(reduced){
  reduced <- data.frame(seqnames=seqnames(reduced),
                        starts=start(reduced),
                        ends=end(reduced),
                        names=c(names(reduced)),
                        scores=0,strand = strand(reduced))
  return(reduced)
}

allAnnotations_plus_ranges_reduced_df = lapply(allAnnotations_plus_ranges_reduced,function(x) reducedToDf(x))

################## minus strand ############################

annotation_custom_negative = rbind(annotation_custom_negative,allEnsembl_negative)

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

### converting back to 0 based annotations : 

allAnnotations$starts = allAnnotations$starts -1

#valid chromosomes
chromosomes_mm = paste0("chr",c(1:25,"M"))
allAnnotations = allAnnotations[is.element(allAnnotations$seqnames, chromosomes_mm),]

write.table(allAnnotations,paste0(BOut, "/final90percent/allAnnotations.bed"),sep="\t",quote = F,row.names = F,col.names = F)

#########################################################################
## counting windows for transcriptional output and multimapping - for this we need the ensembl 3' utr annotations, refSeq 3' utr annotations
## and intergenic peaks that pass the 90% threshold
#########################################################################

refSeq = read.table(paste0(ucscDir, "/processed/refSeq_mrna_utrsPresent.bed"),stringsAsFactors = F)
ensembl = read.delim(paste0(ensemblDir, "/proteinCoding_annotatedUTRs.bed"),stringsAsFactors = F,header = F)

allAnnotations <- cbind(allAnnotations, "250CountWindow")
colnames(allAnnotations) <- paste0("V", 1:7)
refSeq_ensembl = rbind(refSeq,ensembl,allAnnotations)
refSeq_ensembl_positive = refSeq_ensembl %>% filter(V6=="+")
refSeq_ensembl_negative = refSeq_ensembl %>% filter(V6=="-")


#annotation_custom = read.table("/clustertmp/pooja/mESCinput/final90percent//ends_greater90percent_intergenic_n100.txt",stringsAsFactors = F,sep="\t",header = T)
#annotation_custom_intergenic = annotation_custom[is.element(annotation_custom$peakKind,"intergenic"),]
#annotation_custom_intergenic = annotation_custom_intergenic[,c(1:7)]
#colnames(annotation_custom_intergenic) = paste0("V",c(1:ncol(annotation_custom_intergenic)))
#annotation_custom_intergenic_positive = annotation_custom_intergenic %>% filter(V6 == "+") %>% mutate(V2 = V3-250)
#annotation_custom_intergenic_negative = annotation_custom_intergenic %>% filter(V6 == "-") %>% mutate(V3 = V2+250)

total_positive = refSeq_ensembl_positive
#total_positive = total_positive[!duplicated(total_positive[,c(1:5)]),]
total_positive$V2 = total_positive$V2 +1
total_negative = refSeq_ensembl_negative
#total_negative = total_negative[!duplicated(total_negative[,c(1:5)]),]
total_negative$V2 = total_negative$V2 + 1
###### now reducing Ovelrapping annotations per gene



### positive strand 

total_positive_split = split(total_positive,f = total_positive$V4,drop = T )

total_positive_split_ranges = lapply(total_positive_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4))
)


total_positive_reduced = lapply(total_positive_split_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )
#allAnnotations_plus_ranges_reduced = unlist(reduce(split(a[[1]], elementMetadata(a[[1]])$names))) 



total_positive_reduced_df = lapply(total_positive_reduced,function(x) reducedToDf(x))

total_positive_reduced_df = do.call(rbind,total_positive_reduced_df)


####### negative strand 



total_negative_split = split(total_negative,f = total_negative$V4,drop = T )

total_negative_split_ranges = lapply(total_negative_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4))
)


total_negative_reduced = lapply(total_negative_split_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )
#allAnnotations_plus_ranges_reduced = unlist(reduce(split(a[[1]], elementMetadata(a[[1]])$names))) 



total_negative_reduced_df = lapply(total_negative_reduced,function(x) reducedToDf(x))

total_negative_reduced_df = do.call(rbind,total_negative_reduced_df)

countingWindowsTranscriptionalOutput = rbind(total_positive_reduced_df,total_negative_reduced_df)
countingWindowsTranscriptionalOutput$starts = countingWindowsTranscriptionalOutput$starts -1


write.table(countingWindowsTranscriptionalOutput,paste0(BOut, "/final90percent/countingWindows_transcriptionalOutput.bed"),sep="\t",quote = F,row.names = F,col.names = F)



