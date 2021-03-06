

###### this script gets the most distal end per gene. 
### as in input we use all transcript ends from refSeq and ensembl. 


suppressMessages(library(reshape))
suppressMessages(library(checkmate))
suppressMessages(library(Biostrings))
suppressMessages(library(ggplot2))
suppressMessages(library(biomaRt))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))


                                        #library(tidyr)
#### reading in the annotations from refSeq and ensembl

refSeq_unmerged = read.table(paste0(ucscDir, "/processed/refSeq_mrna_utrsPresent.bed"),stringsAsFactors = F)
refFlat <- read.table(paste0(ucscDir, "/refFlat.txt"), stringsAsFactors = F)

refSeqTranscr <- merge(refFlat, refSeq_unmerged, by.x = "V2", by.y = "V7")
nrow(refSeq_unmerged)
nrow(refSeqTranscr)
refSeqTranscr <- unique(refSeqTranscr[,c("V3.x", "V5.x", "V6.x", "V4.y", "V4.x", "V2")])
colnames(refSeqTranscr) <- c("chr", "start", "end", "gid", "strand", "tid")


ensembl = read.delim(paste0(ensemblDir, "/proteinCoding_annotatedUTRs.bed"),stringsAsFactors = F,header = F)
ensemblTrans <- read.delim(paste0(ensemblDir, "/transcriptStartsAndEnds_all.txt"), stringsAsFactors = F, header = T)
head(ensemblTrans)
ensemblTrans <- merge(ensemblTrans, ensembl, by.x = "ensembl_transcript_id", by.y = "V7")
head(ensemblTrans)
ensemblTrans <- unique(ensemblTrans[,c("chromosome_name", "transcript_start", "transcript_end", "external_gene_name", "strand", "ensembl_transcript_id")])
colnames(ensemblTrans) <- c("chr", "start", "end", "gid", "strand", "tid")


allTrans = rbind(refSeqTranscr,ensemblTrans)
allTrans$start <- allTrans$start + 1

gid <- allTrans %>% dplyr::group_by(gid) %>% dplyr::summarize(uid = paste(gid, tid, sep = ",", collapse = '|')) %>% as.data.frame
head(gid)
nrow(gid)
nrow(allTrans)
allTrans <- merge(allTrans, gid, by = "gid")
nrow(allTrans)
head(allTrans)
#allTrans = allTrans[-which(allTrans$gid == ""),]
allTransGR <- with(allTrans, GRanges(chr, IRanges(start = start,end = end),strand = strand,score=0,names=uid))

allTransGRReduced = reduce(split(allTransGR, elementMetadata(allTransGR)$names)) 

allTransGRReducedDF <- do.call(rbind, lapply(allTransGRReduced, as.data.frame))
head(allTransGRReducedDF)

allTransGRReducedDFBed <- cbind(allTransGRReducedDF[,c(1,2,3)], sub(",.*", "", rownames(allTransGRReducedDF)), 0, allTransGRReducedDF$strand, rownames(allTransGRReducedDF))




### in this step, we get the most distal 3' end per gene. 




write.table(allTransGRReducedDFBed,paste0(OutPath, "/toExtend_longestEnsembl_refSeq_n100.txt"),quote = F,sep="\t",row.names = F,col.names = F)
write.table(allTransGRReducedDFBed,paste0(OutPath, "/toExtend_longestEnsembl_refSeq_n100.bed"),quote = F,sep="\t",row.names = F,col.names = F)




# creating a combination of gene annotations from ensembl and refSeq 



############ THE REFSEQ TRANSCRIPT ANNOTATION COULD BE PROCESSED IN THE ANNOTATION SCRIPT i.e  downloading. 

#    refSeqGene = read.table("/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation/annotations_mm10/refSeq_mm10_GRCm38_06-march-2017/refSeqGenes.bed",stringsAsFactors = F)

# refSeqGene = read.table("/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation/intergenicPeaks/annotationsTranscripts/refSeq_ucscTableBrower_mm10_08-march-2017/refSeq_genes.bed",stringsAsFactors = F)



#    ensemblGenes = read.table("/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation/annotations_mm10/ensembl_mm10_Ensembl_Genes_87_06-march-2017/transcriptStartsAndEnds_all.txt",header=T,stringsAsFactors = F)
#    ensemblGenes = ensemblGenes[,c(1:6)]
#    ensemblGenes_rearranged = cbind.data.frame(ensemblGenes$chromosome_name,ensemblGenes$transcript_start,ensemblGenes$transcript_end,ensemblGenes$ensembl_gene_id,0,ensemblGenes$strand)
#    refSeqGene = refSeqGene[,c(1:6)]

#    colnames(refSeqGene) = paste0("V",c(1:ncol(refSeqGene)))
#    colnames(ensemblGenes_rearranged) = paste0("V",c(1:ncol(ensemblGenes_rearranged)))

#    allTranscripts = rbind(refSeqGene,ensemblGenes_rearranged)
#### removing the duplicated entries : 

#    allTranscripts = allTranscripts[!duplicated(allTranscripts[,c("V1","V2","V3","V6"),]),]
    
#    write.table(allTranscripts,"/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation/annotations_mm10/allTranscripts_refSeq_ensembl.bed",sep="\t",quote = F,row.names = F,col.names = F)







    refSeqExon = read.table(paste0("///groups/ameres/bioinformatics/references/danio_rerio/dr10/refSeq_danRer10_exonAnnotations.bed"),stringsAsFactors = F)

# refSeqExon = read.table("/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation/intergenicPeaks/annotationsTranscripts/refSeq_ucscTableBrower_mm10_08-march-2017/refSeq_genes.bed",stringsAsFactors = F)



ensemblExons = read.delim(paste0(ensemblDir, "/exonStartsAndEnds_all.txt"),header=T,stringsAsFactors = F)
    ensemblExons_rearranged = cbind.data.frame(ensemblExons$chromosome_name,ensemblExons$exon_chrom_start,ensemblExons$exon_chrom_end,ensemblExons$ensembl_gene_id,0,ensemblExons$strand)

refSeqExon = refSeqExon[,c(1:6)]

    colnames(refSeqExon) = paste0("V",c(1:ncol(refSeqExon)))
    colnames(ensemblExons_rearranged) = paste0("V",c(1:ncol(ensemblExons_rearranged)))

    allExons = rbind(refSeqExon,ensemblExons_rearranged)
#### removing the duplicated entries : 

    allExons = allExons[!duplicated(allExons[,c("V1","V2","V3","V6"),]),]
    
    write.table(allExons,paste0(OutPath, "/allExons_refSeq_ensembl.bed"),sep="\t",quote = F,row.names = F,col.names = F)




