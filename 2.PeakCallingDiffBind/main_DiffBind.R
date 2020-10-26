# A script to use Diffbind to get Differentially 5hmC Enriched Genes
# Author: Tian

library("DiffBind")

message("Prepare a Sample Sheet for loading all file")

Samples <- dir("../2.PeakCalling/macs2/", pattern="_peaks.narrowPeak")

SampleList <- data.frame(SampleID=substr(Samples,1,8),
                         Tissue=substr(Samples,5,6),
                         Factor="",
                         Condition=substr(Samples,5,6),
                         Treatment="5hmC",
                         Replicate="",
                         bamReads=paste0("../1.Preprocess/myGreyList/", substr(Samples,1,12), ".grey_filtered.bam"),
                         ControlID="",
                         bamControl="../1.Preprocess/MergedInput.bam",
                         Peaks=paste0("../2.PeakCalling/macs2/", substr(Samples, 1, 12) , "_peaks.narrowPeak"),
                         PeakCaller="narrow")

SampleList$Tissue[SampleList$Tissue %in% c("LT", "NL")] <- "Liver"
SampleList$Tissue[SampleList$Tissue %in% c("NC", "TC")] <- "Colon"

SampleList$Condition[SampleList$Condition %in% c("LT", "TC")] <- "Tumour"
SampleList$Condition[SampleList$Condition %in% c("NL", "NC")] <- "Normal"

SampleList$Replicate <- c(1:5, 1:4, 1:5, 1:4)

write.csv(SampleList, file="SampleList.csv",quote=F, row.names=F)

SampleList <- read.csv("./SampleList.csv", header=T)
message("Loading Peaks into DiffBind")
myDBA <- dba(sampleSheet="./SampleList.csv", dir="./")
plot(myDBA)

message("Working on NC/TC comparison")
message("Generate Consensus Peak")
NCTC <- dba(myDBA, mask=myDBA$masks$Colon)

myDBA_consensus <- dba.peakset(NCTC, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2)
consensus <- dba(myDBA_consensus, mask=myDBA_consensus$masks$Consensus, minOverlap=1)
NCTC.OL <- dba.overlap(myDBA_consensus, myDBA_consensus$masks$Consensus)
consensus_peaks <- dba.peakset(consensus, bRetrieve=TRUE)
NCTC <- dba.count(NCTC, peaks=consensus_peaks, summits=250, bParallel=30) # This step will cost a lot of time.

message("Generate Unique/Consensus Peaks")

message("Differential Analysis")
NCTC <- dba.contrast(NCTC, categories=DBA_CONDITION)
NCTC <- dba.analyze(NCTC)
NCTC.peaks <- dba.report(NCTC, bCalled=TRUE, th=1)

message("Annotate Regions")
library(org.Hs.eg.db)
library(ChIPpeakAnno)
data(TSS.human.GRCh38)
NCTC.peaks <- annotatePeakInBatch(NCTC.peaks, AnnotationData=TSS.human.GRCh38)
NCTC.peaks <- addGeneIDs(NCTC.peaks, "org.Hs.eg.db", IDs2Add = c("entrez_id", 'symbol'))

message("Unique Peaks for TC")
UniqueTC <- NCTC.peaks[which(NCTC.peaks$Called1 >=2 & NCTC.peaks$Called2 < 2 & NCTC.peaks$FDR <= 0.1)]

message("Unique Peaks for NC")
UniqueNC <- NCTC.peaks[which(NCTC.peaks$Called1 < 2 & NCTC.peaks$Called2 >= 2 & NCTC.peaks$FDR <= 0.1)]

message("Differential Peaks (both exist but have difference)")
DiffPeak <- NCTC.peaks[which(NCTC.peaks$Called1 >=2 & NCTC.peaks$Called2 >= 2 & NCTC.peaks$FDR <= 0.1)]

#NCUnique <- data.frame(NCTC.OL[[1]])
#TCUnique <- data.frame(NCTC.OL[[2]])
#Consensus <- data.frame(NCTC.OL[[3]])

#NCTC.peaks <- data.frame(NCTC.peaks)


## Now I want to compare LT/TC group
#message("Working on LT/TC comparison")
#message("Generate Consensus Peak")
#LTTC <- dba(myDBA, mask=myDBA$masks$Tumour)
#
#message("Generate Consensus Peak")
#myDBA_consensus <- dba.peakset(LTTC, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=0.66)
#myDBA_consensus <- dba(myDBA_consensus, mask=myDBA_consensus$masks$Consensus, minOverlap=1)
#consensus_peaks <- dba.peakset(myDBA_consensus, bRetrieve=TRUE)
#LTTC <- dba.count(LTTC, peaks=consensus_peaks, summits=250, bParallel=30) # This step will cost a lot of time.
#
#message("Generate Unique/Consensus Peaks")
#NCTC.OL <- dba.overlap(myDBA_consensus, myDBA_consensus$masks$Consensus)
#
#LTTC <- dba.contrast(LTTC, categories=DBA_TISSUE)
#LTTC <- dba.analyze(LTTC)
#LTTC.peaks <- dba.report(LTTC, bCalled=TRUE)
#
#message("Annotate Regions")
#library(ChIPpeakAnno)
#data(TSS.human.GRCh38)
#LTTC.peaks <- annotatePeakInBatch(LTTC.peaks, AnnotationData=TSS.human.GRCh38)
