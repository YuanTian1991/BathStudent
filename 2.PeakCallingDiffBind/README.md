# CRC-2: Peak Calling/Diffbind

After preprocessing, next I want to call peaks for each bound file. After searching seems there are some tools I can use, like Pepr, MACS2. etc, but since MACS2 is kind of popular and troditional, so I will use it here.

## 1. Peak Calling

### 1.1 Merge Inputs files

Firstly, since I have 4 phenotype (LT, TC, NC, NT), but only 2 of them have Input data, so I decided to merge all Input data togather and use all Bound file to call peak against it. I personally think, if we have proper Input data for each phenotype, we should definitly call peaks on each Input, so it's just a compromisation solution here.

```r
samtools merge MergedInput.bam `ls ./myGreyList/*_Inp.grey_filtered.bam`
```

### 1.2 Apply peak calling

Then we can start to do peak calling. The key macs2 code is just online:

```bash
macs2 callpeak -t Bound.bam -c Input.bam -f BAM -g hs --outdir macs2 -n SampleName 2> macs2/SampleName-macs2.log
```

The `gs` parameter is vital as it indicates different species genome length, clearly that human and mouse have different genome length. Here I use hs indicate human. This is written in MACS2 manual. 

Since I want to make things easier and faster, I wrote below R script to do MACS2 calling across all samples automatically and parallel.

```r
# This is script to do MACS2 calling
# Author: Tian

bams <- dir("../1.Preprocess/myGreyList/", pattern="*_bnd.sorted.bam")
bams <- unique(sapply(bams, function(x) strsplit(x,split="[.]")[[1]][1]))

if (!file.exists("./myGreyList")) dir.create("./myGreyList")

commands <- list()
for(i in bams) {

    command <- paste0("macs2 callpeak -t ../1.Preprocess/myGreyList/", i, ".sorted.bam -c ../1.Preprocess/MergedInput.bam -f BAM -g hs --outdir macs2 -n ", i, " -B -q 0.1 2> macs2/",   i, "-macs2.log")
    message(command)
    commands <- c(commands, command)
}

runMacs <- function(i)
{
    system(i)
}

library(doParallel)
detectCores()
cl <- makeCluster(18)
registerDoParallel(cl)
getDoParWorkers()

library(foreach)
x <- foreach(i = commands) %dopar% runMacs(i)

registerDoSEQ()
on.exit(stopCluster(cl))
```

### 1.3 Check Peak Calling Result via IGV

After peak calling, we can directly use IGV to check the peak status, for example, I already know IGR5 is an important gene which show significant 5hmC difference between Stem and Differentiated cell, so I checked a bit on that. Below is an perfect example, as Peaks merely only show up in TC samples.

![CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled.png](CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled.png)

Now in theory, we can check genes who have differential peak enrichment status between phenotypes via IGV. However, there are two problem here: 1, there are too many genes across human genome to allow us to check one by one. 2, if both in one gene, both two phenotype (like TC/NC) are enriched, it's hard to judge if their enrichment status are significantly different or not.

So we continue to use Diffbind R pacakge to solve above question.

I think peaks have many many other ways to do downstream analysis, Diffbind is just one of many of them. At this stage, we can start to searching for noval tools for project-specific task.

## 2. Diffbind Analysis

Diffbind is a popular software to ChIP-seq analysis, it load peaks called by software like MACS, and find **consensus peaks**, then finally calculate if there are differential status between phenotypes (like TC/NC). Apart from that, Diffbind have another kind of analysis called `occupancy analysis` which means it will directly find out regions enriched by peaks, these regions may not be consensus peaks, but uniquely turn up in certain phenotype.

So here is what we want to do, we want to use Diffbind, to do occupancy analysis on all 4 phenotype, so that we would know how these peaks are enriched in each phenotype. Also, for consensus peaks show up in both two phenotype, we want to know if there are any changes happen between phenotypes.

Diffbind is an [R package](https://bioconductor.org/packages/release/bioc/html/DiffBind.html), which means it's something totally reply on R programming environment. So it requires user to know how to code in R. After installation, we can try follow the [vignette](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) written by package author to try analysis the data. 

### 2.1 Prepare Sample Sheet CSV file

According to the vignette of Diffbind, it requires a CSV file to load all sample, I don't know if there are othe ways but I decided to follow it's solution. Firstly I need to generate a `sample_sheet.csv` file for loading.

csv file loading is a very common way in Bioinformatics, many software requires a csv file contains sample, directory, phentypes information .etc. Merely the reason is at the beginning age, sequencing companies return customer with a csv file with these information in...

I wrote below script to organise original files into a SampleList R data frame, and wrote it into current directory, for Diffbind loading

```r
# A script to use Diffbind to get Differentially 5hmC Enriched Genes
# Author: Tian

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
```

### 2.2 Diffbind Loading

Following the phenotype, now I load all samples with below code. After loading, we can use `plot(myDBA)` function to draw a heatmap based on all these peaks.

```r
myDBA <- dba(sampleSheet="./SampleList.csv", dir="./")
plot(myDBA)
```

In the heatmap, we  can roughly see that NL/LT samples are clustered togather. Other samples are not very clear. It's OK here, as it use all peaks across human genome to calculate the correlation value, so differential signals hidden in huge genome may not strong enough to separate them.

![CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%201.png](CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%201.png)

Also, we can have a check on loadded samples, because in peak calling stage, different samples will have different numbers of peak returned.

```r
> myDBA
18 Samples, 523337 sites in matrix (853458 total):
         ID Tissue Condition Treatment Replicate Caller Intervals
1  SUL_LT49  Liver    Tumour      5hmC         1 narrow    163214
2  SUL_LT51  Liver    Tumour      5hmC         2 narrow    340856
3  SUL_LT52  Liver    Tumour      5hmC         3 narrow    299212
4  SUL_LT53  Liver    Tumour      5hmC         4 narrow    125792
5  SUL_LT55  Liver    Tumour      5hmC         5 narrow    172529
6  SUL_NC49  Colon    Normal      5hmC         1 narrow    121753
7  SUL_NC51  Colon    Normal      5hmC         2 narrow    220174
8  SUL_NC52  Colon    Normal      5hmC         3 narrow    170373
9  SUL_NC53  Colon    Normal      5hmC         4 narrow     76938
10 SUL_NL49  Liver    Normal      5hmC         1 narrow    207445
11 SUL_NL51  Liver    Normal      5hmC         2 narrow    382434
```

### 2.3 Select samples (TC/LT here)

Above work have helped to to load all samples into R session, here I want to only extract TC and NC phenotype for analysis. According to the vignette, I can use `mask` parameter to do it. The myDBA mask shows below result:

```r
> names(myDBA_consensus$masks)
 [1] "Colon"             ""                  "Normal"
 [4] "Tumour"            "5hmC"              "narrow"
 [7] "Consensus"         "Replicate.1-2-3-4" "All"
[10] "None"
>
```

which means, if you select `Colon` ,it will only return you Colon tissues, with two phenotypes Normal/Tumour separately.

So I used below command to select Colon samples, which means I now exclude Liver Tumour and Normal Liver, just focused on NCTC analysis.

```r
> NCTC <- dba(myDBA, mask=myDBA$masks$Colon)
> NCTC
8 Samples, 231982 sites in matrix (530488 total):
        ID Tissue Condition Treatment Replicate Caller Intervals
1 SUL_NC49  Colon    Normal      5hmC         1 narrow    121753
2 SUL_NC51  Colon    Normal      5hmC         2 narrow    220174
3 SUL_NC52  Colon    Normal      5hmC         3 narrow    170373
4 SUL_NC53  Colon    Normal      5hmC         4 narrow     76938
5 SUL_TC49  Colon    Tumour      5hmC         1 narrow      8947
6 SUL_TC51  Colon    Tumour      5hmC         2 narrow    165713
7 SUL_TC52  Colon    Tumour      5hmC         3 narrow    125350
8 SUL_TC53  Colon    Tumour      5hmC         4 narrow    216705
```

### 2.4 Consensus Peak

After forming proper data, then an important step here to find consensus peaks. Consensus peaks means peaks show up across samples. However, here there are two problem: one is we should try find consensus peaks for each condition, the other is what if some samples does not have this peak, for example in above code block, there are 4 Normal and 4 Tumour, SUL_TC49 has only 8947 peaks, if we based on this, we can have at most 8947 peaks across all TC samples?

The solution Diffbind proposed is set up `minOverlap` parameter, in below code I set the minOverlap as 2, which means it would select peaks enriched in 2 samples. 

```r
myDBA_consensus <- dba.peakset(NCTC, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2)
myDBA_consensus <- dba(myDBA_consensus, mask=myDBA_consensus$masks$Consensus, minOverlap=1)
consensus_peaks <- dba.peakset(myDBA_consensus, bRetrieve=TRUE)
```

We can see how many consensus peaks are find:

```r
> myDBA_consensus
2 Samples, 73135 sites in matrix:
      ID Tissue Condition Treatment Replicate Caller Intervals
1 Normal  Colon    Normal      5hmC   1-2-3-4 narrow     53624
2 Tumour  Colon    Tumour      5hmC   1-2-3-4 narrow     40588
```

So we have get 53624 peaks in Normal Colon, and 40588 peaks in Tumour Colon, seems Tumour Colon have less 5hmC enriched.  This is core part of occupancy analysis. 

The `consensus_peak` above is calcualted based on these two phenotypes, extract peaks that exist in BOTH of them, at least 2 samples each. Because then we will do differential analysis, and difference can only be called on "consensus peak".

Finally, the last step here is to do count, which is slow, Diffbind would get a score for each peak, then use that score for differential calling.

```bash
NCTC <- dba.count(NCTC, peaks=consensus_peaks, summits=250) # This step will cost a lot of time.
```

### 2.5 Occupancy Analysis

Occupancy analysis is a kind of analysis that DiffBind directly get result from Peaks, and if peaks mapped on one fragment of genome, it would extract this fragment as an "peak". So if you have multiple samples for one phenotype, for example here I have 4-5 samples for LT/TC. etc, it allow you to set a cutoff, that if a peak show up in at least certain number of samples. In above code actually I have done this step, the result is in myDBA_consensus.

We can use below code to see unique peaks between two phenotype

```r
dba.plotVenn(myDBA_consensus, myDBA_consensus$masks$Consensus)
```

![CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%202.png](CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%202.png)

So we can see that there are 33113 unique peaks for Normal Colon, 19955 unique peaks for Tumour Colon, in below section, Diffbind will compare the 20067 peaks, see if there are any difference in them. By the way, there is another function `dba.overlap` helps to find these unique peaks. 

```r
NCTC.OL <- dba.overlap(myDBA_consensus, myDBA_consensus$masks$Consensus)
```

Above result contained the unique peaks in each group.

Note that there are two important things here:

1. Occupancy analysis will not tell you if there are any difference between phenotype, it will only tell you if there are peaks enriched. There is NO WAY to compare which peak is stronger.

2. Based on on my understanding, this Occupancy Analysis is focusing on peaks with a less strict statistic calculation, that's why we normally can get a lot of result here. In below Differential Calling result, strict statistic algorithm would be applied across samples.

### 2.6 Differential Calling

Finally, we reach Differential Calling stage, it's vital and important for Bioinformatic analysis, 80% of Bioinformatic works are focusing on get signals between two status, which means looking for differential status between phenotype, on gene expression, methylation status .etc

In above section, we should have get count result from `dba.count` This is an important step, will convert initial peak information into "comparable" counts, like across 500-bp window, how many differential counts are enriched for each sample. 

Then finally I can do differential calling. The code is simple:

```bash
NCTC <- dba.contrast(NCTC, categories=DBA_CONDITION)
NCTC <- dba.analyze(NCTC)
NCTC.peaks <- dba.report(NCTC, bCalled=TRUE, th=1)
```

The `bCalled` parameter is important here, as it will return previous occupancy analysis result. Parameter `th` means we don't do any p value filtering, but return all peak result.

After running, the result looks like below. It gave use the differential region's between TC and NC. Here there are many columns: The first 5 are chr/position/strand. Their width are 500, so we can roughtly know how this software work, it bins all peaks into 500-bp windows, then use linear regression to calcualte their differential enrichment status, surely it did some normalisation ahead. The rest column we can find annotation in the [manaul](https://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf).

```bash
> NCTC.peaks <- dba.report(NCTC
> knitr::kable(head(NCTC.peaks))

|seqnames |    start|      end| width|strand | Conc| Conc_Normal| Conc_Tumour|  Fold| p.value|      FDR| Called1| Called2|
|:--------|--------:|--------:|-----:|:------|----:|-----------:|-----------:|-----:|-------:|--------:|-------:|-------:|
|chrX     | 46375722| 46376222|   501|*      | 3.63|       -1.27|        4.61| -5.88|       0| 0.000106|       1|       0|
|chr13    | 75607372| 75607872|   501|*      | 3.22|       -1.27|        4.19| -5.46|       0| 0.000106|       1|       0|
|chr6     | 40908421| 40908921|   501|*      | 4.74|        0.66|        5.70| -5.03|       0| 0.000133|       2|       1|
|chr11    | 13013994| 13014494|   501|*      | 3.41|       -0.73|        4.37| -5.10|       0| 0.000133|       3|       0|
|chr10    |  3978236|  3978736|   501|*      | 3.43|       -0.98|        4.39| -5.37|       0| 0.000140|       3|       0|
|chr13    | 51358219| 51358719|   501|*      | 3.25|       -1.27|        4.22| -5.49|       0| 0.000140|       1|       0|
>
```

The last two columns are important, it indicates for each peak (row), how many peak signal show up in Group 1 (TC), and Group 2(NC). We will use this this signal to filter out some peak have less than 2 signal.

### 2.7 Peak Annotation

In above section, we have successfully get the differential result. However, before we look further, we better annotate it abit, here I found a quick an easy way to annotate these peaks with [ChIPpeakAnno](https://www.bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html) R package. The code is below:

```r
library(org.Hs.eg.db)
library(ChIPpeakAnno)
data(TSS.human.GRCh38)
NCTC.peaks <- annotatePeakInBatch(NCTC.peaks, AnnotationData=TSS.human.GRCh38)
NCTC.peaks <- addGeneIDs(NCTC.peaks, "org.Hs.eg.db", IDs2Add = c("entrez_id", 'symbol'))
```

Note that `org.Hs.eg.db` package also need to be installed here.

In short, what ChIPpeakAnno doing is mapping our peaks to genome, then tell us what gene they are related. The addGeneIDs here is only to add entrez_id and symbol to our result, otherwise, we only have Ensembl Gene_ID.

After annotation, our result now looks like: (The data frame is too big, that's how it looks after my exporting)

![CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%203.png](CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%203.png)

Now we can have close look at this data frame. It contains 3 column I would focus on: Called1, Called2 and FDR. By using them, we can find out Unique Peaks for TC, Unique Peaks for NC, and Peaks in both phenotype but showing differential changes duing tumouring.

The code is below, in general, we think if a peak only significantly enriched on one phenotype, but not on another phenotype as `Unique Peaks`.

```r
message("Unique Peaks for TC")
UniqueTC <- NCTC.peaks[which(NCTC.peaks$Called1 >=2 & NCTC.peaks$Called2 < 2 & NCTC.peaks$FDR <= 0.1)]

message("Unique Peaks for NC")
UniqueNC <- NCTC.peaks[which(NCTC.peaks$Called1 < 2 & NCTC.peaks$Called2 >= 2 & NCTC.peaks$FDR <= 0.1)]

message("Differential Peaks (both exist but have difference)")
DiffPeak <- NCTC.peaks[which(NCTC.peaks$Called1 >=2 & NCTC.peaks$Called2 >= 2 & NCTC.peaks$FDR <= 0.1)]
```

Let's have a check on these 3 result, by cross checking Peaks IGV:

FAT1, a gene was enriched by peaks in both NC and TC (In DiffPeak object), but showing higher enriched status in TC.

![CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%204.png](CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%204.png)

RASSFF10, a gene show only enriched in TC, from UniqueTC result.

![CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%205.png](CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%205.png)

So DiffBind succesfully finished it's job. We can now write out these 3 files into Excel with write.csv(XXX, file="XXX.csv").

Until now, we nearly have finished Differential analysis using DiffBind package. DiffBind provided a set of function to visualise these result. I am not going try them now.

## 3. GSEA Analysis

Analysis is normally the last step for a bioinformatic pipeline. Here I found a simply tool to do it, [metascape](https://metascape.org/gp/index.html#/main/step1). This is one of the easiest tool to do this job. Just paste the genes (feature colunmn) into the web tool, it would start working.

Here I inputted the genes I discovered in **UniqueTC** and **DiffPeak** because I think both new-emerged signal and differential changed signal duing tumouring are important, and coule maybe indicate something.

![CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%206.png](CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%206.png)

After a whie, I received the report, in the report I can see that, figure 6 clearly indicated that the genes I inserted are significantly enriched with Colon-Cancer status, which is exactly what I want to see.

![CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%207.png](CRC-2%20Peak%20Calling%20Diffbind%20e52b8863e0fb42abbe102dff621219d9/Untitled%207.png)