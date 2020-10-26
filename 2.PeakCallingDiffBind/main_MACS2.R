# This is script to do MACS2 calling
# Author: Tian

bams <- dir("../1.Preprocess/myGreyList/", pattern="*_bnd.sorted.bam")
bams <- unique(sapply(bams, function(x) strsplit(x,split="[.]")[[1]][1])) 

if (!file.exists("./myGreyList")) dir.create("./myGreyList")

commands <- list()
for(i in bams) {
    
    command <- paste0("macs2 callpeak -t ../1.Preprocess/myGreyList/", i, ".sorted.bam -c ../1.Preprocess/MergedInput.bam -f BAM -g hs --outdir macs2 -n ", i, " -B -q 0.1 2> macs2/", i, "-macs2.log")
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
