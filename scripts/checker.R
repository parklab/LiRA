rm(list=ls())
library(digest)
options(stringsAsFactors = F)

args <- commandArgs(trailingOnly = T)
bam <- normalizePath(args[1])
bed <- normalizePath(args[2])
reference <- normalizePath(args[3])
output <- args[4]

tmp.dir <- paste(".",digest(list(args,date())),sep="")
dir.create(tmp.dir)
setwd(tmp.dir)

system(paste("cat ",bed,
             " | awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3,$1\";\"$3\";\"$4}'",
             " > sites.bed",sep=""))

system(paste("bedtools tag -i ",bam," -files sites.bed -names",
             " | samtools view -h -",
             " | grep -e 'YB:[^\\t]*,' -e '^@'",
             " | samtools view -b -",
             " > reads.bam",
             " && samtools index reads.bam",sep=""))

system(paste("samtools view reads.bam",
             " | grep -o 'YB:[^\\t]*'",
             " | sed 's#YB:Z:##g'",
             " | tr ',' '\n'",
             " | tr ';' '\t'",
             " | awk 'BEGIN{OFS=\"\\t\"}{print $1,$2-1,$2,$3\";\"$4}'",
             " | sort -V",
             " | uniq > sites.bed",sep=""))

local <- function() {
  out <- system(paste("$LIRA_DIR/scripts/run_checker.py --bam reads.bam --fasta ",reference," --bed sites.bed > sites.by.read.txt",sep=""))
  if(as.numeric(system(paste("wc -l sites.by.read.txt | tr ' ' '\t' | cut -f 1",sep=""),intern=T)) == 0) {
    linkage <- data.frame(RR=numeric(0),RV=numeric(0),VR=numeric(0),VV=numeric(0),site1=character(0),site2=character(0))
    return(linkage)
  }
  results <- read.table(paste("sites.by.read.txt",sep=""),comment.char="",quote="",sep="\t",header=F,colClasses=c("character","character","character","numeric","character","character"))
  names(results) <- c("readname","id","chromosome","position","ref","alt")
  results$ref <- toupper(results$ref)
  results$targeted_alt <- gsub("([^;]*);([^;]*)","\\2",results$id)
  results$call <- NA
  results$call[results$alt == results$targeted_alt] <- "V"
  results$call[results$alt == results$ref] <- "R"
  results$call[is.na(results$call)] <- "N"
  keep <- names(which(table(gsub("([^ ]*) ([^ ]*)","\\1",unique(paste(results$read,results$position)))) > 1))
  results <- results[results$readname %in% keep,]
  if(nrow(results) == 0) {
    linkage <- data.frame(RR=numeric(0),RV=numeric(0),VR=numeric(0),VV=numeric(0),site1=character(0),site2=character(0))
    return(linkage)
  }
  
  results$id <- paste(results$chromosome,results$position,results$id,sep=";")
  results <- unique(results)
  results <- split(results,results$readname)
  results <- lapply(results,function(x){
    if(length(unique(x$position)) != nrow(x)) {
      return(data.frame(RR=numeric(0),RV=numeric(0),VR=numeric(0),VV=numeric(0),site1=character(0),site2=character(0)))
    }
    x <- x[order(x$position),]
    combos <- t(combn(x=1:nrow(x),m=2))
    n <- nrow(combos)
    nul <- rep(0,n)
    call <- paste(x$call[combos[,1]],x$call[combos[,2]],sep="")
    ret <- data.frame(RR=as.numeric(call == "RR"),RV=as.numeric(call == "RV"),VR=as.numeric(call == "VR"),VV=as.numeric(call == "VV"),site1=x$id[combos[,1]],site2=x$id[combos[,2]])
    rownames(ret) <- paste(ret$site1,ret$site2)
    return(ret)
  })
  n <- unique(unlist(lapply(results,rownames)))
  nul <- rep(0,length(n))
  tmp <- strsplit(n," ")
  combined <- data.frame(RR=nul,RV=nul,VR=nul,VV=nul,site1=sapply(tmp,function(x){x[1]}),site2=sapply(tmp,function(x){x[2]}))
  rownames(combined) <- n
  for(i in seq_along(results)) {
    combined[rownames(results[[i]]),c("RR","RV","VR","VV")] <- combined[rownames(results[[i]]),c("RR","RV","VR","VV")] + results[[i]][,c("RR","RV","VR","VV")]
  }
  rownames(combined) <- NULL
  combined <- combined[,c("site1","site2","RR","RV","VR","VV")]
  return(combined)
}

linkage <- local()
write.table(linkage,file="output.txt",sep="\t",row.names = F,col.names = F,quote=F)
write.table(linkage[-(1:nrow(linkage)),],file=paste("../",output,sep=""),sep="\t",row.names = F,col.names = T,quote=F)
system(paste("cat output.txt | sort -V >> ../",output,sep=""))
setwd("..")
system(paste("rm -r ",tmp.dir,sep=""))
  
