options(stringsAsFactors = FALSE)
library(stringr)
library(digest)

read.config <- function(config.path) {
  config <- read.table(config.path,sep="\t")
  obj <- list()
  for(i in 1:nrow(config)) {
    obj[[i]] <- config[i,2]
  }
  names(obj) <- config[,1]
  
  if(!grepl("^/",obj$analysis_path)) {
    obj$analysis_path <- paste(getwd(),"/",obj$analysis_path,sep="")
  }
  obj$config <- cbind(names(obj),unlist(obj))
  return(obj)
}
args <- commandArgs(trailingOnly = T)
cmd <- args[2]
if(is.na(cmd)) {
  stop("LiRA error: missing command")
}
config <- read.config(args[1])

#global settings
GAP_REQUIREMENT <- 2000 #gap between SNPs required for regions to be analyzed separately (i.e. it is unlikely there will be any proper reads spanning a pair of positions over this distance)
READS_TARGET <- 100000 #target number of reads per linkage job

log <- function(log,cmd,intern=F,capture.out=F) {
  newcmd <- paste(gsub(" | ",paste(" 2>> ",log," | ",sep=""),cmd,fixed=T)," 2>> ",log,sep="")
  newcmd <- gsub(" && ",paste(" 2>> ",log," && ",sep=""),newcmd,fixed=T)
  if(capture.out) {
    newcmd <- paste(newcmd," >> ",log,sep="")
  }
  write.table(x = c(date(),cmd),file = log,append = T,col.names = F,row.names = F,quote = F)
  tmp <- system(newcmd,intern=intern)
  if(intern) {
    return(tmp)
  } else {
    note <- c(paste("Exit: ",tmp,sep=""),"","")
    write.table(x=note,file=log,append=T,col.names = F,row.names = F,quote=F)
  }
}

out.log <- function(out,line) {
  write.table(x = paste(date(),": ",line,sep=""),file=out,append=T,col.names=F,row.names=F,quote=F)
}

if(cmd == "setup") {
  owd <- getwd()
  dir.create(config$analysis_path)
  setwd(config$analysis_path)
  
  #soft link bams
  res <- system(paste("ln -s ",config$bam," reads.bam",sep=""))
  res <- system(paste("ln -s ",config$bam_index," reads.bam.bai",sep=""))
  
  res <- system(paste("ln -s ",config$vcf," input_calls.vcf.gz",sep=""))
  res <- system(paste("ln -s ",config$vcf_index," input_calls.vcf.gz.tbi",sep=""))
  write.table(x=config$config,file = "config.txt",sep="\t",col.names = F,row.names = F,quote=F)
  
  system("mkdir job_scripts")
  system("touch .setup_done")
  
  setwd(owd)
}

if(cmd == "split_genome") {
  chromosome <- args[3]
  owd <- getwd()
  setwd(config$analysis_path)
  
  #check for setup
  if(is.na(file.info(".setup_done")$size)) {
    stop("LiRA error: setup not done.  Use the setup command prior to using the 'split_genome' command.")
  }
  suppressWarnings(dir.create(chromosome))
  setwd(chromosome)
  
  outfile <- "split_genome.out.txt"
  logfile <- "split_genome.log.txt"
  
  #getting sites
  out.log(outfile,"Getting sites (output: sites.bed)")
  system("rm -rf jobs")
  
  #determine regions that should be analyzed together due to reads spanning multiple SNP sites
  log(logfile,paste("bcftools view -i 'TYPE=\"snp\" & N_ALT=1' ",config$vcf," ",chromosome,
                    " | grep -v '#'",
                    " | cut -f 1,2,4,5",
                    " | awk '{print $1\"\t\"$2-1\"\t\"$2\"\t\"$3\";\"$4}'",
                    " > sites.bed",
                    " && bedtools merge -d ",GAP_REQUIREMENT," -i sites.bed",
                    " > regions.bed",
                    " && bedtools multicov -bams ../reads.bam -bed regions.bed",
                    " > coverage.bed",
                    sep=""))
  coverage <- read.table("coverage.bed",sep="\t",header=FALSE)
  cum.reads <- cumsum(coverage[,4])
  
  #assign regions to jobs targeting READS_TARGET
  run.assignment <- ceiling(cum.reads/READS_TARGET)
  jobs <- split(coverage[,1:3],run.assignment)
  dir.create("jobs")
  for(j in seq_along(jobs)) {
    prefix <- paste("jobs/",j,sep="")
    dir.create(prefix)
    write.table(jobs[[j]],file=paste(prefix,"/regions.bed",sep=""),sep="\t",col.names = F,row.names = F,quote=F)
    
    log(logfile,paste("#Extract relevant sites within the assigned regions.bed"))
    log(logfile,paste("bedtools intersect -wa -a sites.bed -b ",prefix,"/regions.bed > ",prefix,"/sites.bed",sep=""))
    
    js <- paste("../job_scripts/",chromosome,"_",j,".sh",sep="")
    lines <- c("#!/bin/bash",paste("Rscript --vanilla $LIRA_DIR/scripts/main.R ",config$analysis_path,"/config.txt local_region_function ",config$analysis_path,"/",chromosome,"/jobs/",j,sep=""))
    writeLines(text=lines,con=js)
    system(paste("chmod +x ",js,sep=""))
  }
  setwd(owd)
}

if(cmd == "local_region_function") {
  wrapper <- function() {
    work.dir <- args[3]
    owd <- getwd()
    setwd(work.dir)
    outfile <- "local.region.function.out.txt"
    logfile <- "local.region.function.log.txt"
    
    #Check if we're done.
    if(!is.na(file.info(paste("linkage.rda",sep=""))$size)) {
      out.log(outfile,"Job already finished.")
      return(0)
    }
    
    #Check for saved work
    if(class(suppressWarnings(try(load("results.rda"),silent=T))) == "try-error") {
      
      #check if sites by read has been started but not finished
      if(!is.na(file.info("sites.by.read.txt")$size)) {
        out.log(outfile,"Resuming")
        last.pos <- log(logfile,paste("cat sites.by.read.txt | cut -f 4 | uniq | tail -2 | head -1",sep=""),intern=T)
        log(logfile,paste("cat sites.by.read.txt | head -n -1 | awk '{if($4 < ",last.pos,") print}' > tmp.sites.by.read.txt && mv tmp.sites.by.read.txt sites.by.read.txt",sep=""))
        log(logfile,paste("cat sites.bed | awk '{if($2 >= ",last.pos,") print}' > tmp.sites.bed",sep=""))
        out <- log(logfile,paste("$LIRA_DIR/scripts/linkage.py --bam ",config$analysis_path,"/reads.bam --fasta ",config$reference," --bed tmp.sites.bed >> sites.by.read.txt",sep=""))
      } else {
        out.log(outfile,paste("Getting pileup by read (output: ",work.dir,"/sites.by.read.txt)",sep=""))
        out <- log(logfile,paste("$LIRA_DIR/scripts/linkage.py --bam ",config$analysis_path,"/reads.bam --fasta ",config$reference," --bed sites.bed > sites.by.read.txt",sep=""))
      }
      
      if(as.numeric(system(paste("wc -l sites.by.read.txt | tr ' ' '\t' | cut -f 1",sep=""),intern=T)) == 0) {
        linkage <- data.frame(RR=numeric(0),RV=numeric(0),VR=numeric(0),VV=numeric(0))
        save(linkage,file=paste("linkage.rda",sep=""))
        setwd(dir=owd)
        return(0)
      }
      results <- try(read.table(paste("sites.by.read.txt",sep=""),comment.char="",quote="",sep="\t",header=F,colClasses=c("character","character","character","numeric","character","character","numeric")),silent=T)
      if(class(results) == "try-error") {
        system(paste("rm sites.by.read.txt",sep=""))
        stop("Bad sites.by.read -- removed for next run.")
      }
      names(results) <- c("readname","id","chromosome","position","ref","alt","qual")
      results$ref <- toupper(results$ref)
      results$targeted_alt <- gsub("([^;]*);([^;]*)","\\2",results$id)
      results$call <- NA
      results$call[results$alt == results$targeted_alt] <- "V"
      results$call[results$alt == results$ref] <- "R"
      out.log(outfile,paste("total unexpected calls: ",sum(is.na(results$call)),sep=""))
      results$call[is.na(results$call)] <- "N"
      keep <- names(which(table(gsub("([^ ]*) ([^ ]*)","\\1",unique(paste(results$read,results$position)))) > 1))
      results <- results[results$readname %in% keep,]
      
      if(nrow(results) == 0) {
        linkage <- data.frame(RR=numeric(0),RV=numeric(0),VR=numeric(0),VV=numeric(0))
        save(linkage,file=paste("linkage.rda",sep=""))
        setwd(dir=owd)
        return(0)
      }
      
      #Save work
      save(results,file=paste("results.rda",sep=""))
    }
    r <- unique(results$readname)
    batches <- split(r,rep(seq(from=1,to=ceiling(length(r)/1000)),each=1000)[1:length(r)])
    local.local.region.function <- function(read.batch.num) {
      batch <- data.frame(id=character(0),call=character(0),count=numeric(0))
      batch.rda <- paste("batch",read.batch.num,".rda",sep="")
      
      tmp <- results[results$readname %in% batches[[read.batch.num]],]
      tmp$qual <- NULL
      tmp <- aggregate(x = tmp$call,by=list(readname=tmp$readname,id=tmp$id,chromosome=tmp$chromosome,position=tmp$position),function(x) {
        if(length(unique(x)) > 1) {
          return("N")
        } else {
          return(unique(x))
        }
      })
      names(tmp)[5] <- "call"
      tmp <- tmp[!(tmp$call == "N"),]
      if(nrow(tmp) == 0) {
        save(batch,file=batch.rda)
      }
      ind <- aggregate(x = tmp$position, by = list(readname=tmp$readname), function(x) {
        length(unique(x))
      })
      keep <- ind[ind[,2] > 1,1]
      tmp <- tmp[tmp$readname %in% keep,]
      if(nrow(tmp) == 0) {
        save(batch,file=batch.rda)
      }
      tmp$id <- paste(tmp$chromosome,tmp$position,tmp$id,sep=";")
      tmp <- tmp[order(tmp$position),]
      tmp <- split(tmp[,c("id","call")],tmp$readname)
      names(tmp) <- NULL
      tmp <- do.call(rbind,lapply(tmp,function(x){
        ind <- t(combn(1:nrow(x),m=2))
        tp <- data.frame(id=paste(x$id[ind[,1]],x$id[ind[,2]],sep="~"),call=paste(x$call[ind[,1]],x$call[ind[,2]],sep=""))
        return(tp)
      }))
      tmp <- aggregate(x=tmp$call,by=list(id=tmp$id, call=tmp$call),FUN=length)
      names(tmp)[3] <- "count"
      batch <- tmp
      save(batch,file=paste("batch",read.batch.num,".rda",sep=""))
    }
    
    l <- list.files(pattern="batch")
    l <- file.info(l)
    l <- rownames(l)[!(l$size == 0)]
    todo <- seq_along(batches)
    todo <- todo[!(todo %in% as.numeric(gsub("[^0-9]","",l)))]
    if(length(todo) > 0) {
      for(i in todo) {
        local.local.region.function(i)
      }
    }
    
    l <- paste("batch",seq_along(batches),".rda",sep="")
    all.id <- character(0)
    for(i in l) {
      load(i)
      all.id <- unique(c(all.id,batch$id))
    }
    n <- length(all.id)
    linkage <- data.frame(RR=numeric(n),RV=numeric(n),VR=numeric(n),VV=numeric(n))
    rownames(linkage) <- all.id
    for(i in l) {
      load(i)
      for(call in c("RR","RV","VR","VV")) {
        tmp <- batch[batch$call == call,]
        linkage[tmp$id,call] <- linkage[tmp$id,call] + tmp$count
      }
    }
    save(linkage,file=paste("linkage.rda",sep=""))
    out.log(outfile,"Done.")
    setwd(dir = owd)
  }
  wrapper()
}

if(cmd == "get_vcf_info") {
  chromosome <- args[3]
  owd <- getwd()
  setwd(config$analysis_path)
  suppressWarnings(dir.create(chromosome))
  setwd(chromosome)
  outfile <- "get_vcf_info.out.txt"
  logfile <- "get_vcf_info.log.txt"
  cmds <- c(paste("bcftools query -i 'TYPE=\"snp\" & N_ALT=1' -r ",chromosome," -s ",config$sample," -f '",c("%CHROM;%POS;%REF;%ALT\\t%ID",
                                                                                                             "%INFO/DBSNP_CAF",
                                                                                                             "[%GT]",
                                                                                                             "[%AD]\\t[%GQ]",
                                                                                                             "%CHROM\\t%POS"),"\\n' ",config$vcf," 2> /dev/null ",
                  c("> .tmp1",
                    "| tr ',' '\\t' | cut -f 1 | sed 's#^[.]$#1#g' > .tmp2",
                    "> .tmp3",
                    " | tr ',.' '\\t0' > .tmp4",
                    " | awk '{print $1\"\\t\"$2-2\"\\t\"$2+1}' > .tmp5"),sep="",collapse=" && "),
            paste("bedtools getfasta -fi ",config$reference," -bed .tmp5 -fo - | grep -v '>' | tr 'actg' 'ACTG' > .tmp6",sep=""),
            "paste .tmp1 .tmp2 .tmp3 .tmp4 .tmp6 > .tmp7",
            "rm .tmp1",
            "rm .tmp2",
            "rm .tmp3",
            "rm .tmp4",
            "rm .tmp5",
            "rm .tmp6",
            "echo \"site\tid\tpop.ref.freq\tgenotype\tref\talt\tgq\tcontext\" > vcf-info.txt",
            "cat .tmp7 | sed -E 's#^([^\\t]*)\\t([^\\t]*)\\t([^\\t]*)\\t([^\\t]*)\\t([^\\t]*)\\t([^\\t]*)\\t([^\\t]*)\\t*$#\\1\\t\\2\\t\\3\\t\\4\\t\\5\\t\\6\\t0\\t\\7#g' >> vcf-info.txt",
            "rm .tmp7")
  cmd <- paste(cmds,collapse=" && ")
  log(logfile,cmd)
  setwd(owd)
}

if(cmd == "get_alt_counts") {
  chromosome <- args[3]
  owd <- getwd()
  setwd(config$analysis_path)
  suppressWarnings(dir.create(chromosome))
  setwd(chromosome)
  outfile <- "get_alt_counts.out.txt"
  logfile <- "get_alt_counts.log.txt"
  cmd <- paste(c(paste("bcftools query -i 'TYPE=\"snp\" & N_ALT=1' -r ",chromosome," -f '",c("%CHROM\\t%POS","%CHROM;%POS;%REF;%ALT"),"\\n' ",config$vcf," 2> /dev/null ",c(" | awk '{print $1\"\\t\"$2-1\"\\t\"$2}' > .tmp1"," > .tmp2"),sep="",collapse=" && "),
                 "echo \"id\tA\tC\tG\tT\" > alt-counts.txt",
                 "$LIRA_DIR/scripts/bulk_check.py --bam ../reads.bam --bed .tmp1 > .tmp3",
                 "paste .tmp2 .tmp3 >> alt-counts.txt",
                 "rm .tmp1",
                 "rm .tmp2",
                 "rm .tmp3"),sep="",collapse=" && ")
  cmd <- paste(cmd,collapse=" && ")
  log(logfile,cmd)
  setwd(owd)
}

if(cmd == "collect_results"){
  setwd(config$analysis_path)
  scripts <- list.files("job_scripts")
  f <- gsub(".sh","",paste(gsub("_","/jobs/",list.files("job_scripts")),"/linkage.rda",sep=""))
  test <- file.info(f)$size
  if(any(is.na(test))) {
    writeLines(paste(config$analysis_path,"/job_scripts/",scripts[is.na(test)],sep=""))
    stop("These jobs still undone.")
  }
  #linkage <- do.call(rbind,lapply(f,function(x){load(x); return(linkage)}))

}