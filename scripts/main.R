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
  
  for(n in c("analysis_path","reference","bam","bam_index","vcf","vcf_index","phased_vcf")) {
    if(n %in% names(obj)) {
      obj[[n]] <- normalizePath(obj[[n]])
    }
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

if(cmd == "check_results"){
  owd <- getwd()
  setwd(config$analysis_path)
  scripts <- list.files("job_scripts")
  f <- gsub(".sh","",paste(gsub("_","/jobs/",list.files("job_scripts")),"/linkage.rda",sep=""))
  test <- file.info(f)$size
  if(any(is.na(test))) {
    writeLines(paste(config$analysis_path,"/job_scripts/",scripts[is.na(test)],sep=""))
    stop("These jobs still undone.")
  } else {
    print("Jobs finished.")
  }
  setwd(owd)
}

if(cmd == "collect_results") {
  chromosome <- args[3]
  owd <- getwd()
  setwd(config$analysis_path)
  setwd(chromosome)
  f <- paste("jobs/",list.files("jobs"),"/linkage.rda",sep="")
  linkage <- do.call(rbind,lapply(f,function(x){load(x); return(linkage)}))
  site1 <- gsub("([^~]*)~([^~]*)","\\1",rownames(linkage))
  site2 <- gsub("([^~]*)~([^~]*)","\\2",rownames(linkage))
  linkage$site1 <- site1
  linkage$site2 <- site2
  linkage$pair_id <- rownames(linkage)
  rownames(linkage) <- NULL
  pos1 <- as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\2",linkage$site1))
  pos2 <- as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\2",linkage$site2))
  linkage <- linkage[order(pos1,pos2),]
  rownames(linkage) <- 1:nrow(linkage)
  save(linkage,file=paste(config$analysis_path,"/",chromosome,"/linkage.rda",sep=""))
}

if(cmd == "compare_to_bulk") {
  owd <- getwd()
  bulk.config <- read.config(args[3])
  chromosome <- args[4]
  setwd(config$analysis_path)
  setwd(chromosome)
  outfile <- "compare_to_bulk.out.txt"
  logfile <- "compare_to_bulk.log.txt"
  suppressWarnings(dir.create("compare_to_bulk"))
  
  load("linkage.rda")
  single.cell <- linkage
  
  load(paste(bulk.config$analysis_path,"/",chromosome,"/linkage.rda",sep=""))
  bulk <- linkage
  
  #make a site-specific table
  site.frame <- read.table("vcf-info.txt",sep="\t",header=T)
  rownames(site.frame) <- site.frame$site
  site.frame$site <- NULL
  site.frame$chromosome <- chromosome
  site.frame$pos <- as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\2",rownames(site.frame)))
  site.frame$single.cell.ref <- site.frame$ref
  site.frame$single.cell.alt <- site.frame$alt
  site.frame$ref <- gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\3",rownames(site.frame))
  site.frame$alt <- gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\4",rownames(site.frame))
  site.frame$filter <- system(paste("bcftools query -i 'TYPE=\"snp\" & N_ALT=1' -r ",chromosome," -f '%FILTER\\n' ",config$vcf," 2> /dev/null",sep=""),intern=T)
  site.frame <- site.frame[order(site.frame$pos),]
  
  site.frame$context <- paste(site.frame$context,">",site.frame$alt,sep="")
  change <- substr(site.frame$context,2,2) %in% c("G","A")
  site.frame$context[change] <- gsub("(.)(.)(.)","\\3\\2\\1",site.frame$context[change])
  site.frame$context[change] <- gsub("A","W",gsub("C","X",gsub("G","Y",gsub("T","Z",site.frame$context[change]))))
  site.frame$context[change] <- gsub("W","T",gsub("X","G",gsub("Y","C",gsub("Z","A",site.frame$context[change]))))
  
  names(site.frame)[names(site.frame) == "genotype"] <- "single.cell.gt"
  names(site.frame)[names(site.frame) == "gq"] <- "single.cell.gq"
  
  #pull bulk info
  bulk.info <- read.table(paste(bulk.config$analysis_path,"/",chromosome,"/vcf-info.txt",sep=""),sep="\t",header=T)
  rownames(bulk.info) <- bulk.info$site
  bulk.info$site <- NULL
  site.frame[,"bulk.gt"] <- bulk.info$genotype
  site.frame[,"bulk.gq"] <- bulk.info$gq
  site.frame[,"bulk.ref"] <- bulk.info$ref
  site.frame[,"bulk.alt"] <- bulk.info$alt
  
  #correct pop.ref.freq
  ind2 <- is.na(site.frame$pop.ref.freq) & (site.frame$id == ".")
  site.frame$pop.ref.freq[ind2] <- 1
  
  #pull bulk phasing info
  dumb <- system(paste("bcftools query -r ",chromosome," -s ",bulk.config$sample," -f '%CHROM;%POS;%REF;%ALT\t[%GT]\n' ",bulk.config$phased_vcf," 2> /dev/null | grep -e '0|1' -e '1|0' > phasing.txt",sep=""))
  phasing <- read.table("phasing.txt",sep="\t")
  n <- phasing[,1]
  phasing <- phasing[,-1]
  names(phasing) <- n
  phasing <- phasing[names(phasing) %in% rownames(site.frame)]
  site.frame$bulk.phasing <- 0
  site.frame[names(phasing),"bulk.phasing"] <-  as.numeric(as.factor(phasing))
  
  #only look for linkage near bulk het sites in 1KG
  #somatic candidates: bulk.gt == "0/0" or bulk.gt == "./."
  site.frame$in.linkage <- rownames(site.frame) %in% c(single.cell$site1,single.cell$site2)
  site.frame$onek.bulk.het <- (site.frame$pop.ref.freq != 1) & (site.frame$bulk.gt == "0/1")
  site.frame$somatic <- (site.frame$bulk.gt == "0/0")|(site.frame$bulk.gt == "./.")

  onek.bulk.sites <- rownames(site.frame)[site.frame$onek.bulk.het]
  somatic.sites <- rownames(site.frame)[site.frame$somatic]
  single.cell <- single.cell[((single.cell$site1 %in% c(onek.bulk.sites,somatic.sites))&(single.cell$site2 %in% c(onek.bulk.sites,somatic.sites))),]

  single.cell <- single.cell[single.cell$pair_id %in% bulk$pair_id,]
  rownames(single.cell) <- single.cell$pair_id
  single.cell$pair_id <- NULL
  rownames(bulk) <- bulk$pair_id
  bulk <- bulk[rownames(single.cell),]
  bulk$pair_id <- NULL

  bulk$phased.orientation <- NA
  ind <- site.frame[bulk$site1,"bulk.phasing"] == site.frame[bulk$site2,"bulk.phasing"]
  ind1 <- ind2 <- ind
  ind1[is.na(ind1)] <- F
  bulk$phased.orientation[ind1] <- "cis"

  ind2 <- !ind2
  ind2[is.na(ind2)] <- F
  bulk$phased.orientation[ind2] <- "trans"
  
  #onek.bulk to onek.bulk
  tmp <- site.frame[bulk$site1,"onek.bulk.het"] & site.frame[bulk$site2,"onek.bulk.het"]
  site.frame$bulk.onek.bulk.linked <- F
  site.frame$bulk.onek.bulk.linked[rownames(site.frame) %in% c(bulk$site1[tmp],bulk$site2[tmp])] <- T
  
  ########
  bulk.onek.bulk.site1.trans <- bulk[with(bulk,(VR > VV)) & tmp,]
  n <- nrow(bulk.onek.bulk.site1.trans)
  single.cell.onek.bulk.site1.trans <- single.cell[rownames(bulk.onek.bulk.site1.trans),]
  onek.bulk.site1.trans.info <- data.frame(site=bulk.onek.bulk.site1.trans[,"site1"],
                                           hc.bulk=bulk.onek.bulk.site1.trans[,"VR"],
                                           hc.single.cell=single.cell.onek.bulk.site1.trans[,"VR"],
                                           dc.single.cell=single.cell.onek.bulk.site1.trans[,"RR"] + single.cell.onek.bulk.site1.trans[,"VV"],
                                           to=bulk.onek.bulk.site1.trans[,"site2"],
                                           type=rep("1KG",n),
                                           orientation=rep("trans",n),
                                           flags=rep("PASS",n))
  ind <- with(bulk.onek.bulk.site1.trans,(RR != 0) | (VV != 0))
  onek.bulk.site1.trans.info$flags[ind] <- paste(onek.bulk.site1.trans.info$flags[ind],"bad_bulk_haplotypes",sep=";")
  ind <- single.cell[rownames(bulk.onek.bulk.site1.trans),"VR"] == 0
  onek.bulk.site1.trans.info$flags[ind] <- paste(onek.bulk.site1.trans.info$flags[ind],"sc_dropout",sep=";")
  ind <- with(bulk.onek.bulk.site1.trans,phased.orientation == "cis")
  onek.bulk.site1.trans.info$flags[ind] <- paste(onek.bulk.site1.trans.info$flags[ind],"bad_phasing",sep=";")
  onek.bulk.site1.trans.info$flags[onek.bulk.site1.trans.info$flags != "PASS"] <- gsub("PASS;","",onek.bulk.site1.trans.info$flags[onek.bulk.site1.trans.info$flags != "PASS"])
  
  ########
  bulk.onek.bulk.site2.trans <- bulk[with(bulk,(RV > VV)) & tmp,]
  n <- nrow(bulk.onek.bulk.site2.trans)
  single.cell.onek.bulk.site2.trans <- single.cell[rownames(bulk.onek.bulk.site2.trans),]
  
  onek.bulk.site2.trans.info <- data.frame(site=bulk.onek.bulk.site2.trans[,"site2"],
                                           hc.bulk=bulk.onek.bulk.site2.trans[,"RV"],
                                           hc.single.cell=single.cell.onek.bulk.site2.trans[,"RV"],
                                           dc.single.cell=single.cell.onek.bulk.site2.trans[,"RR"] + single.cell.onek.bulk.site2.trans[,"VV"],
                                           to=bulk.onek.bulk.site2.trans[,"site1"],
                                           type=rep("1KG",n),
                                           orientation=rep("trans",n),
                                           flags=rep("PASS",n))
  ind <- with(bulk.onek.bulk.site2.trans,(RR != 0) | (VV != 0))
  onek.bulk.site2.trans.info$flags[ind] <- paste(onek.bulk.site2.trans.info$flags[ind],"bad_bulk_haplotypes",sep=";")
  ind <- single.cell[rownames(bulk.onek.bulk.site2.trans),"RV"] == 0
  onek.bulk.site2.trans.info$flags[ind] <- paste(onek.bulk.site2.trans.info$flags[ind],"sc_dropout",sep=";")
  ind <- with(bulk.onek.bulk.site2.trans,phased.orientation == "cis")
  onek.bulk.site2.trans.info$flags[ind] <- paste(onek.bulk.site2.trans.info$flags[ind],"bad_phasing",sep=";")
  onek.bulk.site2.trans.info$flags[onek.bulk.site2.trans.info$flags != "PASS"] <- gsub("PASS;","",onek.bulk.site2.trans.info$flags[onek.bulk.site2.trans.info$flags != "PASS"])
  
  ########
  bulk.onek.bulk.site1.cis <- bulk[with(bulk,(VV > VR)) & tmp,]
  n <- nrow(bulk.onek.bulk.site1.cis)
  single.cell.onek.bulk.site1.cis <- single.cell[rownames(bulk.onek.bulk.site1.cis),]
  
  onek.bulk.site1.cis.info <- data.frame(site=bulk.onek.bulk.site1.cis[,"site1"],
                                         hc.bulk=bulk.onek.bulk.site1.cis[,"VV"],
                                         hc.single.cell=single.cell.onek.bulk.site1.cis[,"VV"],
                                         dc.single.cell=single.cell.onek.bulk.site1.cis[,"RV"] + single.cell.onek.bulk.site1.cis[,"VR"],
                                         to=bulk.onek.bulk.site1.cis[,"site2"],
                                         type=rep("1KG",n),
                                         orientation=rep("cis",n),
                                         flags=rep("PASS",n))
  
  ind <- with(bulk.onek.bulk.site1.cis,(RV != 0) | (VR != 0))
  onek.bulk.site1.cis.info$flags[ind] <- paste(onek.bulk.site1.cis.info$flags[ind],"bad_bulk_haplotypes",sep=";")
  ind <- single.cell[rownames(bulk.onek.bulk.site1.cis),"VV"] == 0
  onek.bulk.site1.cis.info$flags[ind] <- paste(onek.bulk.site1.cis.info$flags[ind],"sc_dropout",sep=";")
  ind <- with(bulk.onek.bulk.site1.cis,phased.orientation == "trans")
  onek.bulk.site1.cis.info$flags[ind] <- paste(onek.bulk.site1.cis.info$flags[ind],"bad_phasing",sep=";")
  onek.bulk.site1.cis.info$flags[onek.bulk.site1.cis.info$flags != "PASS"] <- gsub("PASS;","",onek.bulk.site1.cis.info$flags[onek.bulk.site1.cis.info$flags != "PASS"])
  
  ########
  bulk.onek.bulk.site2.cis <- bulk[with(bulk,(VV > RV)) & tmp,]
  n <- nrow(bulk.onek.bulk.site2.cis)
  single.cell.onek.bulk.site2.cis <- single.cell[rownames(bulk.onek.bulk.site2.cis),]
  
  onek.bulk.site2.cis.info <- data.frame(site=bulk.onek.bulk.site2.cis[,"site2"],
                                         hc.bulk=bulk.onek.bulk.site2.cis[,"VV"],
                                         hc.single.cell=single.cell.onek.bulk.site2.cis[,"VV"],
                                         dc.single.cell=single.cell.onek.bulk.site2.cis[,"RV"] + single.cell.onek.bulk.site2.cis[,"VR"],
                                         to=bulk.onek.bulk.site2.cis[,"site1"],
                                         type=rep("1KG",n),
                                         orientation=rep("cis",n),
                                         flags=rep("PASS",n))
  
  ind <- with(bulk.onek.bulk.site2.cis,(RV != 0) | (VR != 0))
  onek.bulk.site2.cis.info$flags[ind] <- paste(onek.bulk.site2.cis.info$flags[ind],"bad_bulk_haplotypes",sep=";")
  ind <- single.cell[rownames(bulk.onek.bulk.site2.cis),"VV"] == 0
  onek.bulk.site2.cis.info$flags[ind] <- paste(onek.bulk.site2.cis.info$flags[ind],"sc_dropout",sep=";")
  ind <- with(bulk.onek.bulk.site2.cis,phased.orientation == "trans")
  onek.bulk.site2.cis.info$flags[ind] <- paste(onek.bulk.site2.cis.info$flags[ind],"bad_phasing",sep=";")
  onek.bulk.site2.cis.info$flags[onek.bulk.site2.cis.info$flags != "PASS"] <- gsub("PASS;","",onek.bulk.site2.cis.info$flags[onek.bulk.site2.cis.info$flags != "PASS"])
  
  
  combined.onek <- rbind(onek.bulk.site1.cis.info,onek.bulk.site2.cis.info,onek.bulk.site1.trans.info,onek.bulk.site2.trans.info)
  flags <- sapply(split(combined.onek$flags,combined.onek$site),function(x){paste(x,collapse=", ")})
  site.frame[names(flags),"onek.flags"] <- flags
  combined.onek$flags <- NULL
  
  ########
  tmp <- (site.frame[single.cell$site1,"somatic"]) & (site.frame[single.cell$site2,"onek.bulk.het"]) & with(bulk,(VR == 0) & (VV == 0))
  ########
  bulk.somatic.site1.trans <- bulk[(single.cell$VR > single.cell$VV) & tmp,]
  single.cell.somatic.site1.trans <- single.cell[rownames(bulk.somatic.site1.trans),]
  n <- nrow(single.cell.somatic.site1.trans)
  somatic.site1.trans.info <- data.frame(site=bulk.somatic.site1.trans[,"site1"],
                                         hc.bulk=bulk.somatic.site1.trans[,"RR"],
                                         hc.single.cell=single.cell.somatic.site1.trans[,"VR"],
                                         dc.single.cell=single.cell.somatic.site1.trans[,"RR"],
                                         to=bulk.somatic.site1.trans[,"site2"],
                                         type=rep("somatic",n),
                                         orientation=rep("trans",n))
  ########
  bulk.somatic.site1.cis <- bulk[(single.cell$VV > single.cell$VR) & tmp,]
  single.cell.somatic.site1.cis <- single.cell[rownames(bulk.somatic.site1.cis),]
  n <- nrow(single.cell.somatic.site1.cis)
  somatic.site1.cis.info <- data.frame(site=bulk.somatic.site1.cis[,"site1"],
                                       hc.bulk=bulk.somatic.site1.cis[,"RV"],
                                       hc.single.cell=single.cell.somatic.site1.cis[,"VV"],
                                       dc.single.cell=single.cell.somatic.site1.cis[,"RV"],
                                       to=bulk.somatic.site1.cis[,"site2"],
                                       type=rep("somatic",n),
                                       orientation=rep("cis",n))
  
  ########
  tmp <- (site.frame[single.cell$site1,"onek.bulk.het"]) & (site.frame[single.cell$site2,"somatic"]) & with(bulk,(RV == 0) & (VV == 0))
  ########
  bulk.somatic.site2.trans <- bulk[(single.cell$RV > single.cell$VV) & tmp,]
  single.cell.somatic.site2.trans <- single.cell[rownames(bulk.somatic.site2.trans),]
  n <- nrow(single.cell.somatic.site2.trans)
  somatic.site2.trans.info <- data.frame(site=bulk.somatic.site2.trans[,"site2"],
                                         hc.bulk=bulk.somatic.site2.trans[,"RR"],
                                         hc.single.cell=single.cell.somatic.site2.trans[,"RV"],
                                         dc.single.cell=single.cell.somatic.site2.trans[,"RR"],
                                         to=bulk.somatic.site2.trans[,"site1"],
                                         type=rep("somatic",n),
                                         orientation=rep("trans",n))
  
  bulk.somatic.site2.cis <- bulk[(single.cell$VV > single.cell$RV) & tmp,]
  single.cell.somatic.site2.cis <- single.cell[rownames(bulk.somatic.site2.cis),]
  n <- nrow(single.cell.somatic.site2.cis)
  somatic.site2.cis.info <- data.frame(site=bulk.somatic.site2.cis[,"site2"],
                                       hc.bulk=bulk.somatic.site2.cis[,"VR"],
                                       hc.single.cell=single.cell.somatic.site2.cis[,"VV"],
                                       dc.single.cell=single.cell.somatic.site2.cis[,"VR"],
                                       to=bulk.somatic.site2.cis[,"site1"],
                                       type=rep("somatic",n),
                                       orientation=rep("cis",n))
  
  ########
  tmp <- (site.frame[single.cell$site1,"somatic"]) & (site.frame[single.cell$site2,"somatic"]) & with(bulk,(RV == 0) & (VR == 0) & (VV == 0))
  ########
  bulk.somatic.somatic.trans <- bulk[with(single.cell,((RV + VR)/2 > VV) & (VR > 0) & (RV > 0))  & tmp,]
  single.cell.somatic.somatic.trans <- single.cell[rownames(bulk.somatic.somatic.trans),]
  n <- nrow(single.cell.somatic.somatic.trans)
  som.som.trans.info <- data.frame(site=c(bulk.somatic.somatic.trans[,"site1"],bulk.somatic.somatic.trans[,"site2"]),
                                   hc.bulk=c(bulk.somatic.somatic.trans[,"RR"],bulk.somatic.somatic.trans[,"RR"]),
                                   hc.single.cell=c(single.cell.somatic.somatic.trans[,"VR"],single.cell.somatic.somatic.trans[,"RV"]),
                                   dc.single.cell=c(single.cell.somatic.somatic.trans[,"RR"],single.cell.somatic.somatic.trans[,"RR"]),
                                   to=c(bulk.somatic.somatic.trans[,"site2"],bulk.somatic.somatic.trans[,"site1"]),
                                   type=c(rep("somatic-somatic",n),rep("somatic-somatic",n)),
                                   orientation=rep("trans",n))
  
  bulk.somatic.somatic.cis <- bulk[with(single.cell,((RV + VR)/2 < VV) & (VV > 0)) & tmp,]
  single.cell.somatic.somatic.cis <- single.cell[rownames(bulk.somatic.somatic.cis),]
  n <- nrow(single.cell.somatic.somatic.cis)
  som.som.cis.info <- data.frame(site=c(bulk.somatic.somatic.cis[,"site1"],bulk.somatic.somatic.cis[,"site2"]),
                                 hc.bulk=c(bulk.somatic.somatic.cis[,"RR"],bulk.somatic.somatic.cis[,"RR"]),
                                 hc.single.cell=c(single.cell.somatic.somatic.cis[,"VR"],single.cell.somatic.somatic.cis[,"RV"]),
                                 dc.single.cell=c(single.cell.somatic.somatic.cis[,"RR"],single.cell.somatic.somatic.cis[,"RR"]),
                                 to=c(bulk.somatic.somatic.cis[,"site2"],bulk.somatic.somatic.cis[,"site1"]),
                                 type=c(rep("somatic-somatic",n),rep("somatic-somatic",n)),
                                 orientation=rep("cis",n))
  
  combined.somatic <-  rbind(somatic.site1.cis.info,somatic.site2.cis.info,somatic.site1.trans.info,somatic.site2.trans.info)
  combined.somatic.somatic <- rbind(som.som.trans.info,som.som.cis.info)
  site.frame$somatic.proper <- rownames(site.frame) %in% combined.somatic$site
  
  data <- rbind(combined.onek,combined.somatic,combined.somatic.somatic)
  data$distance <- site.frame[data$site,"pos"] - site.frame[data$to,"pos"]
  data$site.pop.ref.freq <- site.frame[data$site,"pop.ref.freq"]
  data$to.pop.ref.freq <- site.frame[data$to,"pop.ref.freq"]
  
  data <- data[order(data$site,data$to),]
  rownames(data) <- NULL
  save(data,file=paste("data.",bulk.config$name,".rda",sep=""))
  
  tmp <- data[!(data$type == "somatic-somatic"),]
  tmp <- split(tmp,tmp$site)
  stat <- sapply(tmp,function(x){
    tmpstat <- suppressWarnings(max(apply(x[x$dc.single.cell == 0,c("hc.single.cell","hc.bulk")],1,min)))
    if(is.infinite(tmpstat)) {
      tmpstat <- -max(apply(x[,c("dc.single.cell","hc.bulk","hc.single.cell")],1,min))}
    return(tmpstat)})
  
  
  to <- sapply(tmp,function(x){
    tmpto <- suppressWarnings(x$to[which.max(apply(x[x$dc.single.cell == 0,c("hc.single.cell","hc.bulk")],1,min))])
    if(length(tmpto) == 0) {
      tmpto <- x$to[which.max(apply(x[,c("dc.single.cell","hc.bulk","hc.single.cell")],1,min))]
    }
    return(tmpto)})
  
  dist <- sapply(tmp,function(x){
    tmpdist <- suppressWarnings(x$distance[which.max(apply(x[x$dc.single.cell == 0,c("hc.single.cell","hc.bulk")],1,min))])
    if(length(tmpdist) == 0) {
      tmpdist <- x$distance[which.max(apply(x[,c("dc.single.cell","hc.bulk","hc.single.cell")],1,min))]
    }
    return(tmpdist)
  })
  
  orient <- sapply(tmp,function(x){
    tmporient <- suppressWarnings(x$orientation[which.max(apply(x[x$dc.single.cell == 0,c("hc.single.cell","hc.bulk")],1,min))])
    if(length(tmporient) == 0) {
      tmporient <- x$orientation[which.max(apply(x[,c("dc.single.cell","hc.bulk","hc.single.cell")],1,min))]
    }
    return(tmporient)
  })
  s <- sapply(tmp,function(x){x$site[1]})
  
  site.frame[,"stat"] <- NA
  site.frame[s,"stat"] <- stat
  site.frame[,"stat.to"] <- NA
  site.frame[s,"stat.to"] <- to
  site.frame[,"stat.distance"] <- NA
  site.frame[s,"stat.distance"] <- dist
  site.frame[,"stat.orientation"] <- NA
  site.frame[s,"stat.orientation"] <- orient
  
  #get pileup depths
  tmp <- read.table(paste(bulk.config$analysis_path,"/",chromosome,"/alt-counts.txt",sep=""),header=T,sep="\t")
  rownames(tmp) <- tmp$id
  tmp$id <- NULL
  tmp <- tmp[rownames(site.frame),]
  
  germline.stat <- sapply(which(!site.frame$somatic),function(x){
    return(sum(tmp[x,setdiff(c("A","C","G","T"),c(site.frame$alt[x],site.frame$ref[x]))]))
  })
  somatic.stat <- sapply(which(site.frame$somatic),function(x){
    return(sum(tmp[x,site.frame$alt[x]]))
  })
  site.frame$bulk.alt.count <- 0
  site.frame$bulk.alt.count[!site.frame$somatic] <- germline.stat
  site.frame$bulk.alt.count[site.frame$somatic] <- somatic.stat
  
  mutations <- site.frame
  linkage.parsed <- single.cell
  save(mutations,file=paste("mutations.",bulk.config$name,".rda",sep=""))
  save(linkage.parsed,file=paste("linkage.parsed.",bulk.config$name,".rda",sep=""))
  setwd(owd)
}

if(cmd == "collect_compare") {
  owd <- getwd()
  bulk.config <- read.config(args[3])
  setwd(config$analysis_path)
  chromosomes <- as.character(1:22)
  if(config$gender == "female") {
    chromosomes <- c(chromosomes,"X")
  }
  mutations <- do.call(rbind,lapply(chromosomes,function(x){load(paste(x,"/mutations.",bulk.config$name,".rda",sep="")); return(mutations)}))
  tmp <- mutations[mutations$onek.bulk.het & (mutations$bulk.phasing != 0),]
  onek.sites <- data.frame(chromosome=tmp$chromosome,position1=tmp$pos-1,position2=tmp$pos,ref.alt.phase=paste(tmp$ref,tmp$alt,tmp$bulk.phasing,sep=";"))
  write.table(onek.sites,file=paste("onek.",bulk.config$name,".bed",sep=""),sep="\t",row.names = F,col.names = F,quote=F)
  
  linkage <- do.call(rbind,lapply(chromosomes,function(x){load(paste(x,"/data.",bulk.config$name,".rda",sep="")); return(data)}))
  save(linkage,file="linkage.rda")
  rownames(linkage) <- paste(linkage$site,linkage$to)
  
  names(mutations)[names(mutations) == "stat"] <- "composite_coverage"
  all.somatic <- mutations[mutations$somatic.proper & (mutations$bulk.alt.count == 0),]
  all.somatic$pair.id <- paste(rownames(all.somatic),all.somatic$stat.to)
  all.somatic$hc.bulk <- linkage[all.somatic$pair.id,"hc.bulk"]
  all.somatic$dc.single.cell <- linkage[all.somatic$pair.id,"dc.single.cell"]
  all.somatic$hc.single.cell <- linkage[all.somatic$pair.id,"hc.single.cell"]
  
  candidate.somatic <- all.somatic[all.somatic$composite_coverage >= 2,]
  save(mutations,file=paste("mutations.",bulk.config$name,".rda",sep=""))
  save(all.somatic,file=paste("all.somatic.",bulk.config$name,".rda",sep=""))
  save(candidate.somatic,file=paste("candidate.somatic.",bulk.config$name,".rda",sep=""))
  
  #summary
  mutations.tmp <- mutations
  onek <- sum(mutations$onek.bulk)
  mutations.tmp <- mutations.tmp[!is.na(mutations.tmp$composite_coverage),]
  mutations.tmp <- mutations.tmp[mutations.tmp$bulk.onek.bulk.linked | mutations.tmp$somatic.proper,]
  mutations.tmp <- mutations.tmp[!(mutations.tmp$composite_coverage == 0),]
  mutations.tmp$misphased <- mutations.tmp$composite_coverage < 0
  mutations.tmp$cell <- config$name
  mutations.tmp$bulk.analysis <- bulk.config$name
  mutations.tmp <- mutations.tmp[mutations.tmp$single.cell.alt > 0,]
  tmp <- table(mutations.tmp$misphased,mutations.tmp$somatic.proper)
  summary <- data.frame(analysis.name=config$name,bulk.analysis=bulk.config$name,onek=onek,germline.phased=tmp[1,1],germline.misphased=tmp[2,1],somatic.phased=tmp[1,2],somatic.misphased=tmp[2,2])
  save(summary,file=paste("summary.",bulk.config$name,".rda",sep=""))
}

if(cmd == "power_jobs") {
  owd <- getwd()
  bulk.config <- read.config(args[3])
  setwd(config$analysis_path)
  chromosomes <- as.character(1:22)
  if(config$gender == "female") {
    chromosomes <- c(chromosomes,"X")
  }
  jobs <- unlist(lapply(chromosomes,function(x){list.files(paste(x,"/jobs",sep=""),full.names = T)}))
  job.dir <- paste("power.",bulk.config$name,"_job_scripts",sep="")
  suppressWarnings(dir.create(job.dir))
  js <- paste(job.dir,"/",gsub("/jobs/","_",jobs),".sh",sep="")
  for(i in seq_along(jobs)) {
    lines <- c("#!/bin/bash",paste("Rscript --vanilla $LIRA_DIR/scripts/main.R ",config$analysis_path,"/config.txt local_power_function ",bulk.config$analysis_path,"/config.txt ",config$analysis_path,"/",jobs[i],sep=""))
    writeLines(lines,con=js[i])
    system(paste("chmod +x ",js[i],sep=""))
  }
  setwd(owd)
}

if(cmd == "local_power_function") {
  wrapper <- function() {
    bulk.config <- read.config(args[3])
    work.dir <- args[4]
    owd <- getwd()
    setwd(work.dir)
    outfile <- "local.power.function.out.txt"
    logfile <- "local.power.function.log.txt"
    tag <- paste(".",bulk.config$name,sep="")
    bulk.het.bed <- paste(config$analysis_path,"/onek.",bulk.config$name,".bed",sep="")
    
    if(is.na(file.info(paste("sites.by.read.txt",sep=""))$size)) {
      powers <- numeric(0)
      save(powers,file=paste("powers",tag,".rda",sep=""))
      system(paste("touch powers.one",tag,".bed",sep=""))
      system(paste("touch powers.two",tag,".bed",sep=""))
      setwd(owd)
      return(0)
    }
    
    read.results <- function(fn) {
      tmp <- try(read.table(fn,comment.char="",quote="",sep="\t",header=F,colClasses=c("character","character","character","numeric","character","character","numeric")),silent=T)
      if(class(tmp) == "try-error") {
        return(data.frame(readname=character(0),id=character(0),chromosome=character(0),position=numeric(0),ref=character(0),alt=character(0),qual=numeric(0)))
      }
      names(tmp) <- c("readname","id","chromosome","position","ref","alt","qual")
      tmp$ref <- toupper(tmp$ref)
      tmp$targeted_alt <- gsub("([^;]*);([^;]*)","\\2",tmp$id)
      tmp$call <- NA
      tmp$call[tmp$alt == tmp$targeted_alt] <- "V"
      tmp$call[tmp$alt == tmp$ref] <- "R"
      tmp <- tmp[!is.na(tmp$call),]
      return(tmp)
    }
    
    results <- read.results("sites.by.read.txt")
    if(nrow(results) == 0) {
      powers <- numeric(0)
      save(powers,file=paste("powers",tag,".rda",sep=""))
      system(paste("touch powers.one",tag,".bed",sep=""))
      system(paste("touch powers.two",tag,".bed",sep=""))
      setwd(owd)
      return(0)
    }
    results <- results[order(results$chromosome,results$position,results$call),]
    id <- paste(results$chromosome,results$position,results$call,sep=";")
    rg <- split(results$readname,id)
    
    tmp.het.bed <- paste("onek.sites",tag,".bed",sep="")
    log(logfile,paste("bedtools intersect -wa -a ",bulk.het.bed," -b sites.bed > ",tmp.het.bed,sep=""))
    onek.sites <- try(read.table(tmp.het.bed,sep="\t",header=F),silent = T)
    if(class(onek.sites) == "try-error") {
      powers <- numeric(0)
      save(powers,file=paste("powers",tag,".rda",sep=""))
      system(paste("touch powers.one",tag,".bed",sep=""))
      system(paste("touch powers.two",tag,".bed",sep=""))
      setwd(owd)
      return(0)
    }
    
    onek.ids <- paste(onek.sites[,1],onek.sites[,3],sep=";")
    onek.ids <- c(paste(onek.ids,"R",sep=";"),paste(onek.ids,"V",sep=";"))
    var.allele <- gsub("(.);(.);(.)","\\3",onek.sites[,4])
    ref.allele <- rep("",length(var.allele))
    ref.allele[var.allele == "1"] <- "2"
    ref.allele[var.allele == "2"] <- "1"
    allele.key <- data.frame(id=onek.ids,allele=c(ref.allele,var.allele))
    rg <- rg[names(rg) %in% onek.ids]
    allele.key <- allele.key[allele.key$id %in% names(rg),]
    
    tmp <- log(logfile,paste("ls ",bulk.config$analysis_path,"/",onek.sites[1,1],"/jobs/*/sites.bed",sep=""),intern=T)
    select <- sapply(tmp,function(x){
      as.numeric(log(logfile,paste("bedtools intersect -wa -a ",tmp.het.bed," -b ",x," | wc -l ",sep=""),intern=T))
    })
    select <- names(select)[select != 0]
    if(length(select) > 1) {
      print("length > 1 in bulk")
    }
    rg.bulk <- unlist(recursive = F,x = lapply(select,function(x){
      tmp <- read.results(gsub("sites.bed","sites.by.read.txt",x))
      if(nrow(tmp) == 0) {
        return(NULL)
      }
      tmp <- tmp[order(tmp$chromosome,tmp$position,tmp$call),]
      id <- paste(tmp$chromosome,tmp$position,tmp$call,sep=";")
      rg.bulk <- split(tmp$readname,id)
      rg.bulk <- rg.bulk[unique(id)]
      rg.bulk <- rg.bulk[names(rg.bulk) %in% names(rg)]
      return(rg.bulk)
    }))
    rg <- rg[names(rg) %in% names(rg.bulk)]
    rg.bulk <- rg.bulk[names(rg)]
    
    regions <- read.table("regions.bed",sep="\t")
    big.region <- paste(regions[1,1],":",regions[1,2],"-",regions[nrow(regions),3],sep="")
    obj <- list()
    for(i in 1:2) {
      if(i == 1) {
        tmp.tag <- paste(".sc",tag,sep="")
        bam <- paste(config$analysis_path,"/reads.bam",sep="")
        rg.tmp <- rg
      } else {
        tmp.tag <- paste(".bulk",tag,sep="")
        bam <- paste(bulk.config$analysis_path,"/reads.bam",sep="")
        rg.tmp <- rg.bulk
      }
      suppressWarnings(dir.create(tmp.tag))
      setwd(tmp.tag)
      input.readnames <- unique(unlist(rg.tmp))
      write.table(file = "read_names.txt",x = input.readnames,quote = F,row.names = F,col.names = F)
      log(logfile,paste("samtools view -b ",bam," ",big.region,
                        " > tmp.bam",
                        " && picard FilterSamReads I=tmp.bam O=intersection.bam READ_LIST_FILE=read_names.txt FILTER=includeReadList SORT_ORDER=coordinate CREATE_INDEX=true WRITE_READS_FILES=false",
                        " && samtools view -H intersection.bam > header.txt",
                        " && samtools view intersection.bam > intersection.sam",
                        " && samtools view intersection.bam | cut -f 1 > rn.txt",
                        " && samtools view intersection.bam | cut -f 4 > positions.txt",sep=""))
      if(is.na(file.info("intersection.sam")$size)) {
        powers <- numeric(0)
        save(powers,file=paste("../powers",tag,".rda",sep=""))
        system(paste("touch ../powers.one",tag,".bed",sep=""))
        system(paste("touch ../powers.two",tag,".bed",sep=""))
        setwd(owd)
        return(0)
      }
      sam <- readLines(paste("intersection.sam",sep=""))
      header <- readLines(paste("header.txt",sep=""))
      reads <- readLines(paste("rn.txt",sep=""))
      positions <- as.numeric(readLines(paste("positions.txt",sep="")))
      split.sams <- lapply(rg.tmp,function(x){sam[reads %in% x][order(positions[reads %in% x])]})
      suppressWarnings(dir.create("power_bams"))
      res <- numeric(0)
      dumb <- lapply(seq_along(split.sams),function(x) {
        sam.tmp <- paste("power_bams/",x,tmp.tag,".sam",sep="")
        bam.tmp <- paste("power_bams/",x,tmp.tag,".bam",sep="")
        tmp <- c(header,split.sams[[x]])
        write.table(x = tmp,file = sam.tmp,quote=F,row.names = F,col.names = F)
        dumb <- system(paste("samtools view -b ",sam.tmp,
                             " > ",bam.tmp,
                             " && samtools index ",bam.tmp,
                             " && rm ",sam.tmp,sep=""))
      })
      setwd("..")
    }
    dumb <- lapply(seq_along(split.sams),function(x){
      tmp <- try(data.frame(do.call(rbind,strsplit(suppressWarnings(system(paste("samtools depth .sc",tag,"/power_bams/",x,".sc",tag,".bam .bulk",tag,"/power_bams/",x,".bulk",tag,".bam | awk '{if($3 < $4){print $1\"\\t\"$2-1\"\\t\"$2\"\\t\"$3} else {print $1\"\\t\"$2-1\"\\t\"$2\"\\t\"$4}}' | grep -v '\t0$' | grep -v '\t1$'",sep=""),intern=T)),"\t"))),silent=T)
      if(class(tmp) == "try-error") {
        return(data.frame(chromosome=character(0),
                          start=numeric(0),
                          stop=numeric(0),
                          stat=numeric(0)))
      }
      if(nrow(tmp) == 0) {
        return(data.frame(chromosome=character(0),
                          start=numeric(0),
                          stop=numeric(0),
                          stat=numeric(0)))
      }
      tmp[,1] <- as.character(tmp[,1])
      tmp[,2] <- as.numeric(tmp[,2])
      tmp[,3] <- as.numeric(tmp[,3])
      tmp[,4] <- as.numeric(tmp[,4])
      names(tmp) <- c("chromosome","start","stop","stat")
      system(paste("rm .sc",tag,"/power_bams/",x,".sc",tag,".bam &&",
                   "rm .bulk",tag,"/power_bams/",x,".bulk",tag,".bam",sep=""))
      return(tmp)
    })
    ind <- which(names(split.sams) %in% allele.key$id[allele.key$allele == 1])
    one <- do.call(rbind,dumb[ind])
    if(is.null(one)) {
      one <- data.frame(chromosome=character(0),
                        start=numeric(0),
                        stop=numeric(0),
                        stat=numeric(0))
    }
    one <- one[order(one$chromosome,one$start,one$stop),]
    two <- do.call(rbind,dumb[-ind])
    if(is.null(two)) {
      two <- data.frame(chromosome=character(0),
                        start=numeric(0),
                        stop=numeric(0),
                        stat=numeric(0))
    }
    two <- two[order(two$chromosome,two$start,two$stop),]
    
    mergefeature <- function(x) {
      local <- function(i) {
        tmp <- x[x$stat == i,c(1,2,3)]
        write.table(tmp,file="tmp.bed",sep="\t",quote=F,row.names = F,col.names = F)
        system("cat tmp.bed | bedtools merge > tmp.res")
        res <- read.table("tmp.res",sep="\t",header=F,colClasses = c("character","numeric","numeric"))
        names(res) <- c("chromosome","start","stop")
        res$stat <- i
        return(res)
      }
      res <- do.call(rbind,lapply(unique(x$stat),local))
      res <- res[order(res$chromosome,res$start,res$stop),]
      return(res)
    }
    
    if(nrow(one) != 0) {
      one <- aggregate(x=one$stat,by=list(chromosome=one$chromosome,
                                          start=one$start,
                                          stop=one$stop),FUN=max)
      names(one)[4] <- "stat"
      onebed <- mergefeature(one)
      write.table(onebed,paste("powers.one",tag,".bed",sep=""),sep="\t",col.names = F,row.names = F,quote=F)
    } else {
      system(paste("touch powers.one",tag,".bed",sep=""))
    }
    if(nrow(two) != 0) {
      two <- aggregate(x=two$stat,by=list(chromosome=two$chromosome,
                                          start=two$start,
                                          stop=two$stop),FUN=max)
      names(two)[4] <- "stat"
      twobed <- mergefeature(two)
      write.table(twobed,paste("powers.two",tag,".bed",sep=""),sep="\t",col.names = F,row.names = F,quote=F)
    } else {
      system(paste("touch powers.two",tag,".bed",sep=""))
    }
    powers <- table(c(one$stat,two$stat))
    save(powers,file=paste("powers",tag,".rda",sep=""))
    system(paste("rm -r .sc",tag," 2> /dev/null",sep=""))
    system(paste("rm -r .bulk",tag," 2> /dev/null",sep=""))
    system("rm tmp.bed 2> /dev/null && rm tmp.res 2> /dev/null")
    setwd(owd)
  }
  wrapper()
}

if(cmd == "check_power") {
  owd <- getwd()
  bulk.config <- read.config(args[3])
  setwd(config$analysis_path)
  scripts <- list.files(paste("power.",bulk.config$name,"_job_scripts",sep=""))
  prefix <- gsub(".sh","",paste(gsub("_","/jobs/",scripts),"/",sep=""))
  f1 <- paste(prefix,"powers.",bulk.config$name,".rda",sep="")
  f2 <- paste(prefix,"powers.one.",bulk.config$name,".bed",sep="")
  f3 <- paste(prefix,"powers.two.",bulk.config$name,".bed",sep="")
  test <- is.na(file.info(f1)$size) | is.na(file.info(f2)$size) | is.na(file.info(f3)$size)
  if(any(test)) {
    writeLines(paste(config$analysis_path,"/power.",bulk.config$name,"_job_scripts/",scripts[test],sep=""))
    stop("These jobs still undone.")
  } else { 
    print("Jobs finished.")
  }
}