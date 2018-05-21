
out.log.template <- function(line,tag) {
  write.table(x = paste(date(),": [",tag,"]\t",line,sep=""),file="log.txt",append=T,col.names=F,row.names=F,quote=F)
}

out.log.cmd <- function(cmd,tag) {
  out.log.template(cmd,"BASH")
  system(paste("/bin/bash -c ",shQuote(paste("exec 1> >(while read line; do echo \"$(date +'%a %h %d %H:%M:%S %Y'): [STDOUT]\t$line\" >> log.txt; done;); ",
                                             "exec 2> >(while read line; do echo \"$(date +'%a %h %d %H:%M:%S %Y'): [STDERR]\t$line\" >> log.txt; done;); ",
                                             cmd,sep="")),sep=""))
}

read.config.worker <- function(config.path) {
  config <- read.table(config.path,sep="\t",fill=T)
  config <- config[!(config[,1] == ""),]
  obj <- list()
  for(i in 1:nrow(config)) {
    obj[[i]] <- config[i,2]
  }
  names(obj) <- config[,1]
  
  for(n in c("SNPEFF","DBSNP","KGEN","analysis_path","reference","bam","vcf","phased_vcf")) {
    if(n %in% names(obj)) {
      tmp <- system(paste("echo ",obj[[n]],sep=""),intern=T)
      tmp.create <- F
      if(is.na(file.info(tmp)$size)) {
        tmp.create <- T
        system(paste("touch ",tmp,sep=""))
      }
      obj[[n]] <- normalizePath(tmp)
      if(tmp.create) {
        system(paste("rm ",tmp,sep=""))
      }
    }
  }
  
  if("GAP_REQUIREMENT" %in% names(obj)) {
    obj$GAP_REQUIREMENT <- as.numeric(obj$GAP_REQUIREMENT)
  }
  if("READS_TARGET" %in% names(obj)) {
    obj$READS_TARGET <- as.numeric(obj$READS_TARGET)
  }
  if("BATCH_SIZE" %in% names(obj)) {
    obj$BATCH_SIZE <- as.numeric(obj$BATCH_SIZE)
  }
  if("bulk" %in% names(obj)) {
    obj$bulk <- as.logical(obj$bulk)
  }
  if("MAX_DISTANCE_10X" %in% names(obj)) {
    obj$MAX_DISTANCE_10X <- as.numeric(obj$MAX_DISTANCE_10X)
  }
  
  obj$config <- cbind(names(obj),unlist(obj))
  return(obj)
}

read.config <- function(config.path) {
  cwd <- getwd()
  tmp <- strsplit(config.path,"/")[[1]]
  config.path.new <- tmp[length(tmp)]
  tmp <- paste(tmp[-length(tmp)],collapse="/")
  if(tmp != "") {
    setwd(tmp)
  }
  tst <- try(read.config.worker(config.path.new),silent=T)
  if(class(tst) == "try-error") {
    setwd(cwd)
    stop(tst)
  } else {
    setwd(cwd)
    return(tst)
  }
}

reverse.complement <- function(x) {
  x <- gsub("A","W",gsub("C","X",gsub("G","Y",gsub("T","Z",x))))
  x <- gsub("W","T",gsub("X","G",gsub("Y","C",gsub("Z","A",x))))
  x <- sapply(strsplit(x,""),function(x){
    x <- x[seq(from=length(x),by=-1,to=1)]
    return(paste(x,collapse=""))
  })
  return(x)
}

job.loop <- function(bash.commands,job.names) {
  out.log("Submitting jobs...")
  res <- submit.jobs(bash.commands,job.names)
  if(!res) {
    out.log("Parallel job submission failed.")
    stop("Parallel job submission failed.")
  }
  Sys.sleep(30)
  while(check.jobs(job.names)) {
    Sys.sleep(30)
  }
  out.log("All jobs finished.")
  return(T)
}

batcher <- function(bash.commands,batch.size) {
  if(length(bash.commands) == 0) {
    return(character(0))
  } else {
    len <- length(bash.commands)
    batches <- base::split(bash.commands[sample(1:len,size = len,replace = F)],rep(seq(from=1,by=1,to=ceiling(len/batch.size)),each=batch.size)[1:len])
    wrap <- sapply(batches,function(x){paste(x,collapse="; ")})
    names(wrap) <- paste(digest(list(wrap,date())),"_",seq_along(wrap),sep="")
    return(wrap)
  }
}

bam.region.selector <- function(in.bam, regions, out.bam) {
  tmp <- paste(".",digest(c(in.bam,regions,out.bam,date())),sep="")
  local <- function(){
    regions <- base::split(regions,rep(seq(from=1,by=1,to=ceiling(length(regions)/5000)),each=5000)[1:length(regions)])
    dir.create(tmp)
    
    for(i in seq_along(regions)) {
      system(paste("samtools view -b ",in.bam," ",paste(regions[[i]],collapse=" ")," > ",tmp,"/",i,".bam && samtools index ",tmp,"/",i,".bam",sep=""))
    }
    bamlist <- paste(tmp,"/bamfiles.txt",sep="")
    writeLines(paste(tmp,"/",seq_along(regions),".bam",sep=""),con=bamlist)
    system(paste("samtools merge -f -b ",bamlist," ",out.bam,"; samtools index ",out.bam,"; rm -r ",tmp,sep=""))
  }
  result <- try(local(),silent=T)
  if(class(result) == "try-error") {
    system(paste("rm -r ",tmp," 2> /dev/null",sep=""))
    stop(result)
  } else {
    return(0)
  }
}

combine.objects <- function(in.list) {
  name.sort <- function(in.names) {
    if(!any(is.na(as.numeric(in.names)))) {
      in.names <- as.numeric(in.names)
    }
    in.names <- as.character(in.names[order(in.names)])
    return(in.names)
  }
  
  test.obj <- in.list[[1]]
  if(is.null(names(test.obj))) {
    access <- names(dimnames(test.obj))
    n <- lapply(access,function(x) {
      name.sort(unique(unlist(lapply(in.list,function(y){dimnames(y)[x]}))))
    })
  } else {
    n <- name.sort(unique(unlist(lapply(in.list,names))))
    vec <- rep(0,length(n))
    names(vec) <- n
    for(i in seq_along(in.list)) {
      vec[names(in.list[[i]])] <- vec[names(in.list[[i]])] + in.list[[i]]
    }
    return(vec)
  }
}