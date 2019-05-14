
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
  
  for(n in c("SNPEFF","DBSNP","KGEN","EAGLE","EAGLE_HG19_REF","EAGLE_HG38_REF","analysis_path","reference","bam","vcf","phased_vcf")) {
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
  
  if(!("reference_identifier" %in% names(obj))) {
    obj$reference_identifier <- "hg19"
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

get.chromosomes <- function(config) {
  if(config$reference_identifier %in% c("GRCh37","hg38")) {
    chromosomes <- paste("chr",1:22,sep="")
    if(config$gender == "female") {
      chromosomes <- c(chromosomes,"chrX")
    }
  } else if(config$reference_identifier == "hg19") {
    chromosomes <- as.character(1:22)
    if(config$gender == "female") {
      chromosomes <- c(chromosomes,"X")
    }
  } else {
    stop("Cannot get chromosomes.")
  }
  if(!is.null(config$only_chromosomes)) {
    oc <- strsplit(config$only_chromosomes,",")[[1]]
    chromosomes <- chromosomes[chromosomes %in% oc]
  }
  return(chromosomes)
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
  
  test.list <- sapply(in.list,length) == 0
  test.obj <- in.list[!test.list][[1]]
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

annotate.and.phase.vcf <- function(chromosome) {
  out.log(paste("Annotating and phasing calls.vcf.gz",sep=""))
  out.log(paste("Output vcf: calls.id.vcf.gz",sep=""))
  hg19.convert <-  F
  if(config$phasing_software == "shapeit" & config$reference_identifier == "hg19") {
    hg19.convert <- T
  }
  grch37.convert <- F
  if(config$phasing_software == "eagle" & config$reference_identifier == "GRCh37") {
    grch37.convert <- T
  }
  
  if(config$phasing_software == "shapeit") {
    database <- global$DBSNP
  } else if(config$phasing_software == "eagle") {
    if(config$reference_identifier %in% c("hg19","GRCh37")) {
      db.dir <- global$EAGLE_HG19
    } else if(config$reference_identifier == "hg38") {
      db.dir <- global$EAGLE_HG38
    }
    database.tmp <- list.files(db.dir,pattern=paste("ALL.*chr",gsub("chr","",chromosome),"[_.].*vcf.gz$",sep=""),full.names = T)
    out.log.cmd(paste("bcftools view --force-samples  -s . ",database.tmp," -O z -o db.vcf.gz && tabix -f db.vcf.gz",sep=""))
    database <- "db.vcf.gz"
  }
  out.log.cmd(paste("java -jar ",global$SNPEFF,"/SnpSift.jar annotate -exists LIRA_GHET -tabix -name \"DBSNP_\" ",database," calls.vcf.gz > calls.id.vcf && bgzip -f calls.id.vcf &&  tabix -f calls.id.vcf.gz",sep=""))
  
  suppressWarnings(dir.create("phasing"))
  setwd("phasing")
  out.log("Make phasing input from population-polymorphic SNPs")
  out.log.cmd(paste("bcftools view -h ../calls.id.vcf.gz",
                    " | grep -e '##contig' -e '#CHROM' -e '##FORMAT=<ID=GT' -e '##FILTER' -e '##ALT' -e '##fileformat'",
                    " > phasing-input.vcf",sep=""))
  out.log.cmd(paste("bcftools view ../calls.id.vcf.gz",
                    " | awk 'BEGIN{OFS=\"\t\"}{if($8 ~ \"LIRA_GHET\"){$8 = \".\"; $9 = \"GT\"; print $0}}'",
                    " | tr ':' '\t'",
                    " | cut -f 1-10",
                    " | grep -e '#' -e '0/0' -e '0/1' -e '1/0' -e '1/1' -e '0|1' -e '1|1'",
                    " | sed 's#1/0#0/1#g'",
                    " | sed 's#0|1#0/1#g'",
                    " | sed 's#1|0#0/1#g'",
                    ifelse(hg19.convert," | sed 's#^#chr#g'",""),
                    " >> phasing-input.vcf",sep=""))
  if(config$phasing_software == "shapeit") {
    tmp <- paste("[.]chr",gsub("chr","",chromosome),"[.]",sep="")
    if(chromosome == "X") {
      tmp <- ".chrX_(non|NON)PAR."
      sex <- data.frame(sample=config$sample,sex=ifelse(config$gender == "male",1,2))
      write.table(x = sex,file = "sex.ped",quote = F,sep = "\t",row.names = F,col.names = F)
      add.args <- "--chrX --input-sex sex.ped"
    } else {
      add.args <- ""
    }
    l <- list.files(global$KGEN,recursive=T,pattern=tmp,full.names = T)
    legend <- l[grepl("legend",l)]
    hap <- l[grepl("hap",l)]
    map <- l[grepl("genetic_map",l)]
    sample <- list.files(global$KGEN,recursive=T,pattern="sample",full.names=T)
    
    out.log("Run shapeit check")
    out.log.cmd(paste("shapeit -check -V phasing-input.vcf -M ",map,
                      " --input-ref ",hap," ",legend," ",sample,
                      " --output-log shapeit-check.log",sep=""))
    out.log("Run shapeit")
    out.log.cmd(paste("shapeit -V phasing-input.vcf -M ",map,
                      " --input-ref ",hap," ",legend," ",sample,
                      " --exclude-snp shapeit-check.snp.strand.exclude -O phasing-output ",add.args,
                      " && shapeit -convert --input-haps phasing-output --output-vcf phasing-output.vcf",sep=""))
    if(hg19.convert) {
      out.log("Reformat phased vcf")
      out.log.cmd("sed -i 's#^chr##g' phasing-output.vcf && bgzip phasing-output.vcf && tabix -f phasing-output.vcf.gz")
    }
  } else if(config$phasing_software == "eagle") {
    out.log("Run eagle")
    out.log.cmd("bgzip -f phasing-input.vcf && tabix -f phasing-input.vcf.gz")
    vcfRef <- list.files(db.dir,pattern=paste(".chr",gsub("chr","",chromosome),"[_.].*bcf$",sep=""),full.names = T,recursive = T)
    if(config$reference_identifier %in% c("GRCh37","hg19")) {
      map <- list.files(paste(global$EAGLE,"/tables",sep=""),pattern="hg19",full.names = T)
    } else if(config$reference_identifier == "hg38") {
      map <- list.files(paste(global$EAGLE,"/tables",sep=""),pattern="hg38",full.names = T)
    }
    out.log.cmd(paste(global$EAGLE,"/eagle --geneticMapFile ",map," --vcfRef ",vcfRef," --vcfTarget phasing-input.vcf.gz --vcfOutFormat z --outPrefix phasing-output && tabix -f phasing-output.vcf.gz",sep=""))
  }
  
  out.log.cmd("bcftools view phasing-output.vcf.gz | awk '{print $1\";\"$2\";\"$4\";\"$5\"\t\"$10}' > phasing-output.txt")
  tmp <- read.table("phasing-output.txt")
  phasing <- data.frame(phase=tmp[,2])
  rownames(phasing) <- tmp[,1]
  save(phasing,file="phasing.rda")
  setwd("..")
}