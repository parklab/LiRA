options(stringsAsFactors = FALSE,scipen=100)
library(stringr)
library(digest)
owd <- getwd()
lira.dir <- Sys.getenv("LIRA_DIR")

#get utility functions
source(paste(lira.dir,"/scripts/utils.R",sep=""))

<<<<<<< HEAD
=======
read.config <- function(config.path) {
  config <- read.table(config.path,sep="\t")
  obj <- list()
  for(i in 1:nrow(config)) {
    obj[[i]] <- config[i,2]
  }
  names(obj) <- config[,1]
  
  suppressWarnings(dir.create(obj[["analysis_path"]]))
  for(n in c("analysis_path","reference","bam","bam_index","vcf","vcf_index","phased_vcf")) {
    if(n %in% names(obj)) {
      obj[[n]] <- normalizePath(obj[[n]])
    }
  }
  obj$config <- cbind(names(obj),unlist(obj))
  return(obj)
}
>>>>>>> 4e1a51a0d6161c13dda49becf3f4de3c988ea8be
args <- commandArgs(trailingOnly = T)
cmd <- args[2]

#read command
if(is.na(cmd)) {
  cmd <- "null_command"
}

#read config file
if(!is.na(args[1])) {
  if(!args[2] ==  "joint") {
    config <- read.config(args[1])
    if(is.null(config$mode)) {
      config$mode <- "lira"
    }
  } else {
    config <- list()
    config$mode <- "lira"
  }
} else {
  config <- list()
  config$mode <- "lira"
}

#read global settings
global <- read.config(paste(lira.dir,"/global-config.txt",sep=""))

#read parallelization functions
source(paste(lira.dir,"/scripts/",global$PARALLEL_SCRIPT,sep=""))

#read mode functions
if(config$mode == "lira") {
  source(paste(Sys.getenv("LIRA_DIR"),"/scripts/functions.R",sep=""))
} else if(config$mode == "10X") {
  source(paste(Sys.getenv("LIRA_DIR"),"/scripts/functions10x.R",sep=""))
}

if(!(cmd %in% c("null_command","joint"))) {
  suppressWarnings(dir.create(config$analysis_path))
  setwd(config$analysis_path)
} 

if(!(cmd == "null_command")) {
  out.log <- function(line,tag=paste("LiRA: ",cmd,sep="")){out.log.template(line,tag)}
} else {
  out.log <- function(line,tag=""){}
}

if(cmd == "setup") {
<<<<<<< HEAD
  out.log("Start setup")
  setup(config)
  out.log("Setup done")
=======
  owd <- getwd()
  suppressWarnings(dir.create(config$analysis_path))
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
>>>>>>> 4e1a51a0d6161c13dda49becf3f4de3c988ea8be
}

if(cmd == "split") {
  out.log("Start split")
  chromosome <- args[3]
  parallel <- as.logical(args[4])
  overwrite <- as.logical(args[5])
  
  if(!file.exists("progress/.setup")) {
    out.log(paste("Setup not done. Use 'lira setup ...' prior to 'lira split ...'."))
    stop("Error: setup not done.")
  }
  
  if(chromosome != "None" & !parallel) {
    progress.file <- paste("progress/.split_",chromosome,sep="")
    if(file.exists(progress.file)) {
      if(overwrite) {
        out.log(paste("Overwriting chromosome ",chromosome,".",sep=""))
        out.log.cmd(paste("rm ",progress.file,sep=""))
        split(config,chromosome)
      } else {
        message <- paste("Already done with chromosome ",chromosome,". Use -o to overwrite.",sep="")
        out.log(message)
        stop(message)
      }
    } else {
      split(config,chromosome)
      out.log("Split done")
    }
  }
  else if(parallel & chromosome == "None") {
    if(config$gender == "female") {
      chromosomes <- c(1:22,"X")
    } else {
      chromosomes <- as.character(1:22)
    }
    done <- file.exists(paste("progress/.split_",chromosomes,sep=""))
    if(!overwrite) {
      out.log(paste("Chromosomes ",paste(chromosomes[done],collapse=", ")," are already done. Excluding.",sep=""))
      chromosomes <- chromosomes[!done]
    } else {
      out.log(paste("Overwriting chromosomes ",paste(chromosomes[done],collapse=", "),".",sep=""))
      sapply(chromosomes[done],function(x){out.log.cmd(paste("rm progress/.split_",x,sep=""))})
    }
    job.names <- paste(digest(list(config,date(),"split")),"_",chromosomes,sep="")
    cmds <- paste("lira split -c ",normalizePath("config.txt")," -m ",chromosomes,sep="")
    res <- job.loop(cmds,job.names)
  } 
  else {
    stop("ambiguous arguments")
  }
}

if(cmd == "plink") {
  overwrite <- as.logical(args[3])
  batch.size <- suppressWarnings(as.numeric(args[4]))
  if(is.na(batch.size)) {
    batch.size <- global$BATCH_SIZE
  }
  scripts <- list.files("job_scripts")
  done.files <- paste("progress/.",gsub(".sh","",scripts),sep="")
  ind <- file.exists(done.files)
  completed <- done.files[ind]
  
  if(overwrite) {
    out.log("Running regions in parallel")
    out.log("Overwriting previous runs")
    dumb <- sapply(completed,function(x){system(paste("rm ",x,sep=""))})
  } else {
    scripts <- scripts[!ind]
  }
  tot <- length(scripts)
  cmds <- paste(getwd(),"/job_scripts/",scripts,sep="")
  batches <- batcher(cmds,batch.size)
  res <- job.loop(batches,names(batches))
}

<<<<<<< HEAD
if(cmd == "local_region_function") {
  work.dir <- args[3]
  local.region.function(config,work.dir)
=======
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
  
  data <- do.call(rbind,lapply(chromosomes,function(x){load(paste(x,"/data.",bulk.config$name,".rda",sep="")); return(data)}))
  rownames(data) <- paste(data$site,data$to)
  save(data,file=paste("data.",bulk.config$name,".rda",sep=""))
  
  names(mutations)[names(mutations) == "stat"] <- "composite_coverage"
  all.somatic <- mutations[mutations$somatic.proper & (mutations$bulk.alt.count == 0),]
  all.somatic$pair.id <- paste(rownames(all.somatic),all.somatic$stat.to)
  all.somatic$hc.bulk <- data[all.somatic$pair.id,"hc.bulk"]
  all.somatic$dc.single.cell <- data[all.somatic$pair.id,"dc.single.cell"]
  all.somatic$hc.single.cell <- data[all.somatic$pair.id,"hc.single.cell"]
  
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
>>>>>>> 4e1a51a0d6161c13dda49becf3f4de3c988ea8be
}

if(cmd == "local_barcode_function") {
  work.dir <- args[3]
  type <- args[4]
  if(is.na(type)) {
    type <- ""
  }
  local.barcode.function(config,work.dir,type)
}

if(cmd == "compare") {
  bulk.config <- args[3]
  chromosome <- args[4]
  parallel <- as.logical(args[5])
  overwrite <- as.logical(args[6])
  wait <- as.logical(args[7])
  if(config$mode == "lira") {
    if(!grepl("^[/]",bulk.config)) {
      bulk.config <- paste(owd,"/",bulk.config,sep="")
    }
    bulk.config <- read.config(bulk.config)
  } else if(config$mode == "10X") {
    if(bulk.config != "None") {
      warning("10x mode detected - bulk config input ignored") 
    }
  }
  if(chromosome != "None" & !parallel) {
    progress.file <- paste("progress/.compare_",chromosome,sep="")
    if(file.exists(progress.file)) {
      if(overwrite) {
        out.log(paste("Overwriting chromosome ",chromosome,".",sep=""))
        out.log.cmd(paste("rm ",progress.file,sep=""))
      } else {
        message <- paste("Already done with chromosome ",chromosome,". Use -o to overwrite.",sep="")
        out.log(message)
        stop(message)
      }
    }
    if(config$mode ==  "10X") {
      compare(config,chromosome,overwrite,wait)
    } else {
      compare(config,bulk.config,chromosome,overwrite,wait)
    }
    
  } 
  else if(parallel & chromosome == "None") {
    if(config$gender == "female") {
      chromosomes <- c(1:22,"X")
    } else {
      chromosomes <- as.character(1:22)
    }
    done <- file.exists(paste("progress/.compare_",chromosomes,sep=""))
    if(!overwrite) {
      out.log(paste("Chromosomes ",paste(chromosomes[done],collapse=", ")," are already done. Excluding.",sep=""))
      chromosomes <- chromosomes[!done]
    } else {
      out.log(paste("Overwriting chromosomes ",paste(chromosomes[done],collapse=", "),".",sep=""))
      sapply(chromosomes[done],function(x){out.log.cmd(paste("rm progress/.compare_",x,sep=""))})
    }
<<<<<<< HEAD
    job.names <- paste(digest(list(config,date(),"compare")),"_",chromosomes,sep="")
    if(config$mode == "10X") {
      cmds <- paste("lira compare -s ",config$analysis_path,"/config.txt -m ",chromosomes,ifelse(wait," -w",""),sep="")
    } else {
      cmds <- paste("lira compare -s ",config$analysis_path,"/config.txt -b ",bulk.config$analysis_path,"/config.txt -m ",chromosomes,ifelse(overwrite," -o ",""),ifelse(wait," -w",""),sep="") 
=======
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
                        " && java -jar $PICARD FilterSamReads I=tmp.bam O=intersection.bam READ_LIST_FILE=read_names.txt FILTER=includeReadList SORT_ORDER=coordinate CREATE_INDEX=true WRITE_READS_FILES=false",
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
>>>>>>> 4e1a51a0d6161c13dda49becf3f4de3c988ea8be
    }
    res <- job.loop(cmds,job.names)
  } 
  else {
    stop("ambiguous arguments")
  }
}

if(cmd == "local_power_function") {
  bulk.config <- args[3]
  work.dir <- args[4]
  if(config$mode == "lira") {
    if(!grepl("^[/]",bulk.config)) {
      bulk.config <- paste(owd,"/",bulk.config,sep="")
    }
    bulk.config <- read.config(bulk.config)
    local.power.function(config,bulk.config,work.dir)
  } else if(config$mode == "10X") {
    if(bulk.config != "None") {
      out.log("10x mode detected - bulk config input ignored")
    }
  }
}

if(cmd == "ppower") {
  out.log("Start ppower")
  bulk.config <- args[3]
  overwrite <- as.logical(args[4])
  batch.size <- suppressWarnings(as.numeric(args[5]))
  if(config$mode == "lira") {
    if(!grepl("^[/]",bulk.config)) {
      bulk.config <- paste(owd,"/",bulk.config,sep="")
    }
    bulk.config <- read.config(bulk.config)
    job.dir <- paste("power.",bulk.config$name,"_job_scripts",sep="")
    scripts <- list.files(job.dir)
    done.files <- paste("progress/.power_",bulk.config$name,"_",gsub(".sh","",scripts),sep="")
  } else if(config$mode == "10X") {
    job.dir <- "something"
    scripts <- "something"
    done.files <- "something"
  }
  if(is.na(batch.size)) {
    batch.size <- global$BATCH_SIZE
  }
  ind <- file.exists(done.files)
  completed <- done.files[ind]
  
  out.log("Running regions in parallel")
  if(overwrite) {
    out.log("Overwriting previous runs")
    dumb <- sapply(completed,function(x){system(paste("rm ",x,sep=""))})
  } else {
    scripts <- scripts[!ind]
  }
  cmds <- paste(job.dir,"/",scripts,sep="")
  batches <- batcher(cmds,batch.size)
  res <- job.loop(batches,names(batches))
  out.log("ppower done")
}

if(cmd == "varcall") {
  bulk.config <- args[3]
  overwrite <- as.logical(args[4])
  if(config$mode == "lira") {
    if(!grepl("^[/]",bulk.config)) {
      bulk.config <- paste(owd,"/",bulk.config,sep="")
    }
    bulk.config <- read.config(bulk.config)
  } else if(config$mode == "10X") {
    if(bulk.config != "None") {
      warning("10x mode detected - bulk config input ignored") 
    }
  }
  if(overwrite) {
    out.log("Overwriting previous run")
    out.log.cmd(paste("rm progress/.varcall_",bulk.config$name," 2> /dev/null",sep=""))
  }
  varcall(config,bulk.config,overwrite)
}

if(cmd == "joint") {
  config.list <- args[1]
  single.cell.configs <- lapply(readLines(config.list),read.config)
  bulk.config <- args[3]
  out.directory <- args[4]
  if(!grepl("^[/]",bulk.config)) {
    bulk.config <- paste(owd,"/",bulk.config,sep="")
  }
<<<<<<< HEAD
  bulk.config <- read.config(bulk.config)
  joint(single.cell.configs,bulk.config,out.directory)
}
if(cmd != "null_command") {
  setwd(owd)
=======
}

if(cmd == "collect_power") {
  owd <- getwd()
  bulk.config <- read.config(args[3])
  setwd(config$analysis_path)
  chromosomes <- as.character(1:22)
  if(config$gender == "female") {
    chromosomes <- c(chromosomes,"X")
  }
  f1 <- unlist(lapply(chromosomes,function(chr){
    paste(chr,"/jobs/",list.files(paste(chr,"/jobs",sep="")),"/powers.",bulk.config$name,".rda",sep="")
  }))
  f2 <- unlist(lapply(chromosomes,function(chr){
    paste(chr,"/jobs/",list.files(paste(chr,"/jobs",sep="")),"/powers.one.",bulk.config$name,".bed",sep="")
  }))
  f3 <- unlist(lapply(chromosomes,function(chr){
    paste(chr,"/jobs/",list.files(paste(chr,"/jobs",sep="")),"/powers.two.",bulk.config$name,".bed",sep="")
  }))
  powers <- lapply(f1,function(x){load(x); return(powers)})
  powers <- powers[!sapply(powers,length) == 0]
  m <- max(unlist(lapply(powers,function(x){as.numeric(names(x))})))
  aggre <- numeric(m)
  names(aggre) <- 1:m
  for(i in seq_along(powers)) {
    aggre[names(powers[[i]])] <- aggre[names(powers[[i]])] + powers[[i]]
  }
  powers <- aggre
  save(powers,file=paste("powers.",bulk.config$name,".rda",sep=""))
  
  p1 <- paste("powers.one.",bulk.config$name,".bed",sep="")
  p2 <- paste("powers.two.",bulk.config$name,".bed",sep="")
  system(paste("touch ",p1,sep=""))
  system(paste("touch ",p2,sep=""))
  for(i in seq_along(f2)) {
    system(paste("cat ",f2[i]," >> ",p1,".tmp",sep=""))
    system(paste("cat ",f3[i]," >> ",p2,".tmp",sep=""))
  }
  
  system(paste("cat ",p1[1]," | sort -V > ",p1," && rm ",p1,".tmp",sep=""))
  system(paste("cat ",p2[1]," | sort -V > ",p1," && rm ",p2,".tmp",sep=""))
  
  load(paste("data.",bulk.config$name,".rda",sep=""))
  load(paste("mutations.",bulk.config$name,".rda",sep=""))
  mutations <- mutations[mutations$onek.bulk.het,]
  bulk.alt.rate <- (sum(mutations$bulk.alt.count > 0)/nrow(mutations))/2
  data <- data[(data$hc.bulk >= 2),]
  data$pseudostat <- data$hc.single.cell
  data <- data[data$type == "1KG",]
  data$pseudostat[data$hc.single.cell > data$hc.bulk] <- data$hc.bulk[data$hc.single.cell > data$hc.bulk]
  data <- data[data$pseudostat >= 2,]
  data$dc <- data$dc.single.cell > 0
  key <- table(c(data$pseudostat,0),c(data$dc,T))
  key <- key[-1,]
  key.rate <- key[,2]/(key[,1] + key[,2])
  key.rate <- key.rate[names(key.rate) %in% names(powers)]
  
  powers <- powers[names(key.rate)]
  overall.rate.observed <- (1 - bulk.alt.rate) * (1 - key.rate)
  powers <- powers * overall.rate.observed
  save(powers,file=paste("powers.mod.",bulk.config$name,".rda",sep=""))
  setwd(owd)
}

if(cmd == "bootstrap_germline") {
  owd <- getwd()
  bulk.config <- read.config(args[3])
  setwd(config$analysis_path)
  
  load(paste("mutations.",bulk.config$name,".rda",sep=""))
  load(paste("candidate.somatic.",bulk.config$name,".rda",sep=""))
  sample.dir <- paste("samples.",bulk.config$name,sep="")
  suppressWarnings(dir.create(sample.dir))
  mutations <- mutations[mutations$onek.bulk.het,]
  mutations <- mutations[!is.na(mutations$composite_coverage),]
  mutations$site <- rownames(mutations)
  
  thresh <- 2
  candidate.somatic <- candidate.somatic[candidate.somatic$composite_coverage >= thresh,]
  mutations <- mutations[mutations$composite_coverage >= thresh,]
  local <- function(mut,som) {
    mut <- split(mut,mut$stat.distance)
    distances <- table(som$stat.distance)
    new <- list()
    target <- distances * 4
    for(n in names(distances)) {
      if(!(n %in% names(mut))) {
        tmp <- as.numeric(names(mut))
        k <- as.numeric(n)
        n.new <- as.character(tmp[abs(tmp - k) == min(abs(tmp - k))])
        if(length(n.new) > 1) {
          piece <- do.call(rbind,mut[n.new])
        } else {
          piece <- mut[[n.new]]
        }
        new[[n]] <- piece
        mut <- mut[-which(names(mut) %in% n.new)]
      } else {
        new[[n]] <- mut[[n]]
        mut <- mut[-which(names(mut) == n)]
      }
    }
    correct <- names(new)[sapply(new,nrow) < target]
    for(n in correct) {
      while(nrow(new[[n]]) < target[n]) {
        tmp <- as.numeric(names(mut)) - as.numeric(n)
        add <- names(mut)[tmp == min(tmp)]
        if(length(add) > 1) {
          piece <- do.call(rbind,mut[[add]])
        } else {
          piece <- mut[[add]]
        }
        new[[n]] <- rbind(new[[n]],piece)
        mut <- mut[-which(names(mut) %in% add)]
      }
    }
    return(new)
  }
  som <- candidate.somatic
  mut <- local(mutations,som)
  distances <- table(som$stat.distance)
  for(i in 1:100) {
    tmp <- lapply(seq_along(distances),function(x){
      n <- names(distances)[x]
      item <- mut[[n]][sample(1:nrow(mut[[n]]),size=distances[x],replace=F),]
    })
    tmp <- do.call(rbind,tmp)
    booty <- tmp
    save(booty,file=paste(sample.dir,"/",i,".",thresh,".rda",sep=""))
  }
  
  load(paste("candidate.somatic.",bulk.config$name,".rda",sep=""))
  candidate.somatic <- candidate.somatic[candidate.somatic$composite_coverage >= 2,]
  
  m <- max(candidate.somatic$composite_coverage)
  stat.dist <- (table(c(candidate.somatic$composite_coverage,2:m)) - 1)[as.character(2:m)]
  bootstrap <- do.call(rbind,lapply(1:100,function(x){load(paste("samples.",bulk.config$name,"/",x,".2.rda",sep="")); return((table(c(booty$composite_coverage,2:m)) - 1)[as.character(2:m)])}))
  
  obj <- list()
  obj$stat.dist <- stat.dist
  obj$af <- candidate.somatic[,c("composite_coverage","single.cell.alt","single.cell.ref")]
  obj$bootstrap <- colMeans(bootstrap)
  obj$bootaf <- lapply(1:100,function(x){load(paste("samples.",bulk.config$name,"/",x,".2.rda",sep="")); return(booty[,c("composite_coverage","single.cell.alt","single.cell.ref")])})
  save(obj,file=paste("samples.",bulk.config$name,"/obj.rda",sep=""))
  setwd(owd)
}

if(cmd == "call_ssnvs") {
  owd <- getwd()
  bulk.config <- read.config(args[3])
  setwd(config$analysis_path)
  load(paste("samples.",bulk.config$name,"/obj.rda",sep=""))
  obj[["af"]] <- NULL
  obj[["bootaf"]] <-  NULL
  load(paste("powers.mod.",bulk.config$name,".rda",sep=""))
  obj[["powers"]] <- powers
  n <- names(obj$stat.dist)
  obj$stat.dist <- as.numeric(obj$stat.dist)
  names(obj$stat.dist) <- n
  m <- max(as.numeric(names(obj[["stat.dist"]])))
  obj[["powers"]] <- obj[["powers"]][as.numeric(names(obj[["powers"]])) <= m &  as.numeric(names(obj[["powers"]])) >= 2]
  
  data <- data.frame(stat.val=as.numeric(names(obj$stat.dist)),observed.ssnv.count=obj$stat.dist,bootstrap=obj$bootstrap,power=obj$power)
  data$control.rate <- 1e9*((data$bootstrap + 0.5)/(data$power + 0.5))
  data$observed.ssnv.rate <- 1e9 * ((data$observed.ssnv.count+ 0.5)/(data$power + 0.5))
  
  tmp <- sapply(seq_along(data$observed.ssnv.rate),function(x) {
    qbeta(p = c(0.005,0.995),shape1 = 0.5 + data$observed.ssnv.count[x],shape2 = 0.5 + data$power[x] - data$observed.ssnv.count[x]) * 1e9
  })
  data$observed.ssnv.rate.lb <- tmp[1,]
  data$observed.ssnv.rate.ub <- tmp[2,]
  
  sample.observed.ssnv.rate <- function() {
    return(rbeta(n = nrow(data),shape1 = 0.5 + data$observed.ssnv.count,shape2 = 0.5 + data$power - data$observed.ssnv.count)*1e9)
  }
  
  alpha <- data$observed.ssnv.count + 0.5
  beta <- data$power + 0.5
  v <- (alpha * beta)/((alpha + beta)^2 * (alpha + beta + 1))
  v <- v/min(v)
  v <- 1/v
  data$error.rate <- (1/2)^(data$stat.val - 2)
  
  f <- function(params,rate=data$observed.ssnv.rate){
    a <- params[1]; b <- params[2]
    objective <- sum(v * ((a * data$error.rate + b * data$control.rate) - rate)^2)
    if((a < 0)|(b < 0)) {
      objective <- objective^2 
    }
    return(objective)
  }
  
  p <- c(data$observed.ssnv.rate[1],1)
  minimum <- 1000
  tmp <- list()
  class(tmp) <- "try-error"
  while(class(tmp) == "try-error") {
    tmp <- nlm(f=f,p=p,gradtol=1e-8,iterlim=500,steptol=1e-8,fscale = minimum,typsize=p)
  }
  nlm.fit <- tmp
  nlm.fit.bg <- do.call(rbind,lapply(1:100,function(x){
    tmp <- list()
    class(tmp) <- "try-error"
    while(class(tmp) == "try-error") {
      r <- sample.observed.ssnv.rate()
      tmp <- try(nlm(f=function(p){f(p,rate=r)},p=p,gradtol=1e-8,iterlim=500,steptol=1e-8,fscale=minimum,typsize=p)$estimate,silent=T)
    }
    return(tmp)
  }))
  obj <- list()
  obj$data <- data
  obj$nlm.fit <- nlm.fit
  obj$nlm.fit.bg <- nlm.fit.bg
  obj$error.rate.bg <- t(apply(nlm.fit.bg,1,function(x){x[1] * data$error.rate}))
  colnames(obj$error.rate.bg) <- obj$data$stat.val
  
  obj$control.rate.bg <- t(apply(nlm.fit.bg,1,function(x){x[2] * data$control.rate}))
  colnames(obj$control.rate.bg) <- obj$data$stat.val
  
  obj$fit.bg <- obj$error.rate.bg + obj$control.rate.bg
  obj$data$error.rate <- obj$data$error.rate * nlm.fit$estimate[1]
  obj$data$control.rate <- obj$data$control.rate * nlm.fit$estimate[2]
  obj$data$fit.rate <- obj$data$control.rate + obj$data$error.rate
  obj$data$fpr.estimate <- obj$data$error.rate/(obj$data$control.rate + obj$data$error.rate)
  
  cr.max <- max(obj$control.rate.bg[,1])
  cr.min <- min(obj$control.rate.bg[,1])
  
  signal <- round(obj$data$observed.ssnv.count * (1 - obj$data$fpr.estimate))
  signal <- sum(signal) - cumsum(signal) + signal
  noise <- round(obj$data$observed.ssnv.count * obj$data$fpr.estimate)
  noise <- sum(noise) - cumsum(noise) + noise
  names(signal) <- names(noise) <- obj$data$stat.val
  obj$composite_coverage_threshold <- as.numeric(names(which(noise/(signal + noise) <= 0.1)))[1]
  obj$somatic.rate <- obj$data$control.rate[1]
  obj$somatic.rate.lb <- cr.min
  obj$somatic.rate.ub <- cr.max
  obj$fpr <- (noise/(signal + noise))[which(data$stat.val == obj$composite_coverage_threshold)]
  eval(parse(text=paste("obj$rsomatic <- function(n){sample(c(",paste(obj$control.rate.bg[,1],collapse=","),"),replace=F,size=n)}",sep="")))
  print(obj)
  
  pdf(paste("rate-plot.",bulk.config$name,".pdf",sep=""))
  xlim <- c(1,(max(obj$data$stat.val)+0.5)*1.05)
  ylim <- c(0,max(obj$data$observed.ssnv.rate.ub)*1.2)
  plot(x = -10,y = -10,xlim=xlim,ylim=ylim,bty='n',xlab="",ylab="",xaxt='n',yaxt='n')
  points(x=obj$data$stat.val,y=obj$data$observed.ssnv.rate,pch=19)
  bar.w <- 0.2
  for(j in 1:nrow(obj$data)) {
    lines(x = c(obj$data$stat.val[j],obj$data$stat.val[j]),y=c(obj$data$observed.ssnv.rate.lb[j],obj$data$observed.ssnv.rate.ub[j]))
    lines(x = c(obj$data$stat.val[j] - bar.w, obj$data$stat.val[j] + bar.w),y=c(obj$data$observed.ssnv.rate.lb[j],obj$data$observed.ssnv.rate.lb[j]))
    lines(x = c(obj$data$stat.val[j] - bar.w, obj$data$stat.val[j] + bar.w),y=c(obj$data$observed.ssnv.rate.ub[j],obj$data$observed.ssnv.rate.ub[j]))
  }
  error.min <- apply(obj$error.rate.bg,2,min)
  error.max <- apply(obj$error.rate.bg,2,max)
  polygon(x = c(obj$data$stat.val,obj$data$stat.val[seq(from=nrow(obj$data),by=-1,to=1)]),
          y = c(error.min,error.max[seq(from=nrow(obj$data),by=-1,to=1)]),border=NA,col=rgb(1,0,0,alpha=0.25))
  lines(x = obj$data$stat.val,y=obj$data$error.rate,col="red",lwd=1)
  
  cr.min <- apply(obj$control.rate.bg,2,min)
  cr.max <- apply(obj$control.rate.bg,2,max)
  polygon(x = c(obj$data$stat.val,obj$data$stat.val[seq(from=nrow(obj$data),by=-1,to=1)]),
          y = c(cr.min,cr.max[seq(from=nrow(obj$data),by=-1,to=1)]),border=NA,col=rgb(0,0,1,alpha=0.25))
  lines(x = obj$data$stat.val,y=obj$data$control.rate,col="blue",lwd=1)
   
  fit.min <- apply(obj$fit.bg,2,min)
  fit.max <- apply(obj$fit.bg,2,max)
  polygon(x = c(obj$data$stat.val,obj$data$stat.val[seq(from=nrow(obj$data),by=-1,to=1)]),
          y = c(fit.min,fit.max[seq(from=nrow(obj$data),by=-1,to=1)]),border=NA,col=rgb(0,1,0,alpha=0.25))
  lines(x = obj$data$stat.val,y=obj$data$fit.rate,col="green",lwd=1)
   
  ytick.at <- seq(from=0,by=100,to=max(ylim))
  ytick <- as.character(ytick.at)
  axis(side = 2,at = ytick.at,labels = ytick,las=2,cex.axis=1.4)

  xtick.at <- seq(from=max(obj$data$stat.val),by=-2,to=2)[nrow(obj$data):1]
  xtick <- as.character(xtick.at)
  axis(side = 1,at = xtick.at,labels = xtick,cex.axis=1.4)
  axis(side = 1,tick = T,lwd.ticks=0,label=c("",""),at = c(min(obj$data$stat.val),max(obj$data$stat.val)))
   
  mtext(text = paste(config$name," vs. ",bulk.config$name,sep=""),side=3)
  points(x = 1.5,y=obj$somatic.rate,pch=19,col="purple")
  lines(x = c(1.5,1.5),y=c(obj$somatic.rate.lb,obj$somatic.rate.ub),col="purple",lwd=1)
  lines(x = c(1.5 - bar.w,1.5 + bar.w),y = c(obj$somatic.rate.lb,obj$somatic.rate.lb),col="purple",lwd=1)
  lines(x = c(1.5 - bar.w,1.5 + bar.w),y = c(obj$somatic.rate.ub,obj$somatic.rate.ub),col="purple",lwd=1)
  lines(x = c(obj$composite_coverage_threshold,obj$composite_coverage_threshold),y=c(0,ylim[2]),lty="dashed",lwd=1)
  par(xpd=NA)
  mtext(text = "sSNVs/GBp",side = 2,at = -ylim[2]/6,adj = c(0.5,0.5),line = 4)
  par(xpd=F)

  par(xpd=NA)
  mtext(text = "Composite coverage (CC)",side = 1, at = (xlim[1] + xlim[2])/2, adj=c(0.5,0.5),line=3)
  par(xpd=F)
  dev.off()
  
  save(obj,file=paste("call_ssnvs.rate_fit.",bulk.config$name,".rda",sep=""))
  load(paste("all.somatic.",bulk.config$name,".rda",sep=""))
  all.somatic$phred.quality <- NA
  all.somatic$fp.prob <- NA
  all.somatic$phred.quality[all.somatic$composite_coverage < 0] <- 0
  all.somatic$fp.prob[all.somatic$composite_coverage < 0] <- 1
  for(j in 1:nrow(obj$data)) {
    all.somatic$phred.quality[all.somatic$composite_coverage == obj$data$stat.val[j]] <- -10*log10(obj$data$fpr.estimate[j])
    all.somatic$fp.prob[all.somatic$composite_coverage == obj$data$stat.val[j]] <- obj$data$fpr.estimate[j]
  }
  all.somatic$status <- "no_power"
  all.somatic$status[all.somatic$composite_coverage < 0] <- "filtered_FP"
  all.somatic$status[all.somatic$composite_coverage >= 2] <- "uncertain_call"
  all.somatic$status[all.somatic$composite_coverage >= obj$composite_coverage_threshold] <- "call"
  ssnvs <- all.somatic
  save(ssnvs,file=paste("ssnvs.",bulk.config$name,".rda",sep=""))
>>>>>>> 4e1a51a0d6161c13dda49becf3f4de3c988ea8be
}