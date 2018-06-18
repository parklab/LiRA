options(stringsAsFactors = FALSE,scipen=100)
library(stringr)
library(digest)
owd <- getwd()
lira.dir <- Sys.getenv("LIRA_DIR")

#get utility functions
source(paste(lira.dir,"/scripts/utils.R",sep=""))

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

if(!(cmd %in% c("null_command","joint","joint_subset"))) {
  suppressWarnings(dir.create(config$analysis_path))
  setwd(config$analysis_path)
} 

if(!(cmd == "null_command")) {
  out.log <- function(line,tag=paste("LiRA: ",cmd,sep="")){out.log.template(line,tag)}
} else {
  out.log <- function(line,tag=""){}
}

out.log("Start")

if(cmd == "setup") {
  setup(config)
}

if(cmd == "split") {
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
    }
  }
  else if(parallel & chromosome == "None") {
    chromosomes <- get.chromosomes(config)
    done <- file.exists(paste("progress/.split_",chromosomes,sep=""))
    if(!overwrite) {
      if(any(done)) {
        out.log(paste("Chromosome(s) ",paste(chromosomes[done],collapse=", ")," are already done. Excluding.",sep=""))
      }
      chromosomes <- chromosomes[!done]
    } else {
      out.log(paste("Overwriting chromosomes ",paste(chromosomes[done],collapse=", "),".",sep=""))
      sapply(chromosomes[done],function(x){out.log.cmd(paste("rm progress/.split_",x,sep=""))})
    }
    job.names <- paste(digest(list(config,date(),"split")),"_",chromosomes,sep="")
    cmds <- paste("lira split -c ",normalizePath("config.txt")," -m ",chromosomes,sep="")
    if(length(chromosomes) > 0) {
      res <- job.loop(cmds,job.names)
    }
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
  if(tot > 0) {
    res <- job.loop(batches,names(batches))
  }
}

if(cmd == "local_region_function") {
  work.dir <- args[3]
  local.region.function(config,work.dir)
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
    chromosomes <- get.chromosomes(config)
    done <- file.exists(paste("progress/.compare_",chromosomes,sep=""))
    if(!overwrite) {
      if(any(done)) {
        out.log(paste("Chromosome(s) ",paste(chromosomes[done],collapse=", ")," are already done. Excluding.",sep="")) 
      }
      chromosomes <- chromosomes[!done]
    } else {
      out.log(paste("Overwriting chromosomes ",paste(chromosomes[done],collapse=", "),".",sep=""))
      sapply(chromosomes[done],function(x){out.log.cmd(paste("rm progress/.compare_",x,sep=""))})
    }
    job.names <- paste(digest(list(config,date(),"compare")),"_",chromosomes,sep="")
    if(config$mode == "10X") {
      cmds <- paste("lira compare -s ",config$analysis_path,"/config.txt -m ",chromosomes,ifelse(wait," -w",""),sep="")
    } else {
      cmds <- paste("lira compare -s ",config$analysis_path,"/config.txt -b ",bulk.config$analysis_path,"/config.txt -m ",chromosomes,ifelse(overwrite," -o ",""),ifelse(wait," -w",""),sep="") 
    }
    if(length(chromosomes) > 0) {
      res <- job.loop(cmds,job.names)
    }
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
  if(length(scripts) > 0) {
    res <- job.loop(batches,names(batches))
  }
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
  bulk.config <- args[3]
  out.directory <- args[4]
  if(!grepl("^[/]",out.directory)) {
    out.directory  <- paste(owd,"/",out.directory,sep="")
  }
  use.uncertain.call <- as.logical(args[5])
  use.low.power <- as.logical(args[6])
  overwrite <- as.logical(args[7])
  single.cell.configs <- lapply(readLines(config.list),function(x) {
    if(!grepl("^[/]",x)) {
      x <- paste(owd,"/",x,sep="")
    }
    return(read.config(x))
  })
  if(!grepl("^[/]",bulk.config)) {
    bulk.config <- paste(owd,"/",bulk.config,sep="")
  }
  bulk.config <- read.config(bulk.config)
  joint(single.cell.configs,bulk.config,out.directory,use.uncertain.call,use.low.power,overwrite)
}

if(cmd == "joint_subset") {
  bulk.config <- args[3]
  work.dir <- args[4]
  if(!grepl("^[/]",bulk.config)) {
    bulk.config <- paste(owd,"/",bulk.config,sep="")
  }
  bulk.config <- read.config(bulk.config)
  joint.subset(config,bulk.config,work.dir)
}

if(!(cmd %in% c("null_command","joint"))) {
  setwd(config$analysis_path)
  out.log("Finish")
} 


if(cmd != "null_command") {
  setwd(owd)
}