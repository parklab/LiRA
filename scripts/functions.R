setup <- function(config) {
  out.log("Start setup")
  
  #soft link bams
  res <- out.log.cmd(paste("ln -s ",config$bam," reads.bam",sep=""))
  index.possibility.1 <- paste(config$bam,".bai",sep="")
  index.possibility.2 <- gsub(".bam$",".bai",config$bam)
  if(file.exists(index.possibility.1)) {
    res <- out.log.cmd(paste("ln -s ",index.possibility.1," reads.bam.bai",sep=""))
  } else if(file.exists(index.possibility.2)) {
    res <- out.log.cmd(paste("ln -s ",index.possibility.2," reads.bam.bai",sep=""))
  } else {
    stop("Cannot find bam index.")
  }
  
  #soft link vcf
  res <- out.log.cmd(paste("ln -s ",config$vcf," input_calls.vcf.gz",sep=""))
  res <- out.log.cmd(paste("ln -s ",config$vcf,".tbi input_calls.vcf.gz.tbi",sep=""))
  write.table(x=config$config,file = "config.txt",sep="\t",col.names = F,row.names = F,quote=F)
  
  out.log.cmd("mkdir job_scripts")
  out.log.cmd("mkdir progress")
  out.log("Setup done.")
  
  #mark setup as complete
  out.log.cmd("touch progress/.setup")
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
    database.tmp <- list.files(db.dir,pattern=paste("ALL.*chr",gsub("chr","",chromosome),"_.*vcf.gz$",sep=""),full.names = T)
    out.log.cmd(paste("bcftools view --force-samples  -s . ",database.tmp," -O z -o db.vcf.gz && tabix -f db.vcf.gz",sep=""))
    database <- "db.vcf.gz"
  } else if (config$phasing_software == "crossbred_mouse") {
    database <- global$DBSNP
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
                    " | grep -e '#' -e '0/0' -e '0/1' -e '1/0' -e '1/1'",
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
    vcfRef <- list.files(db.dir,pattern=paste(".chr",gsub("chr","",chromosome),"_.*bcf$",sep=""),full.names = T,recursive = T)
    if(config$reference_identifier %in% c("GRCh37","hg19")) {
      map <- list.files(paste(global$EAGLE,"/tables",sep=""),pattern="hg19",full.names = T)
    } else if(config$reference_identifier == "hg38") {
      map <- list.files(paste(global$EAGLE,"/tables",sep=""),pattern="hg38",full.names = T)
    }
    out.log.cmd(paste(global$EAGLE,"/eagle --geneticMapFile ",map," --vcfRef ",vcfRef," --vcfTarget phasing-input.vcf.gz --vcfOutFormat z --outPrefix phasing-output && tabix -f phasing-output.vcf.gz",sep=""))
  } else if (config$phasing_software == "crossbred_mouse") {
    out.log("No phasing software to be run--VCF for heterozygous crossbred mouse is expected")
    # this needs a phasing-output.vcf.gz for the mouse; otherwise, we cannot perform bulk phasing!
    #out.log.cmd(paste(global$EAGLE,"/eagle --geneticMapFile ",map," --vcfRef ",vcfRef," --vcfTarget phasing-input.vcf.gz --vcfOutFormat z --outPrefix phasing-output && tabix -f phasing-output.vcf.gz",sep=""))
    # out.log.cmd(paste("shapeit -V phasing-input.vcf -M ",map,
    #                   " --input-ref ",hap," ",legend," ",sample,
    #                   " --exclude-snp shapeit-check.snp.strand.exclude -O phasing-output ",add.args,
    #                   " && shapeit -convert --input-haps phasing-output --output-vcf phasing-output.vcf",sep=""))
    out.log.cmd("cp phasing-input.vcf phasing-output.temp.vcf")
    out.log.cmd("sed 's#/#|#g' phasing-output.temp.vcf > phasing-output.vcf")
    out.log.cmd("bgzip phasing-output.vcf && tabix -f phasing-output.vcf.gz")
  }
}

split <- function(config,chromosome) {
  out.log(paste("Starting analysis of chromosome ",chromosome,sep=""))
  system(paste("rm -r ",chromosome," 2> /dev/null",sep=""))
  suppressWarnings(dir.create(chromosome))
  setwd(chromosome)
  
  #getting sites
  out.log("Getting sites (output: sites.bed)")
  
  #determine regions that should be analyzed together due to reads spanning multiple SNP sites
  out.log.cmd(paste("bcftools view -s ",config$sample," -i 'TYPE=\"snp\" && N_ALT=1' ../input_calls.vcf.gz ",chromosome," > calls.vcf && bgzip calls.vcf && tabix -f calls.vcf.gz",sep=""))
  out.log.cmd(paste("bcftools view calls.vcf.gz",
                    " | grep -v '#'",
                    " | cut -f 1,2,4,5",
                    " | awk '{print $1\"\t\"$2-1\"\t\"$2\"\t\"$3\";\"$4}'",
                    " > sites.bed",
                    " && bedtools merge -d ",global$GAP_REQUIREMENT," -i sites.bed",
                    " > regions.bed",
                    " && bedtools multicov -bams ../reads.bam -bed regions.bed",
                    " > coverage.bed",
                    sep=""))
  coverage <- read.table("coverage.bed",sep="\t",header=FALSE)
  cum.reads <- cumsum(coverage[,4])
  
  #assign regions to jobs targeting READS_TARGET
  run.assignment <- ceiling(cum.reads/global$READS_TARGET)
  jobs <- base::split(coverage[,1:3],run.assignment)
  dir.create("jobs")
  for(j in seq_along(jobs)) {
    prefix <- paste("jobs/",j,sep="")
    dir.create(prefix)
    write.table(jobs[[j]],file=paste(prefix,"/regions.bed",sep=""),sep="\t",col.names = F,row.names = F,quote=F)
    
    out.log("Extract relevant sites within the assigned regions.bed")
    out.log.cmd(paste("bedtools intersect -wa -a sites.bed -b ",prefix,"/regions.bed > ",prefix,"/sites.bed",sep=""))
    
    js <- paste("../job_scripts/",chromosome,"_",j,".sh",sep="")
    lines <- c("#!/bin/bash",paste("Rscript --vanilla $LIRA_DIR/scripts/main.R ",config$analysis_path,"/config.txt local_region_function ",config$analysis_path,"/",chromosome,"/jobs/",j,sep=""))
    out.log(paste("Write ",js,sep=""))
    writeLines(text=lines,con=js)
    out.log.cmd(paste("chmod +x ",js,sep=""))
  }
  
  #create SnpSift annotated vcf and phased vcf (for bulk)
  annotate.and.phase.vcf(chromosome)
  setwd(config$analysis_path)
  
  out.log(paste("Done with chromosome ",chromosome,sep=""))
  out.log.cmd(paste("touch progress/.split_",chromosome,sep=""))
}

local.region.function <- function(config,work.dir) {
  setwd(work.dir)
  
  wrapper <- function() {
    out.log(paste("Getting pileup by read (output: ",work.dir,"/sites.by.read.txt)",sep=""))
    out <- out.log.cmd(paste("$LIRA_DIR/scripts/linkage.py --bam ",config$analysis_path,"/reads.bam --fasta ",config$reference_file," --bed sites.bed > sites.by.read.txt",sep=""))
    
    if(as.numeric(system(paste("wc -l sites.by.read.txt | tr ' ' '\t' | cut -f 1",sep=""),intern=T)) == 0) {
      linkage <- data.frame(RR=numeric(0),RV=numeric(0),VR=numeric(0),VV=numeric(0))
      save(linkage,file=paste("linkage.rda",sep=""))
      return(0)
    }
    results <- read.table(paste("sites.by.read.txt",sep=""),comment.char="",quote="",sep="\t",header=F,colClasses=c("character","character","character","numeric","character","character","numeric"))
    
    names(results) <- c("readname","id","chromosome","position","ref","alt","qual")
    results$ref <- toupper(results$ref)
    results$targeted_alt <- gsub("([^;]*);([^;]*)","\\2",results$id)
    results$call <- NA
    results$call[results$alt == results$targeted_alt] <- "V"
    results$call[results$alt == results$ref] <- "R"
    out.log(paste("total unexpected calls: ",sum(is.na(results$call))," of ",length(results$call),sep=""))
    results$call[is.na(results$call)] <- "N"
    keep <- names(which(table(gsub("([^ ]*) ([^ ]*)","\\1",unique(paste(results$read,results$position)))) > 1))
    results <- results[results$readname %in% keep,]
    
    if(nrow(results) == 0) {
      linkage <- data.frame(RR=numeric(0),RV=numeric(0),VR=numeric(0),VV=numeric(0))
      save(linkage,file=paste("linkage.rda",sep=""))
      return(0)
    }
    
    #Save work
    save(results,file=paste("results.rda",sep=""))
    
    r <- unique(results$readname)
    batches <- base::split(r,rep(seq(from=1,to=ceiling(length(r)/1000)),each=1000)[1:length(r)])
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
      tmp <- base::split(tmp[,c("id","call")],tmp$readname)
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
    out.log("Done.")
  }
  wrapper()
  
  out.log.cmd("cat sites.bed | cut -f 1,3 > targets")
  
  get_alt_counts <- function() {
    cmd <- paste(c(paste("bcftools query -i 'TYPE=\"snp\" & N_ALT=1' -R targets -f '",c("%CHROM\\t%POS","%CHROM;%POS;%REF;%ALT"),"\\n' ",config$vcf," 2> /dev/null ",c(" | awk '{print $1\"\\t\"$2-1\"\\t\"$2}' > .tmp1"," > .tmp2"),sep="",collapse=" && "),
                   "echo \"id\tA\tC\tG\tT\" > alt-counts.txt",
                   paste("$LIRA_DIR/scripts/bulk_check.py --bam ",config$bam," --bed .tmp1 > .tmp3",sep=""),
                   "paste .tmp2 .tmp3 >> alt-counts.txt",
                   "rm .tmp1",
                   "rm .tmp2",
                   "rm .tmp3"),sep="",collapse=" && ")
    cmd <- paste(cmd,collapse=" && ")
    out.log.cmd(cmd)
  }
  if(config$bulk) {
    out.log("Getting alt counts")
    get_alt_counts()
  }
  
  print(config$phasing_software)
  get_vcf_info <- function() {
    if(config$phasing_software == "eagle") {
      caf.col <- "%INFO/DBSNP_AF"
    } else if(config$phasing_software == "shapeit") {
      caf.col <- "%INFO/DBSNP_CAF"
    } else if (config$phasing_software == "crossbred_mouse") {
      caf.col <- "%INFO/DBSNP_CAF"
    }
    cmds <- c(paste("bcftools query -i 'TYPE=\"snp\" & N_ALT=1' -s ",config$sample," -R targets -f '",c("%CHROM;%POS;%REF;%ALT\\t%ID",
                                                                                                        ifelse(config$bulk,caf.col,"."),
                                                                                                        "[%GT]",
                                                                                                        "[%AD]\\t[%GQ]",
                                                                                                        "%CHROM\\t%POS"),"\\n' ",paste(config$analysis_path,"/",system("cat sites.bed | cut -f 1 | uniq",intern=T),"/",ifelse(config$bulk,"calls.id.vcf.gz","calls.vcf.gz"),sep=""),
                    c("> .tmp1",
                      "| tr ',' '\\t' | cut -f 1 | sed 's#^[.]$#1#g' > .tmp2",
                      "> .tmp3",
                      " | tr ',.' '\\t0' > .tmp4",
                      " | awk '{print $1\"\\t\"$2-2\"\\t\"$2+1}' > .tmp5"),sep="",collapse=" && "),
              paste("bedtools getfasta -fi ",config$reference_file," -bed .tmp5 -fo - | grep -v '>' | tr 'actg' 'ACTG' > .tmp6",sep=""),
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
    # why is pop.ref.freq empty? For this cell, many of the targets do not have genotype information available. 
    cmd <- paste(cmds,collapse=" && ")
    out.log.cmd(cmd)
    site.frame <- read.table("vcf-info.txt",sep="\t",header=T)
    rownames(site.frame) <- site.frame$site
    site.frame$site <- NULL
    site.frame$chromosome <- gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\1",rownames(site.frame))
    site.frame$pos <- as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\2",rownames(site.frame)))
    site.frame$single.cell.ref <- site.frame$ref
    site.frame$single.cell.alt <- site.frame$alt
    site.frame$ref <- gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\3",rownames(site.frame))
    site.frame$alt <- gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\4",rownames(site.frame))
    site.frame$filter <- system(paste("bcftools query -i 'TYPE=\"snp\" & N_ALT=1' -R targets -f '%FILTER\\n' ",config$vcf," 2> /dev/null",sep=""),intern=T)
    site.frame <- site.frame[order(site.frame$pos),]
    names(site.frame)[names(site.frame) == "genotype"] <- "single.cell.gt"
    names(site.frame)[names(site.frame) == "gq"] <- "single.cell.gq"
    
    if(config$bulk) {
      site.frame$context <- paste(site.frame$context,">",site.frame$alt,sep="")
      change <- substr(site.frame$context,2,2) %in% c("G","A")
      site.frame$context[change] <- reverse.complement(site.frame$context[change])
      site.frame$pop.ref.freq[is.na(site.frame$pop.ref.freq)] <- 1
      
      names(site.frame)[names(site.frame) == "single.cell.ref"] <- "bulk.ref"
      names(site.frame)[names(site.frame) == "single.cell.alt"] <- "bulk.alt"
      names(site.frame)[names(site.frame) == "single.cell.gt"] <- "bulk.gt"
      names(site.frame)[names(site.frame) == "single.cell.gq"] <- "bulk.gq"
      #load phasing
      dumb <- system(paste("bcftools query -R targets -f '%CHROM;%POS;%REF;%ALT\t[%GT]\n' ../../phasing/phasing-output.vcf.gz 2> /dev/null | grep -e '0|1' -e '1|0' > phasing.txt",sep=""))
      phasing <- try(read.table("phasing.txt",sep="\t"),silent = T) # is phasing.txt empty? # yes! it appears that phasing.txt has not even been read!
      print(phasing)
      if(class(phasing) != "try-error") {
        vec <- numeric(nrow(phasing))
        names(vec) <- phasing[,1]
        vec[phasing[,2] == "1|0"] <- 1
        vec[phasing[,2] == "0|1"] <- 2
        site.frame$bulk.phasing <- 0
        site.frame[names(vec),"bulk.phasing"] <-  vec
      } else {
        site.frame$bulk.phasing <- 0
      }
      
      
      
      #load alt-counts
      tmp <- read.table("alt-counts.txt",header=T,sep="\t")
      site.frame$raw_count_A <- NA
      site.frame$raw_count_C <- NA
      site.frame$raw_count_G <- NA
      site.frame$raw_count_T <- NA
      site.frame[tmp$id,"raw_count_A"] <- tmp$A
      site.frame[tmp$id,"raw_count_C"] <- tmp$C
      site.frame[tmp$id,"raw_count_G"] <- tmp$G
      site.frame[tmp$id,"raw_count_T"] <- tmp$T
    }
    vcf.info <- site.frame
    save(vcf.info,file="vcf.info.rda")
  }
  out.log("Getting info from VCF")
  get_vcf_info()
  
  tmp <- strsplit(work.dir,"/")[[1]]
  n <- length(tmp)
  progress.file <- paste(config$analysis_path,"/progress/.",paste(tmp[n - 2],"_",tmp[n],sep=""),sep="")
  out.log.cmd(paste("touch ",progress.file,sep=""))
}

compare <- function(config, bulk.config, chromosome, overwrite, wait) {
  setwd(chromosome)
  if(overwrite) {
    system("rm -r compare 2> /dev/null")
  }
  suppressWarnings(dir.create("compare"))
  setwd("compare")
  
  process.local <- function(local.config,undone.message) {
    obj <- list()
    jobs.dir <- paste(local.config$analysis_path,"/",chromosome,"/jobs",sep="")
    j <- list.files(jobs.dir)
    done.files <- paste(local.config$analysis_path,"/progress/.",chromosome,"_",j,sep="")
    done <- file.exists(done.files)
    if(any(!done)) {
      if(!wait){
        out.log("Jobs not done:")
        sapply(done.files[!done],function(x){out.log(x)})
        stop(undone.message)
      } else {
        while(any(!done)) {
          out.log(paste("Waiting for ",bulk.config$name," plink to finish...",sep=""))
          out.log("Undone:")
          sapply(j[!done],function(x){out.log(paste(local.config$analysis_path,"/job_scripts/",chromosome,"_",x,".sh",sep=""))})
          done <- file.exists(done.files)
          Sys.sleep(120)
        }
      }
    }
    print(paste("Getting linkage information from",jobs.dir))
    linkage <- do.call(rbind,lapply(paste(jobs.dir,"/",j,"/linkage.rda",sep=""),function(x){load(x); return(linkage)}))
    linkage$site1 <- gsub("([^~]*)~([^~]*)","\\1",rownames(linkage))
    linkage$site2 <- gsub("([^~]*)~([^~]*)","\\2",rownames(linkage))
    linkage$pair_id <- rownames(linkage)
    rownames(linkage) <- NULL
    pos1 <- as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\2",linkage$site1))
    pos2 <- as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\2",linkage$site2))
    linkage <- linkage[order(pos1,pos2),]
    rownames(linkage) <- 1:nrow(linkage)
    obj$linkage <- linkage
    
    vcf.info <- do.call(rbind,lapply(paste(jobs.dir,"/",j,"/vcf.info.rda",sep=""),function(x){load(x); return(vcf.info)}))
    vcf.info <- vcf.info[order(vcf.info$pos),]
    obj$vcf.info <- vcf.info
    return(obj)
  }
  
  #Load single cell info
  out.log("Loading single cell info 1")
  print("sc 1")
  #out.log("ctrl1")
  #### NEED TO FIGURE OUT WHY THE JOBS EXIT AT THIS POINT
  #print(config)
  
  out.log("Loading bulk info 1")
  #print(config)
  print("bulk 1")
  #print(bulk.config)
  
  obj <- process.local(config,paste("Jobs still undone from plink.  See ",config$analysis_path,"/",chromosome,"/compare/log.txt",sep=""))
  out.log(obj)
  #out.log("ctrl2")
  single.cell.linkage <- obj$linkage
  out.log(single.cell.linkage)
  print("save single cell linkage")
  save(single.cell.linkage,file="single.cell.linkage.rda")
  
  single.cell.vcf.info <- obj$vcf.info
  print("save single cell vcf info")
  save(single.cell.vcf.info,file="single.cell.vcf.info.rda")

  #Load bulk info
  out.log("Loading bulk info")
  print("bulk config")
  obj <- process.local(bulk.config,paste("Bulk jobs still undone from plink.  See ",config$analysis_path,"/",chromosome,"/compare/log.txt",sep=""))
  bulk.linkage <- obj$linkage
  # print("save bulk linkage")
  # print(head(bulk.linkage))
  save(bulk.linkage,file="bulk.linkage.rda")
 
  bulk.vcf.info <- obj$vcf.info
  # print("save bulk vcf")
  save(bulk.vcf.info,file="bulk.vcf.info.rda")
  
  #make a site-specific table
  print("make a site-specific table")
  site.frame <- single.cell.vcf.info
  site.frame$id <- bulk.vcf.info$id
  site.frame$pop.ref.freq <- bulk.vcf.info$pop.ref.freq
  
  #pull bulk info
  print("pull bulk info---------")
  site.frame[,"bulk.gt"] <- bulk.vcf.info$bulk.gt
  site.frame[,"bulk.gq"] <- bulk.vcf.info$bulk.gq
  site.frame[,"bulk.ref"] <- bulk.vcf.info$bulk.ref
  site.frame[,"bulk.alt"] <- bulk.vcf.info$bulk.alt
  site.frame[,"bulk.phasing"] <-  bulk.vcf.info$bulk.phasing
  
  site.frame$bulk.phasing[site.frame$bulk.phasing == 0] <- NA

  #only look for linkage near bulk het sites in 1KG
  #somatic candidates: bulk.gt == "0/0" or bulk.gt == "./."
  print("only look for linkage near bulk het sites in 1KG")
  site.frame$in.linkage <- rownames(site.frame) %in% c(single.cell.linkage$site1,single.cell.linkage$site2)
  site.frame$onek.bulk.het <- (site.frame$pop.ref.freq != 1) & (site.frame$bulk.gt == "0/1")
  print(table(site.frame$pop.ref.freq))
  print(table(site.frame$bulk.gt))
  print(table(site.frame$onek.bulk.het))
  site.frame$somatic <- (site.frame$bulk.gt == "0/0")|(site.frame$bulk.gt == "./.")
  
  print("onekbulk")
  onek.bulk.sites <- rownames(site.frame)[site.frame$onek.bulk.het] # this turns up empty
  print(head(onek.bulk.sites))
  somatic.sites <- rownames(site.frame)[site.frame$somatic]
  print("somatic.sites")
  print(head(somatic.sites))
  single.cell.linkage <- single.cell.linkage[((single.cell.linkage$site1 %in% c(onek.bulk.sites,somatic.sites))&(single.cell.linkage$site2 %in% c(onek.bulk.sites,somatic.sites))),]
  print("single.cell.linkage")
  print(head(single.cell.linkage))
  
  print("sclinkage")
  single.cell.linkage <- single.cell.linkage[single.cell.linkage$pair_id %in% bulk.linkage$pair_id,]
  rownames(single.cell.linkage) <- single.cell.linkage$pair_id
  single.cell.linkage$pair_id <- NULL
  print(head(single.cell.linkage))
  rownames(bulk.linkage) <- bulk.linkage$pair_id
  bulk.linkage<- bulk.linkage[rownames(single.cell.linkage),]
  bulk.linkage$pair_id <- NULL
  print("bl")
  print(head(bulk.linkage))
  
  print("bulklinkage")
  bulk.linkage$phased.orientation <- ""
  ind <- site.frame[bulk.linkage$site1,"bulk.phasing"] == site.frame[bulk.linkage$site2,"bulk.phasing"]
  ind1 <- ind2 <- ind
  ind1[is.na(ind1)] <- F
  bulk.linkage$phased.orientation[ind1] <- "cis"
  print(head(site.frame)) 
  # why are all of the bulk.phasing entries at "NA"?
  # all of the pop.ref.freq are at 1!

  print("trans")
  ind2 <- !ind2
  ind2[is.na(ind2)] <- F
  bulk.linkage$phased.orientation[ind2] <- "trans"
  print(bulk.linkage) # this fills out fine
  
  haplotype.parser <- function(ind, type, site, orientation, hc.bulk.function, hc.single.cell.function, dc.single.cell.function,flag.functions=NULL) {
    bulk.tmp <- bulk.linkage[ind,]
    single.cell.tmp <- single.cell.linkage[ind,]
    n <- nrow(bulk.tmp)
    info <- data.frame(site=bulk.tmp[,ifelse(site == 1,"site1","site2")],
                       hc.bulk=hc.bulk.function(bulk.tmp,single.cell.tmp),
                       hc.single.cell=hc.single.cell.function(bulk.tmp,single.cell.tmp),
                       dc.single.cell=dc.single.cell.function(bulk.tmp,single.cell.tmp),
                       to=bulk.tmp[,ifelse(site == 1,"site2","site1")],
                       type=rep(type,n),
                       orientation=rep(orientation,n))
    if(!is.null(flag.functions)) {
      info$flags <- rep("PASS",n)
      for(i in seq_along(flag.functions)) {
        f.ind <- flag.functions[[i]](bulk.tmp,single.cell.tmp)
        info$flags[f.ind] <- paste(info$flags[f.ind],names(flag.functions)[i],sep=";")
      }
      info$flags[info$flags != "PASS"] <- gsub("PASS;","",info$flags[info$flags != "PASS"])
    }
    return(info)
  }
  
  #onek.bulk to onek.bulk
  print("site.frame")
  tmp <- site.frame[bulk.linkage$site1,"onek.bulk.het"] & site.frame[bulk.linkage$site2,"onek.bulk.het"]
  site.frame$bulk.onek.bulk.linked <- F
  site.frame$bulk.onek.bulk.linked[rownames(site.frame) %in% c(bulk.linkage$site1[tmp],bulk.linkage$site2[tmp])] <- T
  flag.functions <- list()
  # this command is determining if any of the single-cell variants are linked in the bulk sample
  # the site.frame is from the single-cell vcf info ('vcf-info.txt')
  
  
  #trans, site 1
  print("trans 1")
  flag.functions[["bad_bulk_haplotypes"]] <- function(b,s) {(b$RR != 0) | (b$VV != 0)}
  flag.functions[["sc_dropout"]] <- function(b,s){s$VR == 0}
  flag.functions[["bad_phasing"]] <- function(b,s){b$phased.orientation == "cis"}
  onek.bulk.site1.trans.info <- haplotype.parser(ind=with(bulk.linkage,VR > VV) & tmp,
                                                 type="1KG",
                                                 site=1,
                                                 orientation="trans",
                                                 hc.bulk.function=function(b,s){b$VR},
                                                 hc.single.cell.function=function(b,s){s$VR},
                                                 dc.single.cell.function=function(b,s){s$RR + s$VV},#should be just s$RR?
                                                 flag.functions=flag.functions)
  print(onek.bulk.site1.trans.info)
  
  #trans, site 2
  print("trans 2")
  flag.functions[["sc_dropout"]] <- function(b,s){s$RV == 0}
  onek.bulk.site2.trans.info <- haplotype.parser(ind=with(bulk.linkage,RV > VV) & tmp,
                                                 type="1KG",
                                                 site=2,
                                                 orientation="trans",
                                                 hc.bulk.function=function(b,s){b$RV},
                                                 hc.single.cell.function=function(b,s){s$RV},
                                                 dc.single.cell.function=function(b,s){s$RR + s$VV},#should be just s$RR?
                                                 flag.functions=flag.functions)
  
  #cis, site 1
  print("cis 1")
  flag.functions[["bad_bulk_haplotypes"]] <- function(b,s) {(b$RV != 0) | (b$VR != 0)}
  flag.functions[["sc_dropout"]] <- function(b,s) {s$VV == 0}
  flag.functions[["bad_phasing"]] <- function(b,s) {b$phased.orientation == "trans"}
  onek.bulk.site1.cis.info <- haplotype.parser(ind=with(bulk.linkage,VV > VR) & tmp,
                                                 type="1KG",
                                                 site=1,
                                                 orientation="cis",
                                                 hc.bulk.function=function(b,s){b$VV},
                                                 hc.single.cell.function=function(b,s){s$VV},
                                                 dc.single.cell.function=function(b,s){s$RV + s$VR},#should be just s$RV?
                                                 flag.functions=flag.functions)
  
  #cis, site 2
  print("cis 2")
  onek.bulk.site2.cis.info <- haplotype.parser(ind=with(bulk.linkage,VV > RV) & tmp,
                                               type="1KG",
                                               site=2,
                                               orientation="cis",
                                               hc.bulk.function=function(b,s){b$VV},
                                               hc.single.cell.function=function(b,s){s$VV},
                                               dc.single.cell.function=function(b,s){s$RV + s$VR},#should be just s$RV?
                                               flag.functions=flag.functions)
  
  print("onek")
  combined.onek <- rbind(onek.bulk.site1.cis.info,onek.bulk.site2.cis.info,onek.bulk.site1.trans.info,onek.bulk.site2.trans.info)
  #print(combined.onek) # combined.onek is empty
  # print(combined.onek$flags,combined.onek$site) #these are empty
  flags <- sapply(base::split(combined.onek$flags,combined.onek$site),function(x){paste(x,collapse=", ")}) # this list is empty, and it is throwing the error
  print("site.frames")
  ###The compare issue interferes with runtime below here
  print(flags)
  site.frame[names(flags),"onek.flags"] <- flags
  print(site.frame)
  combined.onek$flags <- NULL
  
  #site 1 somatic, site 2 onek.bulk
  print("tmp")
  tmp <- (site.frame[single.cell.linkage$site1,"somatic"]) & (site.frame[single.cell.linkage$site2,"onek.bulk.het"]) & with(bulk.linkage,(VR == 0) & (VV == 0))
  ###The compare issue interferes with runtime above here
  
  #trans
  print("trans site2 somatic")
  somatic.site1.trans.info <- haplotype.parser(ind=with(single.cell.linkage,VR > VV) & tmp,
                                               type="somatic",
                                               site=1,
                                               orientation="trans",
                                               hc.bulk.function=function(b,s){b$RR},
                                               hc.single.cell.function=function(b,s){s$VR},
                                               dc.single.cell.function=function(b,s){s$RR})
  #cis
  print("cis site2 somatic")
  somatic.site1.cis.info <- haplotype.parser(ind=with(single.cell.linkage,VV > VR) & tmp,
                                             type="somatic",
                                             site=1,
                                             orientation="cis",
                                             hc.bulk.function=function(b,s){b$RV},
                                             hc.single.cell.function=function(b,s){s$VV},
                                             dc.single.cell.function=function(b,s){s$RV})
  
  #site 1 onek.bulk, site2 somatic
  tmp <- (site.frame[single.cell.linkage$site1,"onek.bulk.het"]) & (site.frame[single.cell.linkage$site2,"somatic"]) & with(bulk.linkage,(RV == 0) & (VV == 0))
  
  #trans
  print("trans site2 somatic")
  somatic.site2.trans.info <- haplotype.parser(ind=with(single.cell.linkage,RV > VV) & tmp,
                                               type="somatic",
                                               site=2,
                                               orientation="trans",
                                               hc.bulk.function=function(b,s){b$RR},
                                               hc.single.cell.function=function(b,s){s$RV},
                                               dc.single.cell.function=function(b,s){s$RR})
  #cis
  print("cis site2 somatic")
  somatic.site2.cis.info <- haplotype.parser(ind=with(single.cell.linkage,VV > RV) & tmp,
                                             type="somatic",
                                             site=2,
                                             orientation="cis",
                                             hc.bulk.function=function(b,s){b$VR},
                                             hc.single.cell.function=function(b,s){s$VV},
                                             dc.single.cell.function=function(b,s){s$VR})
  #site 1 mosaic, site 2 onek.bulk
  tmp <- !(site.frame[single.cell.linkage$site1,"onek.bulk.het"]) & (site.frame[single.cell.linkage$site2,"onek.bulk.het"])
  
  #trans
  mosaic.site1.trans.info <- haplotype.parser(ind=with(single.cell.linkage,VR > VV) & with(bulk.linkage,(VR > VV) & (RR > VR) & (RV > VR)) & tmp,
                                              type="mosaic",
                                              site=1,
                                              orientation="trans",
                                              hc.bulk.function=function(b,s){b$RR},
                                              hc.single.cell.function=function(b,s){s$VR},
                                              dc.single.cell.function=function(b,s){s$RR})
  #cis
  mosaic.site1.cis.info <- haplotype.parser(ind=with(single.cell.linkage,VV > VR) & with(bulk.linkage,(VV > VR) & (RV > VV) & (RR > VV)) & tmp,
                                              type="mosaic",
                                              site=1,
                                              orientation="cis",
                                              hc.bulk.function=function(b,s){b$RV},
                                              hc.single.cell.function=function(b,s){s$VV},
                                              dc.single.cell.function=function(b,s){s$RV})
  
  #site 2 mosaic, site 1 onek.bulk
  tmp <- (site.frame[single.cell.linkage$site1,"onek.bulk.het"]) & !(site.frame[single.cell.linkage$site2,"onek.bulk.het"])
  
  #trans
  mosaic.site2.trans.info <- haplotype.parser(ind=with(single.cell.linkage,RV > VV) & with(bulk.linkage,(RV > VV) & (RR > RV) & (VR > RV)) & tmp,
                                              type="mosaic",
                                              site=2,
                                              orientation="trans",
                                              hc.bulk.function=function(b,s){b$RR},
                                              hc.single.cell.function=function(b,s){s$RV},
                                              dc.single.cell.function=function(b,s){s$RR})
  #cis
  mosaic.site2.cis.info <- haplotype.parser(ind=with(single.cell.linkage,VV > RV) & with(bulk.linkage,(VV > RV) & (VR > VV) & (RR > VV)) & tmp,
                                              type="mosaic",
                                              site=2,
                                              orientation="trans",
                                              hc.bulk.function=function(b,s){b$RV},
                                              hc.single.cell.function=function(b,s){s$VV},
                                              dc.single.cell.function=function(b,s){s$VR})
  
  combined.somatic <-  rbind(somatic.site1.cis.info,somatic.site2.cis.info,somatic.site1.trans.info,somatic.site2.trans.info)
  rownames(combined.somatic) <- paste(combined.somatic$site,combined.somatic$to)
  combined.mosaic <- rbind(mosaic.site1.cis.info,mosaic.site2.cis.info,mosaic.site1.trans.info,mosaic.site2.trans.info)
  rownames(combined.mosaic) <- paste(combined.mosaic$site,combined.mosaic$to)
  save(combined.somatic,file=paste("combined.somatic.",bulk.config$name,".rda",sep=""))
  save(combined.mosaic,file=paste("combined.mosaic.",bulk.config$name,".rda",sep=""))
  site.frame$somatic.proper <- rownames(site.frame) %in% combined.somatic$site
  site.frame$mosaic.proper <- rownames(site.frame) %in% combined.mosaic$site
  
  process.data <- function(data) {
    data$distance <- site.frame[data$site,"pos"] - site.frame[data$to,"pos"]
    data$stat <- data$hc.bulk
    data$stat[data$hc.single.cell < data$hc.bulk] <- data$hc.single.cell[data$hc.single.cell < data$hc.bulk]
    
    ind <- data$dc.single.cell > 0
    ind2 <- data$dc.single.cell < data$stat
    data$stat[ind & ind2] <- -data$dc.single.cell[ind & ind2]
    data$stat[ind & !ind2] <- -data$stat[ind & !ind2]
    
    tmp <- base::split(data,data$site)
    
    index <- sapply(tmp,function(x){
      if(any(x$stat > 0)) {
        index <- which.max(x$stat)
      } else {
        index <- which.min(x$stat)
      }
      return(index)
    })
    
    extract <- function(column){sapply(seq_along(tmp),function(x){tmp[[x]][index[x],column]})}
    info <- data.frame(stat=extract("stat"),to=extract("to"),distance=extract("distance"),
                       orientation=extract("orientation"),site=extract("site"),
                       pair.id=sapply(seq_along(tmp),function(x){rownames(tmp[[x]])[index[x]]}),
                       hc.bulk=extract("hc.bulk"),hc.single.cell=extract("hc.single.cell"),
                       dc.single.cell=extract("dc.single.cell"))
    if(ncol(info) == 0 & nrow(info) == 0) {
      info <- data.frame(stat=numeric(0),to=character(0),distance=numeric(0),orientation=character(0),site=character(0),pair.id=character(0),hc.bulk=numeric(0),hc.single.cell=numeric(0),dc.single.cell=numeric(0))
    }
    return(info)
  }
  
  data <- rbind(combined.onek,combined.somatic)
  rownames(data) <- paste(data$site,data$to)
  save(data,file=paste("data.",bulk.config$name,".rda",sep=""))
  obj <- process.data(data)
  
  for(c in colnames(obj)[!(colnames(obj) == "site")]) {
    site.frame[,c] <- NA
    site.frame[obj$site,c] <- obj[,c]
  }
  
  data.mosaic <- combined.mosaic
  rownames(data.mosaic) <- paste(data.mosaic$site,data.mosaic$to)
  obj.mosaic <- process.data(data.mosaic)
  colnames(obj.mosaic) <- paste("mosaic.",colnames(obj.mosaic),sep="")
  
  for(c in colnames(obj.mosaic)[!(colnames(obj.mosaic) == "mosaic.site")]) {

    site.frame[,c] <- NA
    site.frame[obj.mosaic$mosaic.site,c] <- obj.mosaic[,c]
  }
  
  #get pileup depths
  tmp <- bulk.vcf.info[,c("raw_count_A","raw_count_C","raw_count_G","raw_count_T")]
  names(tmp) <- c("A","C","G","T")
  site.frame$bulk.alt.count <- 0
  for(alt in names(tmp)) {
    ind <- (site.frame$onek.bulk.het & !(site.frame$alt == alt) & !(site.frame$ref == alt)) | ((!site.frame$onek.bulk.het) & site.frame$alt == alt)
    site.frame[ind,"bulk.alt.count"] <- site.frame[ind,"bulk.alt.count"] + tmp[ind,alt]
  }
  
  mutations <- site.frame
  save(mutations,file=paste("mutations.",bulk.config$name,".rda",sep=""))
  
  tmp <- mutations[mutations$onek.bulk.het & !is.na(mutations$bulk.phasing),]
  onek.sites <- data.frame(chromosome=tmp$chromosome,position1=tmp$pos-1,position2=tmp$pos,ref.alt.phase=paste(tmp$ref,tmp$alt,tmp$bulk.phasing,sep=";"))
  write.table(onek.sites,file=paste("onek.",bulk.config$name,".bed",sep=""),sep="\t",row.names = F,col.names = F,quote=F)

  all.somatic <- mutations[mutations$somatic.proper & (mutations$bulk.alt.count == 0),]
  all.mosaic <- mutations[mutations$mosaic.proper,]
  
  candidate.somatic <- all.somatic[all.somatic$stat >= 2,]
  save(all.somatic,file=paste("all.somatic.",bulk.config$name,".rda",sep=""))
  save(candidate.somatic,file=paste("candidate.somatic.",bulk.config$name,".rda",sep=""))
  save(all.mosaic,file=paste("all.mosaic.",bulk.config$name,".rda",sep=""))
  
  mutations.tmp <- mutations
  onek <- sum(mutations$onek.bulk)
  mutations.tmp <- mutations.tmp[!is.na(mutations.tmp$stat),]
  mutations.tmp <- mutations.tmp[mutations.tmp$bulk.onek.bulk.linked | mutations.tmp$somatic.proper,]
  mutations.tmp <- mutations.tmp[!(mutations.tmp$stat == 0),]
  mutations.tmp$misphased <- mutations.tmp$stat < 0
  mutations.tmp$cell <- config$name
  mutations.tmp$bulk.analysis <- bulk.config$name
  mutations.tmp <- mutations.tmp[mutations.tmp$single.cell.alt > 0,]
  tmp <- table(c(mutations.tmp$misphased,T,F,T,F),c(mutations.tmp$somatic.proper,T,T,F,F)) - 1
  summary <- data.frame(analysis.name=config$name,bulk.analysis=bulk.config$name,onek=onek,germline.phased=tmp[1,1],germline.misphased=tmp[2,1],somatic.phased=tmp[1,2],somatic.misphased=tmp[2,2])
  save(summary,file=paste("summary.",bulk.config$name,".rda",sep=""))
  setwd(config$analysis_path)
  
  jobs <- list.files(paste(chromosome,"/jobs",sep=""),full.names = T)
  job.dir <- paste("power.",bulk.config$name,"_job_scripts",sep="")
  suppressWarnings(dir.create(job.dir))
  js <- paste(job.dir,"/",gsub("/jobs/","_",jobs),".sh",sep="")
  for(i in seq_along(jobs)) {
    lines <- c("#!/bin/bash",paste("Rscript --vanilla $LIRA_DIR/scripts/main.R ",config$analysis_path,"/config.txt local_power_function ",bulk.config$analysis_path,"/config.txt ",config$analysis_path,"/",jobs[i],sep=""))
    writeLines(lines,con=js[i])
    system(paste("chmod +x ",js[i],sep=""))
  }
  
  out.log.cmd(paste("touch progress/.compare_",chromosome,sep=""))
}

local.power.function <- function(config,bulk.config,work.dir) {
  local <- function(){
    setwd(work.dir)
    tag <- paste(".",bulk.config$name,sep="")
    bulk.het.bed <- paste(work.dir,"/../../compare/onek.",bulk.config$name,".bed",sep="")
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
    rg <- base::split(results$readname,id)
    tmp.het.bed <- paste("onek.sites",tag,".bed",sep="")
    out.log.cmd(paste("bedtools intersect -wa -a ",bulk.het.bed," -b sites.bed > ",tmp.het.bed,sep=""))
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
    
    tmp <- system(paste("ls ",bulk.config$analysis_path,"/",onek.sites[1,1],"/jobs/*/sites.bed 2> /dev/null",sep=""),intern=T)
    select <- sapply(tmp,function(x){
      as.numeric(system(paste("bedtools intersect -wa -a ",tmp.het.bed," -b ",x," | wc -l ",sep=""),intern=T))
    })
    select <- names(select)[select != 0]
    if(length(select) > 1) {
      out.log("More than one corresponding job in bulk.")
    }
    rg.bulk <- unlist(recursive = F,x = lapply(select,function(x){
      tmp <- read.results(gsub("sites.bed","sites.by.read.txt",x))
      if(nrow(tmp) == 0) {
        return(NULL)
      }
      tmp <- tmp[order(tmp$chromosome,tmp$position,tmp$call),]
      id <- paste(tmp$chromosome,tmp$position,tmp$call,sep=";")
      rg.bulk <- base::split(tmp$readname,id)
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
      out.log.cmd(paste("samtools view -b ",bam," ",big.region,
                        " > tmp.bam",
                        " && java -jar ",global$PICARD," FilterSamReads I=tmp.bam O=intersection.bam READ_LIST_FILE=read_names.txt FILTER=includeReadList SORT_ORDER=coordinate CREATE_INDEX=true WRITE_READS_FILES=false",
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
      varpos <- as.numeric(gsub("([^;]*);([^;]*);([^;]*)","\\2",names(split.sams)[x]))
      tmp <- try(data.frame(do.call(rbind,strsplit(suppressWarnings(system(paste("samtools depth .sc",tag,"/power_bams/",x,".sc",tag,".bam .bulk",tag,"/power_bams/",x,".bulk",tag,".bam",
                                                                                 " | awk '{if($3 < $4){print $1\"\\t\"$2-1\"\\t\"$2\"\\t\"$3} else {print $1\"\\t\"$2-1\"\\t\"$2\"\\t\"$4}}'",
                                                                                 " | grep -v '\t0$' | grep -v '\t1$'",sep=""),intern=T)),"\t"))),silent=T)
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
                          stat=numeric(0),
                          dist=numeric(0),
                          id=character(0),
                          context=character(0)))
      }
      tmp[,1] <- as.character(tmp[,1])
      tmp[,2] <- as.numeric(tmp[,2])
      tmp[,3] <- as.numeric(tmp[,3])
      tmp[,4] <- as.numeric(tmp[,4])
      names(tmp) <- c("chromosome","start","stop","stat")
      tmp$dist <- tmp$start - (varpos - 1)
      tmp$id <- paste(tmp$chromosome,tmp$start,sep=";")
      system(paste("rm .sc",tag,"/power_bams/",x,".sc",tag,".bam &&",
                   "rm .bulk",tag,"/power_bams/",x,".bulk",tag,".bam",sep=""))
      return(tmp)
    })
    ind <- names(split.sams) %in% allele.key$id[allele.key$allele == 1]
    one <- do.call(rbind,dumb[ind])
    if(is.null(one)) {
      one <- data.frame(chromosome=character(0),
                        start=numeric(0),
                        stop=numeric(0),
                        stat=numeric(0),
                        dist=numeric(0),
                        id=character(0),
                        context=character(0))
    }
    one <- one[order(one$chromosome,one$start,one$stop),]
    two <- do.call(rbind,dumb[!ind])
    if(is.null(two)) {
      two <- data.frame(chromosome=character(0),
                        start=numeric(0),
                        stop=numeric(0),
                        stat=numeric(0),
                        dist=numeric(0),
                        id=character(0),
                        context=character(0))
    }
    two <- two[order(two$chromosome,two$start,two$stop),]
    
    take.max <- function(allele) {
      tst <- table(allele$id)
      unchange <- allele[allele$id %in% names(which(tst == 1)),]
      fix <- allele[allele$id %in% names(which(tst != 1)),]
      fix <- do.call(rbind,lapply(base::split(fix,fix$id),function(x){x[which.max(x$stat),]}))
      fix <- rbind(fix,unchange)
      fix <- fix[order(fix$chromosome,fix$start,fix$stop),]
      return(fix)
    }
    one <- take.max(one)
    two <- take.max(two)
    
    mergefeature <- function(allele,name) {
      local <- function(i) {
        tmp <- allele[allele$stat == i,c(1,2,3)]
        write.table(tmp,file="tmp.bed",sep="\t",quote=F,row.names = F,col.names = F)
        system("cat tmp.bed | bedtools merge > tmp.res")
        res <- read.table("tmp.res",sep="\t",header=F,colClasses = c("character","numeric","numeric"))
        names(res) <- c("chromosome","start","stop")
        res$stat <- i
        return(res)
      }
      res <- do.call(rbind,lapply(unique(allele$stat),local))
      if(is.null(res)) {
        system(paste("touch ",name,sep=""))
        return(0)
      }
      res <- res[order(res$chromosome,res$start,res$stop),]
      if(nrow(res) != 0) {
        res <- aggregate(x=res$stat,by=list(chromosome=res$chromosome,start=res$start,stop=res$stop),FUN=max)
        names(res)[4] <- "stat"
        write.table(res,name,sep="\t",col.names = F,row.names = F,quote = F)
      } else {
        system(paste("touch ",name,sep=""))
      }
    }
    dumb <- mergefeature(one,paste("powers.one",tag,".bed",sep=""))
    dumb <- mergefeature(two,paste("powers.two",tag,".bed",sep=""))
    
    powers <- table(c(one$stat,two$stat))
    save(powers,file=paste("powers",tag,".rda",sep=""))
    
    context.bed <- unique(rbind(one[,1:3],two[,1:3]))
    rownames(context.bed) <- paste(context.bed$chromosome,context.bed$start,sep=";")
    context.bed <- context.bed[order(context.bed[,1],context.bed[,2],context.bed[,3]),]
    context.bed[,2] <- context.bed[,2] - 1
    context.bed[,3] <- context.bed[,3] + 1
    
    write.table(context.bed,file="context.bed",quote=F,row.names = F,col.names = F,sep="\t")
    context <- system(paste("bedtools getfasta -fi ",config$reference_file," -bed context.bed -fo - | grep -v '>'",sep=""),intern=T)
    context[substr(context,2,2) %in% c("A","G")] <- reverse.complement(context[substr(context,2,2) %in% c("A","G")])
    context.bed[,"context"] <- context
    
    
    one$context <- context.bed[one$id,"context"]
    two$context <- context.bed[two$id,"context"]
    combined <- rbind(one[,c("stat","dist","context")],two[,c("stat","dist","context")])
    power.array <- table(combined)
    save(power.array,file="power.array.rda")
    
    system(paste("rm -r .sc",tag," 2> /dev/null",sep=""))
    system(paste("rm -r .bulk",tag," 2> /dev/null",sep=""))
    system("rm tmp.bed 2> /dev/null && rm tmp.res 2> /dev/null")
    setwd(config$analysis_path)
  }
  dumb <- local()
  tmp <- strsplit(work.dir,"/")[[1]]
  n <- length(tmp) 
  progress.file <- paste(config$analysis_path,"/progress/.power_",bulk.config$name,"_",paste(tmp[n - 2],"_",tmp[n],sep=""),sep="")
  out.log.cmd(paste("touch ",progress.file,sep=""))
  setwd(owd)
}

varcall <- function(config,bulk.config,overwrite) {
  dir <- paste("varcall_",bulk.config$name,sep="")
  if(overwrite) {
    out.log("Overwriting...")
    out.log.cmd(paste("rm -r ",dir," 2> /dev/null",sep=""))
  }
  chromosomes <- get.chromosomes(config)
  
  done <- unlist(lapply(chromosomes,function(chromosome){
    tmp <- paste(chromosome,"_",list.files(paste(config$analysis_path,"/",chromosome,"/jobs",sep="")),sep="")
    done <- file.exists(paste(config$analysis_path,"/progress/.power_",bulk.config$name,"_",tmp,sep=""))
    names(done) <- paste(config$analysis_path,"/power.",bulk.config$name,"_job_scripts/",tmp,".sh",sep="")
    return(done)
  }))

  if(any(!done)) {
    out.log("Jobs not done from ppower:")
    sapply(names(done)[!done],function(x){out.log(x)})
    stop("Jobs not done from ppower")
  }
  suppressWarnings(dir.create(dir))
  setwd(dir)
  
  #Load power files

  
  job.dirs <- unlist(lapply(chromosomes,function(x){list.files(paste("../",x,"/jobs",sep=""),full.names = T)}))
  
  powers <- combine.objects(lapply(job.dirs,function(x){
    load(paste(x,"/powers.",bulk.config$name,".rda",sep=""))
    return(powers)
  }))
  save(powers,file=paste("powers.",bulk.config$name,".rda",sep=""))
  
  
  powers.one <- paste(job.dirs,"/powers.one.",bulk.config$name,".bed",sep="")
  powers.two <- paste(job.dirs,"/powers.two.",bulk.config$name,".bed",sep="")
  
  #Create allele specific power bed files
  bed.combine <- function(tag,bed.vec) {
    tmp <- paste("powers.",tag,".bed.tmp",sep="")
    final <- paste("powers.",tag,".bed",sep="")
    out.log.cmd(paste("rm ",tmp," 2> /dev/null",sep=""))
    for(i in seq_along(bed.vec)) {
      out.log.cmd(paste("cat ",bed.vec[i]," >> ",tmp,sep=""))
    }
    out.log.cmd(paste("cat ",tmp," | sort -V > ",final," && rm ",tmp,sep=""))
  }
  dumb <- bed.combine("one",powers.one)
  dumb <- bed.combine("two",powers.two)
  
  #Combine mutation and data files
  mutations <- do.call(rbind,lapply(chromosomes,function(chr){load(paste("../",chr,"/compare/mutations.",bulk.config$name,".rda",sep="")); return(mutations)}))
  save(mutations,file="mutations.rda")
  
  candidate.somatic <- do.call(rbind,lapply(chromosomes,function(chr){load(paste("../",chr,"/compare/candidate.somatic.",bulk.config$name,".rda",sep="")); names(candidate.somatic)[names(candidate.somatic) == "composite_coverage"] <- "stat"; return(candidate.somatic)}))
  save(candidate.somatic,file="candidate.somatic.rda")
  
  all.somatic <- do.call(rbind,lapply(chromosomes,function(chr){load(paste("../",chr,"/compare/all.somatic.",bulk.config$name,".rda",sep="")); names(all.somatic)[names(all.somatic) == "composite_coverage"] <- "stat"; return(all.somatic)}))
  save(all.somatic,file="all.somatic.rda")
  
  data <- do.call(rbind,lapply(chromosomes,function(chr){load(paste("../",chr,"/compare/data.",bulk.config$name,".rda",sep="")); return(data)}))
  save(data,file="data.rda")
  
  all.mosaic <- do.call(rbind,lapply(chromosomes,function(chr){load(paste("../",chr,"/compare/all.mosaic.",bulk.config$name,".rda",sep="")); return(all.mosaic)}))
  save(all.mosaic,file="all.mosaic.rda")
  
  #Make modified powers
  #Rate of spurious observation of somatic alternate allele in bulk
  bulk.alt.rate <- (sum(mutations$bulk.alt.count[mutations$onek.bulk.het] > 0)/sum(mutations$onek.bulk.het))/2
  
  #Rate of spurious disqualification of somatic sites by read discordance based on observed linkage between known polymorphisms
  onekg.linkage <- data[(data$hc.bulk >= 2) & (data$type == "1KG"),]
  onekg.linkage$pseudostat <- onekg.linkage$hc.single.cell
  onekg.linkage$pseudostat[onekg.linkage$hc.single.cell > onekg.linkage$hc.bulk] <- onekg.linkage$hc.bulk[onekg.linkage$hc.single.cell > onekg.linkage$hc.bulk]
  onekg.linkage <- onekg.linkage[onekg.linkage$pseudostat >= 2,]
  onekg.linkage$dc <- onekg.linkage$dc.single.cell > 0
  
  key <- table(c(onekg.linkage$pseudostat,0),c(onekg.linkage$dc,T))
  key <- key[-1,]
  key.rate <- key[,2]/(key[,1] + key[,2])
  key.rate <- key.rate[names(key.rate) %in% names(powers)]
  
  #Adjust powers
  powers.mod <- powers[names(key.rate)]
  overall.rate.observed <- (1 - bulk.alt.rate) * (1 - key.rate)
  powers.mod <- powers.mod * overall.rate.observed
  powers.mod[powers.mod == 0] <- 1
  save(powers.mod,file="powers.mod.rda")
  
  tmp <- seq(from=2,by=1,to=max(candidate.somatic$stat))
  vec <- rep(1,length(tmp))
  names(vec) <- tmp
  powers.mod <- powers.mod[names(powers.mod) %in% names(vec)]
  vec[names(powers.mod)] <- powers.mod
  powers.mod <- vec
  
  #Bootstrap germline
  suppressWarnings(dir.create("bootstrap"))
  suppressWarnings(setwd("bootstrap"))
  onek <- mutations[mutations$onek.bulk.het,]
  onek <- onek[!is.na(onek$stat),]
  onek$site <- rownames(onek)
  thresh <- 2
  candidate.somatic <- candidate.somatic[candidate.somatic$stat >= thresh,]
  onek <- onek[onek$stat >= thresh,]
  local <- function(mut,som) {
    mut <- base::split(mut,mut$distance)
    distances <- table(som$distance)
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
    for(i in 1:global$BOOTSTRAP_REPLICATES) {
      tmp <- lapply(seq_along(distances),function(x){
        n <- names(distances)[x]
        item <- new[[n]][sample(1:nrow(new[[n]]),size=distances[x],replace=F),]
      })
      booty <- do.call(rbind,tmp)
      save(booty,file=paste(i,".",thresh,".rda",sep=""))
    }
  }
  dumb <- local(onek,candidate.somatic)
  m <- max(candidate.somatic$stat)
  bootstrap <- colMeans(do.call(rbind,lapply(1:global$BOOTSTRAP_REPLICATES,function(x){load(paste(x,".",thresh,".rda",sep="")); return((table(c(booty$stat,2:m)) - 1)[as.character(2:m)])})))
  save(bootstrap,file="bootstrap.rda")
  setwd("..")
  
  #Fit two-component model and assign sSNV quality scores
  result <- data.frame(stat.val=as.numeric(names(bootstrap)),
                       observed.ssnv.count=as.numeric(table(c(candidate.somatic$stat,seq(from=2,by=1,to=max(as.numeric(names(bootstrap))))))) -1,
                       bootstrap=bootstrap,
                       power=powers.mod[as.numeric(names(powers.mod)) %in% names(bootstrap)])
  result$control.rate <- 1e9*((result$bootstrap + 0.5)/(result$power + 0.5))
  result$observed.ssnv.rate <- 1e9 * ((result$observed.ssnv.count+ 0.5)/(result$power + 0.5))
  tmp <- sapply(seq_along(result$observed.ssnv.rate),function(x) {
    qbeta(p = c(0.005,0.995),shape1 = 0.5 + result$observed.ssnv.count[x],shape2 = 0.5 + result$power[x] - result$observed.ssnv.count[x]) * 1e9
  })
  result$observed.ssnv.rate.lb <- tmp[1,]
  result$observed.ssnv.rate.ub <- tmp[2,]
  
  sample.observed.ssnv.rate <- function() {
    return(rbeta(n = nrow(result),shape1 = 0.5 + result$observed.ssnv.count,shape2 = 0.5 + result$power - result$observed.ssnv.count)*1e9)
  }
  
  alpha <- result$observed.ssnv.count + 0.5
  beta <- result$power + 0.5
  v <- (alpha * beta)/((alpha + beta)^2 * (alpha + beta + 1))
  v <- v/min(v)
  v <- 1/v
  result$error.rate <- (1/2)^(result$stat.val - 2)
  
  f <- function(params,rate=result$observed.ssnv.rate){
    a <- params[1]; b <- params[2]
    objective <- sum(v * ((a * result$error.rate + b * result$control.rate) - rate)^2)
    if((a < 0)|(b < 0)) {
      objective <- objective^2 
    }
    return(objective)
  }
  
  p <- c(result$observed.ssnv.rate[1],1)
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
  obj$data <- result
  obj$nlm.fit <- nlm.fit
  obj$nlm.fit.bg <- nlm.fit.bg
  obj$error.rate.bg <- t(apply(nlm.fit.bg,1,function(x){x[1] * result$error.rate}))
  colnames(obj$error.rate.bg) <- obj$data$stat.val
  
  obj$control.rate.bg <- t(apply(nlm.fit.bg,1,function(x){x[2] * result$control.rate}))
  colnames(obj$control.rate.bg) <- obj$data$stat.val
  
  obj$fit.bg <- obj$error.rate.bg + obj$control.rate.bg
  obj$data$error.rate <- obj$data$error.rate * nlm.fit$estimate[1]
  obj$data$control.rate <- obj$data$control.rate * nlm.fit$estimate[2]
  obj$data$fit.rate <- obj$data$control.rate + obj$data$error.rate
  obj$data$fpr.estimate <- obj$data$error.rate/(obj$data$control.rate + obj$data$error.rate)
  
  cr.max <- max(obj$control.rate.bg[,1])
  obj$control.rate.bg[,1][obj$control.rate.bg[,1] < 0] <- 0
  cr.min <- min(obj$control.rate.bg[,1])
  
  signal <- round(obj$data$observed.ssnv.count * (1 - obj$data$fpr.estimate))
  signal <- sum(signal) - cumsum(signal) + signal
  noise <- round(obj$data$observed.ssnv.count * obj$data$fpr.estimate)
  noise <- sum(noise) - cumsum(noise) + noise
  names(signal) <- names(noise) <- obj$data$stat.val
  obj$stat_threshold <- as.numeric(names(which(noise/(signal + noise) <= 0.1)))[1]
  if(is.na(obj$stat_threshold)) {
    obj$stat_threshold <- Inf
  }
  obj$somatic.rate <- obj$data$control.rate[1]
  obj$somatic.rate.lb <- cr.min
  obj$somatic.rate.ub <- cr.max
  obj$fpr <- (noise/(signal + noise))[which(data$stat.val == obj$stat_threshold)]
  eval(parse(text=paste("obj$rsomatic <- function(n){sample(c(",paste(obj$control.rate.bg[,1],collapse=","),"),replace=F,size=n)}",sep="")))
  
  #Make summary file
  metrics <- do.call(rbind,lapply(chromosomes,function(x){
    load(paste("../",x,"/compare/summary.",bulk.config$name,".rda",sep="")); return(summary)
  }))[,c("onek","germline.phased","germline.misphased","somatic.phased","somatic.misphased")]
  metrics <- colSums(metrics)
  
  summary <- list()
  summary[["Polymorphic germline SNVs: "]] <- metrics["onek"]
  summary[["Properly phased germline SNVs: "]] <- metrics["germline.phased"]
  summary[["Improperly phased germline SNVs: "]] <- metrics["germline.misphased"]
  summary[["Properly phased sSNV candidates: "]] <- metrics["somatic.phased"]
  summary[["Improperly phased sSNV candidates: "]] <- metrics["somatic.misphased"]
  summary[["sSNV rate (Nuc/Gbp): "]] <- round(obj$somatic.rate)
  summary[["\tLower bound (98%): "]] <- round(obj$somatic.rate.lb)
  summary[["\tUpper bound (98%): "]] <- round(obj$somatic.rate.ub)
  summary[["Composite coverage threshold: "]] <- obj$stat_threshold
  
  
  #Make rate plot
  pdf("rate-plot.pdf")
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
  lines(x = c(obj$stat_threshold,obj$stat_threshold),y=c(0,ylim[2]),lty="dashed",lwd=1)
  par(xpd=NA)
  mtext(text = "sSNVs/GBp",side = 2,at = -ylim[2]/6,adj = c(0.5,0.5),line = 4)
  par(xpd=F)
  
  par(xpd=NA)
  mtext(text = "Composite coverage (CC)",side = 1, at = (xlim[1] + xlim[2])/2, adj=c(0.5,0.5),line=3)
  par(xpd=F)
  dev.off()
  
  save(obj,file="call.rate_fit.rda")
  
  #Annotate final variant list
  all.somatic$phred.quality <- NA
  all.somatic$fp.prob <- NA
  all.somatic$phred.quality[all.somatic$stat < 0] <- 0
  all.somatic$fp.prob[all.somatic$stat < 0] <- 1
  for(j in 1:nrow(obj$data)) {
    all.somatic$phred.quality[all.somatic$stat == obj$data$stat.val[j]] <- -10*log10(obj$data$fpr.estimate[j])
    all.somatic$fp.prob[all.somatic$stat == obj$data$stat.val[j]] <- obj$data$fpr.estimate[j]
  }
  all.somatic$status <- "low_power"
  all.somatic$status[all.somatic$stat < 0] <- "filtered_fp"
  all.somatic$status[all.somatic$stat >= 2] <- "uncertain_call"
  all.somatic$status[all.somatic$stat >= obj$stat_threshold] <- "pass"
  
  tmp <- mutations[all.somatic$to,"bulk.phasing"]
  ind <- all.somatic$orientation == "trans"
  ind2 <- tmp == 2
  tmp[ind2 & ind] <- 1
  tmp[!ind2 & ind] <- 2
  all.somatic$bulk.phasing <- tmp
  ssnvs <- all.somatic
  
  unlinked.somatic <- mutations[mutations$somatic,]
  unlinked.somatic <- unlinked.somatic[!(rownames(unlinked.somatic) %in% rownames(ssnvs)),]
  unlinked.somatic$status <- "unlinked"
  unlinked.somatic[!is.na(unlinked.somatic$stat),"status"] <- "bulk_support"
  unlinked.somatic[!is.na(unlinked.somatic$stat) & (unlinked.somatic$stat < 0),"status"] <- "bulk_support;filtered_fp"
  unlinked.somatic[!is.na(unlinked.somatic$stat) & (unlinked.somatic$stat > 0) & (unlinked.somatic$stat < 2),"status"] <- "bulk_support;low_power"
  fix <- colnames(ssnvs)[!(colnames(ssnvs) %in% colnames(unlinked.somatic))]
  for(f in fix) {
    unlinked.somatic[,f] <- NA
  }
  unlinked.somatic <- unlinked.somatic[,colnames(ssnvs)]
  ssnvs <- rbind(ssnvs,unlinked.somatic)
  
  ssnvs$status[ssnvs$status == "call"] <- "pass"
  mosaic <- ssnvs[!is.na(ssnvs$mosaic.stat),]
  mosaic$status[mosaic$mosaic.stat < 0] <- "mosaic;filtered_fp"
  mosaic$status[mosaic$mosaic.stat > 0] <- "mosaic;uncertain_call"
  mosaic$status[mosaic$mosaic.stat >= obj$stat_threshold] <- "mosaic;pass"
  
  tmp <- mutations[mosaic$mosaic.to,"bulk.phasing"]
  ind <- mosaic$mosaic.orientation == "trans"
  ind2 <- tmp == 2
  tmp[ind2 & ind] <- 1
  tmp[!ind2 & ind] <- 2
  mosaic$bulk.phasing <- tmp
  for(c in colnames(mosaic)[grepl("mosaic",colnames(mosaic))]) {
    if(c != "mosaic.proper") {
      mosaic[,gsub("mosaic.","",c)] <- mosaic[,c]
    }
  }
  onek <- mutations[mutations$onek.bulk.het,]
  onek$phred.quality <- NA
  onek$fp.prob <- NA
  onek$status <- "germline_polymorphic_het"
  
  ssnvs <- ssnvs[is.na(ssnvs$mosaic.stat),]
  ssnvs <- rbind(ssnvs,mosaic,onek)
  
  n <- nrow(ssnvs)
  fmt.gt <- rep("GT",n)
  smpl.gt <- rep("0|1",n)
  ind <- is.na(ssnvs$bulk.phasing)
  fmt.gt[ind] <- gsub("|","/",fmt.gt[ind],fixed = T)
  smpl.gt[ind] <- gsub("|","/",smpl.gt[ind],fixed=T)
  ind <- ssnvs$bulk.phasing == 2
  smpl.gt[ind] <- "1|0"
  ind <- !(ssnvs$status %in% c("mosaic;pass","pass","germline_polymorphic_het"))
  smpl.gt[ind] <- gsub("1",".",smpl.gt[ind],fixed=T)
  
  fmt.hcb <- rep("HCB",n)
  smpl.hcb <- ssnvs$hc.bulk
  ind <- is.na(ssnvs$hc.bulk)
  smpl.hcb[ind] <- "."
  
  fmt.dr <- rep("DR",n)
  smpl.dr <- ssnvs$dc.single.cell
  ind <- is.na(ssnvs$dc.single.cell)
  smpl.dr[ind] <- "."
  
  fmt.hc <- rep("HC",n)
  smpl.hc <- ssnvs$hc.single.cell
  ind <- is.na(ssnvs$hc.single.cell)
  smpl.hc[ind] <- "."
  
  fmt.cc <- rep("CC",n)
  smpl.cc <- ssnvs$stat
  ind <- is.na(ssnvs$stat) | (ssnvs$stat < 0)
  smpl.cc[ind] <- "."
  
  fmt.g <- rep("G",n)
  smpl.g <- ssnvs$to
  ind <- is.na(ssnvs$to)
  smpl.g[ind] <- "."
  
  fmt.o <- rep("O",n)
  smpl.o <- ssnvs$orientation
  ind <- is.na(ssnvs$orientation)
  smpl.o[ind] <- "."
  
  fmt.rb <- rep("RB",n)
  smpl.rb <- ssnvs$bulk.ref
  
  fmt.ab <- rep("AB",n)
  smpl.ab <- ssnvs$bulk.alt
  
  fmt.rs <- rep("RS",n)
  smpl.rs <- ssnvs$single.cell.ref
  
  fmt.as <- rep("AS",n)
  smpl.as <- ssnvs$single.cell.alt
  
  format.lines <- paste(fmt.gt,fmt.hcb,fmt.dr,fmt.hc,fmt.cc,fmt.g,fmt.o,fmt.rb,fmt.ab,fmt.rs,fmt.as,sep=":")
  sample.lines <- paste(smpl.gt,smpl.hcb,smpl.dr,smpl.hc,smpl.cc,smpl.g,smpl.o,smpl.rb,smpl.ab,smpl.rs,smpl.as,sep=":")
  
  ssnvs$format.lines <- format.lines
  ssnvs$sample.lines <- sample.lines
  quality.lines <- as.character(round(ssnvs$phred.quality))
  quality.lines[is.na(ssnvs$phred.quality)] <- "."
  
  vcf.out <- "out.vcf"
  writeLines(c("##fileformat=VCFv4.0",
               paste("##fileDate=",format(Sys.time(),"%Y%m%d"),sep=""),
               paste("##reference=",config$reference_file,sep=""),
               paste("##source=lira varcall -s ",args[1]," -b ",args[3],sep=""),
               sapply(chromosomes,function(x){paste("##contig=<ID=",x,">",sep="")}),
               "##FILTER=<ID=UNLINKED,Description=\"No spaning reads\">",
               "##FILTER=<ID=LOW_POWER,Description=\"Low number of spanning reads in bulk or single-cell (CC = 0 or 1)\">",
               "##FILTER=<ID=FILTERED_FP,Description=\"Artifact inferred from discordant reads\">",
               "##FILTER=<ID=UNCERTAIN_CALL,Description=\"No discordant reads and CC >= 2, but variant failed to pass composite coverage threshold.\">",
               "##FILTER=<ID=PASS,Description=\"All filters passed\">",
               "##FILTER=<ID=MOSAIC,Description=\"sSNV pattern is consistent with mosacism; i.e., sSNV supporting haplotype and unmutated haplotypes are observed in bulk.\">",
               "##FILTER=<ID=BULK_SUPPORT,Description=\"Reads support the sSNV in bulk, but no linked haplotype with a nearby gHet is observed.\">",
               "##FILTER=<ID=GERMLINE_POLYMORPHIC_HET,Description=\"Germline polymorphic heterozygous site.\">",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype based on reference panel based haplotype phasing of sSNV-linked germline SNP.\">",
               "##FORMAT=<ID=HCB,Number=1,Type=Integer,Description=\"Haplotype Coverage in Bulk: Spanning read coverage of sSNV position and highest CC (or for FILTERED_FP, highest spanning read coverage) linked germline variant in bulk.\">",
               "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Discordant Read count in single cell.\">",
               "##FORMAT=<ID=HC,Number=1,Type=Integer,Description=\"Haplotype Coverage: sSNV-supporting spanning read coverage of sSNV position and highest CC (or for FILTERED_FP, highest spanning read coverage) linked germline variant in single cell.\">",
               "##FORMAT=<ID=CC,Number=1,Type=Integer,Description=\"Composite Coverage: minimum of HC and HCB, if DR=0.\">",
               "##FORMAT=<ID=G,Number=1,Type=String,Description=\"Germline variant linked: identity of the germline variant used to generate HCB, DR, and HC\">",
               "##FORMAT=<ID=O,Number=1,Type=String,Description=\"Orientation of linkage: cis/trans\">",
               "##FORMAT=<ID=RB,Number=1,Type=Integer,Description=\"Reference-supporting reads in Bulk (from GATK)\">",
               "##FORMAT=<ID=AB,Number=1,Type=Integer,Description=\"Alt-supporting reads in Bulk (from GATK)\">",
               "##FORMAT=<ID=RS,Number=1,Type=Integer,Description=\"Reference-supporting reads in Single cell (From GATK)\">",
               "##FORMAT=<ID=AS,Number=1,Type=Integer,Description=\"Alt-supporting reads in Single cell (from GATK)\">",
               paste("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",config$sample),
               paste(ssnvs$chromosome,ssnvs$pos,ssnvs$id,ssnvs$ref,ssnvs$alt,quality.lines,toupper(ssnvs$status),".",format.lines,sample.lines,sep="\t")),con=vcf.out)
  out.log.cmd(paste("cat ",vcf.out," | grep '#' > tmp.vcf && cat ",vcf.out," | grep -v '#' | sort -V >> tmp.vcf && mv tmp.vcf ",vcf.out," && bgzip -f ",vcf.out," && tabix -f ",vcf.out,".gz",sep = ""))
  ssnvs <- ssnvs[order(ssnvs$chromosome,ssnvs$pos),]
  save(ssnvs,file=paste("out.rda",sep=""))
  
  summary[["Total passing sSNVs: "]] <- sum(ssnvs$status %in% c("pass","mosaic;pass"))
  if(config$gender == "male") {
    summary[["Estimated total # of sSNVs: "]] <- round((3.122 + 3.227) * obj$somatic.rate)
  } else if (cnofig$gender == "female") {
    summary[["Estimated total # of sSNVs: "]] <- round((3.227  * 2) * obj$somatic.rate)
  }
  summary[["Estimated sensitivity: "]] <- summary[["Total passing sSNVs: "]]/summary[["Estimated total # of sSNVs: "]]
  summary[["Estimated FPR: "]] <- sprintf("%0.2f",mean(ssnvs$fp.prob[ssnvs$status == "pass"]))
  
  summary.lines <- paste(names(summary),summary,sep="")
  writeLines(summary.lines,con="summary.txt")
  
  out.log.cmd(paste("cat powers.one.bed | awk '{if($4 >= ",obj$stat_threshold,") print}' | cut -f 1,2,3 | bedtools merge -d 0 -i - > powered.regions.one.bed",sep=""))
  out.log.cmd(paste("cat powers.two.bed | awk '{if($4 >= ",obj$stat_threshold,") print}' | cut -f 1,2,3 | bedtools merge -d 0 -i - > powered.regions.two.bed",sep=""))
  out.log.cmd("bedtools multiinter -i powered.regions.one.bed powered.regions.two.bed  | cut -f 1,2,3,4 > powered.regions.bed")
  
  
  progress.file <- paste(config$analysis_path,"/progress/.varcall_",bulk.config$name,sep="")
  out.log.cmd(paste("touch ",progress.file,sep=""))
}

joint <- function(single.cell.configs,bulk.config,out.directory,use.uncertain.calls,use.low.power,overwrite) {
  if(overwrite) {
    system(paste("rm -r ",out.directory," 2> /dev/null",sep=""))
  }
  suppressWarnings(dir.create(out.directory))
  setwd(out.directory)
  search <- "PASS"
  if(use.uncertain.calls) {
    search <- c(search,"UNCERTAIN_CALL")
  }
  if(use.low.power) {
    search <- c(search,"LOW_POWER")
  }
  search.string <- paste("grep -v '#' | grep ",paste(" -e ",search,sep="",collapse=" "),sep="")
  system("rm variants.txt 2> /dev/null")
  for(i in seq_along(single.cell.configs)) {
    vcf <- paste(single.cell.configs[[i]]$analysis_path,"/varcall_",bulk.config$name,"/out.vcf.gz",sep="")
    out.log.cmd(paste("bcftools query -i 'ID == \".\"' -f '%CHROM;%POS;%REF;%ALT\\t[%O]\\t[%G]\\t%FILTER\\n' ",vcf," | ",search.string," | cut -f 1,2,3 >> variants.txt",sep=""))
    out.log.cmd(paste("bcftools query -i 'FORMAT/AB == 0 && FORMAT/RB > 0 && FORMAT/AS > 0 && ID == \".\"' -f '%CHROM;%POS;%REF;%ALT\\t[%O]\\t[%G]\\t%FILTER\\n' ",vcf," | cut -f 1,2,3 >> variants.txt",sep=""))
  }
  out.log.cmd("cat variants.txt | sort | uniq > tmp.txt && mv tmp.txt variants.txt")
  variants <- read.table("variants.txt")
  names(variants) <- c("sSNV","orientation","linked_gHet")
  variants <- unique(variants)
  save(variants,file="variants.rda")
  
  bash.commands <- paste("Rscript --vanilla $LIRA_DIR/scripts/main.R ",sapply(single.cell.configs,function(x){x$analysis_path}),"/config.txt joint_subset ",bulk.config$analysis_path,"/config.txt ",getwd(),sep="")
  job.names <- paste(digest(list(bash.commands,date())),seq_along(bash.commands),sep="_")
  job.loop(bash.commands,job.names)
  
  all.ssnvs <- lapply(single.cell.configs,function(x){
    load(paste(x$name,".rda",sep=""))
    return(ssnvs)
  })
  
  data <- lapply(all.ssnvs,function(x){
    obj <- list()
    call <- rep(NA,nrow(x))
    call[x$single.cell.alt > 0] <- 1
    call[x$single.cell.ref > 0 & x$single.cell.alt == 0] <- 0
    obj$call <- call
    
    tmp <- x[x$status %in% tolower(search),]
    obj$varlist <- rownames(tmp)
    return(obj)
  })
  
  rescue.mat <- sapply(data,function(x){
    return(x$call)
  })
  rownames(rescue.mat) <- all.ssnvs[[1]]$var_id
  colnames(rescue.mat) <- sapply(single.cell.configs,function(x){x$name})
  rescue.mat <- t(rescue.mat)
  
  varlist <- unique(unlist(lapply(data,function(x){
    return(x$varlist)
  })))
  combined <- do.call(rbind,lapply(single.cell.configs,function(x) {
    load(paste(x$name,".rda",sep=""))
    cell <- x$name
    
    x <- ssnvs
    x <- x[varlist,]
    return(x)
  }))
  combined <- base::split(combined,combined$var_id)
  save(combined,file="combined.rda")
  save(rescue.mat,file="rescue.mat.rda")
  
  # load("combined.rda")
  combined <- combined[sapply(combined,function(x){any(x$hc.bulk > 0)})]
  combined <- combined[sapply(combined,function(x){sum(x$null_haplotype > 0) > 1})]
  combined <- combined[sapply(combined,function(x){sum(x$var_haplotype > 0 & x$null_haplotype == 0) > 1})]
  combined <- combined[sapply(combined,function(x){sum(x$single.cell.alt >= 1) > 1})]
  combined <- combined[sapply(combined,function(x){sum(grepl(paste(tolower(search),collapse="|"),x$status)) > sum(grepl("filtered_fp",x$status))})]

  combined <- lapply(combined,function(x){x$joint_call <- NA; x$joint_call[x$single.cell.alt > 0] <- 1; x$joint_call[x$null_haplotype > 0] <- 0; x$joint_call[x$null_haplotype > 0 & x$var_hapotype > 0] <- NA; return(x)})

  mat <- sapply(combined,function(x){x$joint_call})
  ind <- apply(mat,2,function(x){var(x,na.rm=T) != 0})
  mat <- mat[,ind]
  rownames(mat) <- combined[[1]]$cell
  
  rescue.mat <- rescue.mat[,(!colnames(rescue.mat) %in% colnames(mat))]
  ind <- apply(rescue.mat,2,function(x){sum(x == 1,na.rm=T) > 1 & sum(x == 0,na.rm = T) > 1})
  rescue.mat <- rescue.mat[,ind]
  
  mat.to.plot <- rbind(mat,rep(0,ncol(mat)))
  mat.to.plot <- cbind(mat.to.plot,rep(0,nrow(mat.to.plot)))
  heatmap(-mat.to.plot,scale="none")
  
  n <- ncol(rescue.mat)
  getbg <- function(n) {
    tmp.mat <- mat
    tmp.rescue.mat <- rescue.mat
    
    tmp.mat[is.na(tmp.mat)] <- -2
    tmp.mat[tmp.mat == 0] <- -1
    tmp.mat[tmp.mat == -2] <- 0
    
    tmp.rescue.mat[is.na(tmp.rescue.mat)] <- -2
    tmp.rescue.mat[tmp.rescue.mat == 0] <- -1
    tmp.rescue.mat[tmp.rescue.mat == -2] <- 0
    
    true.statistic <- t(tmp.mat) %*% tmp.rescue.mat
    for(i in 2:n) {
      tst <- t(tmp.mat) %*%  tmp.rescue.mat[sample(1:nrow(tmp.rescue.mat),size=nrow(tmp.rescue.mat),replace = F),]
      if(i == 1) {
        lim <- tst
      } else {
        lim[tst > lim] <- tst[tst > lim]
      }
      keep.tst <- !apply(true.statistic <= lim,2,all)
      tmp.rescue.mat <- tmp.rescue.mat[,keep.tst]
      lim <- lim[,keep.tst]
      true.statistic <- true.statistic[,keep.tst]
      print(ncol(lim))
    }
    

    n.keep <- unlist(lapply(1:nrow(true.statistic),function(x){names(which(true.statistic[x,] > lim[x]))}))
    return(n.keep)
  }
  
  mat <- cbind(mat,rescue.mat[,getbg(n)])
  mat.to.plot <- rbind(mat,rep(0,ncol(mat)))
  mat.to.plot <- cbind(mat.to.plot,rep(0,nrow(mat.to.plot)))
  heatmap(-mat.to.plot,scale="none")
  
  dev.off()
}

joint.subset <- function(config,bulk.config,work.dir) {
  setwd(work.dir)
  load("variants.rda")
  chromosomes <- get.chromosomes(config)
  load(paste(config$analysis_path,"/varcall_",bulk.config$name,"/out.rda",sep=""))
  
  ghet <- ssnvs[ssnvs$onek.bulk.het,]
  ado <- sum(ghet$single.cell.alt == 0 | ghet$single.cell.ref == 0)/nrow(ghet)
  save(ado,file=paste(config$name,".ado.rda",sep=""))
  
  sites <- variants$sSNV

  linkage <- do.call(rbind,lapply(chromosomes,function(x){load(paste(config$analysis_path,"/",x,"/compare/single.cell.linkage.rda",sep="")); return(single.cell.linkage)}))
  tmp <- rbind(data.frame(site1=variants$sSNV,site2=variants$linked_gHet,somatic.first=T,cis=variants$orientation == "cis",pair_id=paste(variants$sSNV,"~",variants$linked_gHet,sep=""),somatic.variant=variants$sSNV,row=1:nrow(variants)),
               data.frame(site1=variants$linked_gHet,site2=variants$sSNV,somatic.first=F,cis=variants$orientation == "cis",pair_id=paste(variants$linked_gHet,"~",variants$sSNV,sep=""),somatic.variant=variants$sSNV,row=1:nrow(variants)))
  tmp <- tmp[tmp$pair_id %in% linkage$pair_id,]
  rownames(linkage) <- linkage$pair_id
  for(hap in c("RR","RV","VR","VV")) {
    tmp[,hap] <- linkage[tmp$pair_id,hap]
  }
  
  tmp$null_haplotype <- ""
  tmp$var_haplotype <- ""
  
  ind <- tmp$somatic.first & tmp$cis
  tmp$null_haplotype[ind] <- tmp$RV[ind]
  tmp$var_haplotype[ind] <- tmp$VV[ind]
  
  ind <- tmp$somatic.first & !tmp$cis
  tmp$null_haplotype[ind] <- tmp$RR[ind]
  tmp$var_haplotype[ind] <- tmp$VR[ind]
  
  ind <- !tmp$somatic.first & tmp$cis
  tmp$null_haplotype[ind] <- tmp$VR[ind]
  tmp$var_haplotype[ind] <- tmp$VV[ind]
  
  ind <- !tmp$somatic.first & !tmp$cis
  tmp$null_haplotype[ind] <- tmp$RR[ind]
  tmp$var_haplotype[ind] <- tmp$RV[ind]
  
  tmp$null_haplotype <- as.numeric(tmp$null_haplotype)
  tmp$var_haplotype <- as.numeric(tmp$var_haplotype)
  
  tmp <- tmp[,c("somatic.variant","null_haplotype","var_haplotype")]
  tmp2 <- data.frame(do.call(rbind,lapply(base::split(tmp[,c("null_haplotype","var_haplotype")],tmp$somatic.variant),colSums)))
  
  ssnvs <- ssnvs[rownames(ssnvs) %in% sites,]
  ssnvs$null_haplotype <- 0
  ssnvs$var_haplotype <- 0
  ssnvs[rownames(tmp2),"null_haplotype"] <- tmp2$null_haplotype
  ssnvs[rownames(tmp2),"var_haplotype"] <- tmp2$var_haplotype
  ssnvs$var_id <- rownames(ssnvs)
  ssnvs$cell <- config$name
  rownames(ssnvs) <- NULL
  save(ssnvs,file=paste(config$name,".rda",sep=""))
}