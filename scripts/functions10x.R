setup <- function(config) {
  owd <- getwd()
  setwd(config$analysis_path)
  
  #soft link bams
  res <- out.log.cmd(paste("ln -s ",config$bam," reads.bam",sep=""))
  res <- out.log.cmd(paste("ln -s ",config$bam,".bai reads.bam.bai",sep=""))
  res <- out.log.cmd(paste("ln -s ",config$bam_bulk," reads_bulk.bam",sep=""))
  res <- out.log.cmd(paste("ln -s ",config$bam_bulk,".bai reads_bulk.bam.bai",sep=""))
  res <- out.log.cmd(paste("ln -s ",config$vcf," input_calls.vcf.gz",sep=""))
  res <- out.log.cmd(paste("ln -s ",config$vcf,".tbi input_calls.vcf.gz.tbi",sep=""))
  res <- out.log.cmd(paste("ln -s ",config$vcf_bulk," input_calls_bulk.vcf.gz",sep=""))
  res <- out.log.cmd(paste("ln -s ",config$vcf_bulk,".tbi input_calls_bulk.vcf.gz.tbi",sep=""))
  
  #annotate with snpsift
  write.table(x=config$config,file = "config.txt",sep="\t",col.names = F,row.names = F,quote=F)
  suppressWarnings(dir.create("job_scripts"))
  suppressWarnings(dir.create("progress"))
  
  out.log.cmd("touch progress/.setup")
  setwd(owd)
}

split <- function(config,chromosome) {
  suppressWarnings(dir.create(chromosome))
  setwd(chromosome)
  
  suppressWarnings(dir.create("split"))
  setwd("split")
  
  #getting sites
  
  process.local <- function(bulk=F) {
    if(bulk) {
      newdir <- "bulk"
      sample <- config$sample_bulk
      input.vcf <- "../../../input_calls_bulk.vcf.gz"
      input.bam <- "../../../reads_bulk.bam"
    } else {
      newdir <- "single_cell"
      sample <- config$sample
      input.vcf <- "../../../input_calls.vcf.gz"
      input.bam <- "../../../reads.bam"
    }
    suppressWarnings(dir.create(newdir))
    setwd(newdir)
    
    #Remove -I later
    out.log("Getting variants from vcf (output: calls.vcf.gz)")
    out.log.cmd(paste("bcftools view -I -s ",sample," -i 'TYPE=\"snp\" && N_ALT=1' ",input.vcf," ",chromosome," > calls.vcf && bgzip -f calls.vcf && tabix -f calls.vcf.gz",sep=""))
    
    #Create SnpSift annotated vcf
    out.log("Annotating variants with DBSNP (output: calls.id.vcf.gz)")
    out.log.cmd(paste("java -jar ",global$SNPEFF,"/SnpSift.jar annotate -tabix -name \"DBSNP_\" ",global$DBSNP," calls.vcf.gz > calls.id.vcf && bgzip -f calls.id.vcf && tabix -f calls.id.vcf.gz",sep=""))
    
    #Create sites.bed
    out.log("Getting all SNVs as a bed file (output: sites.bed)")
    out.log.cmd(paste("bcftools query -r ",chromosome," -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' calls.id.vcf.gz | awk '{print $1\"\\t\"$2-1\"\\t\"$2\"\\t\"$1\";\"$2\";\"$3\";\"$4}' > sites.bed",sep=""))
    
    #Get chromosomal bam file with filters
    out.log("Getting filtered read set (output: reads.bam)")
    out.log.cmd(paste("samtools view -b -q 60 -f 2 -F 3852 ",input.bam," ",chromosome," > reads.bam && samtools index reads.bam",sep=""))
    
    #Create barcode indexed bam
    suppressWarnings(dir.create("tmp"))
    out.log("Reformating bam file to allow indexing by barcode3 (output: reads.barcode_indexed.bam)")
    out.log.cmd("python ../../../../TXtools/x_toolbox.py -C -b reads.bam -n 1 -p ./tmp/ -o reads.barcode_indexed.bam")
    
    #out.log("Counting reads by barcode (output: block_info.txt)")
    #out.log.cmd("python ../../../../TXtools/x_toolbox.py -B -b reads.bam -o block_info.txt -n 1 -p ./tmp/ 2>&1 > /dev/null")
    
    out.log("Sorting by barcode (output: reads.sorted.barcode_indexed.bam)")
    out.log.cmd("samtools sort -o reads.sorted.barcode_indexed.bam -T ./tmp/temp -O bam reads.barcode_indexed.bam && samtools index reads.sorted.barcode_indexed.bam")
    out.log.cmd("rm -r tmp")
    setwd("..")
  }
  
  #single cell
  dumb <- process.local()
  
  #bulk
  dumb <- process.local(bulk=T)
  
  #aggregate sites across bulk and single cell
  out.log("Aggregating sites across bulk and single cell (output: sites.bed)")
  out.log.cmd("cat single_cell/sites.bed bulk/sites.bed | sort -V | uniq -u > sites.bed")
  
  process.local2 <- function(bulk=F) {
    if(bulk) {
      newdir <- "bulk"
      jobdir <- "jobs_bulk"
      jobtag <- "bulk_"
    } else {
      newdir <- "single_cell"
      jobdir <- "jobs"
      jobtag <- ""
    }
    setwd(newdir)
    
    out.log("Getting barcodes with variation (output: barcodes.txt)")
    out.log.cmd("bedtools tag -i reads.bam -files ../sites.bed -names | samtools view - | grep 'YB:Z' | grep -o 'BX:Z:[^\t]*' | sed 's#BX:Z:##g' | sort -V | tr '-' '_' | uniq > barcodes.txt")
    
    out.log("Getting reduced bam file (output: reads.reduced.sorted.barcode_indexed.bam)")
    barcodes <- readLines("barcodes.txt")
    bam.region.selector("reads.sorted.barcode_indexed.bam",barcodes,"reads.reduced.sorted.barcode_indexed.bam")
    
    out.log("Counting barcodes (output: read_barcodes.txt)")
    out.log.cmd("samtools view reads.reduced.sorted.barcode_indexed.bam | cut -f 3 > read_barcodes.txt")
    
    out.log("Creating barcode sets...")
    read.barcodes <- readLines("read_barcodes.txt")
    read.barcodes <- table(read.barcodes)
    read.barcodes <- read.barcodes[names(read.barcodes) %in% barcodes]
    assignment <- ceiling(cumsum(read.barcodes)/global$READS_TARGET)
    barcode.sets <- base::split(names(assignment),assignment)
    
    out.log("Writing job scripts...")
    suppressWarnings(dir.create(paste("../../",jobdir,sep="")))
    for(i in seq_along(barcode.sets)) {
      suppressWarnings(dir.create(paste("../../",jobdir,"/",i,sep="")))
      write.table(x = barcode.sets[[i]],file = paste("../../",jobdir,"/",i,"/barcodes.txt",sep=""),quote=F,col.names = F,row.names = F)
      js <- paste("../../../job_scripts/",jobtag,chromosome,"_",i,".sh",sep="")
      lines <- c("#!/bin/bash",paste("Rscript --vanilla $LIRA_DIR/scripts/main.R ",config$analysis_path,"/config.txt local_barcode_function ",config$analysis_path,"/",chromosome,"/",jobdir,"/",i,ifelse(bulk," bulk",""),sep=""))
      writeLines(text=lines,con=js)
      out.log.cmd(paste("chmod +x ",js,sep=""))
    }
    setwd("..")
  }
  
  #single cell
  dumb <- process.local2()
  
  #bulk
  dumb <- process.local2(bulk=T)
  
  
  dir.create("shapeit")
  setwd("shapeit")
  out.log("Make shapeit input from population-polymorphic SNPs")
  out.log.cmd(paste("bcftools view -h ../calls_bulk.id.vcf.gz",
                    " | grep -e '##contig' -e '#CHROM' -e '##FORMAT=<ID=GT' -e '##FILTER' -e '##ALT' -e '##fileformat'",
                    " > shapeit.vcf",sep=""))
  out.log.cmd(paste("bcftools view ../calls_bulk.id.vcf.gz",
                    " | awk 'BEGIN{OFS=\"\t\"}{if($8 ~ \"DBSNP_COMMON\"){$8 = \".\"; $9 = \"GT\"; print $0}}'",
                    " | tr ':' '\t'",
                    " | cut -f 1-10",
                    " | grep -e '#' -e '0/0' -e '0/1' -e '1/0' -e '1/1'",
                    ifelse(grepl("chr",chromosome),""," | sed 's#^#chr#g'"),
                    " >> shapeit.vcf",sep=""))
  pre <- paste(global$KGEN,"/1000GP_Phase3",sep="")
  if(chromosome != "X") {
    leg <- ifelse(is.na(file.info(paste(pre,"/1000GP_Phase3_chr",chromosome,".legend.gz",sep=""))$size),paste(pre,"/1000GP_Phase3_chr",chromosome,".legend",sep=""),paste(pre,"/1000GP_Phase3_chr",chromosome,".legend.gz",sep=""))
    hap <- ifelse(is.na(file.info(paste(pre,"/1000GP_Phase3_chr",chromosome,".hap.gz",sep=""))$size),paste(pre,"/1000GP_Phase3_chr",chromosome,".hap",sep=""),paste(pre,"/1000GP_Phase3_chr",chromosome,".hap.gz",sep=""))
    
    out.log("Run shapeit check")
    out.log.cmd(paste("shapeit -check -V shapeit.vcf -M ",pre,"/genetic_map_chr",chromosome,"_combined_b37.txt",
                      " --input-ref ",hap," ",leg," ",pre,"/1000GP_Phase3.sample ",
                      " --output-log shapeit-check.log",sep=""))
    out.log("Run shapeit")
    out.log.cmd(paste("shapeit -V shapeit.vcf -M ",
                      pre,"/genetic_map_chr",chromosome,"_combined_b37.txt --input-ref ",hap," ",leg," ",pre,"/1000GP_Phase3.sample ",
                      " --exclude-snp shapeit-check.snp.strand.exclude -O shapeit-phased",
                      "&& shapeit -convert --input-haps shapeit-phased --output-vcf shapeit-phased.vcf",sep=""))
  } else {
    leg <- ifelse(is.na(file.info(paste(pre,"/1000GP_Phase3_chrX_NONPAR.legend.gz",sep=""))$size),paste(pre,"/1000GP_Phase3_chrX_NONPAR.legend",sep=""),paste(pre,"/1000GP_Phase3_chrX_NONPAR.legend.gz",sep=""))
    hap <- ifelse(is.na(file.info(paste(pre,"/1000GP_Phase3_chrX_NONPAR.legend.hap.gz",sep=""))$size),paste(pre,"/1000GP_Phase3_chrX_NONPAR.hap",sep=""),paste(pre,"/1000GP_Phase3_chrX_NONPAR.hap.gz",sep=""))
    
    out.log("Run shapeit check")
    out.log.cmd(paste("shapeit -check -V shapeit.vcf -M ",pre,"/genetic_map_chrX_nonPAR_combined_b37.txt",
                      " --input-ref ",hap," ",leg," ",pre,"/1000GP_Phase3.sample",
                      " --output-log shapeit-check.log --chrX",sep=""))
    out.log("Run shapeit")
    sex <- data.frame(sample=config$sample,sex=ifelse(config$gender == "male",1,2))
    write.table(x = sex,file = "sex.ped",quote = F,sep = "\t",row.names = F,col.names = F)
    out.log.cmd(paste("shapeit -V shapeit.vcf -M ",
                      pre,"/genetic_map_chrX_nonPAR_combined_b37.txt",
                      " --input-ref ",hap," ",leg," ",pre,"/1000GP_Phase3.sample ",
                      " --exclude-snp shapeit-check.snp.strand.exclude -O shapeit-phased --chrX --input-sex sex.ped",
                      "&& shapeit -convert --input-haps shapeit-phased --output-vcf shapeit-phased.vcf",sep=""))
  }
  out.log("Reformat phased vcf")
  out.log.cmd("sed -i 's#^chr##g' shapeit-phased.vcf && bgzip shapeit-phased.vcf && tabix -f shapeit-phased.vcf.gz")

  setwd(config$analysis_path)
  out.log(paste("Done with chromosome ",chromosome,sep=""))
  out.log.cmd(paste("touch progress/.split_",chromosome,sep=""))
  setwd(owd)
}

local.barcode.function <- function(config,work.dir,type="") {
    owd <- getwd()
    setwd(work.dir)
    
    barcodes <- readLines("barcodes.txt")
    #generate bam file
    if(type == "bulk") {
      bam.region.selector("../../split/bulk/reads.reduced.sorted.barcode_indexed.bam",barcodes,"reads.barcodes.bam")
      out.log.cmd("python $LIRA_DIR/TXtools/x_toolbox.py -K -b reads.barcodes.bam -i ../../split/bulk/reads.bam -o reads.bam")
      bulk <- T
    } else {
      bam.region.selector("../../split/single_cell/reads.reduced.sorted.barcode_indexed.bam",barcodes,"reads.barcodes.bam")
      out.log.cmd("python $LIRA_DIR/TXtools/x_toolbox.py -K -b reads.barcodes.bam -i ../../split/single_cell/reads.bam -o reads.bam")
      bulk <- F
    }
    
    
    suppressWarnings(dir.create("tmp"))
    out.log.cmd("samtools sort -o reads.pos_sort.bam -T ./tmp/temp -O bam reads.bam && samtools index reads.pos_sort.bam")
    out.log.cmd("bedtools bamtobed -i reads.pos_sort.bam | bedtools intersect -a ../../split/sites.bed -b - | uniq > sites.bed")
    out.log.cmd(paste("$LIRA_DIR/scripts/linkage10x.py --bed sites.bed --bam reads.pos_sort.bam --fasta ",config$reference," > sites.by.barcode.txt",sep=""))
    sites.by.barcode <- read.table("sites.by.barcode.txt",sep="\t")
    names(sites.by.barcode) <- c("barcode","read","chromosome","position","reference_allele","observed_allele","variant_id")
    
    tmp <- paste(sites.by.barcode$read,sites.by.barcode$variant_id)
    read.dup <- names(which(table(tmp) > 1))
    tmp.site <- sites.by.barcode[tmp %in% read.dup,]
    sites.by.barcode <- sites.by.barcode[!(tmp %in% read.dup),]
    
    tmp.site <- base::split(tmp.site,tmp[tmp %in% read.dup])
    tmp.site <- lapply(tmp.site,function(x){unique(x)})
    tmp.site <- tmp.site[sapply(tmp.site,nrow) == 1]
    tmp.site <- do.call(rbind,tmp.site)
    sites.by.barcode <- rbind(sites.by.barcode,tmp.site)
    sites.by.barcode$read <- NULL
    
    tmp <- paste(sites.by.barcode$barcode,sites.by.barcode$variant_id)
    site.dup <- names(which(table(tmp) > 1))
    tmp.site <- sites.by.barcode[tmp %in% site.dup,]
    sites.by.barcode <- sites.by.barcode[!(tmp %in% site.dup),]
    
    tmp.site <- base::split(tmp.site,tmp[tmp %in% site.dup])
    tmp.site <- lapply(tmp.site,function(x){unique(x)})
    bad.bc <- unique(unlist(lapply(tmp.site,function(x){if(nrow(x) > 1){return(unique(x$barcode))}else{return(character(0))}})))
    
    tmp.site <- tmp.site[sapply(tmp.site,nrow) == 1]
    tmp.site <- do.call(rbind,tmp.site)
    sites.by.barcode <- rbind(sites.by.barcode,tmp.site)
    sites.by.barcode <- sites.by.barcode[!(sites.by.barcode$barcode %in% bad.bc),]
    
    sites.by.barcode <- sites.by.barcode[!(sites.by.barcode$observed_allele == "*"),]
    sites.by.barcode[,"targeted_alternate_allele"] <- gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\4",sites.by.barcode[,"variant_id"])
    sites.by.barcode$call <- NA
    sites.by.barcode$call[sites.by.barcode$observed_allele == sites.by.barcode$reference_allele] <- "R"
    sites.by.barcode$call[sites.by.barcode$observed_allele == sites.by.barcode$targeted_alternate_allele] <- "V"
    sites.by.barcode <- sites.by.barcode[!is.na(sites.by.barcode$call),]
    
    tmp <- sites.by.barcode[,c("variant_id","call")]
    
    sites.by.barcode <- base::split(sites.by.barcode,sites.by.barcode$barcode)
    sites.by.barcode <- sites.by.barcode[sapply(sites.by.barcode,nrow) > 1]
    linkage <- lapply(sites.by.barcode,function(x){
      x <- x[order(x$position),]
      pairs <- t(combn(1:nrow(x),m=2))
      obj <- list()
      for(i in 1:nrow(pairs)) {
        n <- paste(x$variant_id[pairs[i,1]],x$variant_id[pairs[i,2]],sep="~")
        obj[[n]] <- c(RR=0,RV=0,VR=0,VV=0)
        obj[[n]][paste(x$call[pairs[i,1]],x$call[pairs[i,2]],sep="")] <- 1
      }
      return(do.call(rbind,obj))
    })
    
    site.counts <- do.call(rbind,lapply(base::split(tmp,tmp$variant_id),function(x){data.frame(variant_id=x$variant_id[1],R=sum(x$call == "R"),V=sum(x$call == "V"))}))
    rownames(site.counts) <- site.counts$variant_id
    site.counts$variant_id <- NULL
    order.frame <- as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\2",rownames(site.counts)))
    o <- order(order.frame)
    site.counts <- site.counts[o,]
    
    pairs <- unique(unlist(lapply(linkage,function(x){rownames(x)})))
    order.frame <- data.frame(p1=as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)[~]([^;]*);([^;]*);([^;]*);([^;]*)","\\2",pairs)),
                              p2=as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)[~]([^;]*);([^;]*);([^;]*);([^;]*)","\\6",pairs)))
    o <- order(order.frame[,1],order.frame[,2])
    pairs <- pairs[o]
    tmp <- matrix(data = rep(0,4*length(pairs)),nrow = length(pairs),ncol=4,dimnames = list(pairs,c("RR","RV","VR","VV")))
    for(i in seq_along(linkage)) {
      tmp[rownames(linkage[[i]]),] <- tmp[rownames(linkage[[i]]),] + linkage[[i]]
    }
    
    linkage <- tmp
    obj <- list()
    obj$linkage <- linkage
    obj$site.counts <- site.counts
    linkage <- obj
    
    chr <- gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\1",rownames(site.counts))
    pos <- as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\2",rownames(site.counts)))
    bed <- data.frame(chromosome=chr,pos1=pos - 1, pos2=pos)
    write.table(x = bed,file = "sites.filtered.bed",quote=F,row.names = F,col.names = F,sep="\t")
    
    context.bed <- bed
    context.bed$pos1 <- context.bed$pos1 - 1
    context.bed$pos2 <- context.bed$pos2 + 1
    write.table(x = context.bed,file = "context.filtered.bed",quote=F,row.names = F,col.names = F,sep="\t")
    
    if(bulk) {
      out.log.cmd(paste("bcftools query -R sites.filtered.bed -s ",config$sample_bulk," -f '%CHROM;%POS;%REF;%ALT\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t%FILTER\\t%INFO/DBSNP_CAF\\t[%GT]\\t[%AO]\\t[%RO]\\n' ../../split/bulk/calls.id.vcf.gz > vcf_info.txt",sep=""))
    } else {
      out.log.cmd(paste("bcftools query -R sites.filtered.bed -s ",config$sample," -f '%CHROM;%POS;%REF;%ALT\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t%FILTER\\t%INFO/DBSNP_CAF\\t[%GT]\\t[%AO]\\t[%RO]\\n' ../../split/single_cell/calls.id.vcf.gz > vcf_info.txt",sep=""))
    }
    
    vcf.info <- read.table("vcf_info.txt",sep="\t",header=F)
    names(vcf.info) <- c("lira_id","chromosome","pos","ref","alt","dbsnp_id","filter","dbsnp_caf","genotype","alternate_allele_observations","reference_allele_observations")
    vcf.info$dbsnp_caf <- lapply(strsplit(vcf.info$dbsnp_caf,","),function(x){x[1]})
    vcf.info$dbsnp_caf[vcf.info$dbsnp_caf == "."] <- "1"
    vcf.info$dbsnp_caf <- as.numeric(vcf.info$dbsnp_caf)
    names(vcf.info)[names(vcf.info) == "dbsnp_caf"] <- "pop_ref_allele_frequency"
    
    tmp <- data.frame(lira_id=rownames(site.counts),
                      chromosome=chr,
                      pos=pos,
                      ref=gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\3",rownames(site.counts)),
                      alt=gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\4",rownames(site.counts)),
                      dbsnp_id=NA,
                      filter=NA,
                      pop_ref_allele_frequency=NA,
                      genotype=NA,
                      alternate_allele_observations=NA,
                      reference_allele_observations=NA)
    rownames(tmp) <- tmp$lira_id
    for(i in 1:ncol(vcf.info)) {
      tmp[vcf.info$lira_id,names(vcf.info)[i]] <- vcf.info[,names(vcf.info)[i]]
    }
    vcf.info <- tmp
    vcf.info <- vcf.info[vcf.info$lira_id %in% rownames(site.counts),]
    
    context <- system(paste("bedtools getfasta -fi ",config$reference," -bed context.filtered.bed -fo - | grep -v '>' | tr 'acgt' 'ACGT'",sep=""),intern=T)
    alt <- vcf.info$alt
    ind <- vcf.info$ref %in% c("A","G")
    
    context[ind] <- reverse.complement(context[ind])
    alt[ind] <- reverse.complement(alt[ind])
    vcf.info$context <- paste(context,alt,sep=">")
    vcf.info$lira_id <- NULL
    linkage$vcf.info <- vcf.info
    
    save(linkage,file="linkage.rda")
    setwd(dir = owd)
    
    tmp <- strsplit(work.dir,"/")[[1]]
    n <- length(tmp)
    progress.file <- paste(config$analysis_path,"/progress/.",ifelse(bulk,"bulk_",""),paste(tmp[n - 2],"_",tmp[n],sep=""),sep="")
    out.log.cmd(paste("touch ",progress.file,sep=""))
}

compare <- function(config, chromosome) {
  setwd(chromosome)
  
  suppressWarnings(dir.create("compare"))
  setwd("compare")
  
  process.local <- function(jobs.dir,done.tag,undone.message) {
    j <- list.files(jobs.dir)
    done.files <- paste(config$analysis_path,"/progress/.",done.tag,"_",j,sep="")
    done <- file.exists(done.files)
    if(any(!done)) {
      out.log("Jobs not done:")
      sapply(done.files[!done],function(x){out.log(x)})
      stop(undone.message)
    }
    f.linkage <- paste(jobs.dir,j,"/linkage.rda",sep="")
    tmp <- lapply(f.linkage,function(x){load(x); return(linkage)})
    linkage <- lapply(tmp,function(x){data.frame(x$linkage)})
    linkage <- lapply(linkage,function(x){
      linkage <- x
      linkage$site1 <- gsub("([^~]*)~([^~]*)","\\1",rownames(linkage))
      linkage$site2 <- gsub("([^~]*)~([^~]*)","\\2",rownames(linkage))
      linkage$pos1 <- as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\2",linkage$site1))
      linkage$pos2 <- as.numeric(gsub("([^;]*);([^;]*);([^;]*);([^;]*)","\\2",linkage$site2))
      linkage$pair_id <- rownames(linkage)
      rownames(linkage) <- NULL
      linkage$distance <- linkage$pos2 - linkage$pos1
      ind <- linkage$distance < 0
      replace <- linkage[ind,]
      replace$site1 <- linkage$site2[ind]
      replace$site2 <- linkage$site1[ind]
      replace$pos1 <- linkage$pos2[ind]
      replace$pos2 <- linkage$pos1[ind]
      replace$VR <- linkage$RV[ind]
      replace$RV <- linkage$VR[ind]
      replace$pair_id <- paste(replace$site1,replace$site2,sep="~")
      replace$distance <- -replace$distance
      linkage <- linkage[!ind,]
      linkage <- linkage[linkage$distance < global$MAX_DISTANCE_10X,]
      replace <- replace[replace$distance < global$MAX_DISTANCE_10X,]
      obj <- list()
      obj[[1]] <- linkage
      obj[[2]] <- replace
      return(obj)
    })
    linkage <- unlist(linkage,recursive=F)
    combined <- do.call(rbind,linkage)
    combined[,c("RR","RV","VR","VV")] <- 0
    combined <- unique(combined)
    rownames(combined) <- combined$pair_id
    for(i in seq_along(linkage)) {
      combined[linkage[[i]]$pair_id,c("RR","RV","VR","VV")] <- combined[linkage[[i]]$pair_id,c("RR","RV","VR","VV")] + linkage[[i]][,c("RR","RV","VR","VV")]
    }
    linkage <- combined
    linkage <- linkage[order(linkage$site1,linkage$site2),]
    rownames(linkage) <- 1:nrow(linkage)
    
    vcf.info <- do.call(rbind,lapply(tmp,function(x){tmp <- data.frame(x$vcf.info); tmp$id <- rownames(tmp); return(tmp)}))
    vcf.info <- unique(vcf.info[order(vcf.info$pos),])
    rownames(vcf.info) <- vcf.info$id
    vcf.info$id <- NULL
    
    site.counts <- lapply(tmp,function(x){tmp <- data.frame(x$site.counts); tmp$id <- rownames(tmp); return(tmp)})
    combined <- do.call(rbind,site.counts)
    combined$R <- 0
    combined$V <- 0
    combined <- unique(combined)
    rownames(combined) <- combined$id
    combined$id <- NULL
    for(i in seq_along(site.counts)) {
      combined[site.counts[[i]]$id,c("R","V")] <- combined[site.counts[[i]]$id,c("R","V")] + site.counts[[i]][,c("R","V")]
    }
    vcf.info$alternate_barcode_observations <- NA
    vcf.info$reference_barcode_observations <- NA
    vcf.info[rownames(combined),"alternate_barcode_observations"] <- combined$V
    vcf.info[rownames(combined),"reference_barcode_observations"] <- combined$R
    vcf.info$lira_id <- rownames(vcf.info)
    
    obj <- list()
    obj$linkage <- linkage
    obj$vcf.info <- vcf.info
    return(obj)
  }
  out.log("Loading single-cell info")
  single.cell <- process.local("../jobs/",chromosome,paste("Jobs still undone from plink.  See ",config$analysis_path,"/",chromosome,"/compare/log.txt",sep=""))
  
  out.log("Loading bulk info")
  bulk <- process.local("../jobs_bulk/",paste("bulk_",chromosome,sep=""),paste("Bulk jobs still undone from plink.  See ",config$analysis_path,"/",chromosome,"/compare/log.txt",sep=""))
  
  
  #make a site-specific table
  site.frame <- unique(rbind(single.cell$vcf.info[,c("chromosome","pos","ref","alt","dbsnp_id","pop_ref_allele_frequency","context","lira_id")],
                             bulk$vcf.info[,c("chromosome","pos","ref","alt","dbsnp_id","pop_ref_allele_frequency","context","lira_id")]))
  site.frame <- site.frame[order(site.frame$pos),]
  site.frame <- site.frame[!(is.na(site.frame$dbsnp_id)|is.na(site.frame$pop_ref_allele_frequency)|is.na(site.frame$context)),]
  rownames(site.frame) <- site.frame$lira_id
  site.frame$lira_id <- NULL
  
  for(col in c("filter","genotype","alternate_allele_observations","reference_allele_observations","alternate_barcode_observations","reference_barcode_observations")) {
    if(col %in% c("filter","genotype")) {
      fn <- as.character
    } else{
      fn <- as.numeric
    }
    sc.col <- paste("single_cell_",col,sep="")
    bulk.col <- paste("bulk_",col,sep="")
    site.frame[,sc.col] <- NA
    site.frame[rownames(single.cell$vcf.info),sc.col] <- single.cell$vcf.info[,col]
    
    site.frame[,bulk.col] <- NA
    site.frame[rownames(bulk$vcf.info),bulk.col] <- bulk$vcf.info[,col]
    
    site.frame[,sc.col] <- suppressWarnings(fn(site.frame[,sc.col]))
    site.frame[,bulk.col] <- suppressWarnings(fn(site.frame[,bulk.col]))
  }
  
  #only look for linkage near bulk het sites in 1KG
  #somatic candidates: bulk_genotype not het.
  site.frame$in.linkage <- rownames(site.frame) %in% c(single.cell$linkage$site1,single.cell$linkage$site2)
  site.frame$onek.bulk.het <- (site.frame$pop_ref_allele_frequency != 1) & (site.frame$bulk_genotype %in% c("0/1","1/0","0|1","1|0"))
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
}