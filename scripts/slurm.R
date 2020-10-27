
#A function to check job status.
#Given a vector of job identifiers (in the slurm case, names) this should return TRUE if any of the jobs are still running (or if the job fetch request fails)
#FALSE otherwise.
check.jobs <- function(job.names) {
  if(length(job.names) == 0) {
    return(F)
  }
  tmp <- suppressWarnings(system(paste("squeue -u $LOGNAME --name=",paste(job.names,collapse=",")," --format=\"%j\" 2>&1",sep=""),intern=T))
  if(is.na(tmp[1])) {
    return(T)
  } else {
    if(tmp[1] == "NAME") {
      tmp <- tmp[-1]
      if(length(tmp) == 0) {
        return(F)
      } else {
        return(T)
      }
    } else {
      return(T)
    }
  }
}

#A function submit jobs
#Given a vector of bash commands and corresponding job IDs (in the slurm case, names), this should submit jobs and return TRUE if the submission was successful, FALSE otherwise.
#fix_cluster_issue: added to address R overloading concerns raised by HMS IT.  This initiaties copying of a local copy of R to the cluster node on which a job is running if it is not already there.

submit.jobs <- function(bash.commands, job.ids, memory, wall.time, fix_cluster_issue=T) {
  if(missing(memory)) {
    memory <- "32"
  }
  if(missing(wall.time)) {
    wall.time <- 60
  }
  if(wall.time > 7200) {
    queue <- "long"
  } else if(wall.time > 720) {
    queue <- "medium"
  } else {
    queue <- "short"
  }
  result <- sapply(seq_along(bash.commands),function(i){
    if(fix_cluster_issue) {
      add <- ""
    } else {
      add <- "./checkR.sh; LD_LIBRARY_PATH=/tmp/R-3.3.1/lib PATH=/tmp/R-3.3.1/bin/Rscript:$PATH"
    }
    cmd <- paste("sbatch -p ",queue," -o sbatch.out -e sbatch.error -c 1 --mem-per-cpu=",memory,"G -t ",wall.time," -J ",job.ids[i]," --wrap='",add,bash.commands[i],"' 2>&1",sep="")
    return(system(cmd,intern=T))
  })
  if(all(grepl("Submitted batch job",result))) {
    return(T)
  } else {
    print(result)
    return(F)
  }
}
