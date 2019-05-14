
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
submit.jobs <- function(bash.commands, job.ids) {
  result <- sapply(seq_along(bash.commands),function(i){
    return(system(paste("sbatch -p short -o sbatch.out -e sbatch.error -n 1 --mem-per-cpu=16G -t 720 -J ",job.ids[i]," --wrap='",bash.commands[i],"' 2>&1",sep=""),intern=T))
  })
  if(all(grepl("Submitted batch job",result))) {
    return(T)
  } else {
    print(result)
    return(F)
  }
}
