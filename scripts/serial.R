#Runs jobs serially
#In LiRA, check.jobs is only called after submission, and submission here actually runs the scripts.
#Thus, check.jobs just returns false, and submit.jobs actually performs the analysis.

#A function to check job status.
#Given a vector of job identifiers (in the slurm case, names) this should return TRUE if any of the jobs are still running (or if the job fetch request fails)
#FALSE otherwise.
check.jobs <- function(job.names) {
  return(F)
}

#A function submit jobs
#Given a vector of bash commands and corresponding job IDs (in the slurm case, names), this should submit jobs and return TRUE if the submission was successful, FALSE otherwise.
submit.jobs <- function(bash.commands, job.ids) {
  for(b in bash.commands) {
    system(b)
  }
  return(T)
}