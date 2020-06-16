is.slurm <- function() {
	# Is this a slurm job?
	slurm <- any(grepl("SLURM", names(Sys.getenv())))
	if (slurm) {
		# SLURM job notes - sent to job_*.info
		nThreads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
		jobID <- as.integer(Sys.getenv("SLURM_JOBID"))
		info <- as.matrix(Sys.getenv())
		idx <- grepl("SLURM", rownames(info))
		df <- info[idx, ]
		return(slurm)
	} else {
		return(slurm)
	}
}
