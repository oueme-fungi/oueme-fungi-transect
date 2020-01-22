fwd <- list.files("sequences/trim/SH-2257_001", "R1[fr].trim.fastq.gz", recursive = TRUE, full.names = TRUE)
rev <- list.files("sequences/trim/SH-2257_001", "R2[fr].trim.fastq.gz", recursive = TRUE, full.names = TRUE)

dada2::plotQualityProfile(fwd, aggregate = TRUE)
dada2::plotQualityProfile(rev, aggregate = TRUE)
