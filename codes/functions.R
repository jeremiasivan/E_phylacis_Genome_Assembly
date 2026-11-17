# functions for running Eucalyptus x phylacis genome assembly pipeline

# function: run NanoPlot
f_nanoplot <- function(fn_fastq, dir_out, thread, title, exe_nanoplot) {
    cmd_nanoplot <- paste(exe_nanoplot, 
                          "-t", thread,
                          "--fastq", fn_fastq,
                          "-o", dir_out,
                          "--title", title)
    system(cmd_nanoplot)
}

# function: run pigz (decompression) -> Chopper (filtering) -> bgzip (compression)
f_chopper <- function(fn_fastq, min_q, min_len, fn_fastq_filtered, thread, exe_pigz, exe_chopper, exe_bgzip) {
    cmd_chopper <- paste(exe_pigz, "-dc -p", thread, fn_fastq, "|",
                         exe_chopper, "-q", min_q, "-l", min_len, "|",
                         exe_bgzip, "-@", thread, ">", fn_fastq_filtered)
    system(cmd_chopper)
}

# function: parse NanoPlot summary file
f_parse_nanostats <- function(file) {
    # read in the summary file and trim whitespace
    lines <- trimws(readLines(file))

    # extract individual metrics
    total_bases         <- trimws(unlist(strsplit(grep("^Total bases:",         lines, value=T), split=":"))[2])
    n_reads             <- trimws(unlist(strsplit(grep("^Number of reads:",     lines, value=T), split=":"))[2])
    read_length_mean    <- trimws(unlist(strsplit(grep("^Mean read length:",    lines, value=T), split=":"))[2])
    read_length_median  <- trimws(unlist(strsplit(grep("^Median read length:",  lines, value=T), split=":"))[2])
    read_length_N50     <- trimws(unlist(strsplit(grep("^Read length N50:",     lines, value=T), split=":"))[2])
    read_length_stdev   <- trimws(unlist(strsplit(grep("^STDEV read length:",   lines, value=T), split=":"))[2])
    read_quality_mean   <- trimws(unlist(strsplit(grep("^Mean read quality:",   lines, value=T), split=":"))[2])
    read_quality_median <- trimws(unlist(strsplit(grep("^Median read quality:", lines, value=T), split=":"))[2])
    q10                 <- trimws(unlist(strsplit(grep("^>Q10:",                lines, value=T), split=":"))[2])
    q15                 <- trimws(unlist(strsplit(grep("^>Q15:",                lines, value=T), split=":"))[2])
    q20                 <- trimws(unlist(strsplit(grep("^>Q20:",                lines, value=T), split=":"))[2])
    q25                 <- trimws(unlist(strsplit(grep("^>Q25:",                lines, value=T), split=":"))[2])
    q30                 <- trimws(unlist(strsplit(grep("^>Q30:",                lines, value=T), split=":"))[2])

    # return the extracted metrics
    ls_output <- data.table::data.table(total_bases         = total_bases,
                                        n_reads             = n_reads,
                                        read_length_mean    = read_length_mean,
                                        read_length_median  = read_length_median,
                                        read_length_N50     = read_length_N50,
                                        read_length_stdev   = read_length_stdev,
                                        read_quality_mean   = read_quality_mean, 
                                        read_quality_median = read_quality_median,
                                        q10                 = q10,
                                        q15                 = q15,
                                        q20                 = q20,
                                        q25                 = q25,
                                        q30                 = q30)

    return(ls_output)
}

# function: run KMC and GenomeScope2
f_kmc_genomescope <- function(fn_fastq, kmer, dir_out, thread, exe_kmc, exe_kmctools, exe_genomescope) {
    # store the path to the filtered reads
    writeLines(fn_fastq, file.path(dir_out, "files.txt"))

    # run KMC
    cmd_kmc <- paste(exe_kmc,
                     paste0("-k", kmer),
                     paste0("-t", thread),
                     "-ci1 -cs10000",
                     paste0("@", file.path(dir_out, "files.txt")),
                     file.path(dir_out, "kmc_db"),
                     dir_out)
    system(cmd_kmc)

    # generate the k-mer histogram
    cmd_kmc_tools <- paste(exe_kmctools, "transform",
                           file.path(dir_out, "kmc_db"),
                           "histogram",
                           file.path(dir_out, "histogram.txt"),
                           "-ci1 -cs10000")
    system(cmd_kmc_tools)

    # run GenomeScope2
    cmd_genomescope <- paste(exe_genomescope,
                             "-i", file.path(dir_out, "histogram.txt"),
                             "-o", dir_out,
                             "-k", kmer,
                             "-p", 2)
    system(cmd_genomescope)
}

# function: run Kraken2
f_kraken2 <- function(fn_fastq, dir_db, confidence, dir_out, thread, exe_kraken2) {
    cmd_kraken2 <- paste(exe_kraken2, 
                         "--db", dir_db,
                         "--confidence", confidence,
                         "--threads", thread,
                         "--gzip-compressed",
                         "--report", file.path(dir_out, "kraken2_report.txt"),
                         "--report-minimizer-data",
                         "--use-names",
                         "--output", file.path(dir_out, "kraken2_output.txt"),
                         fn_fastq)
    system(cmd_kraken2)
}