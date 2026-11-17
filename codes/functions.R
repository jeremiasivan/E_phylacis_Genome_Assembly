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