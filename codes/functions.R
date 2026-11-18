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

# function: run Hifiasm
f_hifiasm <- function(fn_fastq, dir_out, hifiasm_params, thread, exe_hifiasm) {
    cmd_hifiasm <- paste(exe_hifiasm,
                         "-o", file.path(dir_out, "hifiasm_asm"),
                         "-t", thread,
                         "--ont --dual-scaf",
                         hifiasm_params,
                         fn_fastq,
                         "2>&1 | tee", file.path(dir_out, "hifiasm_restart.log"))
    system(cmd_hifiasm)
}

# function: run minimap2 (output: paf)
f_minimap2_paf <- function(fn_ref, minimap_params, fn_query, fn_out, thread, exe_minimap2) {
    cmd_minimap2 <- paste(exe_minimap2,
                          minimap_params,
                          "-t", thread,
                          fn_ref,
                          fn_query,
                          ">", fn_out)
    system(cmd_minimap2)
}

# function: run minimap2 (output: bam)
f_minimap2_bam <- function(fn_ref, minimap_params, fn_query, fn_out, thread, exe_minimap2, exe_samtools) {
    cmd_minimap2 <- paste(exe_minimap2,
                          minimap_params,
                          "-t", thread,
                          fn_ref,
                          fn_query,
                          "|", 
                          exe_samtools, "sort -o", fn_out)
    system(cmd_minimap2)

    # index the BAM file
    system(paste(exe_samtools, "index -c", fn_out))
}

# function: run compleasm
f_compleasm <- function(fn_fasta, dir_out, dir_busco_lineage, thread, exe_compleasm) {
    cmd_compleasm <- paste(exe_compleasm, "run",
                           "-a", fn_fasta,
                           "-o", dir_out,
                           "-t", thread,
                           "-l eudicotyledons --odb odb12",
                           "-L", dir_busco_lineage)
    system(cmd_compleasm)
}

# function: change the headers of compleasm output table
f_update_compleasm_table <- function(fn_input, fn_output) {
    # open the input file
    lines <- readLines(fn_input)

    # change the first line to indicate the BUSCO version
    lines[1] <- sub("^# compleasm version is: .*", "# BUSCO version is: 6.1.0", lines[1])

    # change the headers
    lines <- sub("\tSingle\t",      "\tComplete\t",   lines)
    lines <- sub("\tInterspaced\t", "\tFragmented\t", lines)

    # write the modified file
    writeLines(lines, con=fn_output)
}

# function: create BlobTools output directory
f_blobtools_create <- function(fn_fasta, taxid, dir_taxdump, dir_out, exe_blobtools) {
    cmd_blobtools_create <- paste(exe_blobtools, "create",
                                  "--fasta", fn_fasta,
                                  "--taxid", taxid,
                                  "--taxdump", dir_taxdump,
                                  dir_out)
    system(cmd_blobtools_create)   
}

# function: add to BlobTools directory (BAM files)
f_blobtools_add_bam <- function(dir_blobtools, fn_bam, thread, exe_blobtools) {
    cmd_blobtools_add <- paste(exe_blobtools, "add",
                               "--cov", fn_bam,
                               "--thread", thread,
                               dir_blobtools)
    system(cmd_blobtools_add)
}

# function: add to BlobTools directory (busco)
f_blobtools_add_busco <- function(dir_blobtools, fn_busco, exe_blobtools) {
    cmd_blobtools_add <- paste(exe_blobtools, "add",
                               "--busco", fn_busco,
                               dir_blobtools)
    system(cmd_blobtools_add)
}

# function: run blastn for BlobTools
f_blastn <- function(fn_query, fn_db, fn_out, thread, exe_blastn) {
    cmd_blastn <- paste(exe_blastn, "-task megablast",
                        "-db", fn_db,
                        "-query", fn_query,
                        "-outfmt '6 qseqid staxids bitscore std'",
                        "-max_target_seqs 10 -max_hsps 1 -evalue 1e-25",
                        "-num_threads", thread,
                        "| sed 's/_sliding:[0-9]*-[0-9]*//'", 
                        ">", fn_out)
    system(cmd_blastn)
}

# function: run Diamond for BlobTools
f_diamond <- function(fn_query, fn_db, fn_out, thread, exe_diamond) {
    cmd_diamond <- paste(exe_diamond, "blastx", 
                         "--query", fn_query,
                         "--db", fn_db,
                         "--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                         "--sensitive --max-target-seqs 1 --evalue 1e-25",
                         "--threads", thread,
                         "| sed 's/_sliding:[0-9]*-[0-9]*//'",
                         ">", fn_out)
    system(cmd_diamond)
}

# function: add BLAST/DIAMOND results to BlobTools directory
f_blobtools_add_hits <- function(dir_blobtools, hits, dir_taxdump, exe_blobtools) {
    cmd_blobtools_add <- paste(exe_blobtools, "add",
                               hits,
                               "--taxrule bestsumorder",
                               "--taxdump", dir_taxdump,
                               dir_blobtools)
    system(cmd_blobtools_add)
}

# function: export BlobTools plots
f_blobtools_export <- function(dir_blobtools, dir_out, prefix, exe_blobtools, exe_blobtk) {
    # export table
    cmd_blobtools_table <- paste(exe_blobtools, "filter",
                                 "--table", file.path(dir_out, paste0(prefix, "_btk.tsv")),
                                 dir_blobtools)
    system(cmd_blobtools_table)

    # export plots
    system(paste(exe_blobtk, "plot",
                 "-v blob", 
                 "-d", dir_blobtools,
                 "-o", file.path(dir_out, paste0(prefix, "_blob.png"))))
    system(paste(exe_blobtk, "plot",
                 "-v snail",
                 "-d", dir_blobtools,
                 "-o", file.path(dir_out, paste0(prefix, "_snail.svg"))))
}

# function: run oatK
f_oatk <- function(fn_fastq, prefix, fn_ref_mito, fn_ref_pltd, thread, exe_oatk) {
    cmd_oatk <- paste(exe_oatk,
                      "-k 1001",
                      "-c 30",
                      "-t", thread,
                      "-p", fn_ref_pltd,
                      "-m", fn_ref_mito,
                      "-o", prefix,
                      fn_fastq)
    system(cmd_oatk)
}