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

# function: run minimap2 for short reads (output: bam)
f_minimap2_shortreads <- function(fn_ref, minimap_params, fn_fastq_r1, fn_fastq_r2, fn_out, thread, exe_minimap2, exe_samtools, exe_bedtools) {
    cmd_minimap2 <- paste(exe_minimap2,
                          minimap_params,
                          "-t", thread,
                          fn_ref,
                          fn_fastq_r1, fn_fastq_r2,
                          "|", exe_samtools, "sort -o", fn_out)
    system(cmd_minimap2)

    # index the BAM file
    system(paste(exe_samtools, "index", fn_out))
}

# function: extract regions with low coverage
f_lowcov_bed <- function(fn_bam, fn_bed, min_depth, exe_samtools) {
    # calculate depth for every position in the BAM file
    df_depth <- data.table::fread(cmd = paste(exe_samtools, "depth -a", fn_bam), header = FALSE,
                                  col.names = c("contig", "pos", "depth"), colClasses = list(character = 1))

    # extract regions with low depth
    df_low <- df_depth %>%
                  filter(depth < min_depth) %>%
                  arrange(contig, pos) %>%
                  group_by(contig) %>%
                  mutate(run = cumsum(is.na(lag(pos)) | pos != lag(pos) + 1)) %>% 
                  group_by(contig, run) %>%
                  summarise(start = min(pos) - 1, end = max(pos), .groups = "drop") %>%
                  select(contig, start, end)
    
    # save as BED file
    data.table::fwrite(df_low, fn_bed, sep = "\t", col.names = FALSE)
}

# function: variant calling
f_variant_calling <- function(fn_ref, fn_bam, ploidy, fn_vcf, exe_bcftools) {
    cmd_variant_calling <- paste(exe_bcftools, "mpileup",
                                 "-f", fn_ref,
                                 "-q 20 -Q 20 -a AD,DP",
                                 "-Ou", fn_bam,
                                 "|", exe_bcftools, "call -m --ploidy", ploidy, "-Ou",
                                 "|", exe_bcftools, "norm -f", fn_ref, "-Oz -o", fn_vcf)
    system(cmd_variant_calling)

    # index the VCF file
    system(paste(exe_bcftools, "index", fn_vcf))

    # extract SNPs
    system(paste(exe_bcftools, "view -V indels -i 'QUAL>=30'", fn_vcf, "-Oz -o", gsub(".vcf.gz", ".snps.vcf.gz", fn_vcf)))
    system(paste(exe_bcftools, "index", gsub(".vcf.gz", ".snps.vcf.gz", fn_vcf)))
}

# function: call consensus sequence from VCF file
f_call_consensus <- function(fn_ref, fn_bed, fn_vcf, fn_consensus, exe_bcftools) {
    cmd_call_consensus <- paste(exe_bcftools, "consensus",
                                "-f", fn_ref,
                                "-m", fn_bed,
                                "-o", fn_consensus,
                                fn_vcf)
    system(cmd_call_consensus)
}

# function: run IQ-TREE2
f_iqtree2 <- function(fn_fasta, thread, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_fasta,
                         "-T", thread,
                         "-B 1000")
    system(cmd_iqtree2)
}

# function: read PAF file
f_read_paf <- function(fn_paf) {
    # function: extract optional tags from PAF lines
    f_extract_tag <- function(ls_lines, tag) {
        pattern <- paste0("\t", stringr::str_escape(tag), "([^\t]*)")
        stringr::str_match(ls_lines, pattern)[, 2]   # capture group; NA where absent
    }

    # read the file
    ls_lines <- readLines(fn_paf)

    # keep lines with the 12 mandatory columns
    ls_lines <- ls_lines[grepl("^([^\t]*\t){11}", ls_lines)]

    # extract mandatory columns
    df_paf <- data.table::fread(text=sub("^(([^\t]*\t){11}[^\t]*).*$", "\\1", ls_lines),
                                sep="\t", header=FALSE, quote="",
                                col.names=c("qname", "qlen", "qstart", "qend", "strand",
                                            "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq"))

    # extract optional tags
    df_paf$de <- as.numeric(f_extract_tag(ls_lines, "de:f:"))
    df_paf$tp <- f_extract_tag(ls_lines, "tp:A:")

    return(df_paf)
}

# function: calculate total length of overlapping regions
f_calculate_merged_len <- function(start, end) {
    # index the start and end positions
    idx <- order(start, end)
    start <- start[idx]
    end <- end[idx]

    # create a new block when the start position is greater than the previous maximum end position
    block <- cumsum(start > lag(cummax(end), default = -1L))

    # calculate the total length of the merged regions
    merged_len <- sum(tapply(end, block, max) - tapply(start, block, min))
    return(merged_len)
}

# function: extract a list of longer contigs that contain a smaller contig
f_extract_top_targets <- function(start, end, identity, target, n=5) {
    str_targets <- tibble(start, end, identity, target) %>%
                    arrange(start, end, identity, target) %>%
                    mutate(pos = row_number()) %>%
                    group_by(target) %>%
                    summarise(len = sum(end - start), first_pos = min(pos)) %>%
                    arrange(desc(len), first_pos) %>%
                    slice_head(n = n) %>%
                    pull(target) %>%
                    paste(collapse = ",")
    return(str_targets)
}

# function: calculate mean read depth over the regions of a BED file
f_mosdepth <- function(fn_bam, fn_bed, prefix, thread, exe_mosdepth) {
    cmd_mosdepth <- paste(exe_mosdepth,
                          "-t", thread,
                          "--by", fn_bed,
                          "--no-per-base",
                          prefix, fn_bam)

    system(cmd_mosdepth)
}

# function: run nhmmer
f_nhmmer <- function(fn_hmm, fn_fasta, fn_tblout, thread, exe_nhmmer) {
    cmd_nhmmer <- paste(exe_nhmmer,
                        "--rna --cpu", thread,
                        "--noali -E 1e-6 --w_length 8000",
                        "--tblout", fn_tblout,
                        "-o /dev/null",
                        fn_hmm, fn_fasta)
    system(cmd_nhmmer)
}

# function: read an nhmmer --tblout file
f_read_nhmmer_tblout <- function(fn_tblout) {
    # define the column names
    ls_cols   <- c("target", "target_acc", "query", "query_acc", "hmm_from", "hmm_to", "ali_from", "ali_to",
                   "env_from", "env_to", "sq_len", "strand", "evalue", "score", "bias")

    # open the file
    ls_lines  <- readLines(fn_tblout)

    # skip commented lines
    ls_lines  <- ls_lines[!startsWith(ls_lines, "#")]

    # split each line based on whitespace
    ls_fields <- strsplit(trimws(ls_lines), "\\s+")

    # convert lines into data.frame and update headers
    df_hits <- as_tibble(do.call(rbind, lapply(ls_fields, `[`, seq_along(ls_cols))), .name_repair = ~ ls_cols)
    
    return(df_hits)
}

# function: update headers of E. decipiens and E. virginea reference genomes
f_update_headers <- function(ref_dec, ref_vir, ref_dec_out, ref_vir_out) {
    # map the original headers to new headers
    ls_map_decipiens <- c(setNames(sprintf("Chr%02d", 1:11), sprintf("CM0246%02d.1", 9:19)),  "CM024632.1" = "Pltd")
    ls_map_virginea  <- c(setNames(sprintf("Chr%02d", 1:11), sprintf("CM0245%02d.1", 14:24)), "CM024525.1" = "Pltd")

    # update the headers of E. decipiens
    seq          <- Biostrings::readBStringSet(ref_dec)                 # open the file
    ls_accession <- sub("\\s.*$", "", names(seq))                       # extract the chromosome names
    is_mapped    <- ls_accession %in% names(ls_map_decipiens)
    names(seq)[is_mapped] <- ls_map_decipiens[ls_accession[is_mapped]]  # update the chromosome names
    Biostrings::writeXStringSet(seq, filepath=ref_dec_out)

    # update the headers of E. virginea
    seq          <- Biostrings::readBStringSet(ref_vir)                 # open the file
    ls_accession <- sub("\\s.*$", "", names(seq))                       # extract the chromosome names
    is_mapped    <- ls_accession %in% names(ls_map_virginea)
    names(seq)[is_mapped] <- ls_map_virginea[ls_accession[is_mapped]]   # update the chromosome names
    Biostrings::writeXStringSet(seq, filepath=ref_vir_out)
}

# function: score every contig of a PAF against one reference
f_score_paf <- function(fn_paf, min_aln, min_id) {
    # read the PAF file
    df_paf <- f_read_paf(fn_paf)

    # filter alignments
    df_paf <- df_paf %>%
                filter(is.na(tp) | tp == "P") %>%                                # keep primary alignments only
                filter(alen >= min_aln) %>%                                      # keep alignments >= min_aln
                mutate(identity = ifelse(is.na(de), nmatch / alen, 1 - de)) %>%  # calculate identity
                filter(identity >= min_id)                                       # keep alignments with identity >= min_id

    # summarise per contig
    df_score <- df_paf %>%
                    group_by(contig = qname) %>%
                    summarise(len      = first(qlen),
                              aligned  = f_calculate_merged_len(qstart, qend) / first(qlen),
                              identity = 100 * sum(identity * alen) / sum(alen),
                              chrom    = names(which.max(tapply(qend - qstart, tname, sum))),
                              n        = n())

    return(df_score)
}

# function: find the best partners of every hap1 and hap2 contig
f_partners <- function(fn_paf) {
    # read the PAF file
    df_paf <- f_read_paf(fn_paf)

    # calculate the total alignment length of every pair
    df_sum <- df_paf %>%
        group_by(qname, tname) %>%
        summarise(alen = sum(alen))

    # extract the best partner of every hap1 contig, and of every hap2 contig
    df_best_h1 <- df_sum %>%
        group_by(contig = qname) %>%
        summarise(partner = tname[which.max(alen)])
    df_best_h2 <- df_sum %>%
        group_by(contig = tname) %>%
        summarise(partner = qname[which.max(alen)])

    # extract contig lengths
    df_len_h1 <- df_paf %>% distinct(contig = qname, len = qlen)
    df_len_h2 <- df_paf %>% distinct(contig = tname, len = tlen)

    return(list(best_h1 = df_best_h1, best_h2 = df_best_h2, len_h1 = df_len_h1, len_h2 = df_len_h2))
}

# function: read blobtools stats table
f_read_blob <- function(fn_blob) {
    # create an empty data.table with the expected columns
    df_blob <- data.table::data.table(contig=character(), gc=numeric(), cov=numeric(), phylum=character())
    if (is.null(fn_blob) || !file.exists(fn_blob) || file.size(fn_blob) == 0) {
        return(df_blob)
    }

    # extract individual rows and headers
    ls_rows <- strsplit(readLines(fn_blob), "\t", fixed=TRUE)
    ls_header <- tolower(trimws(ls_rows[[1]]))

    # extract the column indices
    col_id     <- grep("identifiers", ls_header)
    col_gc     <- grep("gc", ls_header)
    col_cov    <- grep("cov", ls_header)
    col_phylum <- grep("phylum", ls_header)

    # remove the header row
    ls_rows <- ls_rows[-1]

    # keep only rows with the expected number of columns and a valid contig name
    ls_rows <- ls_rows[sapply(ls_rows, function(row) {
        length(row) >= max(col_id, col_gc, col_cov, col_phylum) && tolower(trimws(row[1])) != "index"
    })]

    # return an empty data.table if no valid rows are found
    if (length(ls_rows) == 0) {
        return(df_blob)
    }

    # store the entries as data.frame
    df_blob <- data.table::data.table(contig = trimws(sapply(ls_rows, `[`, col_id)),
                                      gc     = as.numeric(sapply(ls_rows, `[`, col_gc)),
                                      cov    = as.numeric(sapply(ls_rows, `[`, col_cov)),
                                      phylum = trimws(sapply(ls_rows, `[`, col_phylum)))

    # skip rows with invalid GC or coverage
    df_blob <- df_blob[!is.na(gc) & !is.na(cov)]

    # convert GC to percentage
    df_blob$gc <- ifelse(df_blob$gc <= 1, df_blob$gc * 100, df_blob$gc)

    return(df_blob)
}

# function: remove contigs from a FASTA file
f_remove_contigs <- function(fn_in, ls_remove, fn_out, exe_seqkit, exe_samtools) {
    # save the list of contigs to remove
    fn_list <- paste0(fn_out, ".removed.txt")
    writeLines(ls_remove, fn_list)

    # remove the contigs from the FASTA file
    if (length(ls_remove) > 0) {
        system(paste(exe_seqkit, "grep -v -f", fn_list, fn_in, ">", fn_out))
    } else {
        file.copy(fn_in, fn_out, overwrite = TRUE)
    }

    # index the output FASTA file
    system(paste(exe_samtools, "faidx", fn_out))
}

# function: read mosdepth regions.bed.gz file
f_read_depth <- function(fn_regions) {
    # open the file
    ls_lines <- readLines(gzfile(fn_regions))
    if (length(ls_lines) == 0) {
        return(numeric(0))
    }

    # convert the lines into a table
    df_depth <- data.table::fread(text = paste0(paste(ls_lines, collapse = "\n"), "\n"),
                                  sep = "\t", header = FALSE, colClasses = list(character = 1))
    return(setNames(df_depth[[4]], df_depth[[1]]))
}

# function: call the parent of each contig
f_call <- function(dec_aligned, vir_aligned, dec_id, vir_id, min_aligned, margin) {
  case_when(is.na(dec_aligned) & is.na(vir_aligned)                                ~ "no alignment",
            pmax(coalesce(dec_aligned, 0), coalesce(vir_aligned, 0)) < min_aligned ~ "unaligned",
            is.na(dec_aligned)                                                     ~ "virginea",
            is.na(vir_aligned)                                                     ~ "decipiens",
            dec_id - vir_id >= margin                                              ~ "decipiens",
            dec_id - vir_id <= -margin                                             ~ "virginea",
            .default                                                               = "ambiguous")
}

# function: interleave hap1 and hap2 contigs
f_interleave_partners <- function(df_rows, hap_name) {
    # extract the rows of the specified haplotype
    df_self <- df_rows %>%
                   filter(hap == hap_name) %>%
                   arrange(desc(length)) %>%
                   mutate(pair = row_number(), role = 1)

    # extract the rows of the partner haplotype
    df_partner <- df_self %>%
                      filter(partner != "-") %>%
                      select(pair, partner_name = partner) %>%
                      inner_join(df_rows %>% filter(hap != hap_name), by = join_by(partner_name == contig)) %>%
                      mutate(contig = paste0("\u2192 ", partner_name), role = 2) %>%   # "→ h2tg000123l"
                      select(-partner_name)

    # bind the two data.frames together
    df_out <- bind_rows(df_self, df_partner) %>%
                  arrange(pair, role) %>%
                  select(-pair, -role)

    return(df_out)
}

# function: count k-mers
f_meryl_count <- function(ls_fn_reads, fn_db, k, thread, exe_meryl) {
    cmd_meryl <- paste(exe_meryl,
                       paste0("k=", k),
                       "count",
                       paste0("threads=", thread),
                       "output", fn_db,
                       paste(ls_fn_reads, collapse = " "))
    system(cmd_meryl)
}

# function: filter k-mers
f_meryl_filter <- function(fn_db, fn_out, filter_params, thread, exe_meryl) {
    cmd_meryl <- paste(exe_meryl,
                       paste0("threads=", thread),
                       filter_params,
                       "output", fn_out,
                       fn_db)
    system(cmd_meryl)
}

# function: compare two k-mer databases
f_meryl_operator <- function(operator, fn_db1, fn_db2, fn_out, thread, exe_meryl) {
    cmd_meryl <- paste(exe_meryl, operator,
                       paste0("threads=", thread),
                       fn_db1, fn_db2,
                       "output", fn_out)
    system(cmd_meryl)
}

# function: run Merqury trio
f_merqury_trio <- function(db_reads, db_ref1, db_ref2, hap1, hap2, prefix, exe_merqury) {
    cmd_merqury <- paste(exe_merqury,
                         db_reads,
                         db_ref1, db_ref2,
                         hap1, hap2,
                         prefix)
    system(cmd_merqury)
}

# function: check every contig's assignment against the evidence that was not used to make it
f_check_all <- function(df_contig, ls_opposite, ls_parent_chrom, min_cov) {
    # check every contig's assignment
    df_contig %>%
        mutate(partner_parent = unname(ls_parent_chrom[partner]),
               check_intrinsic = case_when(
                   !is.na(dv_call) & !is.na(hapmer_call) & dv_call != hapmer_call      ~ "conflicting evidence",
                   is.na(dv_call) | is.na(hapmer_call)                                 ~ "no intrinsic call",
                   .default = "consistent"),
               check_pair = case_when(
                   rule == "unassigned"                                                ~ "-",
                   is.na(partner_parent) | partner_cov < min_cov                       ~ "no partner evidence",
                   partner_parent == unname(ls_opposite[parent])                       ~ "consistent",
                   .default = "CHECK"),
               check_synteny = case_when(
                   is.na(dec_chrom) | is.na(vir_chrom) | dec_chrom == "-" | vir_chrom == "-"        ~ "no synteny evidence",
                   dec_chrom == vir_chrom                                                           ~ "consistent",
                   .default = "CHECK"),
               check = case_when(
                   rule == "unassigned"                              ~ "-",
                   check_pair == "CHECK" | check_synteny == "CHECK"  ~ "CHECK",
                   .default = "consistent")) %>%
        select(-partner_parent)
}

# function: one row per contig with every piece of evidence behind its assignment
f_evidence_table <- function(df_contig, df_dvscore, df_hapmers) {
    # create a data.frame with the DV scores
    df_dv <- df_dvscore %>%
                transmute(contig,
                          dec_id,
                          vir_id,
                          delta = round(dec_id - vir_id, 3),
                          dec_aligned,
                          vir_aligned,
                          gc,
                          cov,
                          phylum)

    # create a data.frame with the k-mer counts
    df_kmer <- df_hapmers %>%
        transmute(contig        = Contig,
                  kmer_vir      = as.numeric(.data[["E_virginea.unique"]]),
                  kmer_dec      = as.numeric(.data[["E_decipiens.unique"]]),
                  kmer_total    = kmer_vir + kmer_dec,
                  kmer_frac_vir = if_else(kmer_total > 0, round(kmer_vir / kmer_total, 3), NA_real_),
                  kmer_frac_dec = if_else(kmer_total > 0, round(kmer_dec / kmer_total, 3), NA_real_))

    # combine all the evidence into one table
    df_contig %>%
        left_join(df_dv,   by = "contig") %>%
        left_join(df_kmer, by = "contig") %>%
        select(contig, hap, length, parent, rule, check, check_intrinsic, check_pair, check_synteny,
               partner, partner_cov, partner_identity,
               dv_call, dec_id, vir_id, delta, dec_aligned, vir_aligned, dec_chrom, vir_chrom,
               hapmer_call, kmer_vir, kmer_dec, kmer_total, kmer_frac_vir, kmer_frac_dec,
               intrinsic, gc, cov, phylum)
}

# function: run RagTag
f_ragtag <- function(fn_ref, fn_query, ragtag_params, dir_out, thread, exe_ragtag, fn_exclude = NULL) {
    cmd_ragtag <- paste(exe_ragtag, "scaffold",
                        "-t", thread,
                        ragtag_params)
    
    if (!is.null(fn_exclude)) {
        cmd_ragtag <- paste(cmd_ragtag, "-j", fn_exclude)
    }

    cmd_ragtag <- paste(cmd_ragtag, "-o", dir_out, fn_ref, fn_query)
    system(cmd_ragtag)
}

# function: chromosome binning
f_chromosome_bins <- function(df_assign, df_dvscore, parent, col_chrom, ls_anchors) {
    # filter invalid contigs
    ls_ref_chrom <- df_dvscore %>%
                        select(contig, ref_chrom = all_of(col_chrom)) %>%
                        filter(!is.na(ref_chrom), ref_chrom != "-") %>%
                        deframe()

    # extract the contigs per chromosome
    df_bins <- df_assign %>%
                   filter(parent == !!parent) %>%
                   mutate(anchor    = if_else(rule != "chromosome" & partner %in% ls_anchors, partner, contig),
                          ref_chrom = unname(ls_ref_chrom[anchor])) %>%
                   filter(!is.na(ref_chrom)) %>%
                   select(ref_chrom, contig, rule, anchor) %>%
                   arrange(ref_chrom, desc(rule == "chromosome"), contig)

    return(df_bins)
}

# function: scaffold contigs per chromosome
f_scaffold_pass <- function(df_bins, fn_asm, fn_ref, dir_iter, ragtag_params, fn_exclude, thread, exe_seqkit, exe_ragtag) {
    # create the output directories
    dir_bins <- file.path(dir_iter, "bins")
    dir_scaf <- file.path(dir_iter, "scaffolds")
    for (dir in c(dir_bins, dir_scaf)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }

    # output files
    fn_chromosomes <- file.path(dir_iter, "chromosomes.fa")
    fn_unplaced    <- file.path(dir_iter, "unplaced.fa")
    fn_agp         <- file.path(dir_iter, "combined.agp")

    # create empty output files
    file.create(fn_chromosomes, fn_unplaced, fn_agp)

    # iterate over chromosomes
    for (chr in sort(unique(df_bins$ref_chrom))) {
        # input files
        fn_ids   <- file.path(dir_bins, paste0(chr, ".txt"))
        fn_query <- file.path(dir_bins, paste0(chr, ".fa"))
        fn_chref <- file.path(dir_bins, paste0(chr, ".ref.fa"))

        # output directory
        dir_chr  <- file.path(dir_scaf, chr)
        if (!dir.exists(dir_chr)) {
            dir.create(dir_chr, recursive = TRUE)
        }

        # extract the contigs of this chromosome and scaffold them
        writeLines(df_bins$contig[df_bins$ref_chrom == chr], con=fn_ids)
        system(paste(exe_seqkit, "grep -f", fn_ids, fn_asm, ">", fn_query))
        system(paste(exe_seqkit, "grep -p", chr,    fn_ref, ">", fn_chref))
        f_ragtag(fn_ref        = fn_chref,
                 fn_query      = fn_query,
                 ragtag_params = ragtag_params,
                 dir_out       = dir_chr,
                 thread        = thread,
                 exe_ragtag    = exe_ragtag,
                 fn_exclude    = fn_exclude)

        # extract the scaffolds of the chromosome
        fn_scaf <- file.path(dir_chr, "ragtag.scaffold.fasta")
        system(paste(exe_seqkit, "grep -p",    paste0(chr, "_RagTag"), fn_scaf, ">>", fn_chromosomes))
        system(paste(exe_seqkit, "grep -v -p", paste0(chr, "_RagTag"), fn_scaf, ">>", fn_unplaced))
        system(paste("grep -v '^#'", file.path(dir_chr, "ragtag.scaffold.agp"), ">>", fn_agp))
    }

    # extract contigs outside every bin
    fn_unbinned <- file.path(dir_iter, "unbinned.txt")
    ls_all      <- sub("\\s.*$", "", names(Biostrings::fasta.seqlengths(fn_asm)))
    writeLines(setdiff(ls_all, df_bins$contig), con=fn_unbinned)
    system(paste(exe_seqkit, "grep -f", fn_unbinned, fn_asm, ">>", fn_unplaced))

    return(list(fn_chromosomes = fn_chromosomes, fn_unplaced = fn_unplaced, fn_agp = fn_agp))
}

# function: run tidk for telomere motifs
f_tidk_search <- function(fn_fasta, prefix, dir_out, telomere_motif, exe_tidk) {
    cmd_tidk <- paste(exe_tidk, "search",
                      "-s", telomere_motif,
                      "-o", prefix,
                      "-d", dir_out,
                      fn_fasta, "> /dev/null 2>&1")
    system(cmd_tidk)

    # return the output file
    return(file.path(dir_out, paste0(prefix, "_telomeric_repeat_windows.tsv")))
}

# function: read an AGP file
f_read_agp <- function(fn_agp) {
    # set the column names
    ls_cols <- c("object", "object_beg", "object_end", "part_number", "component_type",
                 "component_id", "component_beg", "component_end", "orientation")

    # open the file
    ls_lines <- readLines(fn_agp)
    ls_lines <- ls_lines[!startsWith(ls_lines, "#") & nzchar(ls_lines)]

    # convert to data.frame
    df_agp <- data.table::fread(text = c(paste(ls_cols, collapse = "\t"), ls_lines), sep = "\t", header = TRUE,
                                colClasses = list(character = c(1, 5, 6, 7, 8, 9)))

    return(df_agp)
}

# function: check for telomere arrays at the ends of scaffolds
f_telomere_check <- function(fn_windows, fn_agp, fn_out, fn_ends_out, fn_exclude_out, min_copies = 100, min_anchor = 1e7, edge = 1e5, prefix = "") {
    # read the AGP file
    df_comp <- f_read_agp(fn_agp) %>%
                    filter(component_type == "W") %>%
                    transmute(scaffold = object,
                              start    = object_beg,
                              end      = object_end,
                              contig   = component_id,
                              length   = as.numeric(component_end) - as.numeric(component_beg) + 1)
    
    # extract the anchor contig
    df_scaf <- df_comp %>%
                    group_by(scaffold) %>%
                    filter(any(length >= min_anchor)) %>%
                    summarise(L         = max(end),
                              anchor    = first(contig[length >= min_anchor]),
                              anc_start = min(start[length >= min_anchor]),
                              anc_end   = max(end[length >= min_anchor])) %>%
                    arrange(desc(L))

    # extract telomere arrays
    df_win <- data.table::fread(fn_windows, select = 1:4, col.names = c("scaffold", "window", "fwd", "rev"), colClasses = list(character = 1)) %>%
                  mutate(scaffold = sub(paste0("^", prefix), "", scaffold))
    
    # extract the window size used for tidk
    wsize  <- df_win %>%
                group_by(scaffold) %>%
                filter(n() > 1) %>%
                summarise(w = diff(sort(window))[1]) %>%
                pull(w) %>%
                head(1)

    # collapse adjacent windows into arrays
    df_arr <- df_win %>%
                filter(fwd + rev >= min_copies) %>%
                arrange(scaffold, window) %>%
                group_by(scaffold) %>%
                mutate(array = cumsum(is.na(lag(window)) | window > lag(window) + 2 * wsize)) %>%
                group_by(scaffold, array) %>%
                summarise(start = min(window) - wsize + 1, end = max(window), copies = sum(fwd + rev), fwd = sum(fwd))

    # function: find the host contig of a telomere array
    f_host <- function(s, pos) {
        i <- which(df_comp$scaffold == s & df_comp$start <= pos & pos <= df_comp$end)
        if (length(i)) {
            df_comp$contig[i[1]]
        } else {
            NA_character_
        }
    }

    # function: check for telomere arrays at one end of a scaffold
    f_end <- function(df_a, df_c, anc_lo) {
        # extract candidate arrays
        cand <- df_a %>% filter(start < anc_lo + edge)
        if (nrow(cand) == 0) {
            return(tibble(telomere = FALSE, copies = NA, at_scaffold_end = NA, n_beyond = 0, beyond_bp = 0, beyond = NA_character_))
        }

        # find the nearest array to the anchor
        nearest <- cand[which.max(cand$start), ]
        i_host  <- which(df_c$start <= min(nearest$end, anc_lo + edge) & min(nearest$end, anc_lo + edge) <= df_c$end)

        # find the contigs that lie beyond the nearest array
        beyond  <- ifelse (length(i_host), df_c %>% filter(end < df_c$start[i_host[1]]), df_c[0, ])

        # return a tibble with the results
        tibble(telomere = TRUE,
               copies = nearest$copies,
               at_scaffold_end = nearest$start <= edge,
               n_beyond = nrow(beyond),
               beyond_bp = sum(beyond$length),
               beyond = ifelse (nrow(beyond) > 0, paste(beyond$contig, collapse = ","), NA_character_))
    }

    # check every scaffold for telomere arrays at its ends
    ls_ends <- list()
    for (i in seq_len(nrow(df_scaf))) {
        # extract the scaffold and its length
        s <- df_scaf$scaffold[i]
        L <- df_scaf$L[i]

        # extract the arrays and components of this scaffold
        df_a <- df_arr  %>% filter(scaffold == s)
        df_c <- df_comp %>% filter(scaffold == s)

        # function: mirror the coordinates of a scaffold
        mirror <- function(d) {
            d %>% mutate(start0 = start, start = L - end, end = L - start0) %>% select(-start0)
        }

        # check the left end and the right end of this scaffold
        ls_ends[[i]] <- bind_rows(left  = f_end(df_a, df_c, df_scaf$anc_start[i]),
                                  right = f_end(mirror(df_a), mirror(df_c), L - df_scaf$anc_end[i]),
                                  .id   = "side") %>%
                            mutate(scaffold = s,
                                   anchor   = df_scaf$anchor[i],
                                   .before  = 1)
    }
    df_ends <- bind_rows(ls_ends)

    # store results in a data.frame
    df_rows <- df_arr %>%
                    filter(scaffold %in% df_scaf$scaffold) %>%
                    left_join(df_scaf %>% select(scaffold, L, anc_start, anc_end), by = "scaffold") %>%
                    mutate(forward_fraction = sprintf("%.2f", fwd / copies),
                           component        = coalesce(mapply(f_host, scaffold, start), "-"),
                           terminal         = if_else(start <= edge | end >= L - edge, "yes", "no"),
                           internal_anchor  = if_else(anc_start + edge < start & end < anc_end - edge, "yes", "no")) %>%
                    arrange(match(scaffold, df_scaf$scaffold), start) %>%
                    select(scaffold, array_start = start, array_end = end, copies, forward_fraction, component, terminal, internal_anchor)

    # extract the list of contigs that lie beyond the telomere arrays
    ls_exclude <- as.character(unlist(strsplit(df_ends$beyond[!is.na(df_ends$beyond)], ",")))
    writeLines(ls_exclude, fn_exclude_out)

    # write the results to output files
    data.table::fwrite(df_rows, fn_out,      sep = "\t", quote = FALSE)
    data.table::fwrite(df_ends, fn_ends_out, sep = "\t", quote = FALSE, na = "-")
    
    return(invisible(list(arrays = df_rows, ends = df_ends, exclude = ls_exclude)))
}

# function: extract chromosome-scale contigs
f_find_anchors <- function(df_fai, ls_df_paf_ref, min_ref_cov = 0.2) {
    # calculate the reference coverage of every contig
    df_ref_cov <- bind_rows(ls_df_paf_ref, .id = "reference") %>%
                    filter(alen >= 10000, identity >= 0.9, startsWith(tname, "Chr")) %>%
                    group_by(reference, contig = qname, ref_chrom = tname) %>%
                    summarise(ref_cov = f_calculate_merged_len(tstart, tend) / first(tlen)) %>%
                    group_by(contig) %>%
                    slice_max(ref_cov, n = 1, with_ties = FALSE)

    # return a data.frame with the contig lengths and reference coverage
    df_anchors <- df_fai %>%
                    left_join(df_ref_cov, by = "contig") %>%
                    mutate(ref_cov = coalesce(ref_cov, 0),
                        anchor = ref_cov >= min_ref_cov) %>%
                    arrange(desc(ref_cov))
    
    return(df_anchors)
}