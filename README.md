# E_phylacis_Genome_Assembly

This repository stores the whole pipeline of assembling *Eucalyptus x phylacis* genome from ONT reads. It summarises the steps from the original <a href="https://github.com/roblanf/MM_genome">repository</a>, which are mainly developed using `bash`. It covers read QC, haplotype-resolved assembly with `hifiasm`, organelle genome assembly, assignment of every contig to one of the two parental genomes (using *E. decipiens* and *E. virginea* as the parental references), and per-chromosome scaffolding against *E. regnans*. It is mainly developed and tested using Linux, so there might be incompatibilities using MacOS and Windows.

## Table of Content
- <a href="#prereqs">Prerequisites</a>
- <a href="#genpipe">General Pipeline</a>
- <a href="#refs">References</a>

## <a id="prereqs">Prerequisites</a>
This pipeline requires a number of software, R packages, and databases to run. We recommend you to use environment management system (e.g. `conda`) to install the prerequisites, but you can also use `install.packages()` built-in function in R or RStudio for the R packages.

### Software
|    Name        |                             Website                              |                             Anaconda                             |
| -------------- | :--------------------------------------------------------------: | :--------------------------------------------------------------: |
| barrnap        | <a href="https://github.com/tseemann/barrnap">Link</a>           | <a href="https://anaconda.org/bioconda/barrnap">Link</a>         |
| BCFtools       | <a href="https://github.com/samtools/bcftools">Link</a>          | <a href="https://anaconda.org/bioconda/bcftools">Link</a>        |
| BEDTools       | <a href="https://github.com/arq5x/bedtools2">Link</a>            | <a href="https://anaconda.org/bioconda/bedtools">Link</a>        |
| BLAST+         | <a href="https://blast.ncbi.nlm.nih.gov/">Link</a>               | <a href="https://anaconda.org/bioconda/blast">Link</a>           |
| BlobToolKit    | <a href="https://github.com/blobtoolkit/blobtoolkit">Link</a>    | <a href="https://anaconda.org/bioconda/blobtoolkit">Link</a>     |
| blobtk         | <a href="https://github.com/blobtoolkit/blobtk">Link</a>         | <a href="https://anaconda.org/bioconda/blobtk">Link</a>          |
| Chopper        | <a href="https://github.com/wdecoster/chopper">Link</a>          | <a href="https://anaconda.org/bioconda/chopper">Link</a>         |
| compleasm      | <a href="https://github.com/huangnengCSU/compleasm">Link</a>     | <a href="https://anaconda.org/bioconda/compleasm">Link</a>       |
| DIAMOND        | <a href="https://github.com/bbuchfink/diamond">Link</a>          | <a href="https://anaconda.org/bioconda/diamond">Link</a>         |
| GenomeScope 2  | <a href="https://github.com/tbenavi1/genomescope2.0">Link</a>    | <a href="https://anaconda.org/bioconda/genomescope2">Link</a>    |
| gfatools       | <a href="https://github.com/lh3/gfatools">Link</a>               | <a href="https://anaconda.org/bioconda/gfatools">Link</a>        |
| hifiasm        | <a href="https://github.com/chhylp123/hifiasm">Link</a>          | <a href="https://anaconda.org/bioconda/hifiasm">Link</a>         |
| HMMER          | <a href="http://hmmer.org/">Link</a>                             | <a href="https://anaconda.org/bioconda/hmmer">Link</a>           |
| HTSlib         | <a href="https://github.com/samtools/htslib">Link</a>            | <a href="https://anaconda.org/bioconda/htslib">Link</a>          |
| IQ-TREE 2      | <a href="http://www.iqtree.org">Link</a>                         | <a href="https://anaconda.org/bioconda/iqtree">Link</a>          |
| KMC            | <a href="https://github.com/refresh-bio/KMC">Link</a>            | <a href="https://anaconda.org/bioconda/kmc">Link</a>             |
| Kraken 2       | <a href="https://github.com/DerrickWood/kraken2">Link</a>        | <a href="https://anaconda.org/bioconda/kraken2">Link</a>         |
| Merqury        | <a href="https://github.com/marbl/merqury">Link</a>              | <a href="https://anaconda.org/bioconda/merqury">Link</a>         |
| meryl          | <a href="https://github.com/marbl/meryl">Link</a>                | <a href="https://anaconda.org/bioconda/meryl">Link</a>           |
| minimap2       | <a href="https://github.com/lh3/minimap2">Link</a>               | <a href="https://anaconda.org/bioconda/minimap2">Link</a>        |
| mosdepth       | <a href="https://github.com/brentp/mosdepth">Link</a>            | <a href="https://anaconda.org/bioconda/mosdepth">Link</a>        |
| NanoPlot       | <a href="https://github.com/wdecoster/NanoPlot">Link</a>         | <a href="https://anaconda.org/bioconda/nanoplot">Link</a>        |
| Oatk           | <a href="https://github.com/c-zhou/oatk">Link</a>                | <a href="https://anaconda.org/bioconda/oatk">Link</a>            |
| pigz           | <a href="https://zlib.net/pigz/">Link</a>                        | <a href="https://anaconda.org/conda-forge/pigz">Link</a>         |
| RagTag         | <a href="https://github.com/malonge/RagTag">Link</a>             | <a href="https://anaconda.org/bioconda/ragtag">Link</a>          |
| SAMtools       | <a href="https://github.com/samtools/samtools">Link</a>          | <a href="https://anaconda.org/bioconda/samtools">Link</a>        |
| SeqKit         | <a href="https://github.com/shenwei356/seqkit">Link</a>          | <a href="https://anaconda.org/bioconda/seqkit">Link</a>          |
| tidk           | <a href="https://github.com/tolkit/telomeric-identifier">Link</a>| <a href="https://anaconda.org/bioconda/tidk">Link</a>            |

BLAST+, DIAMOND, and Kraken 2 are only needed for the optional steps (`run_blast` and `run_kraken2` in `codes/1_main.Rmd`). The pipeline also uses `wget`, `curl`, and `tar`, which are available on most Linux systems.

### R packages
|    Name     |                               CRAN                               |                                   Anaconda                               |
| ----------- | :--------------------------------------------------------------: | :----------------------------------------------------------------------: |
| Biostrings  | <a href="https://bioconductor.org/packages/Biostrings">Link</a>  | <a href="https://anaconda.org/bioconda/bioconductor-biostrings">Link</a> |
| data.table  | <a href="https://cran.r-project.org/package=data.table">Link</a> | <a href="https://anaconda.org/conda-forge/r-data.table">Link</a>         |
| doSNOW      | <a href="https://cran.r-project.org/package=doSNOW">Link</a>     | <a href="https://anaconda.org/conda-forge/r-dosnow">Link</a>             |
| ggplot2     | <a href="https://cran.r-project.org/package=ggplot2">Link</a>    | <a href="https://anaconda.org/conda-forge/r-ggplot2">Link</a>            |
| htmlwidgets | <a href="https://cran.r-project.org/package=htmlwidgets">Link</a>| <a href="https://anaconda.org/conda-forge/r-htmlwidgets">Link</a>        |
| knitr       | <a href="https://cran.r-project.org/package=knitr">Link</a>      | <a href="https://anaconda.org/conda-forge/r-knitr">Link</a>              |
| log4r       | <a href="https://cran.r-project.org/package=log4r">Link</a>      | <a href="https://anaconda.org/conda-forge/r-log4r">Link</a>              |
| optparse    | <a href="https://cran.r-project.org/package=optparse">Link</a>   | <a href="https://anaconda.org/conda-forge/r-optparse">Link</a>           |
| plotly      | <a href="https://cran.r-project.org/package=plotly">Link</a>     | <a href="https://anaconda.org/conda-forge/r-plotly">Link</a>             |
| rmarkdown   | <a href="https://cran.r-project.org/package=rmarkdown">Link</a>  | <a href="https://anaconda.org/conda-forge/r-rmarkdown">Link</a>          |
| tidyverse   | <a href="https://cran.r-project.org/package=tidyverse">Link</a>  | <a href="https://anaconda.org/conda-forge/r-tidyverse">Link</a>          |

### Databases and reference genomes
|    Name                                 |                             Source                                               | Used by             | Download        |
| --------------------------------------- | :------------------------------------------------------------------------------- | :------------------ | :-------------- |
| BUSCO lineage `eudicotyledons_odb12`    | <a href="https://github.com/huangnengCSU/compleasm">Link</a>                     | `3_hifiasm.Rmd`     | by the pipeline |
| NCBI taxonomy (`new_taxdump`)           | <a href="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/">Link</a>        | `3_hifiasm.Rmd`     | by the pipeline |
| OatkDB (`magnoliopsida` HMMs)           | <a href="https://github.com/c-zhou/OatkDB">Link</a>                              | `4_organelle.Rmd`   | by the pipeline |
| barrnap `euk.hmm`                       | <a href="https://github.com/tseemann/barrnap">Link</a>                           | `5_assignment.Rmd`  | with barrnap    |
| *E. decipiens* genome                   | <a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_014182575.1/">Link</a> | `5_assignment.Rmd`  | manual          |
| *E. virginea* genome                    | <a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_014182375.1/">Link</a> | `5_assignment.Rmd`  | manual          |
| *E. regnans* genome                     | <a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_037349125.1/">Link</a> | `6_scaffolding.Rmd` | manual          |
| Kraken 2 PlusPFP database [optional]    | <a href="https://benlangmead.github.io/aws-indexes/k2">Link</a>                  | `2_qc.Rmd`          | manual          |
| NCBI `core_nt` [optional]               | <a href="https://ftp.ncbi.nlm.nih.gov/blast/db/">Link</a>                        | `3_hifiasm.Rmd`     | manual          |
| UniProt/Swiss-Prot (DIAMOND) [optional] | <a href="https://www.uniprot.org/">Link</a>                                      | `3_hifiasm.Rmd`     | manual          |

## <a id="genpipe">General Pipeline</a>
1. **Clone the Git repository** <br>
    ```
    git clone git@github.com:jeremiasivan/E_phylacis_Genome_Assembly.git
    ```

2. **Install the prerequisites** <br>
    - Create a new conda environment
        ```
        conda create -n ephylacis
        conda activate ephylacis
        ```
    - Install prerequisites
        ```
        conda install --strict-channel-priority -c conda-forge -c bioconda \
            barrnap bcftools bedtools blast blobtk blobtoolkit chopper 'compleasm>=0.2.7' diamond \
            genomescope2 gfatools hifiasm=0.25.0 hmmer htslib 'iqtree=2.*' kmc kraken2 merqury meryl \
            minimap2 mosdepth nanoplot oatk pandoc pigz \
            bioconductor-biostrings r-base r-data.table r-dosnow r-htmlwidgets r-knitr r-log4r r-optparse r-plotly r-rmarkdown r-tidyverse \
            ragtag samtools seqkit tidk
        ```
    - Download the reference genomes and the optional databases listed above

3. **Update the parameters in `codes/1_main.Rmd`** <br>

4. **Run the pipeline** <br>
    ```
    Rscript -e 'rmarkdown::render("codes/1_main.Rmd")'
    ```

---
## <a id="refs">References</a>
1. De Coster, W. & Rademakers, R. (<a href="https://doi.org/10.1093/bioinformatics/btad311">2023</a>). **NanoPack2: population-scale evaluation of long-read sequencing data**. *Bioinformatics*, *39*(5), btad311.

2. Adler, M. (<a href="https://zlib.net/pigz/">2015</a>). **pigz: A parallel implementation of gzip for modern multi-processor, multi-core machines**. *Software*.

3. Kokot, M., Długosz, M., & Deorowicz, S. (<a href="https://doi.org/10.1093/bioinformatics/btx304">2017</a>). **KMC 3: counting and manipulating k-mer statistics**. *Bioinformatics*, *33*(17), 2759–2761.

4. Ranallo-Benavidez, T.R., Jaron, K.S., & Schatz, M.C. (<a href="https://doi.org/10.1038/s41467-020-14998-3">2020</a>). **GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes**. *Nature Communications*, *11*, 1432.

5. Wood, D.E., Lu, J., & Langmead, B. (<a href="https://doi.org/10.1186/s13059-019-1891-0">2019</a>). **Improved metagenomic analysis with Kraken 2**. *Genome Biology*, *20*, 257.

6. Cheng, H., Concepcion, G.T., Feng, X., Zhang, H., & Li, H. (<a href="https://doi.org/10.1038/s41592-020-01056-5">2021</a>). **Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm**. *Nature Methods*, *18*, 170–175.

7. Cheng, H., Asri, M., Lucas, J., Koren, S., & Li, H. (<a href="https://doi.org/10.1038/s41592-024-02269-8">2024</a>). **Scalable telomere-to-telomere assembly for diploid and polyploid genomes with double graph**. *Nature Methods*, *21*, 967–970.

8. Li, H. (<a href="https://github.com/lh3/gfatools">2020</a>). **gfatools: tools for manipulating sequence graphs in the GFA and rGFA formats**. *Software*.

9. Huang, N. & Li, H. (<a href="https://doi.org/10.1093/bioinformatics/btad595">2023</a>). **compleasm: a faster and more accurate reimplementation of BUSCO**. *Bioinformatics*, *39*(10), btad595.

10. Li, H. (<a href="https://doi.org/10.1093/bioinformatics/bty191">2018</a>). **Minimap2: pairwise alignment for nucleotide sequences**. *Bioinformatics*, *34*(18), 3094–3100.

11. Danecek, P., et al. (<a href="https://doi.org/10.1093/gigascience/giab008">2021</a>). **Twelve years of SAMtools and BCFtools**. *GigaScience*, *10*(2), giab008.

12. Bonfield, J.K., et al. (<a href="https://doi.org/10.1093/gigascience/giab007">2021</a>). **HTSlib: C library for reading/writing high-throughput sequencing data**. *GigaScience*, *10*(2), giab007.

13. Shen, W., Sipos, B., & Zhao, L. (<a href="https://doi.org/10.1002/imt2.191">2024</a>). **SeqKit2: A Swiss army knife for sequence and alignment processing**. *iMeta*, *3*(3), e191.

14. Challis, R., Richards, E., Rajan, J., Cochrane, G., & Blaxter, M. (<a href="https://doi.org/10.1534/g3.119.400908">2020</a>). **BlobToolKit – interactive quality assessment of genome assemblies**. *G3: Genes, Genomes, Genetics*, *10*(4), 1361–1374.

15. Camacho, C., et al. (<a href="https://doi.org/10.1186/1471-2105-10-421">2009</a>). **BLAST+: architecture and applications**. *BMC Bioinformatics*, *10*, 421.

16. Buchfink, B., Reuter, K., & Drost, H.-G. (<a href="https://doi.org/10.1038/s41592-021-01101-x">2021</a>). **Sensitive protein alignments at tree-of-life scale using DIAMOND**. *Nature Methods*, *18*, 366–368.

17. Zhou, C., Brown, M., Blaxter, M., The Darwin Tree of Life Project Consortium, McCarthy, S.A., & Durbin, R. (<a href="https://doi.org/10.1186/s13059-025-03676-6">2025</a>). **Oatk: a de novo assembly tool for complex plant organelle genomes**. *Genome Biology*, *26*, 235.

18. Quinlan, A.R. & Hall, I.M. (<a href="https://doi.org/10.1093/bioinformatics/btq033">2010</a>). **BEDTools: a flexible suite of utilities for comparing genomic features**. *Bioinformatics*, *26*(6), 841–842.

19. Minh, B.Q., et al. (<a href="https://doi.org/10.1093/molbev/msaa015">2020</a>). **IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era**. *Molecular Biology and Evolution*, *37*(5), 1530–1534.

20. Pedersen, B.S. & Quinlan, A.R. (<a href="https://doi.org/10.1093/bioinformatics/btx699">2018</a>). **Mosdepth: quick coverage calculation for genomes and exomes**. *Bioinformatics*, *34*(5), 867–868.

21. Wheeler, T.J. & Eddy, S.R. (<a href="https://doi.org/10.1093/bioinformatics/btt403">2013</a>). **nhmmer: DNA homology search with profile HMMs**. *Bioinformatics*, *29*(19), 2487–2489.

22. Seemann, T. (<a href="https://github.com/tseemann/barrnap">2018</a>). **barrnap 0.9: rapid ribosomal RNA prediction**. *Software*.

23. Rhie, A., Walenz, B.P., Koren, S., & Phillippy, A.M. (<a href="https://doi.org/10.1186/s13059-020-02134-9">2020</a>). **Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies**. *Genome Biology*, *21*, 245.

24. Alonge, M., et al. (<a href="https://doi.org/10.1186/s13059-022-02823-7">2022</a>). **Automated assembly scaffolding using RagTag elevates a new tomato system for high-throughput genome editing**. *Genome Biology*, *23*, 258.

25. Brown, M.R., Gonzalez De la Rosa, P.M., & Blaxter, M. (<a href="https://doi.org/10.1093/bioinformatics/btaf049">2025</a>). **tidk: a toolkit to rapidly identify telomeric repeats from genomic datasets**. *Bioinformatics*, *41*(2), btaf049.

26. Pagès, H., et al. (<a href="https://doi.org/10.18129/B9.bioc.Biostrings">2026</a>). **Biostrings: Efficient manipulation of biological strings**. *R package*.

27. Barrett, T., et al. (<a href="https://doi.org/10.32614/CRAN.package.data.table">2026</a>). **data.table: Extension of 'data.frame'**. *R package*.

28. Microsoft Corporation & Weston, S. (<a href="https://cran.r-project.org/package=doSNOW">2022</a>). **doSNOW: Foreach Parallel Adaptor for the 'snow' Package**. *R package*.

29. Wickham, H. (<a href="https://doi.org/10.1007/978-3-319-24277-4">2016</a>). **ggplot2: Elegant Graphics for Data Analysis**. *Springer-Verlag New York*.

30. Vaidyanathan, R., et al. (<a href="https://doi.org/10.32614/CRAN.package.htmlwidgets">2025</a>). **htmlwidgets: HTML Widgets for R**. *R package*.

31. Xie, Y. (<a href="https://doi.org/10.32614/CRAN.package.knitr">2026</a>). **knitr: A General-Purpose Package for Dynamic Report Generation in R**. *R package*.

32. White, J.M. & Jacobs, A. (<a href="https://doi.org/10.32614/CRAN.package.log4r">2024</a>). **log4r: A Fast and Lightweight Logging System for R, Based on 'log4j'**. *R package*.

33. Davis, T.L. (<a href="https://doi.org/10.32614/CRAN.package.optparse">2026</a>). **optparse: Command Line Option Parser**. *R package*.

34. Sievert, C. (<a href="https://plotly-r.com">2020</a>). **Interactive Web-Based Data Visualization with R, plotly, and shiny**. *Chapman and Hall/CRC*.

35. Allaire, J.J., et al. (<a href="https://doi.org/10.32614/CRAN.package.rmarkdown">2026</a>). **rmarkdown: Dynamic Documents for R**. *R package*.

36. Wickham, H., et al. (<a href="https://doi.org/10.21105/joss.01686">2019</a>). **Welcome to the tidyverse**. *Journal of Open Source Software*, *4*(43), 1686.

37. Anthropic. (<a href="https://claude.ai/">2026</a>). Claude 5.1 Fable was used to review the code and to generate several functions in `codes/` and this README.

---
*Last update: 23 November 2026 by Jeremias Ivan*
