# miRge miRNA counting

This document describes our usage of miRge2.0 for counting miRNAs in our small RNA sequencing libraries.

### Install miRge2.0 and download zebrafish libraries

Download zebrafish libraries containing sequences, coordinates and indexes from [miRge2.0](https://github.com/mhalushka/miRge) to a directory named "mirge_libs".
*Libraries used were rebuilt 05-06-2018 using miRBase22.0 according to GitHub and we downloaded them 03-02-2020*

```
wget -O zebrafish.tar.gz https://jh.box.com/shared/static/nwn7jzn5ekgm51k7jlk43a6h75aasgr1.gz
tar xvzf zebrafish.tar.gz
```
Install miRge2.0 via conda following the instructions [here](https://github.com/mhalushka/miRge). 

### Run miRge2.0 in R using the package parallel

Launch miRge alignment/quantification for 14806R small RNA fastq files, using the trimmed fastq files that are described in *FastqProcessing.md*. There are 22 samples and 2 lanes so 44 fastq files total. Create a command vector for each fastq file in R and use a cluster to run multiple commands at once. 

> miRge requires unzipped fastq files.

```
# Load libraries
library(data.table)
library(parallel)

# Load sample info table
fq_info = fread("../../14806R_fastq_info.txt")

# Add file_id column, to provide a readable file-specific label.
fq_info[, file_id:=paste(gnomex_id, flowcell_id, fastq_type, sep="_")]

# Subset the info table, to select only INSERT fastq files (described in FastqProcessing.md)
insert_info = fq_info[fastq_type == "INSERT"]

# How many cpu cores should _each_ miRge instance use?
n_cpu_cores = 8

# Construct vector of bash miRge commands, one command for each fastq file.
command_vec = paste(
    "mkdir -p", insert_info$file_id, ";",

    "/opt/anaconda/anaconda2/bin/miRge2.0 annotate",
    "   -s", file.path("../../processed_fastq_unzip", insert_info$fastq_filename),
    "   -d miRBase",
    "   -o", insert_info$file_id,
    "   -pb /usr/local/bowtie-1.1.1/",
    "   -lib ../mirge_libs_20200302/",
    "   -sp zebrafish",
    "   -ai",
    "   -cpu", n_cpu_cores,
    "   >", file.path(insert_info$file_id, "log.txt"), "2>&1")

cl = makeCluster(7)

parLapply(cl=cl, X=command_vec, fun=system)

stopCluster(cl)
```
There was an issue with one mirbase ID when the argument -gff was used. Issue was published in GitHub [here](https://github.com/mhalushka/miRge/issues/27).


### Combine miRge counts for all files

There are multiple files within the output folder for each fastq file, including the counts abd alignment statistics. We want to take the counts file for each sample and combine into one long form datatable in R for downstream analyses.

```
library(data.table)

# Search recursive through current directory to get a
# vector of full paths of all files named miR.Counts.csv
counts_file_vec = list.files(path=".", pattern="miR.Counts.csv",
                             full.names=TRUE, recursive=TRUE)

# Create empty list to add count data.tables
res_list = list()

# In a loop, load and process each counts file, and add
# to results list.
for (i in seq_along(counts_file_vec)) {
    file_path = counts_file_vec[i]
    tmp = fread(file_path) # Load count file in as data.table named tmp.
    fq_filename = names(tmp)[2] # Grab the file name these counts came from,
                                # it's the name of the second column of tmp.
    setnames(tmp, c("mirge_mir_id", "count")) # Change columns names.
    set(tmp, j="fastq_filename", value=fq_filename) # Add new column indicating which
                                                    # fastq file these counts belong to.
    tmp = tmp[!mirge_mir_id %in% "miRNAtotal"] # Remove the row containing the total count.
    res_list[[i]] = tmp
}

full_tab = rbindlist(res_list)

fwrite(full_tab, file="14806R_mirge_counts_long_20200310.txt", sep="\t")
```

>We have also run miRge using the miRGeneDB v2.0 database. The long form results table will be included in this repository for reference.




