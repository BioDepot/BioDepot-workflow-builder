#!/usr/bin/env Rscript
# This script extracts specified samples from read counts dataset and
# experiment design table.

###################### Define internally used functions ######################

# Determine script invocation method.
getInvokeMethod <- function()
{
    # Get all command arguments including R system arguments.
    cmd_args <- commandArgs(FALSE)

    # Get the value of "--file" option.
    matched_pattern <- regexpr("(?<=^--file=).+", cmd_args, perl=TRUE)
    prog_cmd <- regmatches(cmd_args, matched_pattern)
    n_prog_cmd <- length(prog_cmd)

    # Get the availability of "--args" option.
    args_opt_avail <- "--args" %in% cmd_args

    # Determine invocation method based on n_prog_cmd and args_opt_avail.
    invoke_method <- NULL
    if(n_prog_cmd == 0) invoke_method <- "R"
    else if(n_prog_cmd == 1) invoke_method <- "Rscript"
    else invoke_method <- "unknown"

    # Return invocation method.
    return(invoke_method)
}

# Check the command line of script invocation by Rscript.
checkInvokeCommandLine <- function(invoke_method)
{
    # Initialize checking status.
    cmd_status <- TRUE

    # Only check command line for Rscript invocation.
    if(invoke_method == "Rscript")
    {
        # Get all command arguments including R system arguments.
        cmd_args <- commandArgs(FALSE)
        # Get the value of "--file" option.
        matched_pattern <- regexpr("(?<=^--file=).+", cmd_args, perl=TRUE)
        prog_cmd <- regmatches(cmd_args, matched_pattern)
        prog_name <- basename(prog_cmd)
        # Get the command line of script invocation by Rscript.
        cmd_args <- commandArgs(TRUE)
        if(length(cmd_args) < 1)
        {
            cat(paste("Usage:", prog_cmd, "[Counts Directory] [Function Directory]\n"))
            cat("\n")
            cat("Arguments:\n\n")
            cat("\t[Counts Directory] is a directory containing read counts and experiment design.\n\n")
            cat("\t[Function Directory] is a directory with all rquired R functions (Optional).\n\n")
            cmd_status <- FALSE
        }
    }
    else
    {
        warning("Script command line is only available for Rscript invocation!")
        # Comment out the following line to allow running this script in an interactive R console.
        #cmd_status <- FALSE
    }

    # Return checking status.
    return(cmd_status)
}

# Obtain the path of read counts directory.
getCountsDir <- function(invoke_method="Rscript", default_dir=getwd())
{
    # invoke_method: "Rscript" or "R".

    counts_dir <- NULL

    # Determine the function/config path.
    if(invoke_method == "Rscript")
    {
        # For the invocation by Rscript at system terminal.
        cmd_args <- commandArgs(TRUE)
        if(length(cmd_args) > 0) counts_dir <- cmd_args[1]
        else
        {
            warning("No counts directory is provided!")
            counts_dir <- default_dir
        }
        if(is.null(counts_dir)) warning("counts_dir is NULL when invoke_method is \"Rscript\"!")
    }
    else if(invoke_method == "R")
    {
        # For the invocation by R command or "source" at R terminal.
        counts_dir <- default_dir
        if(is.null(counts_dir)) warning("counts_dir is NULL when invoke_method is \"R\"!")
    }
    else warning("invoke_method must be either \"R\" or \"Rscript\"!")

    # Return the function/config path.
    return(counts_dir)
}

# Obtain the path of R function scripts.
getFuncDir <- function(invoke_method="Rscript", default_dir=getwd())
{
    # invoke_method: "Rscript" or "R".

    func_dir <- NULL

    # Determine the function/config path.
    if(invoke_method == "Rscript")
    {
        # For the invocation by Rscript at system terminal.
        # First of all, check the command line for user-specified function path.
        cmd_args <- commandArgs(TRUE)
        if(length(cmd_args) > 1) func_dir <- cmd_args[2]
        else
        {
            # Secondly, retrieve function path from the path of this main program.
            cmd_args <- commandArgs(FALSE)
            # Get the value of "--file" option.
            matched_pattern <- regexpr("(?<=^--file=).+", cmd_args, perl=TRUE)
            prog_cmd <- regmatches(cmd_args, matched_pattern)
            func_dir <- dirname(prog_cmd)
        }
        if(is.null(func_dir))
        {
            # Thirdly, use the default value for function path.
            warning("No argument is assigned to func_dir when invoke_method is \"Rscript\"!")
            func_dir <- default_dir
        }
    }
    else if(invoke_method == "R")
    {
        # For the invocation by R command or "source" at R terminal.
        func_dir <- default_dir
        if(is.null(func_dir)) warning("func_dir is NULL when invoke_method is \"R\"!")
    }
    else warning("invoke_method must be either \"R\" or \"Rscript\"!")

    # Return the function/config path.
    return(func_dir)
}

# Set default directories for differential comparison analysis.
setDefaultDEMDirs <- function(data_set="default")
{
    # Initialize default directories.
    user_home <- Sys.getenv("HOME")
    if(data_set == "default")
    {
        default_func_dir <- file.path(user_home, "LINCS/DEG/Programs")
        default_counts_dir <- file.path(user_home, "LINCS/DEG/Repo/Counts")
    }
    else
    {
        default_func_dir <- NULL
        default_counts_dir <- NULL
    }

    # Copy all variables from current function's environment to its parent's environment.
    for(obj in ls()) assign(obj, get(obj,environment()), parent.env(environment()))
}

################### Set up basic computational environment ###################

# Initialize default directories.
# Set default directories for differential comparison analysis.
setDefaultDEMDirs()

# Determine script invocation method.
invoke_method <- getInvokeMethod()
if(invoke_method == "unknown") stop("Script must be invoked by either R or Rscript!")

# Obtain invoked R script command.
stopifnot(checkInvokeCommandLine(invoke_method))

# Determine directory path of read counts data.
counts_dir <- getCountsDir(invoke_method=invoke_method, default_dir=default_counts_dir)
stopifnot(!is.null(counts_dir) && length(counts_dir)==1 && is.character(counts_dir))
if(!dir.exists(counts_dir)) stop(paste0("Counts directory \"", counts_dir, "\" doesn't exist!"))
counts_dir <- normalizePath(counts_dir)
# Determine directory path of all required R function scripts.
func_dir <- getFuncDir(invoke_method=invoke_method, default_dir=default_func_dir)
stopifnot(!is.null(func_dir) && length(func_dir)==1 && is.character(func_dir))
if(!dir.exists(func_dir)) stop(paste0("Function directory \"", func_dir, "\" doesn't exist!"))
func_dir <- normalizePath(func_dir)

# Load the function for extracting read counts and experiment design.
library(tools)
source(file.path(func_dir, "convertReadCountsFileName.R"), local=TRUE)
source(file.path(func_dir, "extractGeneExpSamples.R"), local=TRUE)

# Extract specified samples from read counts dataset and experiment design table.

write_file <- TRUE
ext_name <- "tsv"

# Plate-1 dataset (20150409).
sample_id_name <- "RNAseq"
sample_id_series <- "20150409"
sample_id <- paste(sample_id_name, sample_id_series, sep="_")
read_counts_file <- paste(sample_id, "unq.refseq.umi", "dat", sep=".")
read_counts_file <- file.path(counts_dir, read_counts_file)
exprt_design_file <- paste(paste(sample_id_name,"Experiment-Design",sep="-"), sample_id_series, "txt", sep=".")
exprt_design_file <- file.path(counts_dir, exprt_design_file)
if(file.exists(read_counts_file) && file.exists(exprt_design_file))
{
    sample_type <- "LINCS"
    cat("Extracting", paste(sample_id_name,sample_id_series,sample_type,sep="-"), "samples...\n")
    read_counts_name <- "Read-Counts"
    read_counts_file_new <- convertReadCountsFileName(read_counts_file, read_counts_name)
    file.copy(read_counts_file, read_counts_file_new)
    read_counts_datasets_20150409_lincs <- extractGeneExpSamples(read_counts_file=read_counts_file_new, mol_name_col=1, sample_name_filter="T96s4_", exprt_design_file=exprt_design_file, sample_names=c(paste("H",6:8,sep="")), sample_keep=FALSE, write_file=write_file, extra_name=sample_type, ext_name=ext_name, swap_series_type=TRUE, func_dir=func_dir)
    cat("Done\n")
}

# Plate-2 dataset (20150503).
sample_id_name <- "RNAseq"
sample_id_series <- "20150503"
sample_id <- paste(sample_id_name, sample_id_series, sep="_")
read_counts_file <- paste(sample_id, "unq.refseq.umi", "dat", sep=".")
read_counts_file <- file.path(counts_dir, read_counts_file)
exprt_design_file <- paste(paste(sample_id_name,"Experiment-Design",sep="-"), sample_id_series, "txt", sep=".")
exprt_design_file <- file.path(counts_dir, exprt_design_file)
if(file.exists(read_counts_file) && file.exists(exprt_design_file))
{
    sample_type <- "LINCS"
    cat("Extracting", paste(sample_id_name,sample_id_series,sample_type,sep="-"), "samples...\n")
    read_counts_name <- "Read-Counts"
    read_counts_file_new <- convertReadCountsFileName(read_counts_file, read_counts_name)
    file.copy(read_counts_file, read_counts_file_new)
    read_counts_datasets_20150503_lincs <- extractGeneExpSamples(read_counts_file=read_counts_file_new, mol_name_col=1, sample_name_filter="T96s2_", exprt_design_file=exprt_design_file, sample_names=c(paste("H",7:12,sep="")), sample_keep=FALSE, write_file=write_file, extra_name=sample_type, ext_name=ext_name, swap_series_type=TRUE, func_dir=func_dir)
    cat("Done\n")
}

# Plate-3 dataset (20150712).
sample_id_name <- "RNAseq"
sample_id_series <- "20150712"
sample_id_series_set2 <- paste(sample_id_series, "Set2", sep="-")
sample_id_set2 <- paste(sample_id_name, sample_id_series_set2, sep="_")
read_counts_file_set2 <- paste(sample_id_set2, "unq.refseq.umi", "dat", sep=".")
read_counts_file_set2 <- file.path(counts_dir, read_counts_file_set2)
exprt_design_file_set2 <- paste(paste(sample_id_name,"Experiment-Design",sep="-"), sample_id_series_set2, "txt", sep=".")
exprt_design_file_set2 <- file.path(counts_dir, exprt_design_file_set2)
sample_id_series_set3 <- paste(sample_id_series, "Set3", sep="-")
sample_id_set3 <- paste(sample_id_name, sample_id_series_set3, sep="_")
read_counts_file_set3 <- paste(sample_id_set3, "unq.refseq.umi", "dat", sep=".")
read_counts_file_set3 <- file.path(counts_dir, read_counts_file_set3)
exprt_design_file_set3 <- paste(paste(sample_id_name,"Experiment-Design",sep="-"), sample_id_series_set3, "txt", sep=".")
exprt_design_file_set3 <- file.path(counts_dir, exprt_design_file_set3)
if(file.exists(read_counts_file_set2) && file.exists(exprt_design_file_set2) && file.exists(read_counts_file_set3) && file.exists(exprt_design_file_set3))
{
    read_counts_name <- "Read-Counts"
    exprt_design_name <- "Experiment-Design"
    sample_type <- "LINCS"
    cat("Extracting", paste(sample_id_name,sample_id_series,sample_type,sep="-"), "samples...\n")

    # Extract dataset 2.
    read_counts_file_set2_new <- convertReadCountsFileName(read_counts_file_set2, read_counts_name)
    file.copy(read_counts_file_set2, read_counts_file_set2_new)
    read_counts_datasets_20150712_set2_lincs <- extractGeneExpSamples(read_counts_file=read_counts_file_set2_new, mol_name_col=1, sample_name_filter="T384s1_", exprt_design_file=exprt_design_file_set2, sample_names=c(paste("P",19:24,sep="")), sample_keep=FALSE, write_file=write_file, extra_name=sample_type, ext_name=ext_name, swap_series_type=TRUE, func_dir=func_dir)
    # Extract dataset 3.
    read_counts_file_set3_new <- convertReadCountsFileName(read_counts_file_set3, read_counts_name)
    file.copy(read_counts_file_set3, read_counts_file_set3_new)
    read_counts_datasets_20150712_set3_lincs <- extractGeneExpSamples(read_counts_file=read_counts_file_set3_new, mol_name_col=1, sample_name_filter="T384s1_", exprt_design_file=exprt_design_file_set3, sample_names=c(paste("P",19:24,sep="")), sample_keep=FALSE, write_file=write_file, extra_name=sample_type, ext_name=ext_name, swap_series_type=TRUE, func_dir=func_dir)

    # Merge dataset 2 and 3 to a single dataset.
    read_counts_datasets_20150712_lincs <- NULL
    read_counts_datasets_20150712_lincs$exprt_design <- read_counts_datasets_20150712_set2_lincs$exprt_design
    read_counts_datasets_20150712_lincs$read_counts <- read_counts_datasets_20150712_set2_lincs$read_counts + read_counts_datasets_20150712_set3_lincs$read_counts
    read_counts_file <- paste(paste(sample_id_name,read_counts_name,sep="-"), sample_type, sample_id_series, ext_name, sep=".")
    read_counts_file <- file.path(counts_dir, read_counts_file)
    exprt_design_file <- paste(paste(sample_id_name,exprt_design_name,sep="-"), sample_type, sample_id_series, ext_name, sep=".")
    exprt_design_file <- file.path(counts_dir, exprt_design_file)
    write.table(read_counts_datasets_20150712_lincs$read_counts, file=read_counts_file, quote=FALSE, sep="\t", na="")
    write.table(read_counts_datasets_20150712_lincs$exprt_design, file=exprt_design_file, quote=FALSE, sep="\t", na="", row.names=FALSE)

    cat("Done\n")
}

# Plate-4 datasets (20151120).
sample_id_name <- "RNAseq"
sample_id_series <- "20151120"
sample_id_series_set1 <- paste(sample_id_series, "Set1", sep="-")
sample_id_set1 <- paste(sample_id_name, sample_id_series_set1, sep="_")
read_counts_file_set1 <- paste(sample_id_set1, "unq.refseq.umi", "dat", sep=".")
read_counts_file_set1 <- file.path(counts_dir, read_counts_file_set1)
exprt_design_file_set1 <- paste(paste(sample_id_name,"Experiment-Design",sep="-"), sample_id_series_set1, "txt", sep=".")
exprt_design_file_set1 <- file.path(counts_dir, exprt_design_file_set1)
sample_id_series_set2_k107 <- paste(sample_id_series, "Set2-K107", sep="-")
sample_id_set2_k107 <- paste(sample_id_name, sample_id_series_set2_k107, sep="_")
read_counts_file_set2_k107 <- paste(sample_id_set2_k107, "unq.refseq.umi", "dat", sep=".")
read_counts_file_set2_k107 <- file.path(counts_dir, read_counts_file_set2_k107)
exprt_design_file_set2_k107 <- paste(paste(sample_id_name,"Experiment-Design",sep="-"), sample_id_series_set2_k107, "txt", sep=".")
exprt_design_file_set2_k107 <- file.path(counts_dir, exprt_design_file_set2_k107)
sample_id_series_set2_k108 <- paste(sample_id_series, "Set2-K108", sep="-")
sample_id_set2_k108 <- paste(sample_id_name, sample_id_series_set2_k108, sep="_")
read_counts_file_set2_k108 <- paste(sample_id_set2_k108, "unq.refseq.umi", "dat", sep=".")
read_counts_file_set2_k108 <- file.path(counts_dir, read_counts_file_set2_k108)
exprt_design_file_set2_k108 <- paste(paste(sample_id_name,"Experiment-Design",sep="-"), sample_id_series_set2_k108, "txt", sep=".")
exprt_design_file_set2_k108 <- file.path(counts_dir, exprt_design_file_set2_k108)
if(file.exists(read_counts_file_set1) && file.exists(exprt_design_file_set1) && file.exists(read_counts_file_set2_k107) && file.exists(exprt_design_file_set2_k107) && file.exists(read_counts_file_set2_k108) && file.exists(exprt_design_file_set2_k108))
{
    read_counts_name <- "Read-Counts"
    exprt_design_name <- "Experiment-Design"
    sample_type <- "LINCS-PC"
    cat("Extracting", paste(sample_id_name,sample_id_series,sample_type,sep="-"), "samples...\n")

    # Extract dataset 1.
    read_counts_file_set1_new <- convertReadCountsFileName(read_counts_file_set1, read_counts_name)
    file.copy(read_counts_file_set1, read_counts_file_set1_new)
    read_counts_datasets_20151120_set1_lincs <- extractGeneExpSamples(read_counts_file=read_counts_file_set1_new, mol_name_col=1, sample_name_filter="T384s1_", exprt_design_file=exprt_design_file_set1, sample_names=c("A1","A24","P1","P24","F9",paste("P",9:23,sep="")), sample_keep=FALSE, write_file=write_file, extra_name=sample_type, ext_name=ext_name, swap_series_type=TRUE, func_dir=func_dir)
    # Extract dataset 2-K107.
    read_counts_file_set2_k107_new <- convertReadCountsFileName(read_counts_file_set2_k107, read_counts_name)
    file.copy(read_counts_file_set2_k107, read_counts_file_set2_k107_new)
    read_counts_datasets_20151120_set2_k107_lincs <- extractGeneExpSamples(read_counts_file=read_counts_file_set2_k107_new, mol_name_col=1, sample_name_filter="T384s1_", exprt_design_file=exprt_design_file_set2_k107, sample_names=c("A1","A24","P1","P24","F9",paste("P",9:23,sep="")), sample_keep=FALSE, write_file=write_file, extra_name=sample_type, ext_name=ext_name, swap_series_type=TRUE, func_dir=func_dir)
    # Extract dataset 2 K108.
    read_counts_file_set2_k108_new <- convertReadCountsFileName(read_counts_file_set2_k108, read_counts_name)
    file.copy(read_counts_file_set2_k108, read_counts_file_set2_k108_new)
    read_counts_datasets_20151120_set2_k108_lincs <- extractGeneExpSamples(read_counts_file=read_counts_file_set2_k108_new, mol_name_col=1, sample_name_filter="T384s1_", exprt_design_file=exprt_design_file_set2_k108, sample_names=c("A1","A24","P1","P24","F9",paste("P",9:23,sep="")), sample_keep=FALSE, write_file=write_file, extra_name=sample_type, ext_name=ext_name, swap_series_type=TRUE, func_dir=func_dir)

    # Merge dataset 1, 2-K107 and 2-K108 to a single dataset.
    read_counts_datasets_20151120_lincs <- NULL
    read_counts_datasets_20151120_lincs$exprt_design <- read_counts_datasets_20151120_set1_lincs$exprt_design
    read_counts_datasets_20151120_lincs$read_counts <- read_counts_datasets_20151120_set1_lincs$read_counts + read_counts_datasets_20151120_set2_k107_lincs$read_counts + read_counts_datasets_20151120_set2_k108_lincs$read_counts
    read_counts_file <- paste(paste(sample_id_name,read_counts_name,sep="-"), sample_type, sample_id_series, ext_name, sep=".")
    read_counts_file <- file.path(counts_dir, read_counts_file)
    exprt_design_file <- paste(paste(sample_id_name,exprt_design_name,sep="-"), sample_type, sample_id_series, ext_name, sep=".")
    exprt_design_file <- file.path(counts_dir, exprt_design_file)
    write.table(read_counts_datasets_20151120_lincs$read_counts, file=read_counts_file, quote=FALSE, sep="\t", na="")
    write.table(read_counts_datasets_20151120_lincs$exprt_design, file=exprt_design_file, quote=FALSE, sep="\t", na="", row.names=FALSE)

    cat("Done\n")
}

# Conventional dataset (20151223).
sample_id_name <- "RNAseq"
sample_id_series <- "20151223"
sample_id <- paste(sample_id_name, sample_id_series, sep="_")
read_counts_file <- paste(sample_id, "ReadCounts", "tsv", sep=".")
read_counts_file <- file.path(counts_dir, read_counts_file)
exprt_design_file <- paste(paste(sample_id_name,"Experiment-Design",sep="-"), sample_id_series, "txt", sep=".")
exprt_design_file <- file.path(counts_dir, exprt_design_file)
if(file.exists(read_counts_file) && file.exists(exprt_design_file))
{
    sample_type <- "LINCS"
    cat("Extracting", paste(sample_id_name,sample_id_series,sample_type,sep="-"), "samples...\n")
    read_counts_name <- "Read-Counts"
    read_counts_file_new <- convertReadCountsFileName(read_counts_file, read_counts_name)
    file.copy(read_counts_file, read_counts_file_new)
    read_counts_datasets_20151223_lincs <- extractGeneExpSamples(read_counts_file=read_counts_file_new, mol_name_col=1, exprt_design_file=exprt_design_file, write_file=write_file, extra_name=sample_type, ext_name=ext_name, swap_series_type=TRUE, func_dir=func_dir)
    cat("Done\n")
}
