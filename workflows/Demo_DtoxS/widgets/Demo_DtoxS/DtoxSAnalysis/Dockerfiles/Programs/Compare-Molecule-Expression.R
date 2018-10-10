#!/usr/bin/env Rscript
# This script run comparison of differential expressed molecules for
# genes and proteins.

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
        if(length(cmd_args) < 2)
        {
            cat(paste("Usage:", prog_cmd, "[Configuration File] [Repository Directory] [Function Directory]\n"))
            cat("\n")
            cat("Arguments:\n\n")
            cat("\t[Configuration File] is a data file with customizable configuration parameters.\n\n")
            cat("\t[Repository Directory] is a directory with \"Counts\", \"Params\", \"Scripts\", and \"Results\" sub-directories.\n\n")
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

# Obtain the path of R function scripts or configuration files.
getConfigFile <- function(invoke_method="Rscript", default_file=NULL)
{
    # invoke_method: "R" or "Rscript".

    config_file <- NULL

    # Determine the function/config path.
    if(invoke_method == "Rscript")
    {
        # Type 2: for the invocation by Rscript at system terminal.
        cmd_args <- commandArgs(TRUE)
        if(length(cmd_args) > 0) config_file <- cmd_args[1]
        else
        {
            warning("No configuration file is provided!")
            config_file <- default_file
        }
        if(is.null(config_file)) warning("config_file is NULL when invoke_method is \"Rscript\"!")
    }
    else if(invoke_method == "R")
    {
        # Type 1: for the invocation by "R" at R terminal.
        config_file <- default_file
        if(is.null(config_file)) warning("config_file is NULL when invoke_method is \"R\"!")
    }
    else warning("invoke_method must be either \"R\" or \"Rscript\"!")

    # Return the function/config path.
    return(config_file)
}

# Obtain the path of data repository.
getRepoDir <- function(invoke_method="Rscript", default_dir=getwd())
{
    # invoke_method: "Rscript" or "R".

    repo_dir <- NULL

    # Determine the function/config path.
    if(invoke_method == "Rscript")
    {
        # For the invocation by Rscript at system terminal.
        cmd_args <- commandArgs(TRUE)
        if(length(cmd_args) > 1) repo_dir <- cmd_args[2]
        else
        {
            warning("No repository directory is provided!")
            repo_dir <- default_dir
        }
        if(is.null(repo_dir)) warning("repo_dir is NULL when invoke_method is \"Rscript\"!")
    }
    else if(invoke_method == "R")
    {
        # For the invocation by R command or "source" at R terminal.
        repo_dir <- default_dir
        if(is.null(repo_dir)) warning("repo_dir is NULL when invoke_method is \"R\"!")
    }
    else warning("invoke_method must be either \"R\" or \"Rscript\"!")

    # Return the function/config path.
    return(repo_dir)
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
        if(length(cmd_args) > 2) func_dir <- cmd_args[3]
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
setDefaultDEMDirs <- function(dataset_name="default")
{
    # Initialize default directories.
    user_home <- Sys.getenv("HOME")
    lincs_program_dir <- file.path(user_home, "LINCSData/Programs/Library")
    lincs_dataset_dir <- file.path(user_home, "LINCSData/Datasets")
    lincs_diff_dataset_dir <- file.path(lincs_dataset_dir, "Difference/LINCS.Dataset")
    default_func_dir <- lincs_program_dir
    # RNAseq Plate 1
    if(dataset_name == "LINCS.Dataset.Gene.LINCS.20150409")
    {
        # RNA sequencing 20150409
        default_repo_dir <- file.path(lincs_diff_dataset_dir, dataset_name)
        default_config_dir <- file.path(default_repo_dir, "Configs")
        default_config_file <- paste("Configs", dataset_name, "tsv", sep=".")
        default_config_file <- file.path(default_config_dir, default_config_file)
    }
    # RNAseq Plate 2
    else if(dataset_name == "LINCS.Dataset.Gene.LINCS.20150503")
    {
        # RNA sequencing 20150503
        default_repo_dir <- file.path(lincs_diff_dataset_dir, dataset_name)
        default_config_dir <- file.path(default_repo_dir, "Configs")
        default_config_file <- paste("Configs", dataset_name, "tsv", sep=".")
        default_config_file <- file.path(default_config_dir, default_config_file)
    }
    # RNAseq Plate 3
    else if(dataset_name == "LINCS.Dataset.Gene.LINCS.20150712")
    {
        # RNA sequencing 20150712
        default_repo_dir <- file.path(lincs_diff_dataset_dir, dataset_name)
        default_config_dir <- file.path(default_repo_dir, "Configs")
        default_config_file <- paste("Configs", dataset_name, "tsv", sep=".")
        default_config_file <- file.path(default_config_dir, default_config_file)
    }
    # RNAseq Plate 4
    else if(dataset_name == "LINCS.Dataset.Gene.LINCS-PC.20151120")
    {
        # RNA sequencing 20151120
        default_repo_dir <- file.path(lincs_diff_dataset_dir, dataset_name)
        default_config_dir <- file.path(default_repo_dir, "Configs")
        default_config_file <- paste("Configs", dataset_name, "tsv", sep=".")
        default_config_file <- file.path(default_config_dir, default_config_file)
    }
    # Conventional RNAseq
    else if(dataset_name == "LINCS.Dataset.Gene.LINCS.20151223")
    {
        # RNA sequencing 20151223
        default_repo_dir <- file.path(lincs_diff_dataset_dir, dataset_name)
        default_config_dir <- file.path(default_repo_dir, "Configs")
        default_config_file <- paste("Configs", dataset_name, "tsv", sep=".")
        default_config_file <- file.path(default_config_dir, default_config_file)
    }
    # PRTseq Plate 1
    else if(dataset_name == "LINCS.Dataset.Protein.LINCS.20151015")
    {
        # Protein spectrum 20151015
        default_repo_dir <- file.path(lincs_diff_dataset_dir, dataset_name)
        default_config_dir <- file.path(default_repo_dir, "Configs")
        default_config_file <- paste("Configs", dataset_name, "tsv", sep=".")
        default_config_file <- file.path(default_config_dir, default_config_file)
    }
    # PRTseq Plate 2
    else if(dataset_name == "LINCS.Dataset.Protein.LINCS.20160212")
    {
        # Protein spectrum 20160212
        default_repo_dir <- file.path(lincs_diff_dataset_dir, dataset_name)
        default_config_dir <- file.path(default_repo_dir, "Configs")
        default_config_file <- paste("Configs", dataset_name, "tsv", sep=".")
        default_config_file <- file.path(default_config_dir, default_config_file)
    }
    else if(dataset_name == "default")
    {
        default_func_dir <- file.path(user_home, "LINCS/Programs")
        default_repo_dir <- file.path(user_home, "LINCS/Dataset")
        default_config_dir <- file.path(default_repo_dir, "Configs")
        default_config_file <- "Configs.tsv"
        default_config_file <- file.path(default_config_dir, default_config_file)
    }
    else
    {
        default_func_dir <- NULL
        default_repo_dir <- NULL
        default_config_file <- NULL
    }

    # Copy all variables from current function's environment to its parent's environment.
    for(obj in ls()) assign(obj, get(obj,environment()), parent.env(environment()))
}

################### Set up basic computational environment ###################

# Initialize default directories.
# Set default directories for differential comparison analysis.
setDefaultDEMDirs("LINCS.Dataset.Protein.LINCS.20151015")

# Determine script invocation method.
invoke_method <- getInvokeMethod()
if(invoke_method == "unknown") stop("Script must be invoked by either R or Rscript!")

# Obtain invoked R script command.
stopifnot(checkInvokeCommandLine(invoke_method))

# Determine customizable configuration file.
config_file <- getConfigFile(invoke_method=invoke_method, default_file=default_config_file)
stopifnot(!is.null(config_file) && length(config_file)==1 && is.character(config_file))
if(!file.exists(config_file)) stop(paste0("Config file \"", config_file, "\" doesn't exist!"))
config_file_base <- basename(config_file)
config_file_dir <- normalizePath(dirname(config_file))
config_file <- file.path(config_file_dir, config_file_base)
# Determine directory path of data repository.
repo_dir <- getRepoDir(invoke_method=invoke_method, default_dir=default_repo_dir)
stopifnot(!is.null(repo_dir) && length(repo_dir)==1 && is.character(repo_dir))
if(!dir.exists(repo_dir)) stop(paste0("Repository directory \"", repo_dir, "\" doesn't exist!"))
repo_dir <- normalizePath(repo_dir)
# Determine directory path of all required R function scripts.
func_dir <- getFuncDir(invoke_method=invoke_method, default_dir=default_func_dir)
stopifnot(!is.null(func_dir) && length(func_dir)==1 && is.character(func_dir))
if(!dir.exists(func_dir)) stop(paste0("Function directory \"", func_dir, "\" doesn't exist!"))
func_dir <- normalizePath(func_dir)

# Load the function of differential expression comparison.
source(file.path(func_dir, "compareMoleculeExpression.R"), local=TRUE)

# Run the analysis of differential expression comparison.
read_counts_datasets <- compareMoleculeExpression(func_dir=func_dir, repo_dir=repo_dir, config_file=config_file)
