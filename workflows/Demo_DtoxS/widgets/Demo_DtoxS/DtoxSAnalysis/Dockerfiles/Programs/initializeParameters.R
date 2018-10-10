# Initialize-Parameters.R
# Initializing various input parameters for differential expression analysis.

initializeParameters <- function(config_file, counts_dir, params_dir, scripts_dir)
{
    # Parameters for reading read-counts and experiment-design data files.

    # Read counts data files.
    read_counts_files <- NULL
    # Experiment design scheme files.
    exprt_design_files <- NULL
    # A split character is used to separate the names of selected columns in
    # a read counts data file, such that each separated part can be used as
    # a standard field name.
    split_char <- NULL
    # The standard field names corresponding to separated parts of column names.
    field_names <- NULL
    # field_names_positions includes column indices of specified experiment
    # desgin fields in given experiment design table, in a form of:
    # c(Species=1, Subject=2, State=4, Culture=7, Replicate=8, Measure=5, Sample=6, Time=3)
    field_names_positions <- NULL
    # The data types which each field of experiment design table belongs to.
    # This parameter must be specified for type conversion of data field.
    field_types <- c(Species="character", Subject="character", State="character", Culture="character", Replicate="character", Measure="character", Sample="character", Time="character")
    # A series of default field values when no experiment design is given.
    species_default <- "Human"
    subject_default <- "A"
    state_default <- "CTRL"
    culture_default <- 1
    replicate_default <- 1
    measure_default <- 1
    time_default <- 48
    # Specify which data column is used as row names of read counts matrix,
    # and it can be either an integer number or the name of a selected column.
    # Specify this parameter when read counts data file doesn't dedicate its
    # first column as row names.
    row_name_column <- NULL
    # A regular expression string to specify which columns contain raw read
    # counts data in data files. NULL means all columns in data files contain
    # read counts data.
    read_counts_columns <- NULL
    # Number of lines to skip in reading read counts files.
    read_counts_skip <- 0
    # Number of lines to skip in reading experiment design files.
    exprt_design_skip <- 0
    # Method for merging read counts data: intersect or union.
    read_counts_merge_method <- "intersect"
    # Value for replacing missing read counts with: NULL (no replacing) or any other numeric values.
    missing_value <- NULL
    # File name of processReadCountsDatasets R function.
    # Perform specific processing of read counts datasets of gene/protein expression.
    # Specify this function file if customized processing is required for specific
    # gene/protein read counts datasets.
    process_read_counts_datasets_func_file <- NULL
    processReadCountsDatasets <- NULL
    # File name of generateReadCountsColumnNames R function.
    # Specify this parameter when experiment design scheme files are not set,
    # so that scheme files needed to be generated from read counts data file
    # by parsing the column names in header line to create standard field names.
    generate_read_counts_column_names_func_file <- NULL
    generateReadCountsColumnNames <- NULL

    # Parameters for data files.

    # Drug-group definition file.
    drug_groups_file <- NULL
    # Outlier-cutoff threshold file.
    outlier_cutoffs_file <- NULL

    # Parameters for outlier removal.

    # Flag of outlier filtering
    filter_outlier <- TRUE
    # Minimum number of sample replicates for filtering outlier and
    # comparing groups.
    min_samples <- 3
    # "all": all groups must have at least min_samples replicates.
    # "any": any one of the groups must have at least min_samples replicates.
    rep_req <- "all"
    # Distance cutoff for removing outliers of a sample group with a group
    # size greater than min_samples replicates.
    dist_cutoff_outlier <- 0.01
    # Distance cutoff for discarding entire sample group of a group size
    # less than min_samples replicates. Therefore entire such sample group
    # is kept even if its size is greater than dist_cutoff_outlier.
    dist_cutoff_group <- 0.015
    # Plot cluster.
    plot_clust <- TRUE
    # Merge sample replicates across different plates.
    merge_plates <- FALSE
    merged_plate_number <- 0

    # Parameters for distance heat map.

    # Plot heat map
    plot_heatmap <- TRUE
    # Lower and upper limits for heatmap plot
    lower_limit <- 0
    upper_limit <- 0.02

    # Parameters for differential comparison.

    # Name of CTRL condition.
    state_name_ctrl <- "CTRL"
    # Method of calculating dispersion: "common", "trended", "tagwise", or "auto".
    dispersion <- "trended"
    # Top-n genes/proteins selected
    n_top <- 100
    # Generate statistic results for only DEGs or all genes/proteins.
    deg_only <- FALSE
    # Base of read counts normalization: 1e+6 for RNA data and 1e+4 for protein data.
    norm_base <- 1e+6
    # Remove genes/proteins with small read counts.
    remove_small <- TRUE
    # Plot BCV and Mean-Variance.
    plot_bcv <- TRUE
    plot_smear <- TRUE
    # Characters for plotting points: small dot (.) for plotBCV and big dot (•) for plotSmear.
    pt_chars <- c(16,19)
    # Flag of replacing "+" character in state/drug names.
    sub_char <- TRUE

    # Parameters for writing various data files.

    # Write read counts files
    write_read_counts_file <- TRUE
    # Write experiment design files.
    write_exp_configs_file <- TRUE
    # Write DEG summary and list files.
    write_deg_file <- TRUE
    # Write TOP summary and list files.
    write_top_file <- TRUE
    # Write detailed calculation result files.
    write_calc_file <- TRUE
    # Write figure plot files.
    write_plot_file <- TRUE
    # RData file for read_counts_datasets results.
    read_counts_datasets_rdata_file <- NULL

    # Parameters for screen output.

    # "verbose1" switch controls the verbose output of all functions, so it overrides "verbose".
    verbose1 <- TRUE
    # "verbose" switch controls the verbose output of compareGeneExpDiffs function.
    verbose2 <- FALSE

    # Load customized configuration parameters.
    if(!is.null(config_file) && length(config_file)==1 && is.character(config_file) && file.exists(config_file))
    {
        # Import config parameters into a local environment.
        local({
            # Read in config file.
            config_vars <- read.delim(file=config_file, quote="", comment.char="#", stringsAsFactors=FALSE)
            # Remove the parameters with zero-length values.
            config_vars <- config_vars[nchar(config_vars[,"Value"])>0,]

            for(idx in 1:nrow(config_vars)) assign(config_vars[idx,"Parameter"], config_vars[idx,"Value"])

            # Convert the variables other than character string to their appropriate data types.

            # Parameters for reading read-counts and experiment-design data files.
            if(exists("read_counts_files",inherits=FALSE)) read_counts_files <- file.path(counts_dir, strsplit(read_counts_files, split=" *, *")[[1]])
            if(exists("exprt_design_files",inherits=FALSE)) exprt_design_files <- file.path(counts_dir, strsplit(exprt_design_files, split=" *, *")[[1]])
            if(exists("field_names",inherits=FALSE)) field_names <- strsplit(field_names, split=" *, *")[[1]]
            if(exists("field_names_positions_key",inherits=FALSE) || exists("field_names_positions_value",inherits=FALSE))
            {
                stopifnot(exists("field_names_positions_key",inherits=FALSE) && exists("field_names_positions_value",inherits=FALSE))
                field_names_positions_key <- strsplit(field_names_positions_key, split=" *, *")[[1]]
                field_names_positions_value <- as.numeric(strsplit(field_names_positions_value, split=" *, *")[[1]])
                stopifnot(length(field_names_positions_key) == length(field_names_positions_value))
                field_names_positions <- field_names_positions_value
                names(field_names_positions) <- field_names_positions_key
            }
            if(exists("field_types_key",inherits=FALSE) || exists("field_types_value",inherits=FALSE))
            {
                stopifnot(exists("field_types_key",inherits=FALSE) && exists("field_types_value",inherits=FALSE))
                field_types_key <- strsplit(field_types_key, split=" *, *")[[1]]
                field_types_value <- strsplit(field_types_value, split=" *, *")[[1]]
                stopifnot(length(field_types_key) == length(field_types_value))
                field_types <- field_types_value
                names(field_types) <- field_types_key
            }
            if(exists("culture_default",inherits=FALSE)) culture_default <- as.numeric(culture_default)
            if(exists("replicate_default",inherits=FALSE)) replicate_default <- as.numeric(replicate_default)
            if(exists("measure_default",inherits=FALSE)) measure_default <- as.numeric(measure_default)
            if(exists("time_default",inherits=FALSE)) time_default <- as.numeric(time_default)
            if(exists("row_name_column",inherits=FALSE))
            {
                if(length(row_name_column) > 1) row_name_column <- row_name_column[1]
                if(suppressWarnings(!is.na(as.numeric(row_name_column)))) row_name_column <- as.numeric(row_name_column)
            }
            if(exists("read_counts_skip",inherits=FALSE)) read_counts_skip <- as.numeric(read_counts_skip)
            if(exists("exprt_design_skip",inherits=FALSE)) exprt_design_skip <- as.numeric(exprt_design_skip)
            if(exists("missing_value",inherits=FALSE)) missing_value <- as.numeric(missing_value)
            if(exists("process_read_counts_datasets_func_file",inherits=FALSE))
            {
                process_read_counts_datasets_func_file <- file.path(scripts_dir, process_read_counts_datasets_func_file)
                if(file.exists(process_read_counts_datasets_func_file)) source(process_read_counts_datasets_func_file, local=TRUE)
                else warning(paste("Customizable R function", process_read_counts_datasets_func_file, "doesn't exist!"))
            }
            if(exists("generate_read_counts_column_names_func_file",inherits=FALSE))
            {
                generate_read_counts_column_names_func_file <- file.path(scripts_dir, generate_read_counts_column_names_func_file)
                if(file.exists(generate_read_counts_column_names_func_file)) source(generate_read_counts_column_names_func_file, local=TRUE)
                else warning(paste("Customizable R function", generate_read_counts_column_names_func_file, "doesn't exist!"))
            }

            # Parameters for drug-group definition and outlier-cutoff threshold files.

            if(exists("drug_groups_file",inherits=FALSE)) drug_groups_file <- file.path(params_dir, drug_groups_file)
            if(exists("outlier_cutoffs_file",inherits=FALSE)) outlier_cutoffs_file <- file.path(params_dir, outlier_cutoffs_file)

            # Parameters for outlier removal.
            if(exists("filter_outlier",inherits=FALSE)) filter_outlier <- as.logical(filter_outlier)
            if(exists("min_samples",inherits=FALSE)) min_samples <- as.numeric(min_samples)
            if(exists("dist_cutoff_outlier",inherits=FALSE)) dist_cutoff_outlier <- as.numeric(dist_cutoff_outlier)
            if(exists("dist_cutoff_group",inherits=FALSE)) dist_cutoff_group <- as.numeric(dist_cutoff_group)
            if(exists("plot_clust",inherits=FALSE)) plot_clust <- as.logical(plot_clust)
            if(exists("merge_plates",inherits=FALSE)) merge_plates <- as.logical(merge_plates)
            if(exists("merged_plate_number",inherits=FALSE)) merged_plate_number <- as.numeric(merged_plate_number)

            # Parameters for distance heat map.
            if(exists("plot_heatmap",inherits=FALSE)) plot_heatmap <- as.logical(plot_heatmap)
            if(exists("lower_limit",inherits=FALSE)) lower_limit <- as.numeric(lower_limit)
            if(exists("upper_limit",inherits=FALSE)) upper_limit <- as.numeric(upper_limit)

            # Parameters for differential comparison.
            if(exists("n_top",inherits=FALSE)) n_top <- as.numeric(n_top)
            if(exists("deg_only",inherits=FALSE)) deg_only <- as.logical(deg_only)
            if(exists("norm_base",inherits=FALSE)) norm_base <- as.numeric(norm_base)
            if(exists("remove_small",inherits=FALSE)) remove_small <- as.logical(remove_small)
            if(exists("plot_bcv",inherits=FALSE)) plot_bcv <- as.logical(plot_bcv)
            if(exists("plot_smear",inherits=FALSE)) plot_smear <- as.logical(plot_smear)
            if(exists("pt_chars",inherits=FALSE)) pt_chars <- as.numeric(strsplit(pt_chars, split=" *, *")[[1]])
            if(exists("sub_char",inherits=FALSE)) sub_char <- as.logical(sub_char)

            # Parameters for writing various data files.
            if(exists("write_read_counts_file",inherits=FALSE)) write_read_counts_file <- as.logical(write_read_counts_file)
            if(exists("write_exp_configs_file",inherits=FALSE)) write_exp_configs_file <- as.logical(write_exp_configs_file)
            if(exists("write_deg_file",inherits=FALSE)) write_deg_file <- as.logical(write_deg_file)
            if(exists("write_top_file",inherits=FALSE)) write_top_file <- as.logical(write_top_file)
            if(exists("write_calc_file",inherits=FALSE)) write_calc_file <- as.logical(write_calc_file)
            if(exists("write_plot_file",inherits=FALSE)) write_plot_file <- as.logical(write_plot_file)

            # Parameters for screen output.
            if(exists("verbose1",inherits=FALSE)) verbose1 <- as.logical(verbose1)
            if(exists("verbose2",inherits=FALSE)) verbose2 <- as.logical(verbose2)

            # Copy all variables from current local environment to parent function's environment.
            for(obj in ls()) assign(obj, get(obj,environment()), parent.env(environment()))
        })
    }
    else warning("Specified config file doesn't exist, so built-in parameter values are used!")

    # Copy all variables from current function's environment to its parent's environment.
    for(obj in ls()) assign(obj, get(obj,environment()), parent.env(environment()))
}
