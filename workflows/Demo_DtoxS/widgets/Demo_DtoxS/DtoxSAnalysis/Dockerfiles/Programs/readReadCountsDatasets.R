# This script fulfills the following tasks:
# Read in multiple gene/protein read counts datasets and their experiment
# design schemes, but can generate new schemes from read counts datasets
# if no scheme is given.
#
# Customizable function generateReadCountsColumnNames for pre-processing the
# column names of read-counts data before generating experiment design table,
# such that the pre-processed column names consists of a set of standard
# fields (specified by field_names) concatenated by split_char.
#
# For example, output column names such as Gel_3.Experiment_39.PAZ consists of
# 3 standard fields, namely Measure (Gel_3), Culture (Experiment_39), and State
# (PAZ), which are concatenated by the dot character (.) specified by split_char.

readReadCountsDatasets <- function(read_counts_files, generateReadCountsColumnNames=NULL, split_char=NULL, field_names=NULL, field_types=c(Species="character", Subject="character", State="character", Culture="character", Replicate="character", Measure="character", Sample="character", Time="character"), species_default=species_default, subject_default=subject_default, state_default=state_default, culture_default=culture_default, measure_default=measure_default, time_default=time_default, row_name_column=NULL, read_counts_columns=NULL, read_counts_skip=0, exprt_design_files=NULL, field_names_positions=NULL, exprt_design_skip=0, func_dir=NULL)
{
    # generateReadCountsColumnNames function is to transformed the names of
    # read counts columns to the form of multiple field names specified
    # by field_names and separated by split_char.
    # For example, a column name "Gel_1_Experiment_30_REG" will be converted
    # to "1.30.REG" that contains Measure, Culture, and State, as specified
    # by field_names and separated by a split_char of ".".

    # row_name_column:
    # a character string or an integer number for a column.

    # read_counts_columns:
    # 1) a vector of character strings or integer numbers for multiple columns.
    # 2) a single character string as a regular expression pattern.
    # 3) NULL for all columns.

    # Load required library.
    library(tools)

    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "generateExperimentDesignTable.R"), local=TRUE)
    source(file.path(func_dir, "convertExperimentDesignTable.R"), local=TRUE)

    # Verify the validity of generateReadCountsColumnNames function.
    if(is.null(exprt_design_files))
    {
        stopifnot(!is.null(generateReadCountsColumnNames))
        generateReadCountsColumnNames <- match.fun(generateReadCountsColumnNames)
        stopifnot(is.function(generateReadCountsColumnNames))
    }

    # Verify data validity when exprt_design_files is provided.
    if(!is.null(exprt_design_files))
    {
        # Verify length equality of read_counts_files and exprt_design_files.
        stopifnot(length(read_counts_files) == length(exprt_design_files))
        # Verify the availability of field_names_positions.
        stopifnot(!is.null(field_names_positions))
    }

    # Verify the validity of read_counts_columns.
    if(!is.null(read_counts_columns))
    {
        stopifnot(is.vector(read_counts_columns))
        stopifnot(is.character(read_counts_columns)||is.numeric(read_counts_columns))
        stopifnot(length(read_counts_columns)>0)
    }

    # Read in and merge multiple datasets.
    read_counts_list <- NULL
    exprt_design_list <- NULL
    for(idx in 1:length(read_counts_files))
    {
        read_counts_file <- read_counts_files[idx]
        read_counts_file_name <- basename(read_counts_file)

        # Name conventions for read counts data files and experiment design files:
        #
        # [Name].[Type].[Series].[Suffix]
        #
        # [Name]: e.g. Read-Counts, Experiment-Design, RNA-Seq, Protein-Seq.
        # [Type]: e.g. LINCS, PromoCell, IPSC.
        # [Series]: e.g. series/version number, date.
        # [Suffix]: e.g. tsv, R.tsv.

        # Current dataset name.
        dataset_name <- gsub("-", "_", file_path_sans_ext(read_counts_file_name))
        # Read in read-counts data file.
        if(!is.null(row_name_column)) read_counts_dataset <- read.delim(read_counts_file, quote="", row.names=row_name_column, comment.char="#", stringsAsFactors=FALSE, skip=read_counts_skip)
        else read_counts_dataset <- read.delim(read_counts_file, quote="", comment.char="#", stringsAsFactors=FALSE, skip=read_counts_skip)
        dataset_col_names <- colnames(read_counts_dataset)
        # Generate read counts matrix.
        if(is.vector(read_counts_columns))
        {
            if(is.character(read_counts_columns) && length(read_counts_columns)==1)
            {
                # read_counts_columns is a regular expression pattern.
                read_counts_cols_matched <- grep(read_counts_columns, dataset_col_names, value=TRUE)
                read_counts <- as.matrix(read_counts_dataset[,read_counts_cols_matched])
            }
            else read_counts <- as.matrix(read_counts_dataset[,read_counts_columns])
        }
        else read_counts <- as.matrix(read_counts_dataset)
        # Generate a new experiment design table or convert an existing one.
        if(is.null(exprt_design_files))
        {
            # Generate experiment design table if it doesn't exist.
            read_counts_cols_prep <- generateReadCountsColumnNames(colnames(read_counts))
            exprt_design <- generateExperimentDesignTable(read_counts_column_names=read_counts_cols_prep, split_char=split_char, field_names=field_names, field_types=field_types, species_default=species_default, subject_default=subject_default, state_default=state_default, culture_default=culture_default, replicate_default=replicate_default, measure_default=measure_default, sample_default=1:length(read_counts_cols_prep), time_default=time_default)
        }
        else
        {
            # Convert an existing experiment design table to a standard one.
            exprt_design_file <- exprt_design_files[idx]
            # Read in experiment-design data file.
            exprt_design_dataset <- read.delim(exprt_design_file, quote="", comment.char="#", stringsAsFactors=FALSE, skip=exprt_design_skip)
            exprt_design <- convertExperimentDesignTable(exprt_design=exprt_design_dataset, field_names_positions=field_names_positions, field_types=field_types, species_default=species_default, subject_default=subject_default, state_default=state_default, culture_default=culture_default, replicate_default=replicate_default, measure_default=measure_default, sample_default=1:nrow(exprt_design_dataset), time_default=time_default)
        }
        # Set column names to a combination of multiple fields:
        #
        # Standardized nomenclature:
        # [State].[Subject].[Culture].[Replicate].[Measure].[Sample]
        #
        # Conventional nomenclature:
        # [State].[Cell].[Experiment].[Dish].[Plate].[Well]
        #
        # Therefore each column has a unique name across all datasets.
        #
        # Note: all these 6 fields must contain non-empty values for all rows.
        # Check the N/A status of all 6 fields.
        stopifnot(all(!is.na(exprt_design$State)) && all(!is.na(exprt_design$Subject)) && all(!is.na(exprt_design$Culture)) && all(!is.na(exprt_design$Replicate)) && all(!is.na(exprt_design$Measure)) && all(!is.na(exprt_design$Sample)))
        # Assemble the ID field from all 6 fields.
        id_name <- paste(exprt_design$State, exprt_design$Subject, exprt_design$Culture, exprt_design$Replicate, exprt_design$Measure, exprt_design$Sample, sep=".")
        colnames(read_counts) <- id_name
        # Also add a new field ID to experiment design table for indexing
        # the columns of read counts matrix.
        exprt_design <- data.frame(exprt_design, ID=id_name, stringsAsFactors=FALSE)
        # Save generated read-counts matrix and experiment design table.
        read_counts_list[[dataset_name]] <- read_counts
        exprt_design_list[[dataset_name]] <- exprt_design
    }

    # Return read counts list and experiment design list.
    read_counts_datasets <- list(read_counts_list=read_counts_list, exprt_design_list=exprt_design_list)
    return(read_counts_datasets)
}
