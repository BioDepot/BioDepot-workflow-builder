# Extract specified samples from read counts data file.

extractGeneExpSamples <- function(read_counts_file, mol_name_col=NULL, sample_name_filter=NULL, read_counts_skip=0, exprt_design_file, sample_names=NULL, sample_keep=TRUE, exprt_design_skip=0, write_file=FALSE, extra_name="extra", ext_name=NULL, swap_series_type=FALSE, func_dir=NULL)
{
    # Load required library
    require(tools)

    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "swapSeriesAndTypeInFileName.R"), local=TRUE)

    # extra_name must not be NULL for writing data files.
    stopifnot(!is.null(extra_name))

    # Read in read-count data file.
    if(!is.null(mol_name_col)) read_counts <- read.delim(file=read_counts_file, quote="", row.names=mol_name_col, comment.char="#", stringsAsFactors=FALSE, skip=read_counts_skip)
    else read_counts <- read.delim(file=read_counts_file, quote="", comment.char="#", stringsAsFactors=FALSE, skip=read_counts_skip)
    # Remove specified sub-string contained in sample names of read-count data.
    if(!is.null(sample_name_filter)) colnames(read_counts) <- gsub(sample_name_filter, "", colnames(read_counts))
    # Read in experiment-design data file.
    exprt_design <- read.delim(file=exprt_design_file, quote="", comment.char="#", row.names=NULL, stringsAsFactors=FALSE, skip=exprt_design_skip)

    # Select the specified samples.
    if(!is.null(sample_names))
    {
        if(sample_keep)
        {
            if(is.character(sample_names)) exprt_design <- exprt_design[(exprt_design$Well %in% sample_names),]
            else if(is.numeric(sample_names)) exprt_design <- exprt_design[sample_names,]
        }
        else
        {
            if(is.character(sample_names)) exprt_design <- exprt_design[!(exprt_design$Well %in% sample_names),]
            else if(is.numeric(sample_names)) exprt_design <- exprt_design[!(rownames(exprt_design) %in% sample_names),]
        }
    }
    read_counts <- read_counts[,exprt_design$Well]

    # Write experiment design table and read counts dataset.
    if(write_file)
    {
        # Write extracted read counts data file.
        read_counts_file_base <- basename(read_counts_file)
        read_counts_file_dir <- dirname(read_counts_file)
        read_counts_file_main_name <- file_path_sans_ext(read_counts_file_base)
        read_counts_file_ext_name <- file_ext(read_counts_file_base)
        if(!is.null(ext_name)) read_counts_file_extra <- paste(read_counts_file_main_name, extra_name, ext_name, sep=".")
        else read_counts_file_extra <- paste(read_counts_file_main_name, extra_name, read_counts_file_ext_name, sep=".")
        if(swap_series_type) read_counts_file_extra <- swapSeriesAndTypeInFileName(read_counts_file_extra)
        read_counts_file_extra <- file.path(read_counts_file_dir, read_counts_file_extra)
        write.table(read_counts, file=read_counts_file_extra, quote=FALSE, sep="\t", na="")

        # Write extracted experiment design table file.
        exprt_design_file_base <- basename(exprt_design_file)
        exprt_design_file_dir <- dirname(exprt_design_file)
        exprt_design_file_main_name <- file_path_sans_ext(exprt_design_file_base)
        exprt_design_file_ext_name <- file_ext(exprt_design_file_base)
        if(!is.null(ext_name)) exprt_design_file_extra <- paste(exprt_design_file_main_name, extra_name, ext_name, sep=".")
        else exprt_design_file_extra <- paste(exprt_design_file_main_name, extra_name, exprt_design_file_ext_name, sep=".")
        if(swap_series_type) exprt_design_file_extra <- swapSeriesAndTypeInFileName(exprt_design_file_extra)
        exprt_design_file_extra <- file.path(exprt_design_file_dir, exprt_design_file_extra)
        write.table(exprt_design, file=exprt_design_file_extra, quote=FALSE, sep="\t", na="", row.names=FALSE)
    }

    # Return read counts list and experiment design list.
    read_counts_datasets <- list(read_counts=read_counts, exprt_design=exprt_design)
    return(read_counts_datasets)
}
