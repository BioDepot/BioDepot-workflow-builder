# Merge multiple read-counts datasets and experiment design schemes.

mergeReadCountsDatasets <- function(read_counts_datasets, merge_method="intersect", missing_value=NULL, func_dir=NULL)
{
    # All experiment design tables have the same set of standard field
    # columns, such that they can be concatenate by rows directly.

    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "cbindTable.R"), local=TRUE)

    # Check merge_method.
    if(merge_method!="intersect" && merge_method!="union") stop("merge_method must be either \"intersect\" or \"union\"!")

    # Merge all experiment design datasets.
    exprt_design_list <- read_counts_datasets$exprt_design_list
    exprt_design_merged <- NULL
    for(dataset_name in names(exprt_design_list))
    {
        exprt_design <- exprt_design_list[[dataset_name]]
        # Merge experiment design table.
        exprt_design_merged <- rbind(exprt_design_merged, exprt_design)
    }

    # Merge all read counts datasets.
    # Multiple read counts matrices are concatenated by columns over a
    # common set of rows across all read counts matrices.

    # Find intersect rows across all read counts matrices.
    read_counts_list <- read_counts_datasets$read_counts_list
    row_names_merged <- NULL
    for(dataset_name in names(read_counts_list))
    {
        read_counts <- read_counts_list[[dataset_name]]
        row_names <- rownames(read_counts)
        if(is.null(row_names_merged)) row_names_merged <- row_names
        else
        {
            if(merge_method == "intersect") row_names_merged <- intersect(row_names_merged, row_names)
            else row_names_merged <- union(row_names_merged, row_names)
        }
    }
    # Merge all read counts matrices for intersect rows.
    if(merge_method == "intersect")
    {
        # The initial read counts matrix for intersect merging is set to
        # an empty variable because the read counts matrix to be merged
        # at each column merging step has the same number of rows.
        read_counts_merged <- NULL
    }
    else
    {
        # The initial read counts matrix for union merging is set to an
        # empty matrix with a complete set of rows from all samples, so
        # that the read counts matrix at each column merging step can be
        # integrated with already merged read counts matrix.
        read_counts_merged <- matrix(nrow=length(row_names_merged), ncol=0, dimnames=list(row_names_merged))
    }
    for(dataset_name in names(read_counts_list))
    {
        read_counts <- read_counts_list[[dataset_name]]
        if(merge_method == "intersect")
        {
            # Merge fixed common rows of read counts data across all samples.
            read_counts_merged <- cbind(read_counts_merged, read_counts[row_names_merged,])
        }
        else
        {
            # Merge all rows of read counts data across all samples.
            read_counts_merged <- cbindTable(read_counts_merged, read_counts)
        }
    }

    # Replace missing values of reaad counts matrix.
    if(!is.null(missing_value) && !is.na(missing_value)) read_counts_merged[is.na(read_counts_merged)] <- missing_value
    # Sort merged read counts matrix according to its row names.
    read_counts_merged <- read_counts_merged[order(rownames(read_counts_merged)),]

    # Return merged read counts matrix and experiment design table.
    read_counts_datasets_merged <- list(read_counts=read_counts_merged, exprt_design=exprt_design_merged)
    return(read_counts_datasets_merged)
}
