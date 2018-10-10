# Convert the name of a read counts data file from the one generated
# by merge_and_count.py to the one required by extractGeneExpSamples.R,
# in the following way:
#
# Before: [SampleID Name]_[SampleID Series].[Remain Name].[Ext Name]
# After: [SampleID Name]-[Remain Name].[SampleID Series].[Ext Name]
#
# For example:
# Before:RNAseq_20150409.unq.refseq.umi.dat
# After: RNAseq-unq-refseq-umi.20150409.dat
#
# If remain_name is given, then [Remain Name] will be replaced by it.
#
# For example, if remain_name is "Read-Counts", then:
# Before:RNAseq_20150409.unq.refseq.umi.dat
# After: RNAseq-Read-Counts.20150409.dat

convertReadCountsFileName <- function(read_counts_file, remain_name=NULL)
{
    # Load required library.
    require(tools)

    # Extract basic elements of a file name.
    file_base <- basename(read_counts_file)
    file_dir <- dirname(read_counts_file)
    file_main_name <- file_path_sans_ext(file_base)
    file_ext_name <- file_ext(file_base)

    # Split main file name with "."
    split_parts_1 <- strsplit(file_main_name, "\\.")[[1]]
    n_split_parts_1 <- length(split_parts_1)
    sample_id <- split_parts_1[1]
    if(n_split_parts_1 > 1) main_name_remain <- paste0(split_parts_1[2:n_split_parts_1], collapse="-")
    else main_name_remain <- NULL

    # Split sample id with "_"
    split_parts_2 <- strsplit(sample_id, "_")[[1]]
    n_split_parts_2 <- length(split_parts_2)
    sample_id_name <- split_parts_2[1]
    if(n_split_parts_2 > 1) sample_id_series <- split_parts_2[2]
    else sample_id_series <- NULL

    # Reorder and splice name elements.
    if(!is.null(remain_name))
    {
        if(!is.null(main_name_remain)) main_name_new <- paste(sample_id_name, remain_name, sep="-")
        else main_name_new <- sample_id_name
    }
    else
    {
        if(!is.null(main_name_remain)) main_name_new <- paste(sample_id_name, main_name_remain, sep="-")
        else main_name_new <- sample_id_name
    }
    if(!is.null(sample_id_series)) main_name_new <- paste(main_name_new, sample_id_series, sep=".")
    file_name_new <- paste(main_name_new, file_ext_name, sep=".")
    if(file_dir != ".") file_name_new <- file.path(file_dir, file_name_new)

    # Return the new file name.
    return(file_name_new)
}
