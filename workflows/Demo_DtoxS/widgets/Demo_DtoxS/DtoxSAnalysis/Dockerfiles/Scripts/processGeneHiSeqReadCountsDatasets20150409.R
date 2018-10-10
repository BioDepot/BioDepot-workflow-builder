# Perform specific processing of read counts datasets of gene/protein expression.

processReadCountsDatasets <- function(read_counts_datasets)
{
    # Adjust the names of promo cell lines.
    exprt_design_list <- read_counts_datasets$exprt_design_list
    read_counts_list <- read_counts_datasets$read_counts_list
    for(dataset_name in names(exprt_design_list))
    {
        exprt_design <- exprt_design_list[[dataset_name]]
        read_counts <- read_counts_list[[dataset_name]]
        exprt_design$State <- gsub("CTRL-.*", "CTRL", exprt_design$State)
        exprt_design$ID <- paste(exprt_design$State, exprt_design$Subject, exprt_design$Culture, exprt_design$Replicate, exprt_design$Measure, exprt_design$Sample, sep=".")
        colnames(read_counts) <- exprt_design$ID
        exprt_design_list[[dataset_name]] <- exprt_design
        read_counts_list[[dataset_name]] <- read_counts
    }

    # Save modified experiment design and read counts.
    read_counts_datasets$exprt_design_list <- exprt_design_list
    read_counts_datasets$read_counts_list <- read_counts_list

    # Return processed datasets.
    return(read_counts_datasets)
}
