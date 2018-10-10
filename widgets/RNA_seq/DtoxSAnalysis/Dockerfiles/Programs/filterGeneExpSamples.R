# Filter outlier samples by removing outlier samples.

filterGeneExpSamples <- function(targets_merged, read_counts_merged, drug_names, dist_cutoffs, dist_cutoff_outlier=0.01, dist_cutoff_group=0.015, min_samples=3, filter_outlier=TRUE, keep_under_samples=FALSE, plot_clust=FALSE, verbose=FALSE, func_dir=NULL)
{
    # Load required library
    require(stats)

    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()

    source(file.path(func_dir, "getCutoff.R"), local=TRUE)
    source(file.path(func_dir, "clusterGeneExpSamples.R"), local=TRUE)
    source(file.path(func_dir, "selectGeneExpSamples.R"), local=TRUE)

    # Sort and filter sample replicates in each drug-treated condition for
    # all drug groups and cell lines.
    targets_merged_sorted <- NULL
    cell_lines <- levels(factor(targets_merged$Cell))
    for(drug_name in drug_names)
    {
        targets_merged_drug <- targets_merged[targets_merged$State==drug_name,]
        for(cell_line in cell_lines)
        {
            targets_merged_drug_cell <- targets_merged_drug[targets_merged_drug$Cell==cell_line,]
            plates <- levels(factor(targets_merged_drug_cell$Plate))
            for(plate in plates)
            {
                # Prepar the dataset of current condition.
                targets_merged_drug_cell_plate <- targets_merged_drug_cell[targets_merged_drug_cell$Plate==plate,]
                sample_names_drug_cell_plate <- targets_merged_drug_cell_plate$ID
                read_counts_drug_cell_plate <- read_counts_merged[,sample_names_drug_cell_plate]
                if(!is.vector(read_counts_drug_cell_plate)) read_counts_drug_cell_plate <- read_counts_drug_cell_plate[rowSums(is.na(read_counts_drug_cell_plate))==0,]
                else read_counts_drug_cell_plate <- read_counts_drug_cell_plate[!is.na(read_counts_drug_cell_plate)]

                # Calculate divisive clusters for the samples in current cell line and drug-treated condition.
                title_text <- paste(drug_name, "on Cell", cell_line, "in Plate", plate)
                divclust_results <- clusterGeneExpSamples(read_counts=read_counts_drug_cell_plate, plot_clust=plot_clust, drug_name=drug_name, cell_line=cell_line, plate_number=plate, dist_cutoffs=dist_cutoffs, dist_cutoff_outlier=dist_cutoff_outlier, dist_cutoff_group=dist_cutoff_group, min_samples=min_samples, title_text=title_text, func_dir=func_dir)

                # Remove outlier samples.
                if(filter_outlier)
                {
                    # Filter outlier sample replicates according to given cutoffs.
                    # Special handling of the CTRL in plate 3 of cell line B.
                    # Set cutoff values for outlier samples.
                    cutoff <- getCutoff(state=drug_name, cell=cell_line, plate=plate, cutoffs=dist_cutoffs, single=dist_cutoff_outlier, group=dist_cutoff_group)
                    sample_names_drug_cell_plate_sel <- selectGeneExpSamples(read_counts_drug_cell_plate, cutoff_outlier=cutoff[1], cutoff_group=cutoff[2], min_samples=min_samples, keep_under_samples=keep_under_samples, verbose=verbose, func_dir=func_dir)
                    sample_names_drug_cell_plate_idx <- sample_names_drug_cell_plate %in% sample_names_drug_cell_plate_sel
                    targets_merged_drug_cell_plate <- targets_merged_drug_cell_plate[sample_names_drug_cell_plate_idx,]
                    # Plot divisive clusters of filtered samples.
                    if(nrow(targets_merged_drug_cell_plate) > 0)
                    {
                        sample_names_drug_cell_plate <- targets_merged_drug_cell_plate$ID
                        read_counts_drug_cell_plate <- read_counts_merged[,sample_names_drug_cell_plate]
                        read_counts_drug_cell_plate <- read_counts_drug_cell_plate[rowSums(is.na(read_counts_drug_cell_plate))==0,]

                        # Calculate divisive clusters for filtered samples after outlier removal.
                        title_text <- paste(drug_name, "on Cell", cell_line, "in Plate", plate, "(Filtered)")
                        clusterGeneExpSamples(read_counts=read_counts_drug_cell_plate, plot_clust=plot_clust, drug_name=drug_name, cell_line=cell_line, plate_number=plate, dist_cutoffs=dist_cutoffs, dist_cutoff_outlier=dist_cutoff_outlier, dist_cutoff_group=dist_cutoff_group, min_samples=min_samples, ylim=divclust_results$ylim, hline=divclust_results$hline, title_text=title_text, func_dir=func_dir)
                    }
                }
                # Save filtered samples.
                if(nrow(targets_merged_drug_cell_plate)>0) targets_merged_sorted <- rbind(targets_merged_sorted, targets_merged_drug_cell_plate)
            }
        }
    }

    # Return filtered targets.
    return(targets_merged_sorted)
}
