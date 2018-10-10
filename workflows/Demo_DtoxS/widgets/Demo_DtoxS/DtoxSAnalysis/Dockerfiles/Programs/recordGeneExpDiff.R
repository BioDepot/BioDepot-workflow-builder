# Record multiple data files generated from DEG comparison between two conditions.

recordGeneExpDiff <- function(group_pair, deg_read_counts, targets_state1, group_names_state1, read_counts_state1, targets_state2, group_names_state2, read_counts_state2, cell_line, plate, n_top, read_counts_file=NULL, exp_configs_file=NULL, top_stats_file=NULL, deg_stats_file=NULL, verbose=FALSE)
{
    # Load required library
    require(tools)

    # Prepare file names.
    if(!is.null(read_counts_file))
    {
        read_counts_file_base <- basename(read_counts_file)
        read_counts_file_dir <- dirname(read_counts_file)
        read_counts_file_main_name <- file_path_sans_ext(read_counts_file_base)
        read_counts_file_ext_name <- file_ext(read_counts_file_base)
    }
    if(!is.null(exp_configs_file))
    {
        exp_configs_file_base <- basename(exp_configs_file)
        exp_configs_file_dir <- dirname(exp_configs_file)
        exp_configs_file_main_name <- file_path_sans_ext(exp_configs_file_base)
        exp_configs_file_ext_name <- file_ext(exp_configs_file_base)
    }
    if(!is.null(top_stats_file))
    {
        top_stats_file_base <- basename(top_stats_file)
        top_stats_file_dir <- dirname(top_stats_file)
        top_stats_file_main_name <- file_path_sans_ext(top_stats_file_base)
        top_stats_file_ext_name <- file_ext(top_stats_file_base)
    }
    if(!is.null(deg_stats_file))
    {
        deg_stats_file_base <- basename(deg_stats_file)
        deg_stats_file_dir <- dirname(deg_stats_file)
        deg_stats_file_main_name <- file_path_sans_ext(deg_stats_file_base)
        deg_stats_file_ext_name <- file_ext(deg_stats_file_base)
    }

    ######## Data file 1: raw read counts of state 1 ########
    # Prepare the dataset of state 1.
    sample_names_state1_compare <- colnames(read_counts_state1)
    # Add common statistics of raw read counts of state 1.
    read_counts_state1_stats_compare <- cbind(Min=rowMins(read_counts_state1), Max=rowMaxs(read_counts_state1), Median=rowMedians(read_counts_state1), Mean=rowMeans(read_counts_state1), Std=rowSds(read_counts_state1))
    # Save raw read counts of state 1.
    if(!is.null(read_counts_file))
    {
        read_counts_state1_file <- NULL
        read_counts_state1_file <- paste(read_counts_file_main_name, group_pair[1], sep="-")
        read_counts_state1_file <- paste(read_counts_state1_file, read_counts_file_ext_name, sep=".")
        read_counts_state1_file <- file.path(read_counts_file_dir, read_counts_state1_file)
        write.table(cbind(Gene=rownames(read_counts_state1), read_counts_state1, read_counts_state1_stats_compare), file=read_counts_state1_file, sep="\t", quote=FALSE, row.names=FALSE)
    }

    ######## Data file 2: experimental configuraitons of state 1 ########
    # Prepare the dataset of state 1.
    # Save experimental configuraitons of state 1.
    if(!is.null(exp_configs_file))
    {
        targets_state1_file <- NULL
        targets_state1_file <- paste(exp_configs_file_main_name, group_pair[1], sep="-")
        targets_state1_file <- paste(targets_state1_file, exp_configs_file_ext_name, sep=".")
        targets_state1_file <- file.path(exp_configs_file_dir, targets_state1_file)
        write.table(targets_state1[,colnames(targets_state1)!="Use"], file=targets_state1_file, sep="\t", quote=FALSE, row.names=FALSE)
    }

    ######## Data file 3: raw read counts of state 2 ########
    # Prepare the dataset of state 2.
    sample_names_state2_compare <- colnames(read_counts_state2)
    # Add common statistics of raw read counts of state 2.
    read_counts_state2_stats_compare <- cbind(Min=rowMins(read_counts_state2), Max=rowMaxs(read_counts_state2), Median=rowMedians(read_counts_state2), Mean=rowMeans(read_counts_state2), Std=rowSds(read_counts_state2))
    # Save raw read counts of state 2.
    if(!is.null(read_counts_file))
    {
        read_counts_state2_file <- NULL
        read_counts_state2_file <- paste(read_counts_file_main_name, group_pair[2], sep="-")
        read_counts_state2_file <- paste(read_counts_state2_file, read_counts_file_ext_name, sep=".")
        read_counts_state2_file <- file.path(read_counts_file_dir, read_counts_state2_file)
        write.table(cbind(Gene=rownames(read_counts_state2), read_counts_state2, read_counts_state2_stats_compare), file=read_counts_state2_file, sep="\t", quote=FALSE, row.names=FALSE)
    }

    ######## Data file 4: experimental configuraitons of state 2 ########
    # Prepare the dataset of state 2.
    # Save experimental configuraitons of state 2.
    if(!is.null(exp_configs_file))
    {
        targets_state2_file <- NULL
        targets_state2_file <- paste(exp_configs_file_main_name, group_pair[2], sep="-")
        targets_state2_file <- paste(targets_state2_file, exp_configs_file_ext_name, sep=".")
        targets_state2_file <- file.path(exp_configs_file_dir, targets_state2_file)
        write.table(targets_state2[,colnames(targets_state2)!="Use"], file=targets_state2_file, sep="\t", quote=FALSE, row.names=FALSE)
    }

    ######## Data file 5: filtered and normalized read counts of both states ########
    # Prepare the dataset of filtered and normalized read counts.
    sample_names_state1_norm_compare <- paste(sample_names_state1_compare, "Norm", sep=".")
    sample_names_state2_norm_compare <- paste(sample_names_state2_compare, "Norm", sep=".")
    sample_names_norm_compare <- c(sample_names_state1_compare, sample_names_state2_compare, sample_names_state1_norm_compare, sample_names_state2_norm_compare)
    read_counts_norm_compare <- deg_read_counts[,sample_names_norm_compare]
    # Save raw and normalized read counts.
    if(!is.null(read_counts_file))
    {
        read_counts_norm_file <- NULL
        read_counts_norm_file <- paste(read_counts_file_main_name, "Norm", sep="-")
        read_counts_norm_file <- paste(read_counts_norm_file, paste(group_pair[1],group_pair[2],sep="."), sep="-")
        read_counts_norm_file <- paste(read_counts_norm_file, read_counts_file_ext_name, sep=".")
        read_counts_norm_file <- file.path(read_counts_file_dir, read_counts_norm_file)
        write.table(cbind(Gene=rownames(read_counts_norm_compare), read_counts_norm_compare), file=read_counts_norm_file, sep="\t", quote=FALSE, row.names=FALSE)
    }

    ######## Data file 6: Statistics of top-ranked genes ########
    # Sort genes according to their p-values.
    deg_read_counts_group_pair_sorted <- deg_read_counts[sort.list(deg_read_counts[,"PValue"], decreasing=FALSE),]
    # Extract the top-n genes from the sorted gene list.
    top_read_counts <- deg_read_counts_group_pair_sorted[1:n_top,]
    # Remove FDR and Regulation columns.
    field_names <- colnames(top_read_counts)
    field_names <- field_names[field_names!="FDR"&field_names!="Regulation"]
    top_read_counts <- top_read_counts[,field_names]
    top_names <- rownames(top_read_counts)
    # Save the DEG statistics of current plate of current cell line of current drug group.
    top_states_names <- data.frame(Cell=rep(cell_line,n_top), Plate=rep(plate,n_top), Condition1=rep(group_pair[1],n_top), Condition2=rep(group_pair[2],n_top), Gene=top_names, logFC=top_read_counts[top_names,"logFC"], logCPM=top_read_counts[top_names,"logCPM"], PValue=top_read_counts[top_names,"PValue"], stringsAsFactors=FALSE)
    # Save DEG statistics of top-ranked genes.
    if(!is.null(top_stats_file))
    {
        top_stats_file <- NULL
        top_stats_file <- paste(top_stats_file_main_name, sep="-")
        top_stats_file <- paste(top_stats_file, paste(group_pair[1],group_pair[2],sep="."), sep="-")
        top_stats_file <- paste(top_stats_file, top_stats_file_ext_name, sep=".")
        top_stats_file <- file.path(top_stats_file_dir, top_stats_file)
        write.table(top_states_names, file=top_stats_file, sep="\t", quote=FALSE, row.names=FALSE)
    }

    ######## Data file 7: Statistics of DEGs ########
    # Extract DEG statistics of DEGs.
    deg_read_counts <- deg_read_counts[deg_read_counts[,"Regulation"]!=0,]
    deg_names <- rownames(deg_read_counts)
    n_deg <- length(deg_names)
    if(verbose) print(paste(n_deg, "DEGs are found between", group_pair[1], "and", group_pair[2]), sep=" ")
    # Save the names of DEGs of current pair conditions.
    deg_states_names <- data.frame(Cell=rep(cell_line,n_deg), Plate=rep(plate,n_deg), Condition1=rep(group_pair[1],n_deg), Condition2=rep(group_pair[2],n_deg), Gene=deg_names, logFC=deg_read_counts[deg_names,"logFC"], logCPM=deg_read_counts[deg_names,"logCPM"], PValue=deg_read_counts[deg_names,"PValue"], FDR=deg_read_counts[deg_names,"FDR"], stringsAsFactors=FALSE)
    # Save DEG statistics of DEGs.
    if(!is.null(deg_stats_file))
    {
        deg_stats_file <- NULL
        deg_stats_file <- paste(deg_stats_file_main_name, sep="-")
        deg_stats_file <- paste(deg_stats_file, paste(group_pair[1],group_pair[2],sep="."), sep="-")
        deg_stats_file <- paste(deg_stats_file, deg_stats_file_ext_name, sep=".")
        deg_stats_file <- file.path(deg_stats_file_dir, deg_stats_file)
        write.table(deg_states_names, file=deg_stats_file, sep="\t", quote=FALSE, row.names=FALSE)
    }

    # Return some datasets.
    stats_dataset <- list(top_read_counts=top_read_counts, deg_states_names=deg_states_names, top_states_names=top_states_names)
    return(stats_dataset)
}
