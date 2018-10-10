# Summarize the counts from DEG statistics.

sumDEGCounts <- function(deg_stats, deg_counts_file=NULL, verbose=FALSE, func_dir=NULL)
{
    # Load required library
    require(tools)

    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "sortTableRows.R"), local=TRUE)

    # Create a new field to indicate the comparison condition.
    deg_stats_table <- data.frame(deg_stats, Condition=paste(deg_stats$Condition1,deg_stats$Condition2,sep="."), stringsAsFactors=FALSE)
    # Summarize DEG counts of each comparison condition for each combination of cell and plate.
    deg_counts_sum <- NULL
    cell_lines <- levels(factor(deg_stats_table$Cell))
    for(cell_line in cell_lines)
    {
        deg_stats_table_cell <- deg_stats_table[deg_stats_table$Cell==cell_line,]
        plates <- levels(factor(deg_stats_table_cell$Plate))
        for(plate in plates)
        {
            deg_stats_table_cell_plate <- deg_stats_table_cell[deg_stats_table_cell$Plate==plate,]
            # Get the counts of DEGs at each drug-treated condition.
            drug_table <- as.data.frame(table(deg_stats_table_cell_plate$Condition), stringsAsFactors=FALSE)
            colnames(drug_table) <- c("Condition", "DEG")
            # Get the names of grouped drugs at their original order.
            cond_names <- unique(deg_stats_table_cell_plate$Condition)
            # Reorder DEG count table according to the original order of grouped drugs.
            drug_table <- sortTableRows(cond_names, drug_table, "Condition")
            deg_counts_sum <- rbind(deg_counts_sum, data.frame(Condition=drug_table$Condition, Cell=rep(cell_line, nrow(drug_table)), Plate=rep(plate, nrow(drug_table)), DEG=drug_table$DEG, stringsAsFactors=FALSE))
        }
    }
    # Remove the Condition field.
    deg_stats_table$Condition <- NULL
    # Sort DEG summary in the order of Drug, Cell, Plate, and DEG.
    deg_counts_sum <- deg_counts_sum[order(xtfrm(deg_counts_sum[,"Condition"]),xtfrm(deg_counts_sum[,"Cell"]),xtfrm(deg_counts_sum[,"Plate"]),xtfrm(deg_counts_sum[,"DEG"])),]
    # Restore the "Condition" field to original "Condition1" and "Condition2" fields.
    conds12 <- data.frame(t(sapply(deg_counts_sum$Condition,function(x){unlist(strsplit(x, " *\\. *"))})), row.names=NULL, stringsAsFactors=FALSE)
    colnames(conds12) <- c("Condition1", "Condition2")
    deg_counts_sum$Condition <- NULL
    deg_counts_sum <- data.frame(conds12, deg_counts_sum, row.names=NULL, stringsAsFactors=FALSE)

    # Print DEG summary.
    if(verbose) print("DEG Summary:")
    if(verbose) print(deg_counts_sum)

    # Prepare file names.
    if(!is.null(deg_counts_file))
    {
        deg_counts_file_base <- basename(deg_counts_file)
        deg_counts_file_dir <- dirname(deg_counts_file)
        deg_counts_file_main_name <- file_path_sans_ext(deg_counts_file_base)
        deg_counts_file_ext_name <- file_ext(deg_counts_file_base)
        # Save the DEG summary file.
        deg_counts_file_name <- paste(deg_counts_file_main_name, deg_counts_file_ext_name, sep=".")
        deg_counts_sum_file <- file.path(deg_counts_file_dir, deg_counts_file_name)
        write.table(deg_counts_sum, deg_counts_sum_file, sep="\t", quote=FALSE, row.names=FALSE)
    }

    # Return the summary
    return(deg_counts_sum)
}
