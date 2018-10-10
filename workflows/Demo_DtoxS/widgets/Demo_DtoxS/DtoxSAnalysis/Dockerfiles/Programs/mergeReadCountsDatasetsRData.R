# Extract and combine UMI read counts, p-values, and log fold changes from
# multiple read-counts RData files and Calc text files.

# User home directory
user_home <- Sys.getenv("HOME")

# The path of the functions.
func_dir <- file.path(user_home, "Documents/Repos/works/MSSM/LINCS/Programs/Library")
utils_dir <- file.path(func_dir, "Utils")

# Load required library.
library(tools)

# Set data diretories.
lincs_dataset_dir <- file.path(user_home, "LINCSData/Datasets")
lincs_diff_dataset_dir <- file.path(lincs_dataset_dir, "Difference/LINCS.Dataset")

# Set up input RData files.
date_tags <- c("20150409", "20150503", "20150712", "LINCS-PC.20151120")
dataset_names <- paste("LINCS.Dataset.Gene", date_tags, sep=".")
n_dataset_names <- length(dataset_names)
dataset_rdata_dirs <- file.path(dataset_names, "Results")
dataset_rdata_dirs <- file.path(lincs_diff_dataset_dir, dataset_rdata_dirs)
dataset_calc_dirs <- file.path(dataset_rdata_dirs, "FDR-0.1")
rdata_file_ext_name <- "RData"
dataset_rdata_files <- paste(dataset_names, rdata_file_ext_name, sep=".")
dataset_rdata_files <- file.path(dataset_rdata_dirs, dataset_rdata_files)

# Set up output table files.
table_file_ext_name <- "tsv"
date_tag <- "20150409-20150503-20150712-20151120"
dataset_table_dir <- file.path(lincs_diff_dataset_dir, paste("LINCS.Dataset.Gene",date_tag,sep="."))
dataset_table_sub_dirs <- file.path(dataset_table_dir, dataset_names)
# UMI read counts table files.
dataset_table_files <- paste(dataset_names, table_file_ext_name, sep=".")
dataset_table_files <- file.path(dataset_table_dir, dataset_table_files)
dataset_table_file <- paste(basename(dataset_table_dir), table_file_ext_name, sep=".")
dataset_table_file <- file.path(dataset_table_dir, dataset_table_file)
dataset_table_file_dir <- dirname(dataset_table_file)
dataset_table_file_base <- basename(dataset_table_file)
dataset_table_file_main_name <- file_path_sans_ext(dataset_table_file_base)
dataset_table_with_multi_headers_file <- paste(dataset_table_file_main_name, "with-Multi-Header", table_file_ext_name, sep=".")
dataset_table_with_multi_headers_file <- file.path(dataset_table_file_dir, dataset_table_with_multi_headers_file)
# p-value table files.
pvalue_table_file <- paste(dataset_table_file_main_name, "pvalue", table_file_ext_name, sep=".")
pvalue_table_file <- file.path(dataset_table_file_dir, pvalue_table_file)
pvalue_rank_table_file <- paste(dataset_table_file_main_name, "pvalue", "rank", table_file_ext_name, sep=".")
pvalue_rank_table_file <- file.path(dataset_table_file_dir, pvalue_rank_table_file)
pvalue_table_with_multi_headers_file <- paste(dataset_table_file_main_name, "pvalue", "with-Multi-Header", table_file_ext_name, sep=".")
pvalue_table_with_multi_headers_file <- file.path(dataset_table_file_dir, pvalue_table_with_multi_headers_file)
pvalue_rank_table_with_multi_headers_file <- paste(dataset_table_file_main_name, "pvalue", "rank", "with-Multi-Header", table_file_ext_name, sep=".")
pvalue_rank_table_with_multi_headers_file <- file.path(dataset_table_file_dir, pvalue_rank_table_with_multi_headers_file)
# LogFC table files.
logfc_table_file <- paste(dataset_table_file_main_name, "logfc", table_file_ext_name, sep=".")
logfc_table_file <- file.path(dataset_table_file_dir, logfc_table_file)
logfc_rank_table_file <- paste(dataset_table_file_main_name, "logfc", "rank", table_file_ext_name, sep=".")
logfc_rank_table_file <- file.path(dataset_table_file_dir, logfc_rank_table_file)
logfc_table_with_multi_headers_file <- paste(dataset_table_file_main_name, "logfc", "with-Multi-Header", table_file_ext_name, sep=".")
logfc_table_with_multi_headers_file <- file.path(dataset_table_file_dir, logfc_table_with_multi_headers_file)
logfc_rank_table_with_multi_headers_file <- paste(dataset_table_file_main_name, "logfc", "rank", "with-Multi-Header", table_file_ext_name, sep=".")
logfc_rank_table_with_multi_headers_file <- file.path(dataset_table_file_dir, logfc_rank_table_with_multi_headers_file)
# Comparison design files.
comparison_design_file <- paste(dataset_table_file_main_name, "Comparison-Design", table_file_ext_name, sep=".")
comparison_design_file <- file.path(dataset_table_file_dir, comparison_design_file)
# RData file including read counts, p-values (and their ranks), and logFC (and their ranks).
readcounts_pvalue_logfc_rdata_file <- paste(dataset_table_file_main_name, rdata_file_ext_name, sep=".")
readcounts_pvalue_logfc_rdata_file <- file.path(dataset_table_file_dir, readcounts_pvalue_logfc_rdata_file)

# Extract the UMI read counts from multiple dataset directories,
# and save the merged read counts of each dataset to output directory.

# Load datasets and extract UMI read counts.
read_counts_datasets_object_name <- "read_counts_datasets"

# Merge read_counts_datasets objects from multiple RData files.

mergeReadCountsDatasetsRData <- function(dataset_rdata_files, read_counts_datasets_object_name="read_counts_datasets")
{
    feature_names_all <- NULL
    original_exprt_design_list_all <- NULL
    original_read_counts_list_all <- NULL
    merged_exprt_design_all <- NULL
    merged_read_counts_all <- NULL
    filtered_exprt_design_all <- NULL
    filtered_read_counts_all <- NULL
    degstats_all <- NULL
    n_dataset_names <- length(dataset_rdata_files)
    for(idx in 1:n_dataset_names)
    {
        if(exists(read_counts_datasets_object_name)) remove(list=read_counts_datasets_object_name)
        dataset_rdata_file <- dataset_rdata_files[idx]
        load(dataset_rdata_file)

        # Original datasets

        # Combine original experiment design list.
        exprt_design_list <- read_counts_datasets$original$exprt_design_list
        for(exprt_design_name in names(exprt_design_list))
        {
            original_exprt_design_list_all[[exprt_design_name]] <- exprt_design_list[[exprt_design_name]]
        }
        # Combine original read counts list.
        read_counts_list <- read_counts_datasets$original$read_counts_list
        for(read_counts_name in names(read_counts_list))
        {
            original_read_counts_list_all[[read_counts_name]] <- read_counts_list[[read_counts_name]]
        }

        # Merged datasets

        # Combine merged experiment design.
        merged_exprt_design_all <- rbind(merged_exprt_design_all, read_counts_datasets$merged$exprt_design)
        # Combine merged read counts.
        read_counts <- read_counts_datasets$merged$read_counts
        # Obtain all feature names.
        if(idx == 1) feature_names_all <- rownames(read_counts)
        # Combine datasets.
        merged_read_counts_all <- cbind(merged_read_counts_all, read_counts)
        write.table(data.frame(Gene=rownames(read_counts), read_counts, stringsAsFactors=FALSE, check.names=FALSE), file=dataset_table_files[idx], quote=FALSE, sep="\t", na="", row.names=FALSE)

        # Filtered datasets

        # Combine filtered experiment design.
        filtered_exprt_design_all <- rbind(filtered_exprt_design_all, read_counts_datasets$filtered$exprt_design)
        # Combine filtered read counts.
        read_counts <- read_counts_datasets$filtered$read_counts
        # Create a dataset that contains all feature names and includes all missing values.
        read_counts_refill <- matrix(nrow=length(feature_names_all), ncol=ncol(read_counts), dimnames=list(feature_names_all,colnames(read_counts)))
        read_counts_refill[rownames(read_counts),] <- read_counts
        # Combine datasets.
        filtered_read_counts_all <- cbind(filtered_read_counts_all, read_counts_refill)

        # DEG statistics datasets

        # Combine DEG statistics.
        for(drug_group in names(read_counts_datasets$degstats))
        {
            for(cell_line in names(read_counts_datasets$degstats[[drug_group]]))
            {
                for(plate_num in names(read_counts_datasets$degstats[[drug_group]][[cell_line]]))
                {
                    for(drug_cond in names(read_counts_datasets$degstats[[drug_group]][[cell_line]][[plate_num]]))
                    {
                        degstats_all[[drug_group]][[cell_line]][[plate_num]][[drug_cond]] <- read_counts_datasets$degstats[[drug_group]][[cell_line]][[plate_num]][[drug_cond]]
                    }
                }
            }
        }
    }
    # Save the combined read_counts_datasets object.
    read_counts_datasets_all <- NULL
    read_counts_datasets_all$original$exprt_design_list <- original_exprt_design_list_all
    read_counts_datasets_all$original$read_counts_list <- original_read_counts_list_all
    read_counts_datasets_all$merged$exprt_design <- merged_exprt_design_all
    read_counts_datasets_all$merged$read_counts <- merged_read_counts_all
    read_counts_datasets_all$filtered$exprt_design <- filtered_exprt_design_all
    read_counts_datasets_all$filtered$read_counts <- filtered_read_counts_all
    read_counts_datasets_all$degstats <- degstats_all

    # Save the read counts table combined from multiple datasets.
    read_counts_all_dataframe <- data.frame(Gene=rownames(merged_read_counts_all), merged_read_counts_all, stringsAsFactors=FALSE, check.names=FALSE)
    write.table(read_counts_all_dataframe, file=dataset_table_file, quote=FALSE, sep="\t", na="", row.names=FALSE)

    # Retrun the merged dataset
    return(read_counts_datasets_all)
}
