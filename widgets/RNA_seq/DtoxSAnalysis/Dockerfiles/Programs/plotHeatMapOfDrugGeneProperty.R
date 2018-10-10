# Plot the heat map of specified property (LogFC/p-Value) of top-ranked genes
# across multiple drugs.

dist.euclidean <- function(x)
{
    return(dist(x,method="euclidean"))
}

dist.pearson <- function(x)
{
    return(as.dist((1-cor(t(x),method="pearson"))/2))
}

add_prop_trans <- function(prop_type, prop_trans_type, read_counts_calc, transfunc)
{
    # Transform p-value.
    if(prop_type=="PValue" || prop_type=="FDR")
    {
        n_args_trans <- length(formals(transfunc))
        #stopifnot(n_args_trans==2 || n_args_trans==3)
        if(n_args_trans < 2) stop(paste("Transformation function for", prop_type, "needs 2 input arguments at least!"))
        else if(n_args_trans == 2) add_prop_trans_expr <- paste0("cbind(read_counts_calc,", prop_trans_type, "=transfunc(read_counts_calc[,prop_type]))")
        else if(n_args_trans == 3) add_prop_trans_expr <- paste0("cbind(read_counts_calc,", prop_trans_type, "=transfunc(read_counts_calc[,prop_type],read_counts_calc[,\"logFC\"]))")
        else stop(paste("Transformation function for", prop_type, "allows 3 input arguments at most!"))
    }
    # Transform other property.
    else add_prop_trans_expr <- paste0("cbind(read_counts_calc,", prop_trans_type, "=transfunc(read_counts_calc[,prop_type]))")
    return(eval(parse(text=add_prop_trans_expr)))
}

# Transformation function for p-value/FDR or other properties:
# Regular logarithmic transformation of p-value/FDR.
# log10_neg_unsigned <- function(values, low_lim=1e-3)
# {
#     values[values<low_lim] <- low_lim
#     return(-log10(values))
# }
#
# Logarithmic transformation of p-value/FDR signed with regulation direction.
# log10_neg_signed <- function(values, logfc=NULL, low_lim=1e-3)
# {
#     values[values<low_lim] <- low_lim
#     values_trans <- -log10(values)
#     if(!is.null(logfc)) values_trans <- values_trans * sign(logfc)
#     return(values_trans)
# }

plotHeatMapOfDrugGeneProperty <- function(drug_table, deg_read_counts, comn_genes=NULL, comn_drugs=NULL, prop_type="logfc", top_gene_type="pvalue", n_genes=100, pvalue_fdr=0.05, merge_type="intersect", mis_val_rep=NULL, gene_clust=TRUE, drug_clust=FALSE, dendrogram=TRUE, dist_type="euclidean", prop_name=NULL, prop_lab=prop_name, gene_name=NULL, drug_name=NULL, low_val=NULL, high_val=NULL, color_pal=c("blue","white","red"), plot_margin=8, state_name_ctrl="CTRL", transfunc=NULL, verbose=FALSE, func_dir=NULL)
{
    # Drug table: Group, Type, Drug, Cell, Plate, and Replicate.
    # drug_table <- drug_table_nibs
    #
    # A set of specified common drugs and genes.
    # comn_drugs <- NULL
    # comn_genes <- NULL
    #
    # 2 types of common genes:
    # merge_type <- "intersect"
    # merge_type <- "union"
    #
    # 2 types of thresholds:
    # p-value/FDR threshold
    # pvalue_fdr <- 0.1
    # Number of top-ranked genes in the sorted p-value/FDR list of genes.
    # n_genes <- 5000
    # Notes:
    # These thresholds are used only when comn_genes is NULL.
    #
    # 3 types of properties plotted in heat map:
    # p-Value
    # prop_type <- "pvalue"
    # False discovery rate (FDR)
    # prop_type <- "fdr"
    # Logrithmic fold change
    # prop_type <- "logfc"
    #
    # 4 ways of selecting top-ranked genes:
    # Top-ranked genes whose p-values are under specified p-value threshold.
    # top_gene_type <- "pvalue"
    # Top-ranked genes whose FDRs are under specified FDR threshold.
    # top_gene_type <- "fdr"
    # The first n_genes number of genes in the sorted p-value list of genes.
    # top_gene_type <- "top.pvalue"
    # The first n_genes number of genes in the sorted FDR list of genes.
    # top_gene_type <- "top.fdr"
    # (Only used when comn_genes is NULL)
    #
    # Transformation function for p-value/FDR or other properties to be plotted.
    # transfunc <- NULL
    # transfunc <- log10_neg_signed
    #
    # 2 types of distance metrics:
    # "euclidean" works for any types of numerical matrix.
    # dist_type <- "euclidean"
    # "pearson" doesn't work for sparse matrix with abundant zero values.
    # dist_type <- "pearson"
    #
    # Replacement of missing values for 2 types of distance metrics.
    # Euclidean:
    # mis_val_rep can be any non-NA value for correct clustering analysis.
    # 1) If merge_type is "intersect", then mis_val_rep can be either NULL
    # or non NULL, because the calculation process of intersect will remove
    # the data containing NAs, and therefore heatmap.2 can work properly.
    # 2) If merge_type is "union", then mis_val_rep needs to be set to a
    # non-NA value, so that heatmap.2 can generate clusters for genes and
    # drugs.
    # Pearson:
    # 1) If merge_type is "intersect", then mis_val_rep can be NULL only,
    # meaning that missing values are kept as NA, so that no clustering analysis
    # will be carried out, because inappropriate replacement for missing values
    # may cause Pearson-based distance metrics to fail. The data containing
    # missing value will be removed by the calculation process of intersect.
    # 2) If merge_type is "union", then mis_val_rep must be set to any non-NA
    # value, otherwise the subsequent Pearson calculation will fail, because
    # the union process doesn't remove the data containing NAs which then will
    # be used by Pearson calculation.
    # Notes:
    # 1) If a transformation function is specified, the replacement value
    # will be applied to the transformed prop_type values.
    # 2) Such missing-value replacment is only used for making heat map plot,
    # not for searching for a common set of genes if needed.
    # No replacement
    # mis_val_rep <- NULL
    # Replace missing values by zero
    # mis_val_rep <- 0
    #
    # Lower and upper limits of plotted values in heat map
    # low_lim_log10_neg <- -log10(1e-3)
    # low_val <- -low_lim_log10_neg
    # high_val <- low_lim_log10_neg
    #
    # Clustering options
    # Cluster genes
    # gene_clust <- TRUE
    # Cluster drugs
    # drug_clust <- TRUE
    # Draw dendrograms for clustered dimensions
    # dendrogram <- TRUE
    #
    # Color palette of heat map
    # color_pal <- c("blue","white","red")
    #
    # Plot layout
    # plot_margin <- 8
    #
    # Text labels in heat map
    # Property name in title
    # prop_name <- "Log10 p-Value"
    # prop_name <- "Log10 FDR"
    # prop_name <- "Log2 Fold Change"
    # Property lable in the histogram
    # prop_lab <- "Log10(p)"
    # prop_lab <- "Log10(FDR)"
    # prop_lab <- "Log2(FC)"
    # Name of gene type in title
    # gene_name <- NULL
    # Name of drug type in title
    # drug_name <- "NIBS"

    # Load required library
    require(matrixStats)
    require(gplots)

    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "generateColorGradient.R"), local=TRUE)

    # Check prop_type.
    if(prop_type!="logfc" && prop_type!="pvalue" && prop_type!="fdr")
    {
        warning("prop_type must be one of logfc, pvalue, or fdr!")
        return(NULL)
    }
    if(prop_type == "logfc") prop_type <- "logFC"
    else if(prop_type == "pvalue") prop_type <- "PValue"
    else if(prop_type == "fdr") prop_type <- "FDR"
    else prop_type <- NULL

    # Check top_gene_type.
    if(top_gene_type!="pvalue" && top_gene_type!="fdr" && top_gene_type!="top.pvalue" && top_gene_type!="top.fdr")
    {
        warning("top_gene_type must be one of pvalue, fdr, top.pvalue, or top.fdr!")
        return(NULL)
    }

    # Check merge_type.
    if(merge_type!="intersect" && merge_type!="union")
    {
        warning("merge_type must be either intersect or union!")
        return(NULL)
    }

    # Check dist_type
    if(dist_type!="euclidean" && dist_type!="pearson")
    {
        warning("dist_type must be either euclidean or pearson!")
        return(NULL)
    }
    # Set distance function.
    if(dist_type == "euclidean") distfun <- dist.euclidean
    else if(dist_type == "pearson") distfun <- dist.pearson
    else
    {
        warning("dist_type must be either euclidean or pearson!")
        return(NULL)
    }

    # Check transformation function.
    if(!is.null(transfunc))
    {
        transfunc <- match.fun(transfunc)
        prop_trans_type <- paste(prop_type, "trans", sep=".")
    }
    else prop_trans_type <- NULL

    # Set comn_drugs.
    drug_names <- unique(drug_table$Drug)
    if(is.null(comn_drugs)) comn_drugs <- drug_names
    comn_drugs <- gsub("\\+", "_", comn_drugs)

    # Determine a common/union set of genes from all drug-treated conditions.
    # Pairwise comparisons of all drugs included in drug_table.
    #drug_combo_idx <- t(combn(length(drug_names), 2))
    #drug_combo <- paste(drug_names[drug_combo_idx[,1]], drug_names[drug_combo_idx[,2]], sep=".")
    drug_combo <- paste(state_name_ctrl, drug_table$Drug, sep=".")
    read_counts_calc_exprs <- paste("deg_read_counts", drug_table$Group, drug_table$Cell, paste0("P",drug_table$Plate), drug_combo, sep="$")
    read_counts_calc_exprs <- gsub("\\+", "_", read_counts_calc_exprs)
    if(is.null(comn_genes))
    {
        for(read_counts_calc_expr in read_counts_calc_exprs)
        {
            state_name <- unlist(strsplit(read_counts_calc_expr, " *\\. *"))[2]
            if(state_name %in% comn_drugs)
            {
                read_counts_calc <- eval(parse(text=read_counts_calc_expr))
                # Transform property value.
                if(!is.null(prop_trans_type)) read_counts_calc <- add_prop_trans(prop_type=prop_type, prop_trans_type=prop_trans_type, read_counts_calc=read_counts_calc, transfunc=transfunc)
                # Sort genes by p-value.
                if(top_gene_type=="pvalue" || top_gene_type=="top.pvalue") calc_cond_sorted <- read_counts_calc[sort.list(read_counts_calc[,"PValue"], decreasing=FALSE), ]
                else if(top_gene_type == "fdr" || top_gene_type=="top.fdr") calc_cond_sorted <- read_counts_calc[sort.list(read_counts_calc[,"FDR"], decreasing=FALSE), ]
                if(merge_type == "intersect")
                {
                    # Intersect set of genes.
                    if(top_gene_type == "pvalue")
                    {
                        if(is.null(comn_genes)) comn_genes <- rownames(calc_cond_sorted[calc_cond_sorted[,"PValue"]<=pvalue_fdr,])
                        else comn_genes <- intersect(comn_genes, rownames(calc_cond_sorted[calc_cond_sorted[,"PValue"]<=pvalue_fdr,]))
                    }
                    if(top_gene_type == "fdr")
                    {
                        if(is.null(comn_genes)) comn_genes <- rownames(calc_cond_sorted[calc_cond_sorted[,"FDR"]<=pvalue_fdr,])
                        else comn_genes <- intersect(comn_genes, rownames(calc_cond_sorted[calc_cond_sorted[,"FDR"]<=pvalue_fdr,]))
                    }
                    else
                    {
                        if(is.null(comn_genes)) comn_genes <- rownames(calc_cond_sorted[1:n_genes,])
                        else comn_genes <- intersect(comn_genes, rownames(calc_cond_sorted[1:n_genes,]))
                    }
                    if(verbose)
                    {
                        if(top_gene_type=="pvalue" || top_gene_type=="top.pvalue") print(paste(paste0(state_name, ":"), length(comn_genes), "common genes and max p-value is", calc_cond_sorted[comn_genes[length(comn_genes)],"PValue"]))
                        else if(top_gene_type=="fdr" || top_gene_type=="top.fdr") print(paste(paste0(state_name, ":"), length(comn_genes), "common genes and max FDR is", calc_cond_sorted[comn_genes[length(comn_genes)],"FDR"]))
                    }
                }
                else if(merge_type == "union")
                {
                    # Union set of genes.
                    if(top_gene_type == "pvalue") comn_genes <- union(comn_genes, rownames(calc_cond_sorted[calc_cond_sorted[,"PValue"]<=pvalue_fdr,]))
                    else if(top_gene_type == "fdr") comn_genes <- union(comn_genes, rownames(calc_cond_sorted[calc_cond_sorted[,"FDR"]<=pvalue_fdr,]))
                    else comn_genes <- union(comn_genes, rownames(calc_cond_sorted[1:n_genes,]))
                    if(verbose) print(paste(paste0(state_name, ":"), length(comn_genes), "common genes"))
                }
                else
                {
                    warning("merge_type must be either intersect or union!")
                    break
                }
            }
        }
    }
    else
    {
        if(verbose) print(paste(length(comn_genes), "common genes are specified"))
    }

    # Obtain and plot a matrix of property values of every selected gene for each drug.
    prop_drug_gene <- NULL
    if(!is.null(comn_genes) && length(comn_genes)>1)
    {
        prop_drug_gene_orig <- NULL

        # Obtain a matrix of property values of every selected gene for each drug.
        cond_pairs <- NULL
        for(read_counts_calc_expr in read_counts_calc_exprs)
        {
            read_counts_calc <- eval(parse(text=read_counts_calc_expr))
            # Transform property value.
            if(!is.null(prop_trans_type)) read_counts_calc <- add_prop_trans(prop_type=prop_type, prop_trans_type=prop_trans_type, read_counts_calc=read_counts_calc, transfunc=transfunc)
            if(is.null(prop_trans_type)) prop_drug_gene <- cbind(prop_drug_gene, read_counts_calc[comn_genes,prop_type])
            else prop_drug_gene <- cbind(prop_drug_gene, read_counts_calc[comn_genes,prop_trans_type])
            read_counts_calc_expr_parts <- unlist(strsplit(read_counts_calc_expr, split="\\$"))
            drug_cond <- unlist(strsplit(read_counts_calc_expr_parts[5], split="\\."))[2]
            cell_line <- read_counts_calc_expr_parts[3]
            plate_number <- unlist(strsplit(read_counts_calc_expr_parts[4], split="P"))[2]
            cond_pair <- paste(drug_cond, cell_line, plate_number, sep=".")
            cond_pairs <- c(cond_pairs, cond_pair)
        }
        rownames(prop_drug_gene) <- comn_genes

        # Save the prop_drug_gene before replacing missing values and removing
        # the genes with zero stds.
        prop_drug_gene_orig <- prop_drug_gene

        # Finding a common set of genes between the intersect genes of comn_drugs and the available genes of the rest drugs.
        if(merge_type == "intersect")
        {
            # Replace the NA property values with supplied value if available.
            if(!is.null(mis_val_rep)) prop_drug_gene[is.na(prop_drug_gene)] <- mis_val_rep
            # Intersect: remove the genes whose property values contain NAs to get an intersect of non-NA values.
            na_flags <- is.na(prop_drug_gene)
            prop_drug_gene <- prop_drug_gene[rowSums(na_flags)==0,]
            comn_genes <- rownames(prop_drug_gene)
        }
        else if(merge_type == "union")
        {
            # Union: replace the NA property values with supplied value if available.
            if(!is.null(mis_val_rep)) prop_drug_gene[is.na(prop_drug_gene)] <- mis_val_rep
        }
        else
        {
            warning("merge_type must be either intersect or union!")
            return(NULL)
        }
        if(verbose) print(paste(length(comn_genes), "common genes are used for plot"))

        # Remove any gene rows and drug columns whose standard deviation is zero,
        # as required by the output of Pearson's correlation fed into the hclust
        # function of heatmap.2.
        if(dist_type == "pearson")
        {
            if(gene_clust)
            {
                prop_drug_gene_nonzero_std_rows <- rowSds(prop_drug_gene)!=0
                if(any(is.na(prop_drug_gene_nonzero_std_rows)))
                {
                    warning("Combined dataset contains N/As!")
                    return(NULL)
                }
                prop_drug_gene <- prop_drug_gene[prop_drug_gene_nonzero_std_rows,]
                comn_genes <- comn_genes[prop_drug_gene_nonzero_std_rows]
            }
            if(drug_clust)
            {
                prop_drug_gene_nonzero_std_cols <- colSds(prop_drug_gene)!=0
                if(any(is.na(prop_drug_gene_nonzero_std_cols)))
                {
                    warning("Combined dataset contains N/As!")
                    return(NULL)
                }
                prop_drug_gene <- prop_drug_gene[,prop_drug_gene_nonzero_std_cols]
                cond_pairs <- cond_pairs[prop_drug_gene_nonzero_std_cols]
            }
        }

        # Plot a heat map ofproperty values of every selected gene for all drugs.
        if(nrow(prop_drug_gene) > 1)
        {
            # Set the names of rows and columns.
            colnames(prop_drug_gene) <- cond_pairs
            rownames(prop_drug_gene) <- comn_genes

            # Set dendrogram parameters.
            row_reord <- NULL
            col_reord <- NULL
            if(!any(is.na(prop_drug_gene)) && !any(is.infinite(prop_drug_gene)) && !any(is.nan(prop_drug_gene)))
            {
                if(any(rowSds(prop_drug_gene)==0)) row_reord <- FALSE
                else row_reord <- gene_clust
                if(any(colSds(prop_drug_gene)==0)) col_reord <- FALSE
                else col_reord <- drug_clust
            }
            else
            {
                row_reord <- FALSE
                col_reord <- FALSE
            }
            if(dendrogram)
            {
                if(gene_clust && drug_clust) dendrogram <- "both"
                else if(gene_clust) dendrogram <- "row"
                else if(drug_clust) dendrogram <- "column"
                else dendrogram <- "none"
            }
            else dendrogram <- "none"

            # Set text parameters.
            prop_title <- if(!is.null(prop_name)) prop_name else {if(prop_type=="logFC") "Log Fold Change" else if(prop_type=="pvalue") "p-Value" else "FDR"}
            title_text <- paste(prop_title, "of", nrow(prop_drug_gene), if(!is.null(gene_name)) gene_name, "Genes at", ncol(prop_drug_gene), if(!is.null(drug_name)) drug_name, "Drug Conditions")
            title_text <- gsub(" {1,}", " ", title_text)
            key.xlab <- if(!is.null(prop_lab)) prop_lab else {if(prop_type=="logFC") "LogFC" else if(prop_type=="pvalue") "p-Value" else "FDR"}
            key.title <- paste(key.xlab, "Histogram")

            # Set margin parameters.
            margin_bottom <- plot_margin
            margin_right <- margin_bottom

            # Set color map parameters.
            n_colors <- 500
            color_map <- generateColorGradient(pal=color_pal, n=n_colors)

            # Plot the heat map of gene expression across multiple cell-plates.
            if(is.null(low_val)||is.null(high_val))
            {
                heatmap.2(prop_drug_gene, Rowv=row_reord, Colv=col_reord, distfun=distfun, dendrogram=dendrogram, scale="none", trace="none", col=color_map, key.title=key.title, key.xlab=key.xlab, keysize=1, margin=c(margin_bottom,margin_right), main=title_text)
            }
            else
            {
                breaks <- seq(low_val, high_val, length.out=n_colors+1)
                heatmap.2(prop_drug_gene, Rowv=row_reord, Colv=col_reord, distfun=distfun, dendrogram=dendrogram, scale="none", trace="none", col=color_map, breaks=breaks, key.title=key.title, key.xlab=key.xlab, keysize=1, margin=c(margin_bottom,margin_right), main=title_text)
            }

            # Save the original and filtered drug-gene property values.
            prop_drug_gene <- list(Original=prop_drug_gene_orig, Filtered=prop_drug_gene)
        }
        else
        {
            warning(paste("No", merge_type, "gene expression matrix can be generated for plot!"))
            # Save the original drug-gene property values.
            prop_drug_gene <- list(Original=prop_drug_gene_orig)
        }
    }
    else warning(paste("No common genes can be found for plot!"))

    # Return the drug-gene matrix of logFC values.
    invisible(prop_drug_gene)
}
