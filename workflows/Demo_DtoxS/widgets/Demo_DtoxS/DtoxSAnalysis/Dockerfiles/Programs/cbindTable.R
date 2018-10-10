# Combine the columns of multiple matrix or data frames in different sizes
# according to row names, meanwhile keeping all rows:
# 1) All columns are organized side by side.
# 2) All rows are kept and missing values are reppresented as NA.

# Reference: gdata package
cbindX <- function (...)
{
    x <- list(...)
    test <- sapply(x, function(z){return(is.matrix(z)|is.data.frame(z))})
    if(any(!test)) stop("Only matrices and data.frames can be used!")
    tmp <- sapply(x, nrow)
    maxi <- which.max(tmp)
    test <- tmp < tmp[maxi]
    for(i in 1:length(tmp))
    {
        if(test[i])
        {
            add <- matrix(nrow=tmp[maxi]-tmp[i], ncol=ncol(x[[i]]))
            if(is.data.frame(x[[i]])) add <- as.data.frame(add)
            colnames(add) <- colnames(x[[i]])
            x[[i]] <- rbind(x[[i]], add)
        }
    }
    ret <- x[[1]]
    for(i in 2:length(tmp)) ret <- cbind(ret, x[[i]])
    return(ret)
}

# Combine all columns of multiple matrix or data frame tables.
cbindTable <- function(..., table.name=FALSE)
{
    # tables contain multiple input matrices or data frames.
    tables <- NULL

    # Check the validity of input arguments as matrix or data.frame.
    n_args <- nargs()

    if(n_args == 0) stop("At least a list argument is required!")
    else if(n_args == 1)
    {
        # Only one argument at command line.
        tables <- list(...)
        if(length(tables) > 0)
        {
            tables <- tables[[1]]
            if(!is.list(tables))
            {
                # The only argument must be a list.
                stop("A single argument must be a list!")
            }
            else
            {
                # The list argument must have more than one element.
                if(length(tables) == 1) stop("List argument must have more than one element!")
                # The elements of the list argument must be matrix and/or data.frame.
                flags <- sapply(tables, function(x){return(is.matrix(x)|is.data.frame(x))})
                if(!all(flags)) stop("Input arguments must be matrix and/or data.frame!")
            }
        }
        else stop("At least a list argument is required!")
    }
    else
    {
        # More than one argument at command line.
        tables <- list(...)
        if(length(tables) == 1)
        {
            tables <- tables[[1]]
            # The only argument must be a list.
            if(is.list(tables))
            {
                # The list argument must have more than one element.
                if(length(tables) == 1) stop("List argument must have more than one element!")
                # The elements of the list argument must be matrix and/or data.frame.
                flags <- sapply(tables, function(x){return(is.matrix(x)|is.data.frame(x))})
                if(!all(flags)) stop("Input arguments must be matrix and/or data.frame!")
            }
            else stop("A single argument must be a list!")
        }
        else
        {
            # Multiple arguments must be matrix and/or data.frame.
            flags <- sapply(tables, function(x){return(is.matrix(x)|is.data.frame(x))})
            if(!all(flags)) stop("Input arguments must be matrix or data.frame!")
        }
    }

    # Check the availability of argument names.
    tab_name <- names(tables)
    if(!is.null(tab_name))
    {
        flags <- sapply(names(tables), function(x){return(nchar(x)==0)})
        if(any(flags))
        {
            warning("Ignore partial argument names!")
            names(tables) <- 1:length(tables)
        }
    }
    else names(tables) <- as.character(1:length(tables))

    # Check the type of input tables.
    flags <- sapply(tables, function(x){return(is.data.frame(x))})
    if(any(flags)) table_type <- "data.frame"
    else table_type <- "matrix"

    # Check the availability of the names of rows and columns of input tables.
    flags <- sapply(tables, function(x){ rlens<-nchar(rownames(x)); clens<-nchar(colnames(x)); return(any(rlens==0)|any(clens==0)) })
    if(any(flags)) stop("Input tables must have non-empty names for all rows and columns!")

    # Disable adding table name to column name if any table has no columns.
    flags <- sapply(tables, function(x){ return(ncol(x)==0) })
    if(any(flags)) table.name <- FALSE

    # Combine multiple input tables according to their row names.

    # Create an empty table that contains a complet set of row names from
    # all tables.
    row_names <- NULL
    for(tab_name in names(tables)) row_names <- c(row_names, unique(rownames(tables[[tab_name]])))
    row_names <- unique(row_names)
    #ctable <- matrix(nrow=length(row_names), ncol=0, dimnames=list(row_names))
    ctable <- data.frame(row.names=row_names, stringsAsFactors=FALSE)

    # Combine all tables into a single table by including all columns and rows:
    for(tab_name in names(tables))
    {
        # Retrieve each input table.
        tab <- tables[[tab_name]]
        # Append current table name to its column names.
        if(table.name) colnames(tab) <- paste(colnames(tab), tab_name, sep=".")

        # Update row names according to row order of current combined table.
        row_names <- rownames(ctable)

        # Get common rows between combined table and current table to be combined.
        common_row_flags <- row_names %in% rownames(tab)

        # Reorder the rows of combined table so that common rows is above different
        # rows.
        #ctable_comm <- matrix(ctable[common_row_flags,], nrow=sum(common_row_flags), ncol=ncol(ctable), dimnames=list(row_names[common_row_flags],colnames(ctable)))
        #ctable_diff <- matrix(ctable[!common_row_flags,], nrow=sum(!common_row_flags), ncol=ncol(ctable), dimnames=list(row_names[!common_row_flags],colnames(ctable)))
        ctable_comm <- ctable[common_row_flags,]
        ctable_diff <- ctable[!common_row_flags,]
        # Sort the sub table with common rows according to the row order of current table.
        ctable_comm <- ctable_comm[order(match(rownames(ctable_comm), rownames(tab))),]
        # Splice two sub tables by rows:
        # The sub table with common rows followed by that with different rows.
        if(ncol(ctable) > 0) ctable <- rbind(ctable_comm, ctable_diff)
        else ctable <- data.frame(row.names=c(rownames(ctable_comm),rownames(ctable_diff)), stringsAsFactors=FALSE)
        # Combine reordered combined ctable table with current table and
        # fill the non-common rows with NA.
        # Note: current ctable table is always combined with the upper
        # subset of already combined table with common rows, and the lower
        # subset of combined table with non-common rows can then be filled with
        # NA.
        ctable <- cbindX(ctable, tab)
    }
    # Sort combined ctable table according to its row names.
    ctable <- ctable[order(rownames(ctable)),]

    # Convert combined table to matrix type if needed.
    if(table_type == "matrix") ctable <- as.matrix(ctable)

    # Return the combined table.
    return(ctable)
}
