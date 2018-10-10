# Sort table rows according to given keys.

sortTableRows <- function(keys, valtab, keyidx)
{
    # Arguments
    #   keys: a list of keys to which valtab will be sorted.
    # valtab: a table of keys and values.
    # keyidx: the column index of keys in valtab.
    #
    # Return
    # A table of keys and values sorted according to given keys.
    
    # Get the logical indicator of valtab rows matched with keys.
    row_sel <- keys %in% valtab[,keyidx]
    # Get the row index of matched keys in valtab.
    key_row_idx <- match(keys,valtab[,keyidx])[row_sel]
    # Get the rwo-reordered valtab according to the keys.
    valtab_sorted <- valtab[key_row_idx,]
    # Return the reordered valtab.
    return(valtab_sorted)
}
