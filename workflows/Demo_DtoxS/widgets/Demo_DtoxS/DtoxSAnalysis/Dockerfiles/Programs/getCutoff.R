# Find outlier cutoff values for matched sample conditions.

getCutoff <- function(state, cell, plate, cutoffs, single=0.01, group=0.015)
{
    # state: sample state.
    # cell: sample cell line.
    # plate: sample plate.
    # cutoffs: a table of cutoff values.
    # single: cutoff value for removing single sample.
    # group: cutoff value for removing a group of samples.
    
    if(nrow(cutoffs) > 0)
    {
        state_idx <- cutoffs$State==state
        cell_idx <- cutoffs$Cell==cell
        plate_idx <- cutoffs$Plate==plate
        row_idx <- state_idx&cell_idx&plate_idx
        matched_row <- cutoffs[row_idx,]
        if(nrow(matched_row) > 0)
        {
            cutoff_single <- matched_row$Single
            if(is.na(cutoff_single)) cutoff_single <- single
            cutoff_group <- matched_row$Group
            if(is.na(cutoff_group)) cutoff_group <- group
        }
        else
        {
            cutoff_single <- single
            cutoff_group <- group
        }
    }
    else
    {
        cutoff_single <- single
        cutoff_group <- group
    }
    # The first is the cutoff of single sample and the second is the one of group samples.
    cutoff <- c(cutoff_single, cutoff_group)
    return(cutoff)
}
