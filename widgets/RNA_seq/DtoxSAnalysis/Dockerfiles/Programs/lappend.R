# Append new elements to an existing list (like c() for vector)

lappend <- function (lst, ...)
{
    lst <- c(lst, list(...))
    return(lst)
}
