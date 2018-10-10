# Format the digits after decimal point.

format_decimal <- function(x, k=2)
{
    # x is a numerical number.
    # k is the number of digits after decimal point.
    
    return (format(round(x, k), nsmall=k))
}
