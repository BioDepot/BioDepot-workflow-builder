# Find the inverse 10-based exponent of a decimal number.
# For example:
# The exponent for 0.27 is 1.
# The exponent for 0.0027 is 3.
# The exponent for 2.7 is 0.

calcInvExp10ofDecimal <- function(x)
{
    x <- abs(x)
    if(x < 1)
    {
        n <- 0
        while(floor(x * 10^n) < 1) n <- n + 1
    }
    else n <- 0
    return(n)
}
