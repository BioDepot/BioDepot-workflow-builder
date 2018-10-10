library(devtools)
install.packages(c('IRdisplay', 'repr', 'evaluate', 'crayon','pbdZMQ', 'uuid', 'digest' ),repos = 'http://cran.us.r-project.org')
devtools::install_github('IRkernel/IRkernel')
IRkernel::installspec()
