FROM biodepot/bioc-r:3.6__r-3.4.3__ubuntu-16.04__818d4296__c361c682__d890a38f
RUN apt-get update && apt-get -y install libgit2-dev libboost-dev
# Install required R Packages
RUN echo 'library("BiocInstaller",lib="/usr/local/rhome/lib/")' >> /tmp/script.R && \
    echo 'install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.2.tar.gz", repos=NULL, type="source")' >> /tmp/script.R && \
    echo 'install.packages(c("bitops","RSQLite"), repos="cran.wustl.edu/")'>> /tmp/script.R && \
    echo 'install.packages("https://cran.r-project.org/src/contrib/Archive/caTools/caTools_1.17.1.1.tar.gz", repos=NULL, type="source")' >>/tmp/script.R && \
    echo 'biocLite(c("edgeR", "matrixStats", "gplots", "compare"), dependencies=TRUE)' >> /tmp/script.R 
    
RUN Rscript /tmp/script.R && \
    rm /tmp/script.R
    
ENV HOME /home/user
WORKDIR $HOME

# Add input directories/files
ADD Programs/ $HOME/Programs/
ADD Scripts/ $HOME/Scripts/
