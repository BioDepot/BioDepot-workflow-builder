FROM biodepot/bioconductor:3.7__ubuntu-18.04__R-3.5.1__081318
RUN apt-get update && apt-get -y install firefox \ 
    build-essential python3-all python3-pip libncurses5-dev libncursesw5-dev libzmq3-dev \
    && pip3 install --upgrade pip \
    && pip install jupyter \
    && apt-get -y remove build-essential \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages(c('IRdisplay', 'repr', 'devtools', 'evaluate', 'crayon','pbdZMQ', 'uuid', 'digest' ),repos = 'http://cran.us.r-project.org'); devtools::install_github('IRkernel/IRkernel'); IRkernel::installspec()"
#set mozilla preferences to not launch their homepage
COPY .mozilla /root/.mozilla
