FROM biodepot/jupyter:6.5.2__r-4.2.2__bookworm-slim__8ce1efff__83930df4__2efaa883
RUN echo R_UNZIPCMD="'"internal"'"  > /root/.Renviron
RUN Rscript -e 'install.packages("BiocManager",repos="http://cran.us.r-project.org",lib="/usr/local/rhome/lib/")'
