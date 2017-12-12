#!/bin/bash
cd /setup
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98431/suppl/GSE98431_Reference-Library.tar.gz
tar -zvxf GSE98431_Reference-Library.tar.gz
rm GSE98431_Reference-Library.tar.gz
mv GSE98431_Referemce-Library References