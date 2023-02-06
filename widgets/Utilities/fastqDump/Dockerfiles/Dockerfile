ARG SRATOOLS_VER="3.0.3"
ARG NCBI_VDB_VER="3.0.2"

FROM debian:bookworm-slim as build
ARG SRATOOLS_VER
ARG NCBI_VDB_VER

RUN apt-get update && apt-get install -y curl make cmake git gcc g++ flex bison uuid-runtime  \
 && curl -L https://github.com/ncbi/ncbi-vdb/archive/refs/tags/${NCBI_VDB_VER}.tar.gz | tar -zvxf - \
 && curl -L https://github.com/ncbi/sra-tools/archive/refs/tags/${SRATOOLS_VER}.tar.gz | tar -zxvf - \
 && cd ncbi* && ./configure && make -j 8 && make install \
 && cd ../sra* && ./configure && make -j 8 && make install

FROM debian:bookworm-slim
ARG SRATOOLS_VER
ARG NCBI_VDB_VER

COPY --from=build /etc /etc
COPY --from=build /usr/local/ncbi/sra-tools /usr/local/ncbi/sra-tools
ENV PATH="$PATH:/usr/local/ncbi/sra-tools/bin" 

LABEL "tag.0.0"="$SRATOOLS_VER"
LABEL "tag.0.1"="$NCBI_VDB_VER"
LABEL "latest"="True"