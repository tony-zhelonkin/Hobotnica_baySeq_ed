FROM docker.io/bioconductor/bioconductor_docker:RELEASE_3_13

COPY install.R .
COPY run.R .
RUN mkdir source
COPY source/ /source
RUN mkdir hobotnica
COPY hobotnica/ /hobotnica
RUN mkdir data

RUN Rscript install.R




