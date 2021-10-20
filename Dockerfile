FROM docker.io/bioconductor/bioconductor_docker:RELEASE_3_13

COPY install.R .
RUN Rscript install.R

COPY source/ .
COPY hobotnica/ .
COPY run.R .


#RUN cd hobotnica && R CMD INSTALL --no-multiarch --with-keep.source .


