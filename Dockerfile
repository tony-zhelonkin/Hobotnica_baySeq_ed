FROM docker.io/bioconductor/bioconductor_docker:RELEASE_3_13

COPY install.R .
RUN Rscript install.R

COPY hobotnica/ .
RUN cd hobotnica && R CMD INSTALL --no-multiarch --with-keep.source .

ENTRYPOINT ["/bin/bash"]