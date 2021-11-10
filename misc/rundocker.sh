#!/bin/bash

docker cp input/countmatrix.txt de_container:data/countmatrix.txt
docker cp input/annotation.txt de_container:data/annotation.txt
docker container start de_container
docker exec -it de_container Rscript run.R data/countmatrix.txt data/annotation.txt output -a stdout -a stdin > hobotnica-debug.txt
docker cp de_container:output .
docker container stop de_container
