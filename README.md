Comparing differential expression analysis tools with Hobotnica
--------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------


Version: 1.0
------------

General usage notes
-------------------
- Program computes and visualizes gene-level RNA-seq differential expression
between two experimental conditions by 5 tools: DESeq, EBSeq, edgeR, 
limma voom and NOISeq. 
- Program also uses Hobotnica to compare results and calculate the best 
differential expression tool for input data. It also visualizes intersection 
of subsets.

----------------------------------------------------------------
Input data format
-----------------
Program uses matrix of un-normalized counts, table of annotation for 
each sample and name of directory for output. Optional flag is "-s", 
which allows saving interim results.

The value in the i-th row and the j-th column of the count matrix tells 
how many reads (or fragments, for paired-end RNA-seq) can be assigned to 
gene i in sample j. Separator for values of matrix elements is ','. First 
column has no name. You can see an example in 
**'data/TCGA_prostate_countmatrix.txt'** file.

Table of annotation contains the sample name in the first column and 
the condition in the second. Separator for values of tables elements is ','.
Columns have names '0' and '1'. You can see an example in 
**'data/annotation_TCGA_prostate.txt'** file.


-----------------------------------------------------------------------------
Output data
-----------
Program gives several files for each tool as output. They all are in 
requested directory.

**'\*_plot.pdf'** files contain visualisation of differential expression analysis.
Expressed genes are signed. Files names are **'DESeq_plot.pdf'**, 
**'EBSeq_plot.pdf'**, **'edgeR_plot.pdf'**, **'NOISeq_plot.pdf'**
and **'voom_plot.pdf'**.

**'\*_sig.txt'** files contain lists of subsets for each tool corresponding 
chosen signature. Default signature is top-30 most expressed genes. Files 
names are **'DESeq_sig.txt'**, **'EBSeq_sig.txt'**, **'edgeR_sig.txt'**,
**'NOISeq_sig.txt'** and **'voom_sig.txt'**.

Program also plots Venn diagram that shows intersection of different tools 
subsets. File name is **'venn_diagram.pdf'**.

File **'hobotnica_scores.txt'** contains results of Hobotnica computing 
for each tool.

Optionally, using "-s" flag you can save this files:

- **'\*_res.rds'** files contain differential expression analysis results. They 
are different for each tool and might be of interest only for users who 
want to write their own R code with that results. Files names are 
**'DESeq_res.rds'**, **'EBSeq_res.rds'**, **'edgeR_res.rds'**, 
**'NOISeq_res.rds'** and **'voom_res.rds'**.

- **'\*_sig.txt.distmatrix'** files contain distance matrixes for each 
**'\*_sig.txt'** subset. Matrixes are ready to be used in Hobotnica. Files names are 
**'DESeq_sig.txt.distmatrix'**, **'EBSeq_sig.txt.distmatrix'**, 
**'edgeR_sig.txt.distmatrix'**, **'NOISeq_sig.txt.distmatrix'** and 
**'voom_sig.txt.distmatrix'**.



---------------------------------------------------------------------------

Quick start
------------
Download repository and change to the directory where this repository 
is stored. There are to ways to run program.
- **From your computer**

Type in command line this commands:

> **$ Rscript install.R**
> 
> **$ Rscript run.R countmatrix annotation output** 

Where **countmatrix** is matrix of un-normalized counts file name and 
**annotation** is table of annotation file name. **output** is a name
of directory for results.

- **From Docker**

Create an image **diffexprimage** and start container **de_container**
(it automatically uses Dockerfile from repository):
> **$ docker build -t diffexprimage .**
> 
> **$ docker run --name de_container diffexprimage echo "hello from docker"**
>
Copy your input data to container. **countmatrix.txt** is matrix of 
un-normalized counts file name and **annotation.txt** is table of 
annotation file name.
> 
> **$ docker cp countmatrix.txt de_container:data/countmatrix.txt**
> 
> **$ docker cp annotation.txt de_container:data/annotation.txt**
>
Start work in container.
> **$ docker start de_container**
> 
> **$ docker exec -it de_container /bin/bash**
> 
Start work with program from Docker container.
> **/# Rscript run.R data/countmatrix.txt data/annotation.txt output**
>
where **output** is a name of directory for results.

-------------------------------------

Source files description
-----------------

- **install.R** \
Installs all required R packages.\
**Example**: 
`Rscript install.R`


- **run.R** \
Contains all the pipeline. It needs three command arguments - count 
matrix file name, table annotation file name and output directory name. 
See output of run.R in output data paragraph of readme.\
Imports functions from **source/DESeq.R**, **source/EBSeq.R**, 
**source/edgeR.R**, **source/NOISeq.R**, **source/voom.R**, 
**source/signatures_utils.R**, **source/calculate_distmatrix_utils.R** \
**Example**: \
`Rscript run.R data/TCGA_prostate_countmatrix.txt 
data/annotation_TCGA_prostate.txt ouput`


- **source/pipeline_parts/de_analysis.R**\
Computes differential expression for all tools. Returns **'\*_res.rds'**
files and **'\*_plot.pdf'** files. Imports functions from **source/DESeq.R**,
**source/EBSeq.R**, **source/edgeR.R**, **source/NOISeq.R**, 
**source/voom.R**.  It needs three command arguments - count 
matrix file name, table annotation file name and output directory name. \
**Example**: 
`Rscript source/pipeline_parts/de_analysis.R data/TCGA_prostate_countmatrix.txt 
data/annotation_TCGA_prostate.txt ouput`


- **source/pipeline_parts/signatures_making.R**\
Makes **'\*_sig.txt'** files using **'\*_res.rds'** files and chosen 
signature. Imports functions from **source/signatures_utils.R**.
It needs three command arguments - count matrix file name, table 
annotation file name and output directory name. \
**Example**: \
`Rscript source/pipeline_parts/signatures_making.R data/TCGA_prostate_countmatrix.txt 
data/annotation_TCGA_prostate.txt ouput`


- **source/pipeline_parts/calculate_distmatrix.R**\
Makes **'\*_sig.txt.distmatrix'** files using **'\*_sig.txt'** files. 
Imports functions from **source/calculate_distmatrix_utils.R**
It needs three command arguments - count matrix file name, table 
annotation file name and output directory name.\
**Example**: \
`Rscript source/pipeline_parts/calculate_distmatrix.R data/TCGA_prostate_countmatrix.txt 
data/annotation_TCGA_prostate.txt ouput`

- **source/pipeline_parts/calculate_hobotnica.R**\
Makes **'hobotnica_scores.txt'** file using **'\*_sig.txt.distmatrix'** 
files. 
Imports functions from **hobotnica/R/Hobotnica.R**
It needs three command arguments - count matrix file name, table 
annotation file name and output directory name.
**Example**: 
`Rscript source/pipeline_parts/calculate_distmatrix.R data/TCGA_prostate_countmatrix.txt 
data/annotation_TCGA_prostate.txt ouput`


Other files were not made for executing and contain functions for different parts of the pipelines.

-------------------------------------------------------------------------

Current problems
----------------
- ***filtered_signature*** signature is unrepresentative. Perhaps every 
tool has their own sensitiveness for logFC, p-value etc.
- Probably need to add some other features.





