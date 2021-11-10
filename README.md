Comparing differential expression analysis tools with Hobotnica
--------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------


Version: 1.2
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
each sample and name of directory for output.

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

**'de_plots.pdf'** file contains visualisation of differential expression analysis.
Expressed genes are signed.

**'\*_sig.txt'** files contain lists of subsets for each tool corresponding 
chosen signature. Default signature is top-30 most expressed genes. You can manually change it by set "-n" option in run.R. Files 
names are **'DESeq_sig.txt'**, **'EBSeq_sig.txt'**, **'edgeR_sig.txt'**,
**'NOISeq_sig.txt'** and **'voom_sig.txt'**.

Program also plots Venn diagram that shows intersection of different tools 
subsets. File name is **'venn_diagram.pdf'**.

File **'crossing.csv'** contains table with information which tool count this gene as expressed.
 
File **'hobotnica_scores.txt'** contains results of Hobotnica computing 
for each tool.

File **'crossing_with_best.csv'** contains results of the best tool and its intersection with other tools results.  

Using "-s" option in run.R you can save this files:

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
> **$ docker run -it --name de_container diffexprimage /bin/bash"**
>
>**# exit**
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

Copy output to your computer.
>**docker cp de_container:output computer_output**
>
where **computer_output** is a directory for results.

-------------------------------------

Source files description
-----------------

- **install.R** \
Installs all required R packages.\
**Example**: 
`Rscript install.R`


- **run.R countmatrix annotation output [-s] [-n top]** \
Contains all the pipeline. It needs three command arguments - count 
matrix file name **countmatrix**, table **annotation** file name and output directory name **output**.\
Options:\
**-s** allows to save extra information (see output)\
**-n top** allows to set **top** number of genes in signature.\ 
See output of run.R in output data paragraph of readme.\
Imports functions from **source/DESeq.R**, **source/EBSeq.R**, 
**source/edgeR.R**, **source/NOISeq.R**, **source/voom.R**, 
**source/signatures_utils.R**, **source/calculate_distmatrix_utils.R** \
**Example**: \
`Rscript run.R data/TCGA_prostate_countmatrix.txt 
data/annotation_TCGA_prostate.txt ouput`

Other files were not made for executing and contain functions for different parts of the pipelines.

-------------------------------------------------------------------------

Misc directory
--------------
Here one can find some convenience scripts:
- **rundocker.sh**
- **checkhobotnica.sh**

The first executes R scripts in isolated docker container. Requires **input/annotation.txt** and **input/countmatrix.txt** in directory from which is called. Outputs folder **output** and writes logs to **Hobotnica.log**. Most recommended way is to call from **$HOME**

The second one shows information about Hobotnica computation process in percent, according to logs in **$HOME/Hobotnica.log**



Current problems
----------------
- ***filtered_signature*** signature is unrepresentative. Perhaps every 
tool has their own sensitiveness for logFC, p-value etc.
- Probably need to add some other features.





