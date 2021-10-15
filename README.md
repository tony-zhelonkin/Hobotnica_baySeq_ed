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
- Program also computes distance matrixes using different types of signatures,
that can be used in Hobotnica to compare results and calculate the best 
differential expression tool for input data. It also visualizes intersection 
of subsets.

----------------------------------------------------------------
Input data format
-----------------
Program uses matrix of un-normalized counts and table of annotation for 
each sample. 

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
Program gives several files for each tool as output. They all are in data 
directory.

**'\*_res.rds'** files contain differential expression analysis results. They 
are different for each tool and might be of interest only for users who 
want to write their own R code with that results. Files names are 
**'DESeq_res.rds'**, **'EBSeq_res.rds'**, **'edgeR_res.rds'**, 
**'NOISeq_res.rds'** and **'voom_res.rds'**.

**'\*_plot.pdf'** files contain visualisation of differential expression analysis.
Expressed genes are signed. Files names are **'DESeq_plot.pdf'**, 
**'EBSeq_plot.pdf'**, **'edgeR_plot.pdf'**, **'NOISeq_plot.pdf'**
and **'voom_plot.pdf'**.

**'\*_sig.txt'** files contain lists of subsets for each tool corresponding 
chosen signature. Default signature is top-30 most expressed genes. Files 
names are **'DESeq_sig.txt'**, **'EBSeq_sig.txt'**, **'edgeR_sig.txt'**,
**'NOISeq_sig.txt'** and **'voom_sig.txt'**.

**'\*_sig.txt.distmatrix'** files contain distance matrixes for each 
**'\*_sig.txt'** subset. Matrixes are ready to be used in Hobotnica. Files names are 
**'DESeq_sig.txt.distmatrix'**, **'EBSeq_sig.txt.distmatrix'**, 
**'edgeR_sig.txt.distmatrix'**, **'NOISeq_sig.txt.distmatrix'** and 
**'voom_sig.txt.distmatrix'**.

Program also plots Venn diagram that shows intersection of different tools 
subsets. File name is **'venn_diagram.pdf'**.

---------------------------------------------------------------------------

Quick start
------------
Type in command line this two commands:

>**Rscript install.R**
> 
> **Rscript run.R countmatrix annotation** 

Where **countmatrix** is matrix of un-normalized counts file name and 
**annotation** is table of annotation file name

-------------------------------------

Source files description
-----------------

- **install.R**

Installs all required R packages.

**Example**: Rscript install.R

- **run.R** 

Contains all the pipeline. It needs two command arguments - count 
matrix file name and table annotation file name. See output in output data.
Imports functions from **source/DESeq.R**, **source/EBSeq.R**, 
**source/edgeR.R**, **source/NOISeq.R**, **source/voom.R**, 
**source/signatures_utils.R**, **source/calculate_distmatrix_utils.R**

**Example**: Rscript run.R data/TCGA_prostate_countmatrix.txt 
data/annotation_TCGA_prostate.txt

- **source/de_analysis.R**

Computes differential expression for all tools. Returns **'\*_res.rds'**
files and **'\*_plot.pdf'** files. Imports functions from **source/DESeq.R**,
**source/EBSeq.R**, **source/edgeR.R**, **source/NOISeq.R**, 
**source/voom.R**

**Example**: Rscript source/de_analysis.R

- **source/signatures_making.R**

Makes **'\*_sig.txt'** files using **'\*_res.rds'** files and chosen 
signature. Imports functions from **source/signatures_utils.R**

**Example**: Rscript source/signatures_making.R

- **calculate_distmatrix.R**

Makes **'\*_sig.txt.distmatrix'** files using **'\*_sig.txt'** files. 
It needs one command arguments - count matrix file name. Imports 
functions from **source/calculate_distmatrix_utils.R**

Example: Rscript source/calculate_distmatrix.R 
data/TCGA_prostate_countmatrix.txt


- **source/DESeq.R**, **source/EBSeq.R**, **source/edgeR.R**, 
**source/NOISeq.R**, **source/voom.R**

Contain functions related to a specific tool. 

***\*\_f(counts, coldata)*** - computes differential expression, where *count*
is matrix of count and *coldata* is prepared for usage information about
conditions.

***\*\_v(\*\_res)*** - visualizes differential expression. *\*\_res* - is a 
result of differential expression analysis.

***\*\_top(results, n)*** - makes subset of top-n expressed genes. *results* 
- is a result of differential expression analysis, *n* - is a number 
of genes in top.

***\*\_filtered(results)*** - makes subset of filtered genes. *results* - is a
result of differential expression analysis. *Have a strange results, better
not to use.*

- source/signatures_utils.R

Contains functions for subset making.

***top_signature(results_deseq2, results_ebseq, results_edger, results_voom,
results_noiseq, n)*** - makes **'\*_sig.txt'** files using top-n signature. 
Arguments are results of differential expression analysis and *n* - 
number of genes in top.

***filtered_signature(results_deseq2, results_ebseq, results_edger, 
results_voom, results_noiseq)*** - makes **'\*_sig.txt'** files using 
filtering. Arguments are results of differential expression analysis. *Have
a strange results, better not to use.*

***draw_venn_diag()*** - plots Venn diagram. Uses as constant that names:
**'DESeq_sig.txt.distmatrix'**, **'EBSeq_sig.txt.distmatrix'**, 
**'edgeR_sig.txt.distmatrix'**, **'NOISeq_sig.txt.distmatrix'** and 
**'voom_sig.txt.distmatrix'**. If one of the subset does not exist, 
skip it.

- **source/calculate_distmatrix_utils.R** 

Contains functions for distance matrix making.

***calculate_distmatrix(countMatrixFile, signatureFile)*** - makes 
**'\*_sig.txt.distmatrix'** from count matrix *countMatrixFile* and 
subset *signatureFile*.

*calculate_distmatrixes (countMatrixFile)* - makes all 
**'\*_sig.txt.distmatrix'** files from from count matrix 
***countMatrixFile***. Uses as constant that names: 
**'DESeq_sig.txt.distmatrix'**, **'EBSeq_sig.txt.distmatrix'** , 
**'edgeR_sig.txt.distmatrix'**, **'NOISeq_sig.txt.distmatrix'** and 
**'voom_sig.txt.distmatrix'**.

-------------------------------------------------------------------------

Current problems
----------------
- ***draw_venn_diag()*** creates unuseful log file for every call.
- A lot of unuseful information in output.
- Informative value of **'NOISeq_plot.pdf'** and **'ebseq_plot.pdf'** is 
questionable.
- ***filtered_signature*** signature is unrepresentative. Perhaps every 
tool has their own sensitiveness for logFC, p-value etc. 
- ***noiseq_filtered*** need to be checked. It has some problems 
with adding filters on Log2FC.
- ***noiseq_top*** need to be checked. It doesn't have any log2FC checks. 
- Dockerfile need to be updated





