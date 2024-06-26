Differential expression analysis of macrophage RNA sequencing data using the Hobotnica tool
--------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------


Version: 1.2.baySeq_edition
------------

This is not the original Hobotnica
-------------------
- This is not the originial Hobotnica tool repository. For the original repo 
and all the instructions on installation and use please refer to 
https://gitlab.at.ispras.ru/mirachirkova/differential-expression 
- The original paper explaining Hobotnica`s maths guts can be found here 
https://f1000research.com/articles/10-1260 

This version of Hobotnica is an edition of the original tool that was created as a 
semester project in [Bioinformatics Institute](https://bioinf.me/en), Saint-Petersburg, Russia 

Supervisors: 
- Evgeny Karpulevich

Students:
- Anton Zhelonkin [@tony-zhelonkin](https://github.com/tony-zhelonkin)
- Alexandra Belyaeva [@BelyaevaAlex](https://github.com/BelyaevaAlex)


What`s the original Hobotnica used for?
-------------------
- Original Hobotnica computes and visualizes gene-level RNA-seq differential expression 
between two experimental conditions by 5 tools: DESeq, EBSeq, edgeR, 
limma voom and NOISeq
- Hobotnica also compares results and calculates the best 
differential expression tool for input data. It also visualizes intersection 
of subsets 

----------------------------------------------------------------
Aim of the project
-----------------
Find the best tool for DE analysis of RNA-seq macrophages data

----------------------------------------------------------------
Objectives
-----------------
- Add additional DE tools into Hobotnica docker container: baySeq, DEGseq
- Compare results from different differential expression analysis tools inside Hobotnica
 
----------------------------------------------------------------

What does this version add?
-------------------
- This eddition adds bayesian differential expression tool baySeq to the 
already implemented 5 tools: DESeq, EBSeq, edgeR, limma voom and NOISeq
- Visualises the volcano-plot, heatmap based on the results of the baySeq 
- Original Hobotnica guts compare baySeq results with the other 5 differential 
expression tools, outputs all the visualisations for the rest of them and calculates 
the best tool from the results of all 6 
- All the instructions on the installation and use from the original repo 
at https://gitlab.at.ispras.ru/mirachirkova/differential-expression hold true
for this version of Hobotnica entirely

What does this version break?
-------------------
- The original Hobotnica draws a Venn diagramm of intersections between 
the 5 implemented DE tools 
- This version breaks the Venn diagramm. For the purpose of smooth run 
the code for Venn diagram was silenced 

----------------------------------------------------------------
Pipeline of the project
-----------------

![](https://github.com/tony-zhelonkin/Hobotnica_baySeq_ed/blob/main/Hobotnica_baySeq_ed.png)


Input data format for Hobotnica
-----------------
Hobotnica uses matrix of un-normalized counts, table of annotation for 
each sample and name of directory for output. An example of the input 
format may be found in the data folder. The data folder contains 
prepared un-normalised reads from RNA-seq experiment on M0 and M1 macrophages (**'RNA_macrophages_reads.txt'**) 
and annotation table (**'annotation_for_macrophages.txt'**) 
The unprepaired reads are available at 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163165 
The R script written to prepare the aforementioned reads may be found in the data folder 
under the name **'reads_preparation.Rmd'** 

-----------------------------------------------------------------------------
Output data
-----------
Program gives several files for each tool as output. They all are in 
requested directory.

**'de_plots.pdf'** and **'heatmaps.pdf'** file contains visualisation of 6 differential expression analysis.
Expressed genes are signed.

**'\*_sig.txt'** files contain lists of subsets for each tool corresponding 
chosen signature. Default signature is top-30 most expressed genes. You can manually change it by set "-n" option in run.R. Files 
names are **'DESeq_sig.txt'**, **'EBSeq_sig.txt'**, **'edgeR_sig.txt'**,
**'NOISeq_sig.txt'**, **'voom_sig.txt'** and **'baySeq_sig.txt'** (added in this version) 

File **'crossing.csv'** contains table with information which tool count this gene as expressed.
 
File **'hobotnica_scores.txt'** contains results of Hobotnica computing 
for each tool.

File **'crossing_with_best.csv'** contains results of the best tool and its intersection with other tools results.  

The Venn diagramm functionality is disabled in this version.


Using "-s" option in run.R you can save this files:

- **'\*_res.rds'** files contain differential expression analysis results. They 
are different for each tool and might be of interest only for users who 
want to write their own R code with that results. Files names are 
**'DESeq_res.rds'**, **'EBSeq_res.rds'**, **'edgeR_res.rds'**, 
**'NOISeq_res.rds'**, **'voom_res.rds'** and **'baySeq_res.rds'**.

- **'\*_sig.txt.distmatrix'** files contain distance matrixes for each 
**'\*_sig.txt'** subset. Matrixes are ready to be used in Hobotnica. Files names are 
**'DESeq_sig.txt.distmatrix'**, **'EBSeq_sig.txt.distmatrix'**, 
**'edgeR_sig.txt.distmatrix'**, **'NOISeq_sig.txt.distmatrix'**, 
**'voom_sig.txt.distmatrix'** and **'baySeq_sig.txt.distmatrix'**.


---------------------------------------------------------------------------

Hobotnica quick start guide
------------
The easiest pathway is described here 
1. Download repository 
2. Move via terminal into the directory of the downloaded repo 
3. Type in commandline the following lines to install all the necessary files and libraries
> **$ Rscript install.R**
4. To run the analysis type in the commandline
> **$ Rscript run.R countmatrix annotation output** 

Where **countmatrix** is matrix of un-normalized counts file name and 
**annotation** is table of annotation file name. **output** is a name
of directory for results. 
The **countmatrix** and **annotation** do not have to be in the same folder with the Hobotnica insides 
The analysis runs smoothly if you explicitly configure the path to the necessary files and the output folder
> **$ Rscript run.R /*your_path*/countmatrix.txt /*your_path*/annotation.txt /*your_path*/output** 

Detailed instructions on use and run from Docker are desribed at https://gitlab.at.ispras.ru/mirachirkova/differential-expression

-------------------------------------

Source files description
-----------------

- **install.R** \
Installs all required R packages, including `baySeq` and `snow`, 
which is required for parallel computations during bayesian analysis.\
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
**source/signatures_utils.R**, **source/calculate_distmatrix_utils.R**, 
**source/baySeq.R** (added in this version) \

Other files were not made for executing and contain functions for different parts of the pipelines.


----------------------------------------------------------------
Contribution to the original Hobotnica source code
-----------------
- **'source/baySeq.R'** is an entirely new script with all baySeq analysis and volcano-plot visualisations (up & running) 
- **'run.R'** several lines of code added to implement baySeq under Hobotnica`s pipeline 
- **source/signatures_utils.R** added code to extract top differentially expressed genes from baySeq, 
created a duplicate heatmap function for the visualisation of baySeq results 
- Venn diagram is disabled until further notice

----------------------------------------------------------------

Results of RNA-seq macrophages data using Hobotnica
-----------------
![](https://github.com/tony-zhelonkin/Hobotnica_baySeq_ed/blob/main/hobotnica_mf)

Although the phenotypic variation under the experiment is obvious, the pattern itself is variable. 
In our case there were almost no intersections with the best tool 
on top differentially expressed genes identified with other DE tools 

![](https://github.com/tony-zhelonkin/Hobotnica_baySeq_ed/blob/main/signature_vis) 

According to the result of the Hobotnica analysis limma-voom seems to be the tool of choice 
in DE analysis of RNA-seq macrophage data 

![](https://github.com/tony-zhelonkin/Hobotnica_baySeq_ed/blob/main/h_scores.png) 

----------------------------------------------------------------

Conclusion on the analysis
-----------------

Results of DE tools comparison using Hobotnica expose current problems of existing software in DE analysis.  
Each tool has its own specificities in the assumptions on the statistical properties inherent to RNA-seq data. 
No consensus as of today exists on the subject. Hobotnica is intended to put some clarity and may be a 
helpful instrument in the hands of a thoughtful bioinformatician

----------------------------------------------------------------

Perspectives
----------------
- Add more DE tools into Hobotnica
- Add new features: comparing DE tools run on simulated counts data, etc.
- Validate the results with RT-PCR, proteomics analysis data





