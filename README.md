# CTdeconv

## Description

CTdeconv wraps three agrithom-signature combinations, which has shown high performance in our benchmark datasets, to deconvolute the relative proportions of six major immune cell types, including B cell, CD4 T cell, CD8 T cell, NK cell, Mono/Macro and Neutrophil. 


## Usage

The main function in this package is CTdeconv. It needs as input a matrix of the gene expression (from either RNA-seq or macroarray) of the samples for which to estimate cell proportions. Each column of expression matrix represents a sample and each row represents a gene.  Therefore, the rownames of the matrix are gene HGNC symbols while the colnames are sample IDs. 

Some parts of CTdeconv analysis are based on CIBERSORT to estimate cell type proportions. However, we can not provide the source code for this algorithm due to the license requirements. So you need to download it from its website http://cibersort.stanford.edu seperately, which is released under the Stanford Non-Commercial License.

In order to use CTdeconv, you will need to:

  1. Got to http://cibersort.stanford.edu
  2. Register and log in
  3. Download the latest R source code from the Download section. It would be a R file named as **'CIBERSORT.R'**
  4. Provide the Path to the R source code to the function `CTdeconv()` 


```R
# library(CTdeconv) ## Load the package or use CTdeconv::CTdeconv() directly.
ciberPath <- 'D:/User/xiergo/Documents/CIBERSORT.R'
res <- CTdeconv(mix = bulkSamplesMatrix, cibersortPath=ciberPath)
```
The bulk expression matrix can also be given as a path to a `.txt` file. In this file, the first row are sample IDs and the first column are gene symbols, and the seperator of each column is `'\t'`.
```R
bulkSampleFile <- 'D:/User/xiergo/Documents/bulkSampleExp.txt'
res <- CTdeconv(mix = bulkSampleFile, cibersortPath=ciberPath)
```
If the platform of your expression data is RNA-Seq rather than macroarray, you need to set RNAseq=T. In this case, the quantile normalization will be skipped in Cibersort analysis process, which is recommended to disabled by the author of Cibersort (see CIBERSORT website http://cibersort.stanford.edu).

```R
res <- EPIC(mix = bulkSamplesMatrix, cibersortPath=ciberPath, RNAseq=T)
```

The output `res` is a matrix with six columns. Each row represents a sample and each column represents a celltype. It provides the proprotions of six cell types in the bulk samples. Note that the proportion is relative since we conduct a normalization to make the sum of proportions in each sample be one. You can save the results as a file conveniently just by providing `CTdeconv()` with an output file path:

```R
outF <- 'D:/User/xiergo/Documents/CTdeconv_result.txt'
res <- EPIC(mix = bulkSamplesMatrix, cibersortPath=ciberPath, RNAseq=T, filename=outF)
```

Various other options are available and are well documented in the help pages from EPIC:
```R
?CTdeconv
```
## Installation

CTdeconv is implemented as an R package, which can be installed from GitHub by:
```R
# install devtools if necessary
install.packages('devtools')

# install CTdeconv
devtools::install_github('xiergo/CTdeconv')

# load
library(CTdeconv)
```

## License
CTdeconv can be used freely by academic groups for non-commercial purposes. 

Contact information
Xie Yuhao (xyhao@bjmu.edu.cn)

## FAQ

Who should I contact in case of a technical or other issue?
Xie Yuhao (xyhao@bjmu.edu.cn). Please provide as much details as possible and ideally send also an example input file (and/or reference profiles) that is causing the issue.
