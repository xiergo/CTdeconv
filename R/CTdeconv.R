#' CTdeconv: An immune cell-type deconvolution method.
#'
#' @docType package
#' @name CTdeconv
#' @import preprocessCore
#' @import parallel
NULL


#' Cell types deconvolution
#'
#' This function integrates deconvolution results of three agrithom-signature
#' combinations which show best performance in our benchmark datasets and provides
#' the relative proportions of six immune cell types in mixture samples.
#'
#' @param mix A matrix (n Genes * m Samples) of genes expression from bulk samples
#' . It can be a string which specifies the path of a expression 'tsv' file. The
#' rownames must be gene symbols and the colnames are the sample id.
#'
#' @param cibersortPath Path to cibersort R script. CIBERSORT is freely available
#' to academic users. Cibersort source script can be obtained
#' from \url{https://cibersort.stanford.edu}.
#
#' @param RNAseq Logical. If TURE, the expression data is RNA-Seq data rather than
#' Macroarray data, and the quantile normalization will not be done during Cibersort
#' analysis process. The default is FALSE.
#'
#' @param cellFrac Logical. If FALSE (default), the proportion of mRNA coming from each
#'  cell type will be returned. If TRUE, a renormalization process will be performed on
#'  the mRNA proportion to account for different mRNA / cell values in different cell types,
#'  as is recommended in EPIC. In this way, the proportion of cells from each cell type
#'  will be returned
#'
#' @param filename A string indecates the path of the file to be saved. If NULL (default),
#' no file will be saved.
#'
#' @return A matrix (m Samples * 6 Celltypes). It provides the proprotions of
#' six cell types in the bulk samples. Note that the proportion is relative since we
#' conduct a normalization to make the sum of proportions in each sample be one.
#'
#' @export
#'
#' @examples
#' # You need to provide path to CIBERSORT.R
#' path='D:/Users/xiergo/Documents/CIBERSORT.R'
#' res=CTdeconv(mix,cibersortPath=path)
CTdeconv <- function(mix,cibersortPath,RNAseq=F,cellFrac=F,filename=NULL) {
  # lm6='.gg/signature_rnaseq_geo60424_LM6.txt'
  # lm22='.gg/LM22.txt'
  # mix='.gg/examplemixture.TXT'
  # RNAseq=F
  if(isSingleString(mix)) mix=read.delim(mix,row.names = 1,check.names = F)
  mix=data.matrix(mix)
  ref=EPIC::BRef
  # require('EPIC')
  epicRes=EPIC::EPIC(mix,ref)
  epicRes=epicRes[['mRNAProportions']]
  qn=!RNAseq
  source(cibersortPath,local = T)

  tmp_mix = tempfile()
  tmp_lm6 = tempfile()
  tmp_lm22 = tempfile()
  write.mat(mix,tmp_mix)
  write.mat(lm6,tmp_lm6)
  write.mat(lm22,tmp_lm22)

  cib_lm6_res=CIBERSORT(tmp_lm6,tmp_mix,QN=qn)
  cib_lm22_res=CIBERSORT(tmp_lm22,tmp_mix,QN=qn)
  cename=list(B=c('B cells naive','B cells memory',"Bcells","B cells"),
              CD4=c('CD4_Tcells',"CD4 T cell",'CD4 T cells',"CD4.T.cells",
                           'T cells CD4 naive','T cells CD4 memory resting'),
              CD8=c('T cells CD8','CD8_Tcells','CD8 T cells',"CD8 T cell",'CD8.T.cells'),
              NK=c('NKcells','NK cells',"NK cell","NK.cells",'NK cells activated'),
              Mono_Macro=c('Monocytes','Macrophages M0','Macrophages M1','Macrophages M2'),
              Neutro=c('Neutrophil','Neutrophils'))
  resls=list(epicRes,cib_lm6_res,cib_lm22_res)
  resls1=lapply(resls, function(x){
    y=sapply(cename, function(i){
      x1=x[,colnames(x) %in% i,drop=F]
      x1=rowSums(x1)
    })
    colnames(y)=names(cename)
    y=y/rowSums(y)
    y
  })
  res=Reduce('+',resls1)/length(resls1)
  # print(res)
  renormalize_vector=c(
    'B'= 0.4016,
    'CD4'= 0.3952,
    'CD8'= 0.3952,
    'NK'= 0.4396,
    'Mono_Macro'= 1.4196,
    'Neutro'= 0.1300)
  if(cellFrac){
    res=t(t(res)/renormalize_vector)
    res=res/rowSums(res)
  }
  # print(res)
  if(!is.null(filename)){
    res1=data.frame(SampleId=rownames(res),res,check.names = F)
    write.table(res1,filename,row.names = F,sep='\t',quote = F)
  }
  return(res)
}

write.mat=function(mat,file){
  df=data.frame(ID=rownames(mat),mat,check.names = F)
  write.table(df,file = file,row.names = F,sep='\t',quote = F)
  return(NULL)
}

#single string detector
isSingleString <- function(input) {
  is.character(input) & length(input) == 1
}
