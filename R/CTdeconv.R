
#' Cell types deconvolution
#'
#' This function integrates deconvolution results of three agrithom-signature
#' combinations which show best performance in our benchmark datasets and provides
#' the relative proportions of six immune cell types in mixture samples.
#'
#' @param mix A matrix (n Genes * m Samples) of genes expression from bulk samples
#' . It can be a string which specifies the path of a expression 'tsv' file. The rownames
#' must be gene symbols and the colnames are the sample id.
#
#' @param RNAseq Logical. if TURE, the expression data is RNA-Seq data rather than
#' Macroarray data, and the quantile normalization will not be done during Cibersort
#' analysis process. The default is FALSE.
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
#' res=CTdeconv(mix)
CTdeconv <- function(mix,RNAseq=F,filename=NULL) {
  # lm6='.gg/signature_rnaseq_geo60424_LM6.txt'
  # lm22='.gg/LM22.txt'
  # mix='.gg/examplemixture.TXT'
  # RNASeq=F
  if(isSingleString(mix)) mix=read.delim(mix,row.names = 1,check.names = F)
  mix=data.matrix(mix)
  ref='BRef'
  epicRes=EPIC(mix,ref)
  epicRes=epicRes[['mRNAProportions']]
  qn=!RNASeq
  cib_lm6_res=CIBERSORT(lm6,mix,QN=qn)
  cib_lm22_res=CIBERSORT(lm22,mix,QN=qn)
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
  if(!is.null(filename)){
    res1=data.frame(SampleId=rownames(res),res,check.names = F)
    write.table(res1,filename,row.names = F,sep='\t',quote = F)
  }
  return(res)
}


