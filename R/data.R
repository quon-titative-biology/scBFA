#' scRNA-seq dataset(exprdata)
#'
#' A matrix contains 950 cells and 500 genes. The source of this dataset is cDC/
#' pre-DC cells(see supplementary files)
#' We subset most variant 100 genes  as  example scRNA-seq dataset(exprdata)
#'
#' @name exprdata
#' @docType data
#' @references \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89232}
#' @keywords data
#' @usage data(exprdata)
NULL

#' Cell types as labels of example scRNA-seq dataset(exprdata)
#'
#' A vector contains the cell types  as labels for cells in
#' example scRNA-seq dataset(exprdata)
#' @name celltype
#' @docType data
#' @references \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89232}
#' @keywords data
#' @usage data(celltype)
NULL

#' example zinb object after fitting a toy dataset with 5 cells and 10 genes
#'
#' The toy dataset is generated from the following code
#' require(zinbwave)
#' GeneExpr = matrix(rpois(50,1),nrow = 10,ncol = 5)
#' rownames(GeneExpr) = paste0("gene",seq_len(nrow(GeneExpr)))
#' colnames(GeneExpr) = paste0("cell",seq_len(ncol(GeneExpr)))
#' celltype = as.factor(sample(c(1,2,3),5,replace = TRUE))
#' zinb = zinbFit(Y = GeneExpr,K=2)
#' @name zinb_toy
#' @docType data
#' @keywords data
#' @usage data(zinb_toy)
NULL



#' toy cell type vector with 3 cell types generated for 5 cells in toy dataset
#'
#' The cell type vector is generated from the following code
#'
#' celltype = as.factor(sample(c(1,2,3),5,replace = TRUE))
#'
#' @name celltype_toy
#' @docType data
#' @keywords data
#' @usage data(celltype_toy)
NULL



#' Reference dataset(disperPlot)
#'
#' A dataframe contains all the gene-wise dispersion estimates loess curve for 14 datasets we benchmarked in Figure 2.a
#'
#' The variable in the columns are:
#' fitted_dispersion: the log value of gene-wise dispersion after fitting a loess curve with respect to TPM value.
#' Note that the genes at the top 2.5% and buttom 2.5% of the TPM value has been removed
#' meantpm is average tpm value calculated per gene
#' dataset are nams for datasets
#' variance is  gene selection method, here is HEG vs HVG
#' @name disperPlot
#' @docType data
#' @keywords data
#' @usage data(disperPlot)
NULL

