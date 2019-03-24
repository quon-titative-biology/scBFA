#' scRNA-seq dataset(exprdata)
#'
#' A matrix contains 950 cells and 500 genes. The source of this dataset is cDC/
#' pre-DC cells(see supplementary files)
#' We subset most variant 100 genes  as  example scRNA-seq dataset(exprdata)
#'
#' @format A 100 x 950 matrix with gene counts value
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89232}
"exprdata"

#' Cell types as labels of example scRNA-seq dataset(exprdata)
#'
#' A vector contains the cell types  as labels for cells in
#' example scRNA-seq dataset(exprdata)
#'
#' @format A vector contains 4 cell types with its length 950
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89232}
"celltype"

#' example zinb object after fitting a toy dataset with 5 cells and 10 genes
#'
#' The toy dataset is generated from the following code
#' require(zinbwave)
#' GeneExpr = matrix(rpois(50,1),nrow = 10,ncol = 5)
#' rownames(GeneExpr) = paste0("gene",seq_len(nrow(GeneExpr)))
#' colnames(GeneExpr) = paste0("cell",seq_len(ncol(GeneExpr)))
#' celltype = as.factor(sample(c(1,2,3),5,replace = TRUE))
#' zinb = zinbFit(Y = GeneExpr,K=2)
#'
#' @format A ZinbModel class contains all of its parameters after fitting
#' a toy dataset
"zinb_toy"


#' toy cell type vector with 3 cell types generated for 5 cells in toy dataset
#'
#' The cell type vector is generated from the following code
#'
#' celltype = as.factor(sample(c(1,2,3),5,replace = TRUE))
#'
#' @format A vector of length 5 contains 3 cell types
"celltype_toy"
