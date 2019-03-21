#' Cell types as labels of example scRNA-seq dataset(exprdata)
#'
#' A vector contains the cell types  as labels for cells in
#' example scRNA-seq dataset(exprdata)
#'
#' @format A vector contains 15 cell types with its length 1059
#' @source \url{http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE48968-GPL13112.rds}
"celltype"


#' Processed scRNA-seq expression profile of dendritic dataset
#'
#' A gene count matrix with its rows are genes and columns are cells. The entry
#' of the matrix stands for the abundance of mRNA  in certain cell.
#'
#' @format A matrix with interger entries and have 2000 rows and 1059 columns
#' @source \url{http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE48968-GPL13112.rds}
"exprdata"
