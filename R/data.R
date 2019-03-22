#' Cell types as labels of example scRNA-seq dataset(exprdata)
#'
#' A vector contains the cell types  as labels for cells in
#' example scRNA-seq dataset(exprdata)
#'
#' @format A vector contains 4 cell types with its length 950
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89232}
"celltype"


#' A zinb object after fitting a ZINB-WaVE on cDC/pre-DC dataset
#'
#' A zinb model object contains all the parameters of ZINB-WaVE model after
#' fitting dentritic dataset on
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89232
#
#'
#' @format A ZinbModel object contains all the parameters of ZINB-WaVE
#' @source \url{https://github.com/drisso/zinbwave}
"zinb"
