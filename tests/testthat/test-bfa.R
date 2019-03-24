context("scBFA")

test_that("bfa halts on numeric matrix input", {

    library(Seurat)
    library(SingleCellExperiment)

    ## Input expression profile, 5 genes x 3 cells

    GeneExpr = matrix(rpois(15,1),nrow = 5,ncol = 3)
    rownames(GeneExpr) = paste0("gene",seq_len(nrow(GeneExpr)))
    colnames(GeneExpr) = paste0("cell",seq_len(ncol(GeneExpr)))
    celltype = as.factor(sample(c(1,2,3),3,replace = TRUE))
    #'
    #'## Create cell level technical batches
    #'
    #'batch = sample(c("replicate 1","replicate 2","replicate 2"))

    ## Create cell level technical batches

    batch = sample(c("replicate 1","replicate 2","replicate 2"))
    X = matrix(NA,nrow = length(batch),ncol = 1)
    X[which(batch =="replicate 1"), ] = 0
    X[which(batch =="replicate 2"), ] = 1
    rownames(X) = colnames(GeneExpr)

    ## run BFA with raw count matrix
    GeneExpr = matrix(as.character(GeneExpr),nrow = nrow(GeneExpr),ncol = ncol(GeneExpr))


    ## Build the scAlign class object and compute PCs
    expect_error(scbfa(scData = GeneExpr,X = scale(X),numFactors =2))

})
