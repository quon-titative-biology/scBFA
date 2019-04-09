#' This function should be called to initialize input parameters into the
#' main scbfa function
#'
#' @return A model environment containing the following parameters:
#' {A,Z,V,U,\eqn{\beta},\eqn{\gamma},\eqn{\epsilon}}.
#'
#' @param modelEnv Empty R environment variable to contain following parameters:
#' {A,Z,V,U,\eqn{\beta},\eqn{\gamma},\eqn{\epsilon}}
#' @param GeneExpr G by N rawcount matrix,
#' in which rows are genes and columns are cells
#' @param X N by C cell-specific covariate matrix(e.g batch effect),
#' in which rows are cells,columns are number of covariates.
#'  If no such covariates are available, then X = NULL
#' @param Q G by T gene-specific covariate matrix(e.g quality control measures),
#' in which rows are genes columns are number of covariates,
#' If no such covariates are available, then Q = NULL
#' @param numFactors Numeric value, number of latent dimensions
#' @param epsilon Numeric value, parameter to control the strength of
#' regularization
#' @param initCellcoef Initialization of C by G gene-specific coefficient matrix
#' as user-defined coefficient \eqn{\beta}.
#' Such user defined coefficient can be applied to address confounding batch
#' effect
#' @param updateCellcoef Logical value, parameter to decide whether to
#' update C by G gene-specific coefficient matrix.
#' Again, when the cell types are confounded with technical batches or
#' there is no cell level covariate matrix,
#' the user can keep the initialization of coefficients as known estimate.
#' @param updateGenecoef Logical value, parameter to decide whether to update
#' N by T gene-specific coefficient matrix.
#' Again, when there is no gene level covariate matrix,
#' this value should be FALSE by default.
#'
#' @keywords internal

InitBinaryFA <- function(modelEnv,
                        GeneExpr,
                        numFactors,
                        epsilon,
                        X=NULL,
                        Q = NULL,
                        initCellcoef,
                        updateCellcoef,
                        updateGenecoef){

    modelEnv$numCells <- ncol(GeneExpr);

    modelEnv$numGenes <- nrow(GeneExpr);

    modelEnv$numFactors <- numFactors


    # Note that X doesn't contain intercept since intercept is not regularized
    #Cell wise Intercept is parameterized as U in the later initialization
    if(!is.null(X)){
        modelEnv$X <- X

        modelEnv$numCoef_X <- ncol(modelEnv$X)

    }else if(is.null(X)){

        modelEnv$X <- matrix(0,nrow = modelEnv$numCells,ncol = 1)

        modelEnv$numCoef_X <- ncol(modelEnv$X)

    }

    # Note that Q doesn't contain intercept since intercept is not regularized
    #Gene wise Intercept is parameterized as V in the later initialization
    if(!is.null(Q)){

        modelEnv$Q <- Q

        modelEnv$numCoef_Q<- ncol(modelEnv$Q)

    }else if(is.null(Q)){

        modelEnv$Q <- matrix(0,nrow = modelEnv$numGenes,ncol = 1)

        modelEnv$numCoef_Q <- ncol(modelEnv$Q)

    }




    # Initialization of N by K low dimensional embedding matrix
    modelEnv$ZZ <- matrix(rnorm(modelEnv$numCells * modelEnv$numFactors),
                        nrow = modelEnv$numCells, ncol = modelEnv$numFactors)
    # Initialization of G by K compressed feature space matrix
    modelEnv$AA <- matrix(rnorm(modelEnv$numGenes * modelEnv$numFactors),
                        nrow = modelEnv$numGenes, ncol = modelEnv$numFactors)
    # Initialization of N by T gene specific coefficient matrix
    modelEnv$gamma <- matrix(rnorm(modelEnv$numCells* modelEnv$numCoef_Q),
                        nrow = modelEnv$numCells,ncol = modelEnv$numCoef_Q)
    # Initialization of N by 1 cellwise intercept
    modelEnv$UU <- matrix(1,nrow = modelEnv$numCells, ncol = 1)
    # Initialization of G by 1 genewise offset
    modelEnv$VV <- matrix(1,nrow = modelEnv$numGenes, ncol = 1)
    # Initialization of C by G coefficient matrix
    if(is.null(initCellcoef)){
        # if the user doesn't provide initialization of coefficient beta
        # then randomly initialize coefficient beta
        modelEnv$beta <- matrix(rnorm(modelEnv$numCoef_X *modelEnv$numGenes ),
                            nrow = modelEnv$numCoef_X,ncol = modelEnv$numGenes)
    }else if(!is.null(initCellcoef)){
        # if the user provides initialization of coefficient beta,
        # then initialize coefficient beta with user-provided initialization
        modelEnv$beta <- initCellcoef
    }

    modelEnv$updateCellcoef = updateCellcoef

    # if there is no cell covariate matrix,
    # then we don't update coefficient matrix to avoid extra computation
    if(is.null(X)){
        print("No cell covariate presents,
                not updating cell coefficient matrix")
        modelEnv$updateCellcoef = FALSE
    }
    # if there is no gene covariate matrix
    # then we don't update coefficient matrix to avoid extra computation
    if(is.null(Q)){
        print("No gene covariate presents,
                not updating gene coefficient matrix")
        modelEnv$updateGenecoef = FALSE
    }
    # Matrix B: Denoting whether a gene is detected
    modelEnv$BB <- t((GeneExpr+0)!=0) + 0
    # Epsilon/number of cells (epsilon_1 in online methods)
    modelEnv$regularize_per_cell = (epsilon/modelEnv$numCells)
    # Epsilon/number of genes (epsilon_2 and epsilon_3 in online methods)
    modelEnv$regularize_per_gene = (epsilon/modelEnv$numGenes)
    # Construc a long vector since the input of optim()
    # needs parameters in the form of vector instead of list.
    modelEnv$parameters <- c(AA = modelEnv$AA, ZZ = modelEnv$ZZ,
                            beta = modelEnv$beta, gamma = modelEnv$gamma,
                        UU = modelEnv$UU, VV = modelEnv$VV,epsilon = epsilon)

    modelEnv$epsilon = epsilon

    return(modelEnv)

}

#' Restore the vector of parameter space into their seperated parameterization
#'
#' @return A list parameters containing the following parameters:
#' {A,Z,U,V,beta,gamma,epsilon}
#'
#' @param parameters: Vectorized parameter space.
#' @param modelEnv: Environment variable contains parameter space ,
#' and global variables such as N,G,C,T,detection matrix B etc
#'
#' @keywords internal
restore <- function(parameters,modelEnv){

    para_names <- names(parameters)
    param = list()
    param$ZZ <- matrix(parameters[grepl("ZZ",para_names)],
                        nrow = modelEnv$numCells, ncol = modelEnv$numFactors)
    param$AA <- matrix(parameters[grepl("AA",para_names)],
                        nrow = modelEnv$numGenes,ncol = modelEnv$numFactors)
    param$beta <- matrix(parameters[grepl("beta",para_names)],
                        nrow = modelEnv$numCoef_X,ncol = modelEnv$numGenes)
    param$gamma <- matrix(parameters[grepl("gamma",para_names)],
                        nrow = modelEnv$numCells,ncol = modelEnv$numCoef_Q)
    param$UU <- parameters[grepl("UU",para_names)]
    param$VV <- parameters[grepl("VV",para_names)]
    param$epsilon <- parameters[grepl("epsilon",para_names)]
    return(param)

}


#' Calculate negative penalized likelihood, used for calls
#' to the optim() function.
#'
#' The penalized likelihood function:
#' \eqn{f(A,Z,\beta,O,U) = \sum [lnP(B;A,Z,U,V,\beta,\gamma)]_ij -
#' \epsilon_1 * ||A||_2^2 - \epsilon_2 * ||Z||_2^2 -
#' \epsilon_3*||\beta||_2^2 - \epsilon_2 * ||\gamma||_2^2}
#'
#' @return Scalar penalized likelihood
#' @param parameters Vectorized parameter space.
#' @param modelEnv Environment variable contains parameter space,
#' and global variables such as N,G,C,detection matrix B, etc
#'
#' @keywords internal

neg_loglikelihood <- function(parameters,modelEnv){

    param <- restore(parameters,modelEnv)
    # Matrix product Z * A
    WW <- tcrossprod(param$ZZ,param$AA)
    # linear predictor
    linearpredictor <-t(t(WW+(modelEnv$X %*% param$beta)+(param$UU))+param$VV) +
        param$gamma %*% t(modelEnv$Q)

    Log1pexp <- log(1+exp(linearpredictor))

    infix = is.infinite(Log1pexp)
    # Fix numberic issues if exp(linearpredictor) is too large
    if(sum(infix)>0){

        Log1pexp[infix] <- linearpredictor[infix]

    }

    LL <- sum(modelEnv$BB * linearpredictor - Log1pexp)
    penalty_cell=0.5*(modelEnv$regularize_per_cell)*norm(param$ZZ,type = "F")^2+
        0.5*(modelEnv$regularize_per_cell)*norm(param$gamma,type = "F")^2
    penalty_gene=0.5*(modelEnv$regularize_per_gene)*norm(param$AA,type = "F")^2+
        0.5*(modelEnv$regularize_per_gene)*norm(param$beta,type = "F")^2
    penalizedLL <- LL - penalty_cell - penalty_gene

    if(is.na(sum(penalizedLL))){stop("NA generated in likelihood calculation")}

    return(-penalizedLL)

}

#' Calculate gradient of the negative log likelihood,
#' used for calls to the optim() function.
#'
#' @return Vectorized gradient
#'
#' @param parameters Vectorized parameter space.
#' @param modelEnv Environment variable contains parameter space,
#' and global variables such as N,G,C,detection matrix B, etc
#'
#' @keywords internal

gradient <- function(parameters,modelEnv){

    param <- restore(parameters,modelEnv)
    # Matrix product Z * A
    WW = tcrossprod(param$ZZ,param$AA)
    # linearpredictors
    linearpredictor <-t(t(WW+(modelEnv$X %*% param$beta)+(param$UU))+param$VV) +
        param$gamma %*% t(modelEnv$Q)
    #Common calulation value needed in the calculation of gradient
    common <- modelEnv$BB - (1/(1 + exp(-linearpredictor)))
    # Gradient of A
    gr_AA <- t(common) %*% param$ZZ - param$AA * (modelEnv$regularize_per_gene)
    # Gradient of Z
    gr_ZZ <- common %*% param$AA- param$ZZ * (modelEnv$regularize_per_cell)
    # Gradient of beta
    if(modelEnv$updateCellcoef == TRUE){
        gr_beta<-t(t(common)%*%modelEnv$X)-
            param$beta*(modelEnv$regularize_per_gene)
    }else if(modelEnv$updateCellcoef == FALSE){
        gr_beta = 0 * param$beta
    }
    # Gradient of gamma
    if(modelEnv$updateGenecoef == TRUE){
        gr_gamma<-common %*% modelEnv$Q-
            param$gamma * (modelEnv$regularize_per_cell)
    }else if(modelEnv$updateGenecoef == FALSE){
        gr_gamma = 0 * param$gamma
    }

    # Gradient of U
    gr_UU <- rowSums(common)
    # Gradient of O
    gr_VV <- colSums(common)
    # Gradient vector
    grad<- c(gr_AA, gr_ZZ, gr_beta,gr_gamma,gr_UU,gr_VV, 0)

    if(is.na(sum(grad))){stop("NA generated in gradient calculation")}

    return(-grad)

}

#' Optimize parameters of BFA's likelihood function
#'
#' @return The entire model environment
#'
#' @param modelEnv Environment variable contains parameter space,
#' and global variables such as N,G,C,detection matrix B, etc
#' @param maxit Maximum number of iteration with respect to objective function,
#' default is 300 iterations
#' @param method Optimization method, default is the conjugate gradient approach
#' L-BFGS-B is recommended for smaller dataset less than 10k cells
#'
#' @keywords internal
OptimBFA <- function(modelEnv,maxit,method){

    opt <- optim(par = modelEnv$parameters,
                fn = neg_loglikelihood,
                gr = gradient,
                modelEnv = modelEnv,
                method = method,
                control = list(maxit =maxit))

    param <- restore(opt$par,modelEnv)

    modelEnv$AA <- param$AA;
    modelEnv$ZZ <- param$ZZ;
    modelEnv$beta <- param$beta;
    modelEnv$gamma <- param$gamma;
    modelEnv$UU <- matrix(param$UU,ncol = 1);
    modelEnv$VV = matrix(param$VV,ncol = 1);


    modelEnv$parameters <- c( AA = modelEnv$AA,
                            ZZ = modelEnv$ZZ,
                            beta = modelEnv$beta,
                            gamma = modelEnv$gamma,
                            UU = modelEnv$UU,
                            VV = modelEnv$VV,
                            epsilon = modelEnv$epsilon)

    return(modelEnv)

}



#' Perform BFA model on the expression profile
#'
#' @return A model environment containing all parameter space of a BFA model
#' as well as global variables needed for calculation:
#' @return \eqn{A}: \eqn{G} by \eqn{K} compressed feature space matrix
#' @return \eqn{Z}: \eqn{N} by \eqn{K} low dimensional embedding matrix
#' @return \eqn{\beta}: \eqn{C} by \eqn{G} cell level coefficient matrix
#' @return \eqn{\gamma}: \eqn{N} by \eqn{T} gene level coefficient matrix
#' @return \eqn{V}: \eqn{G} by \eqn{1} offset matrix
#' @return \eqn{U}: \eqn{N} by \eqn{1} offset matrix
#'
#' @param scData can be a raw count matrix,
#' in which rows are genes and columns are cells;
#' can be a seurat object; can be a SingleCellExperiment object.
#' @param X \eqn{N} by \eqn{C} covariate matrix,e.g batch effect,
#' in which rows are cells,columns are number of covariates.Default is NULL
#' @param Q G by T gene-specific covariate matrix(e.g quality control measures),
#' in which rows are genes columns are number of covariates,
#' If no such covariates are available, then Q = NULL
#' @param numFactors  Numeric value, number of latent dimensions
#' @param maxit Numeric value, parameter to control the
#' Maximum number of iterations in the optimization, default is 300.
#' @param method Method of optimization,default is conjugate gradient approach.
#' @param initCellcoef Initialization of C by G gene-specific coefficient matrix
#' as user-defined coefficient \eqn{\beta}.
#' Such user defined coefficient can be
#'  applied to address confounding batch effect
#' @param updateCellcoef Logical value, parameter to decide whether to
#' update C by G gene-specific coefficient matrix.
#' Again, when the cell types are confounded with technical batches or
#' there is no cell level covariate matrix,
#' the user can keep the initialization of coefficients as known estimate.
#' @param updateGenecoef Logical value, parameter to decide whether to update
#' N by T gene-specific coefficient matrix.
#' Again, when there is no gene level covariate matrix,
#' this value should be FALSE by default.
#'
#' @importFrom zinbwave orthogonalizeTraceNorm
#' @importFrom SummarizedExperiment assay
#' @examples
#'
#' ## Working with Seurat or SingleCellExperiment object
#'
#'library(Seurat)
#'library(SingleCellExperiment)
#'
#'
#' ## Input expression profile, 5 genes x 3 cells
#'
#'GeneExpr = matrix(rpois(15,1),nrow = 5,ncol = 3)
#'rownames(GeneExpr) = paste0("gene",seq_len(nrow(GeneExpr)))
#'colnames(GeneExpr) = paste0("cell",seq_len(ncol(GeneExpr)))
#'celltype = as.factor(sample(c(1,2,3),3,replace = TRUE))
#'
#'## Create cell level technical batches
#'
#'batch = sample(c("replicate 1","replicate 2","replicate 2"))
#'X = matrix(NA,nrow = length(batch),ncol = 1)
#'X[which(batch =="replicate 1"), ] = 0
#'X[which(batch =="replicate 2"), ] = 1
#'rownames(X) = colnames(GeneExpr)
#'
#'## run BFA with raw count matrix
#'
#'bfa_model = scbfa(scData = GeneExpr,X = scale(X),numFactors =2)
#'
#'## Create Seurat object for input to BFA
#'
#'scData = CreateSeuratObject(raw.data = GeneExpr,project="sc",min.cells = 0)
#'
#'## Standardize the covariate matrix should be a default operation
#'
#'bfa_model = scbfa(scData = scData, X = scale(X), numFactors = 2)
#'
#'## Build the SingleCellExperiment object for input to BFA
#'
#'## Set up SingleCellExperiment class
#'
#'sce <- SingleCellExperiment(assay = list(counts = GeneExpr))
#'
#'## Standardize the covariate matrix should be a default operation
#'
#'bfa_model = scbfa(scData = sce, X = scale(X), numFactors = 2)
#'
#' @keywords export
#' @export
scbfa <- function(scData,
                numFactors,
                X=NULL,
                Q = NULL,
                maxit = 300,
                method = "L-BFGS-B",
                initCellcoef = NULL,
                updateCellcoef = TRUE,
                updateGenecoef = TRUE) {
    # extract raw count matrix
    # extract the class scData object
    objClass = class(scData)

    if(objClass == "matrix"){
        GeneExpr = scData
    }else if(objClass == "seurat"){
        sce = Convert(from = scData, to = "sce")
        GeneExpr = as.matrix(assay(sce))
    }else if(objClass == "SingleCellExperiment"){
        GeneExpr = as.matrix(assay(scData))
    }


    # initialization
    modelEnv = new.env();

    InitBinaryFA( modelEnv,
                GeneExpr = GeneExpr,
                X= X,
                Q = Q,
                numFactors = numFactors,
                epsilon =max(dim(GeneExpr)),
                initCellcoef = initCellcoef,
                updateCellcoef = updateCellcoef,
                updateGenecoef = updateGenecoef);
    # optimization
    modelEnv = OptimBFA(modelEnv,maxit = maxit,method = method)
    # Orthogonalization
    orthogonalizefactors<-orthogonalizeTraceNorm(U = modelEnv$ZZ,
                                        V = t(modelEnv$AA),
                                        a=0.5*(modelEnv$regularize_per_cell),
                                        b=0.5*(modelEnv$regularize_per_gene))

    modelEnv$AA <- t(orthogonalizefactors$V) ;
    modelEnv$ZZ <- orthogonalizefactors$U

    return(modelEnv)
}


#' Function to get low dimensional embedding matrix
#'
#' @return Z: \eqn{N} by \eqn{K} low dimensional embedding
#'
#' @param modelEnv output environment variable
#'
#' @examples
#' GeneExpr = matrix(rpois(15,1),3,5)
#' bfa_model = scbfa(scData = GeneExpr,X = NULL,numFactors =2)
#' Z = getScore(bfa_model)
#' @keywords export
#' @export
getScore <- function(modelEnv){return(modelEnv$ZZ)}

#' Function to get low dimensional loading matrix
#'
#' @return \eqn{A}: \eqn{G} by \eqn{K} compressed feature space
#'
#' @param modelEnv output environment variable
#'
#' @examples
#' GeneExpr = matrix(rpois(15,1),3,5)
#' bfa_model = scbfa(scData = GeneExpr,X = NULL,numFactors =2)
#' A = getLoading(bfa_model)
#' @keywords export
#' @export
getLoading <- function(modelEnv){return(modelEnv$AA)}





#' Performs Binary PCA (as outlined in our paper).
#' This function take the input of gene expression profile and
#' perform PCA on gene detection pattern
#'
#' @return A list with class "prcomp",containing the following components:
#' @return sdev: the standard deviations of the principal components
#' (i.e., the square roots of the eigenvalues of
#' the covariance/correlation matrix, though the calculation is
#' actually done with the singular values of the data matrix).
#' @return rotation: the matrix of variable loadings
#' (i.e., a matrix whose columns contain the eigenvectors).
#' The function princomp returns this in the element loadings.
#' @return x: the rotated data (the centred (and scaled if requested)
#' data multiplied by the rotation matrix) is returned.
#' Hence, cov(x) is the diagonal matrix diag(sdev^2).
#' @return center, scale. centering and scaling used, or FALSE.
#'
#' @param scData can be a raw count matrix,
#' in which rows are genes and columns are cells;
#' can be a seurat object; can be a SingleCellExperiment object.
#' @param X \eqn{N} by \eqn{C} covariate matrix,e.g batch effect,
#' in which rows are cells,columns are number of covariates.
#' If no such covariates available X = NULL
#' @param scale. Logical value isndicating whether the variables should be
#' scaled to have unit variance before the analysis takes place.
#' In general scaling is not advisable, since we think the variance in the
#' gene detection space is potentially associated with celltypes
#' (e.g cell type specific markers)
#' @param center Logical value indicating whether the variables
#' should be shifted to be zero centered
#'
#' @import SingleCellExperiment
#' @import Seurat
#' @importFrom SummarizedExperiment assay
#' @importFrom stats lm optim prcomp residuals var
#'
#' @examples
#'
#' ## Working with Seurat or SingleCellExperiment object
#'
#'library(Seurat)
#'library(SingleCellExperiment)
#'
#' ## Input expression profile, 5 genes x 3 cells
#'
#'GeneExpr = matrix(rpois(15,1),nrow = 5,ncol = 3)
#'rownames(GeneExpr) = paste0("gene",seq_len(nrow(GeneExpr)))
#'colnames(GeneExpr) = paste0("cell",seq_len(ncol(GeneExpr)))
#'celltype = as.factor(sample(c(1,2,3),3,replace = TRUE))
#'
#'## Create cell level technical batches
#'
#'batch = sample(c("replicate 1","replicate 2","replicate 2"))
#'X = matrix(NA,nrow = length(batch),ncol = 1)
#'X[which(batch =="replicate 1"), ] = 0
#'X[which(batch =="replicate 2"), ] = 1
#'rownames(X) = colnames(GeneExpr)
#'
#'##run BFA with raw count matrix
#'
#'bpca_model = BinaryPCA(scData = GeneExpr,X = scale(X))
#'
#'## Create Seurat object for input to BFA
#'
#'scData = CreateSeuratObject(raw.data = GeneExpr,project = "sc",min.cells = 0)
#'
#'## Standardize the covariate matrix should be a default operation
#'
#'bpca_model = BinaryPCA(scData = scData, X = scale(X))
#'
#'## Build the SingleCellExperiment object for input to BFA
#'
#'## Set up SingleCellExperiment class
#'
#'sce <- SingleCellExperiment(assay = list(counts = GeneExpr))
#'
#'## Standardize the covariate matrix should be a default operation
#'
#'bpca_model = BinaryPCA(scData = sce, X = scale(X))
#'
#' @keywords export
#' @export
BinaryPCA = function(scData,X=NULL,scale. = FALSE,center = TRUE){

    # extract raw count matrix
    # extract the class scData object
    objClass = class(scData)

    if(objClass == "matrix"){
        GeneExpr = scData
    }else if(objClass == "seurat"){
        sce = Convert(from = scData, to = "sce")
        GeneExpr = as.matrix(assay(sce))
    }else if(objClass == "SingleCellExperiment"){
        GeneExpr = as.matrix(assay(scData))
    }

    binaryEntry = t(GeneExpr!=0)+0
    # calculating the variance across gene detection pattern.
    binaryVariance = apply(binaryEntry,2,var)

    # construct covariate matrix under the condition
    #whether cell-level covariates exists
    if(!is.null(X)){
        X <- scale(X)
    }else if(is.null(X)){
        X <- matrix(0,nrow =nrow(binaryEntry), ncol = 1)
    }
    # obtain residuals of linear regression of binarized expression profile.
    # Note that genes with gene detection pattern doesn't vary across cells
    # gets removed at the beginning
    res <- residuals(lm(binaryEntry[,binaryVariance!=0] ~X))
    pca <- prcomp(res,scale. = scale.,center = center)
    return(pca)
}
