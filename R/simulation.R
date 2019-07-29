#' simulation to generate scRNA-seq data with varying level of
#' gene detection noise versus gene count noise
#'
#' @return GeneExpr,a count matrix with rows number of genes
#' and columns number of cells
#' @return celltype,a vector specify the corresponding
#' celltype after QC measures.
#'
#' @param zinb a ZINB-WaVE object representing ZINB-WaVE fit
#' to real data to get realistic simulation parameters
#' @param celltype a factor to specify the ground-truth cell types in the
#' original dataset that the parameter of zinb object is fit to.
#' Since we filter out some simulated cells due to low amount of genes detected
#' in that cell,
#' we subset the ground truth cell types correspondingly
#' @param disper numeric value, parameter to control the size factor
#' \eqn{r} in \eqn{NB(\mu, r)}.
#' r is varied in the set {0.5,1,5} in our simulation(as outlined in our paper)
#' @param var_count numeric value, parameter to control the noise level
#' added to a common embedding space to generate gene count matrix.
#' This parameter is formulated as \eqn{\sigma_\mu}
#' and and in the paper is selected from the set {0.1, 0.5, 1, 2, 3}
#' @param var_dropout numeric value, parameter to control the noise level added
#' to a common embedding space for to generate gene detection matrix.
#' This parameter is formulated as \eqn{\sigma_\pi}
#' and in the paper is selected from the set {0.1, 0.5, 1, 2, 3}
#' @param delta intercept to control the overall gene detection rate.
#' and in the paper is selected from the set {-2, -0.5, 1,2.5,4}
#'
#' @import zinbwave
#' @import MASS
#' @import SingleCellExperiment
#' @importFrom stats rbinom rnbinom rnorm
#' 
#' @keywords export
#'
#' @export
#' @examples
#'
#' ## raw counts matrix with rows are genes and columns are cells
#' data("zinb_toy",package = "scBFA", envir = environment())

#' ## a vector specify the ground truth of cell types provided by conquer database
#' data("celltype_toy",package = "scBFA",envir = environment())
#'
#' scData = scNoiseSim(zinb = zinb_toy,
#'          celltype = celltype_toy,
#'          disper = 1,
#'          var_dropout =1,
#'          var_count = 1,
#'          delta = 1)
#'

scNoiseSim = function(zinb,
                      celltype,
                      disper,
                      var_dropout = 1,
                      var_count = 1,
                      delta){
    # number of cells
    numCells = nrow(zinb@W)
    # number of genes
    numGenes = ncol(zinb@alpha_mu)
    # number of numFactors
    numFactors =ncol(zinb@W)
    # simulate a N by K low dimensional embedding for mean of gene count,
    # the mean of every row of the embedding space is the corresponding row of
    # embedding in zinb-wave
    score_count=t(apply(zinb@W,1,function(x){
        mvrnorm(1,mu = x,Sigma = diag(rep(var_count,numFactors),numFactors))}
    ))
    # compute a N by G mean expression level for gene count
    mu = exp(zinb@X %*% zinb@beta_mu +
                 t(zinb@gamma_mu) %*%  t(zinb@V) +
                 score_count %*% (zinb@alpha_mu) + zinb@O_mu)
    # simulate a N by K low dimensional embedding for gene detection space,
    # the mean of every row of the embedding space is the corresponding
    # row of embedding in zinb-wave
    score_dropout = t(apply(zinb@W,1,function(x){
        mvrnorm(1,mu = x,Sigma = diag(rep(var_dropout,numFactors),numFactors))}))
    # linear term to parameterize the  probability matrix to determine whether
    # a gene is zero or sampled from a NB distirbution
    linear = zinb@X %*% zinb@beta_pi +t(zinb@gamma_pi) %*%  t(zinb@V) +
        score_dropout %*% (zinb@alpha_pi) - delta
    # probability for a gene to sample zero
    Pi_sample_zero = 1/(1 + exp(-linear))
    # probability for a gene to sample from NB distribution
    Pi_sample_count = 1- Pi_sample_zero
    # observed count matrix O
    pseudoCounts = matrix(0,nrow = numCells, ncol = numGenes)

    for(ii in seq_len(numCells)){
        for(jj in seq_len(numGenes)){
            # whether to sample a zero or a negative binomial distribution
            pi = rbinom(1, size = 1, prob =Pi_sample_zero[ii,jj])

            # if pi = 1, the gene has zero count
            # if pi = 0, the gene is sampled from NB distribution
            if(pi == 0){
                pseudoCounts[ii,jj] = rnbinom(n = 1, size = disper , mu = mu[ii,jj])
            }
        }
    }
    # check if there is NA generated due to numeric reason.
    message(paste("the amount of NA: ",sum(is.na(pseudoCounts))))

    message(paste0("total detection rate:",sum(pseudoCounts != 0)/length(mu)))
    # calculate gene detection rate
    gdr = colSums(pseudoCounts != 0)/numCells
    # calculate cell detection rate
    cdr = rowSums(pseudoCounts!=0)/numGenes

    message("summary of gene detection rate")
    message(summary(gdr))

    message("summary of cell detection rate")
    message(summary(cdr))
    # quality control: filter out genes with gene detection rate less than 0.01
    useGene = gdr> 0.01
    # quality control fileter out cells with cell detection rate less than 0.01
    useCell = cdr >0.01

    pseudoCounts =  pseudoCounts[useCell,useGene]
    
    sce <- SingleCellExperiment(assay = list(counts = t(pseudoCounts)),colData = data.frame(celltype =celltype[useCell] ))
    

    return(sce)

}
