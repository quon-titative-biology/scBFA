#' Perform diagnoisis of dispersion on the expression profile to check whether
#' scBFA works on specific dataset
#'
#' @return A Figure to tell the where the input data's dispersion ~ tpm curve
#' align to the 14 benchmark datasets in Figure 2.a or Gene detection rate
#'
#' @param scData can be a raw count matrix,
#' in which rows are genes and columns are cells;
#' can be a seurat object; can be a SingleCellExperiment object.
#' @param sampleInfo sample level feature matrix,e.g batch effect,experimental conditions
#' in which rows are cells,columns are number of covariates.Default is NULL
#' @param disperType a parameter to tell which dispersion estimate the user can plot
#' DESeq2 offers stepwise dispersion estimate, a gene wise dispersion estimate using
#' "GeneEst", dispersion estimate from fitted disperions ~ mean curve (using "Fitted")
#' And final MAP estimate,using "Map". Default value is "Fitted"
#' @param diagnose_feature a parameter to determin whether the user want to check GDR or dispersion.
#' @import ggplot2
#' @import DESeq2
#' @importFrom SummarizedExperiment mcols
#' @importFrom stats loess predict quantile
#' @importFrom utils data
#' @keywords export
#'
#' @export
#' @examples
#'
#' data(exprdata)
#' diagnose(scData = exprdata)

diagnose <- function(scData,sampleInfo = NULL,disperType = "Fitted",diagnose_feature="dispersion"){
        # colData should be data.frame that consists of feature matrix
        objClass = class(scData)
        # read data based on which class the input scData is
        if(objClass == "matrix"){
            GeneExpr = scData
        }else if(objClass == "seurat"){
            sce = as.SingleCellExperiment(scData)
            GeneExpr = as.matrix(assay(sce))
        }else if(objClass == "SingleCellExperiment"){
            GeneExpr = as.matrix(assay(scData))
        }
        # if sampleInfo matrix is NULL, make them to be column 1
        if(is.null(sampleInfo)){

            sampleInfo = data.frame(rep(1,ncol(GeneExpr)))

        }else{sampleInfo = as.data.frame(sampleInfo)}
        # sampleInfo matrix and expression matrix should share the same sample name
        rownames(sampleInfo) = colnames(GeneExpr)

        if(!identical(rownames(sampleInfo), colnames(GeneExpr))){
                print("Notes: Feature data and expression matrix have different sample names")
        }

    if(diagnose_feature == "dispersion"){

        data(disperPlot)
        # Note DESeq2 package replace design to be ~ 1 when estimating dispersion.
        ddsFullCountTable <- DESeqDataSetFromMatrix(
                countData = GeneExpr,
                colData = sampleInfo,
                design = ~ 1)

        counts<-counts(ddsFullCountTable)
        # obtain geometric means
        geoMeans = apply(counts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row + 1))))
        # estimate size factors
        dds = estimateSizeFactors(ddsFullCountTable, geoMeans=geoMeans)
        # estimate dispersions
        dds = estimateDispersions(dds)
        # which dispersion estimate to use
        if(disperType == "Map"){
            disper = dispersions(dds)
        }else if(disperType == "Fitted"){
            disper = mcols(dds)$dispFit
        }else if(disperType == "GeneEst"){
            disper = mcols(dds)$dispGeneEst
        }

        exprplus1 = t(GeneExpr) + 1
        # compute mean of tpm value
        meantpm = colMeans(log((exprplus1/rowSums(exprplus1)) * (10^6)))

        df = data.frame(disper = disper,meantpm = meantpm)
        # local linear regression to fit the curve
        loessMod <- loess(log(disper) ~ (meantpm), data=df, span=0.3)
        # fitted dispersion value after fitting mean of log(tpm)
        df$fitted_disper = predict(loessMod)

        quant = quantile(df$meantpm, probs = seq(0, 1, 0.025))
        # address border effect of LOESS
        df$plotOutlier <- !((df$meantpm > quant["2.5%"] & df$meantpm< quant["97.5%"]))

        df = df[!df$plotOutlier,]

        df = df[,!names(df) %in% c("disper","plotOutlier")]

        df$dataset = "new"

        df$variance = "new data"

        plotdf = rbind(disperPlot,df)

        plotdf$dataset_selection = paste(plotdf$dataset,plotdf$variance)

        alldatasets = c("mESCs","DC","Intestinal",
          "Pancreatic","HSCs","HSPC","H7-ESC",
          "LPS","LSK","MGE","Myeloid",
          "Memory T","PBMC","Dendritic")

        ll = paste(rep(alldatasets,2),rep(c("HVG","HEG"),each = length(alldatasets)))
        ll = c(ll,"new new data")
        plotdf$dataset_selection = factor(plotdf$dataset_selection,levels =ll)

        colorvalue = c("deepskyblue1","deepskyblue1","deepskyblue1","deepskyblue1",
                       "deepskyblue1","deepskyblue1","deepskyblue1","deepskyblue1",
                       "red3","deepskyblue1","deepskyblue1","deepskyblue1",
                       "deepskyblue1","deepskyblue1","red3","deepskyblue1",
                       "deepskyblue1","deepskyblue1","deepskyblue1","red3",
                       "red3","red3","red3","deepskyblue1", "deepskyblue1",
                       "deepskyblue1","deepskyblue1","deepskyblue1","black")

        p1 <- ggplot(plotdf,aes(x = meantpm,y = fitted_disper,colour =dataset_selection,linetype =variance)) + geom_line(size = 1)

        p1 <- p1 + scale_colour_manual(values = colorvalue,guide = FALSE) + scale_linetype_manual(values = c("solid","dashed","twodash"))
        #p1 <- p1 + ylim(-2.5,5)
        p1 <- p1 + ylab("log(dispersion)") + xlab("log(TPM)")

        p1 <- p1 + theme(legend.position = "none",
                         panel.background = element_blank(),
                         legend.key = element_blank(),
                         axis.text=element_text(size=20),
                         #axis.text = element_blank(),
                         axis.title = element_text(size=20),
                         axis.line.x = element_line(color="black"),
                         axis.line.y = element_line(color="black"),
                         legend.title = element_text(size=20),
                         legend.text = element_text(size=16))




        p1 + ggtitle("black: input    red: don't use scBFA (Group II)    blue: use scBFA (Group I)")

        p1

    }else if(diagnose_feature=="GDR"){

        GDR = sum(GeneExpr!=0)/length(c(GeneExpr))

        genewiseGDR= colMeans(t(GeneExpr!=0))

        return(list(GDR = GDR,genewiseGDR= genewiseGDR))
    }

}
