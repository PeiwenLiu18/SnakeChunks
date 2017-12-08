
edger <- function(){

    library(SARTools)

    baseDir <- getwd()

    targetFile <- snakemake@params[["targetFile"]]
    setwd(paste(snakemake@wildcards[["diffexpr_dir"]], "/", snakemake@wildcards[["test"]], "_vs_", snakemake@wildcards[["ref"]], sep=""))
    new.loc <- paste("edgeR/", targetFile, sep="")
    file.copy(targetFile, new.loc)

    setwd("edgeR")

    colors <- c("green", "red", "blue", "pink")

    # setting params

    projectName <- snakemake@params[["projectName"]]
    author <- snakemake@params[["author"]]

    rawDir <- snakemake@params[["rawDir"]]
    featuresToRemove <- snakemake@params[["featuresToRemove"]]
    varInt <- snakemake@params[["varInt"]]
    condRef <- snakemake@params[["condRef"]]
    batch <- snakemake@params[["batch"]]
    alpha <- snakemake@params[["alpha"]]
    pAdjustMethod <- snakemake@params[["pAdjustMethod"]]
    cpmCutoff <- 1 #snakemake@params[["cpmCutoff"]]
    gene.selection <- snakemake@params[["gene_selection"]]
    normalizationMethod <- snakemake@params[["normalizationMethod"]]
    workDir <- snakemake@params[["wd"]]


# checking parameters
checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
                      rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# edgeR analysis
out.edgeR <- run.edgeR(counts=counts, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

# MDS + clustering
exploreCounts(object=out.edgeR$dge, group=target[,varInt], gene.selection=gene.selection, col=colors)

# summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts=counts, alpha=alpha, col=colors)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))


writeReport.edgeR_modif <- function(target,counts,out.edgeR,summaryResults,majSequences,
                                    workDir,projectName,author,targetFile,rawDir,
                                    featuresToRemove,varInt,condRef,batch,
                                    alpha,pAdjustMethod,colors,gene.selection,
                                    normalizationMethod, cpmCutoff){
  rmarkdown::render(input=system.file("report_edgeR.rmd", package="SARTools"),
                    output_file=paste0(projectName, "_report.html"),
                    output_dir=workDir,
                    intermediates_dir=workDir,
                    knit_root_dir=workDir,
                    run_pandoc=TRUE,
                    quiet=TRUE,
                    clean=TRUE)
  cat("HTML report created\n")
}

# generating HTML report
writeReport.edgeR_modif(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
                  majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                  targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                  condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, colors=colors,
                  gene.selection=gene.selection, normalizationMethod=normalizationMethod,
                  cpmCutoff=cpmCutoff)


}

edger()


