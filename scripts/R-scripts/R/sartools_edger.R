
edger <- function(){

    rm(list=ls())

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

    print("\n\nTEST1\n\")
    # checking parameters
    checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
          rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
          condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
          cpmCutoff=cpmCutoff,gene.selection=gene.selection,
          normalizationMethod=normalizationMethod,colors=colors)

    print("\n\nTEST2\n\")
    print(projectName)
    print(targetFile)

    # loading target file
    target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

    # loading counts
    counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

    # description plots
    majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

    print("\n\nTEST3\n\")

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

    # generating HTML report
    writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
      majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
      targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
      condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, colors=colors,
      gene.selection=gene.selection, normalizationMethod=normalizationMethod)

    # get list of gene_ids of up/down genes
    #up <- as.vector(read.table(paste("tables/", list.files(path = "tables", pattern = "snakemake@wildcards[["test}vssnakemake@wildcards[["ref}.up.txt$")[1], sep=""))[,1])
    #down <- as.vector(read.table(paste("tables/", list.files(path = "tables", pattern = "snakemake@wildcards[["test}vssnakemake@wildcards[["ref}.down.txt$")[1], sep=""))[,1])

    up <- as.vector(read.table(paste("tables/", snakemake@wildcards[["test"]], "vs", snakemake@wildcards[["ref"]], ".up.txt", sep="")[1])[,1])
    down <- as.vector(read.table(paste("tables/", snakemake@wildcards[["test"]], "vs", snakemake@wildcards[["ref"]], ".down.txt", sep="")[1])[,1])


    setwd(baseDir)

    gene_list <- c(up[2:length(up)], down[2:length(down)])
    print(gene_list)
    write.table(gene_list, file=snakemake@output[["gene_list"]], row.names=F, col.names=F, quote=F)
}

edger()


