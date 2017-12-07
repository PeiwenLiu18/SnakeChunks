
#f <- function(x) {
#  print(x)
#}
#targetFile <- snakemake@params[["targetFile"]]
#f(targetFile)

deseq2 <- function(){

    rm(list=ls())

    #sink("{log}")

    library(SARTools)

    baseDir <- getwd()

    targetFile <- snakemake@params[["targetFile"]]
    setwd(paste(snakemake@wildcards[["diffexpr_dir"]], "/", snakemake@wildcards[["test"]], "_vs_", snakemake@wildcards[["ref"]], sep=""))
    new.loc <- paste("DESeq2/", targetFile, sep="")
    file.copy(targetFile, new.loc)

    setwd("DESeq2")

    colors <- c("green", "red", "blue", "pink")

    # setting params

    projectName <- snakemake@params[["projectName"]]
    author <- snakemake@params[["author"]]

    rawDir <- snakemake@params[["rawDir"]]
    featuresToRemove <- snakemake@params[["featuresToRemove"]]
    varInt <- snakemake@params[["varInt"]]
    condRef <- snakemake@params[["condRef"]]
    batch <- snakemake@params[["batch"]]
    fitType <- snakemake@params[["fitType"]]
    cooksCutoff <- snakemake@params[["cooksCutoff"]]
    independentFiltering <- snakemake@params[["independentFiltering"]]
    typeTrans <- snakemake@params[["typeTrans"]]
    locfunc <- snakemake@params[["locfunc"]]
    alpha <- snakemake@params[["alpha"]]
    pAdjustMethod <- snakemake@params[["pAdjustMethod"]]
    workDir <- snakemake@params[["wd"]]

    # checking parameters
    checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
           rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
           condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
           independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
           typeTrans=typeTrans,locfunc=locfunc,colors=colors)

    # loading target file
    target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

    # loading counts
    counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

    # description plots
    majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

    # analysis with DESeq2
    out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
     locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
     cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

    # PCA + clustering
    exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

    # summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
    summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
      independentFiltering=independentFiltering,
      cooksCutoff=cooksCutoff, alpha=alpha)

    # save image of the R session
    save.image(file=paste0(projectName, ".RData"))

    # generating HTML report
    writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
       majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
       targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
       condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
       independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
       typeTrans=typeTrans, locfunc=locfunc, colors=colors)

    # get list of gene_ids of up/down genes
    #up <- as.vector(read.table(paste("tables/", list.files(path = "tables", pattern = "up.txt$")[1], sep=""))[,1])
    #down <- as.vector(read.table(paste("tables/", list.files(path = "tables", pattern = "down.txt$")[1], sep=""))[,1])

    up <- as.vector(read.table(paste("tables/", snakemake@wildcards[["test"]], "vs", snakemake@wildcards[["ref"]], ".up.txt", sep="")[1])[,1])
    down <- as.vector(read.table(paste("tables/", snakemake@wildcards[["test"]], "vs", snakemake@wildcards[["ref"]], ".down.txt", sep="")[1])[,1])

    setwd(baseDir)


    gene_list <- c(up[2:length(up)], down[2:length(down)])
    print(gene_list)
    write.table(gene_list, file=snakemake@output[["gene_list"]], row.names=F, col.names=F, quote=F)

    #sink()
}

deseq2()
