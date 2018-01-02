library("mosaics")

# Fetch Snakemake params, inputs & outputs
input.treatment <- snakemake@input[["treatment"]]
input.control <- snakemake@input[["control"]]


output.peaks_broadPeak <- snakemake@output[["peaks_broadPeak"]]
output.peaks_narrowPeak <- snakemake@output[["peaks_narrowPeak"]]
output.peaks_bed <- snakemake@output[["peaks_bed"]]

params.fdr <- snakemake@params[["fdr"]]

params.frag_len = snakemake@params[["frag_len"]]
params.bin_size = snakemake@params[["bin_size"]]
params.bin_treatment = snakemake@params[["bin_treatment"]]
params.bin_control = snakemake@params[["bin_control"]]

# Run peak-calling with Mosaics
constructBins(  input.treatment,  fileFormat="bed", outfileLoc=params.outdir,
                byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL,
                PET=FALSE, fragLen={params.frag_len}, binSize={params.bin_size}, capping=0 )

constructBins(  input.control,  fileFormat="bed", outfileLoc=params.outdir,
                byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL,
                PET=FALSE, fragLen={params.frag_len}, binSize={params.bin_size}, capping=0 )

fileName <- c(params.bin_treatment, params.bin_control)

bin <- readBins( type=c("chip","input"), fileName=fileName )
fit <- mosaicsFit( bin, analysisType="IO", bgEst="automatic" )
peak <- mosaicsPeak( fit, signalModel="2S", FDR={params.FDR}, maxgap=400, minsize=36, thres=10 )  ###
peak <- extractReads( peak,
    input.treatment, chipFileFormat="bed", chipPET=FALSE, chipFragLen={params.frag_len},
    input.control, controlFileFormat="bed", controlPET=FALSE, controlFragLen={params.frag_len}, parallel=FALSE, nCore=1 )

peak  <-  findSummit(peak, parallel=FALSE,  nCore=1  )
export(  peak,  type="bed",  filename=output.peaks_bed  )
export(  peak,  type="narrowPeak",  filename=output.peaks)

#export(  peak,  type="txt",  filename=".txt"  )
#export(  peak,  type="gff",  filename=".gff"  )

