sink(log)

library(caTools)
library(Rsamtools)
library(spp)
print(sessionInfo())

getwd()

## Fetch snakemake params, inputs and output
input.treatment <- snakemake@input[["treatment"]]
input.control <- snakemake@input[["control"]]

output.pdf <- snakemake@output[["pdf"]]
output.sm_density_wig <- snakemake@output[["sm_density_wig"]]
output.enrich_est_wig <- snakemake@output[["enrich_est_wig"]]
output.peaks_broadPeak <- snakemake@output[["peaks_broadPeak"]]
output.bdg_positions <- snakemake@output[["bdg_positions"]]
output.peaks_narrowPeak <- snakemake@output[["peaks_narrowPeak"]]
output.peaks_bed <- snakemake@output[["peaks_bed"]]

## Analyse data
treatment.data <- read.bam.tags(input.treatment)
control.data <- read.bam.tags(input.control)

binding.characteristics <- get.binding.characteristics(treatment.data,srange=c(50,500),bin=5, accept.all.tags = T) 

pdf(file=output.pdf,width=5,height=5)
par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
plot(binding.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation")
abline(v=binding.characteristics$peak$x,lty=2,col=2)
dev.off()

treatment.data.info <- select.informative.tags(treatment.data)
control.data.info <- select.informative.tags(control.data)

treatment.data.qua <- remove.local.tag.anomalies(treatment.data.info)
control.data.qua <- remove.local.tag.anomalies(control.data.info)

tag.shift <- round(binding.characteristics$peak$x/2)
smoothed.density <- get.smoothed.tag.density(treatment.data.qua,control.tags=control.data.qua,bandwidth=200,step=100,tag.shift=tag.shift)
writewig(smoothed.density,output.sm_density_wig,"s_1 smoothed, background-subtracted tag density")

enrichment.estimates <- get.conservative.fold.enrichment.profile(treatment.data.qua,control.data.qua,fws=500,step=100,alpha=0.01)
writewig(enrichment.estimates,output.enrich_est_wig,"conservative fold-enrichment/depletion estimates shown on log2 scale")

broad.clusters <- get.broad.enrichment.clusters(treatment.data.qua, control.data.qua, window.size=1e3, z.thr=3,tag.shift=round(binding.characteristics$peak$x/2))
write.broadpeak.info(broad.clusters,output.peaks_broadPeak)

fdr <- {params.fdr}
detection.window.halfsize <- binding.characteristics$whs;
bp <- find.binding.positions(signal.data=treatment.data.qua,control.data=control.data.qua,fdr=fdr,whs=detection.window.halfsize)
print(paste("detected",sum(unlist(lapply(bp$npl,function(d) length(d$x)))),"peaks"))
output.binding.results(bp,output.bdg_positions);

#bp <- add.broad.peak.regions(treatment.data.qua, control.data.qua, bp, window.size=1000,z.thr=3)
write.narrowpeak.binding(bp,output.peaks_narrowPeak)

# check if peak file is empty or not
bed <- matrix()
if (file.size(output.peaks_narrowPeak) > 0){
    bed <- read.table(output.peaks_narrowPeak)
    # replace negative values with 0
    for(i in 1:nrow(bed)){
if(bed[i,2] < 0){
    bed[i,2] <- 0
}
    }
write.table(bed, output.peaks_bed, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
} else {
    file.create(output.peaks_narrowPeak, output.peaks_bed, output.pdf, output.sm_density_wig, output.enrich_est_wig, output.peaks_broadPeak, output.bdg_positions)
}




sink()
