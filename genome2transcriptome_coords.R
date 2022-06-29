library(ensembldb)
library(genomation)
library(data.table)

args <- commandArgs(trailingOnly=T)

bed <- readBed(sprintf("%s.bed",args[1]))
gtf <- "reference_files/gencode.vM25.annotation.gtf"
DB <-  ensDbFromGtf(gtf=gtf,organism = "Mus_musculus",genomeVersion="GRCm38",version=100)
edb <- EnsDb(DB)
edb <- EnsDb("reference_files/Mus_musculus.GRCm38.100.sqlite")

bed_tx <- genomeToTranscript(bed, edb)

bed_dt <- as.data.table(bed_tx)

write.table(bed_dt,sprintf("%s.tx_coords.tab",args[1]),sep="\t",row.names = FALSE,quote=FALSE)
