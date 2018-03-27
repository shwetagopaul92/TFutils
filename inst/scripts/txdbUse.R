
library(TFutils)
data("named_tf")
fimo_bs = getBS("VDR", named_tf, GenomicRanges::GRanges("chr1", IRanges::IRanges(1,25000)),importFIMO_local)
library(feather)
library(dplyr)
library(magrittr)
library(data.table)
tfhash = read_feather(system.file("feather/tfhash2.feather", package = "TFutils"))
tfhash_vdr = tfhash %>% filter(tfname %like% 'VDR')
#Grab some target genes of VDR
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
seqlevels(txdb) <- "chr1" #looking at chromosome 1 
keys = as.character(tfhash_vdr$targets) # these are geneids
#select(txdb, keys = keys, columns="TXNAME", keytype="GENEID")
#GR <- transcripts(txdb, filter=list(tx_chrom = "chr1"))
GR <- genes(txdb, filter=list(tx_chrom = "chr1"))
library(GenomicRanges)
subsetByOverlaps(GR, fimo_bs)
