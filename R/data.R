#datalist		gwascat_hg19.rda	named_tf.rda
#fimoMap.rda		motif_id_map.rda	symbol_map.rda
#gwascat_GRCh38.rda	named_metadata_tf.rda	tftColl.rda

#' gwascat_hg19: GRanges of march 21 2018 EBI gwascat
#' @importFrom utils data
#' @docType data
#' @format GenomicRanges GRanges instance
#' @source gwascat::makeCurrentGwascat, with gwascat:::lo38to19 applied
#' @examples
#' data(gwascat_hg19)
#' gwascat_hg19[,1:5]
"gwascat_hg19"

#' gwascat_GRCh38: GRanges of march 21 2018 EBI gwascat
#' @importFrom utils data
#' @docType data
#' @format GenomicRanges GRanges instance
#' @source gwascat::makeCurrentGwascat, with a little extra metadata
#' @examples
#' data(gwascat_GRCh38)
#' gwascat_GRCh38[,1:5]
"gwascat_GRCh38"

#' fimoMap: table with Mnnnn (motif PWM tags) and HGNC symbols for TFs
#' @importFrom utils data
#' @docType data
#' @format data.frame
#' @source Kimberly Glass (rekrg@channing.harvard.edu) 
#' @examples
#' data(fimoMap)
#' head(fimoMap)
"fimoMap"

#' tftColl: GSEABase GeneSetCollection for transcription factor targets
#' @importFrom utils data
#' @docType data
#' @format GSEABase GeneSetCollection instance
#' @source broad institute
#' @examples
#' data(tftColl)
#' tftColl
"tftColl"
