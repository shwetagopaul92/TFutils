
#
# task 1: given a conventional name for a TF, return
# addresses of putative binding sites (as recorded using FIMO)
# in a GRanges.  you must define a source and a genomic region
# for 
#

#' utility to generate link to biocfound bucket for FIMO TFBS scores
#' @param tag character(1) token identifying TF, can be an HGNC gene name or Mnnnn PWM tag.
#' It must be findable in TFutils::fimoMap table.
#' @return character(1) URL
#' @examples
#' URL_s3_tf
#' @export
URL_s3_tf = function(tag="M3433") {
 tab = TFutils::fimoMap
 pwmind = which(tab[,"PWMid"] == tag)
 if (length(pwmind)==1)
  return(sprintf("http://s3.amazonaws.com/biocfound-tfbs/%s.02sort.bed.gz",
    tag))
 else gind = which(tab[,"HGNC"] == tag)
 if (length(gind)==1)
  return(sprintf("http://s3.amazonaws.com/biocfound-tfbs/%s.02sort.bed.gz",
    tab[gind,"PWMid"]))
 else stop("could not locate tag in TFutils::fimoMap")
 }
