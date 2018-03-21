
#
# task 1: given a conventional name for a TF, return
# addresses of putative binding sites (as recorded using FIMO)
# in a GRanges.  you must define a source and a genomic region
# for 
#

#' utility to generate link to biocfound bucket for FIMO TFBS scores
#' @param tag character(1) token identifying TF
#' @examples
#' URL_s3_tf
#' @export
URL_s3_tf = function(tag="M3433") {
 sprintf("http://s3.amazonaws.com/biocfound-tfbs/%s.02sort.bed.gz",
    tag)
 }

#' utility to read FIMO outputs (after sorting and tabix indexing applied as in contents of biocfound-tfbs S3 bucket)
#' @importFrom Rsamtools scanTabix
#' @importFrom utils read.delim
#' @param tf character(1) file id
#' @param which a GRanges delimiting the extract
#' @examples
#' requireNamespace("GenomicRanges")
#' requireNamespace("IRanges")
#' importFIMO_S3("M5946_1", 
#'     which=GenomicRanges::GRanges("chr1", IRanges::IRanges(1,15000)))
#' @export
importFIMO_S3 = function( tf, which ) {
 stopifnot(length(tf)==1, is(tf, "character"))
 url = URL_s3_tf(tag=tf)
 tmp = scanTabix( url, param = which )
 con = textConnection(paste(tmp[[1]], collapse="\n"))
 on.exit(close(con))
 read.delim(con, h=FALSE)
}

#' acquire a GRanges representing binding sites for a TF
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors "mcols<-"
#' @importFrom S4Vectors "mcols"
#' @importFrom GenomeInfoDb "genome<-"
#' @param conv_name character(1) conventional name for TF
#' @param map data.frame with columns PWMid (in form Mnnnn_n,
#' denoting the PWM defining the motif used for FIMO scoring) and
#' HGNC (the conventional name for the TF, as a 'gene symbol')
#' @param which a GRanges instance limiting the request
#' @param meth a function with two parameters, tf (a character(1) Mnnnn
#' string) and which (a GRanges instance) that retrieves TF 'bed'-like
#' records from a resource
#' @param genomeTag character(1) tag for genome build in use
#' @note We assume that the retrieval yields a data.frame in the format
#' ```
#'     V1    V2    V3               V4        V5 V6       V7 
#'  chr1 11496 11511 chr1:11496-11511 -1.555560  - 5.92e-04 
#'  chr1 11520 11535 chr1:11520-11535 -3.161620  + 9.12e-04
#' ```
#' As of 21 Mar 2018 there are only two TFs in S3: VDR, and
#' @examples
#' requireNamespace("GenomicRanges")
#' requireNamespace("IRanges")
#' tfbs_as_GRanges(conv_name="VDR", map=TFutils::fimoMap, 
#'      which=GenomicRanges::GRanges("chr1", IRanges::IRanges(1,25000), 
#'      genomeTag="hg19"), meth=importFIMO_S3)
#' @export
tfbs_as_GRanges = function( conv_name, map, which, meth, genomeTag ) {
    stopifnot(!missing(which), is(which, "GRanges"))
    stopifnot(is(meth, "function"))
    mname = map$PWMid[which(map$HGNC == conv_name)][1]
    if (is.na(mname)) stop( paste("cannot find PWMid for", conv_name) )
    text = meth( mname, which )
    ans = GRanges(text$V1, IRanges(text$V2, text$V3), strand=text$V6)
    mcols(ans) = text[,c("V5", "V7")]
    names(mcols(ans)) = c("score", "p")
    if (!missing(genomeTag)) genome(ans) = genomeTag
    ans
}
