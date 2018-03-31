#' define a structure to hold information about TFs from diverse reference sources
#' @importFrom methods new show
#' @importFrom stats na.omit
#' @aliases "TFCatalog-class"
#' @export
setClass("TFCatalog", representation(name="character",
 nativeIds="character", HGNCmap="data.frame", metadata="ANY"))
#' Constructor for TFCatalog
#' @param name informative character(1) for collection
#' @param nativeIds character() vector of identifiers used by collection creators
#' @param HGNCmap data.frame with column 1 nativeIds, column 2 HGNC or hgnc.heur for MSigDb
#' and any other columns of use
#' @param metadata a list of metadata elements
#' @examples
#' TFCatalog
#' @export
TFCatalog = function(name, nativeIds, HGNCmap, metadata) {
  if (missing(metadata)) metadata=list()
  new("TFCatalog", name=name, nativeIds=nativeIds,
      HGNCmap=HGNCmap, metadata=metadata)
}
#' produce a concise report on TFCatalog instnace
#' @aliases show,TFCatalog-method
#' @param object instance of TFCatalog
#' @export
setMethod("show", "TFCatalog", function(object) {
 cat("TFutils TFCatalog instance", object@name, "\n")
 cat(sprintf(" %d native Ids, including\n", length(object@nativeIds)))
 cat(Biobase::selectSome(object@nativeIds, max=2), "\n")
 cat(sprintf(" %d unique HGNC tags, including\n", length(unique(object@HGNCmap[,2]))))
 cat(Biobase::selectSome(na.omit(object@HGNCmap[,2])), "\n")
})
