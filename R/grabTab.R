#' create table of TF targets and related metadata
#' @import dplyr
#' @import magrittr
#' @importFrom methods as is
#' @importFrom S4Vectors mcols
#' @param tfstub character(1) gene-like symbol for TF; will be grepped in names(gscoll)
#' @param gscoll a GSEABase GeneSetCollection
#' @param orgdb an instance of OrgDb as defined in AnnotationDbi
#' @param gwrngs a GRanges representing EBI gwascat, must have `DISEASE/TRAIT`, `MAPPED_GENE`
#' @examples
#' gt = grabTab("VDR", gscoll=TFutils::tftColl,
#'    orgdb=org.Hs.eg.db::org.Hs.eg.db, gwrngs=TFutils::gwascat_hg19_chr17)
#' dim(gt)
#' head(gt)
#' @export
grabTab = function(tfstub="STAT1", gscoll=TFutils::tftColl, 
     orgdb=org.Hs.eg.db::org.Hs.eg.db, gwrngs=TFutils::gwascat_hg19_chr17) {
MAPPED_GENE <- CHR_ID <- `DISEASE/TRAIT` <- NULL
CHR_POS <- REGION <- NULL
allst1 = unlist(lapply(gscoll[ grep(tfstub, names(gscoll)) ], geneIds))
st1syms = AnnotationDbi::mapIds(orgdb, keys=allst1, keytype="ENTREZID", column="SYMBOL")
chk = as(mcols(gwrngs), "data.frame") %>% 
     filter(MAPPED_GENE %in% st1syms) %>% select(`DISEASE/TRAIT`, MAPPED_GENE, CHR_ID, CHR_POS, REGION)
cbind(TF=tfstub, chk)
}

grabTab.old = function(tfHGNC, pgmap, gscoll, gwrngs) {
    ns = names(gscoll)
    nspl = strsplit(ns, "_")
    base = sapply(nspl, "[", 1)
    ans = pgmap[pgmap$HGNC == tfHGNC, ]
    if (length(hit <- which(tfHGNC == base))>0) {
     # direct hit -- there is a gene set whose prefix matches tfHGNC exactly
           gs = gscoll[[hit[1]]]  # FIXME should iterate if length(hit)>1
           geneIdType(gs) = SymbolIdentifier("org.Hs.eg.db")
           ans = data.frame(
              query=tfHGNC, setname=setName(gs), targs=geneIds(gs), stringsAsFactors=FALSE)
           targs = ans$targs
           hits = intersect(targs, gwrngs$MAPPED_GENE)
           if (length(hits)>0) {
             allt = rep(" ", length(targs))
             names(allt) = targs
             allt[hits] = "yes"
             ans$mappedInGwascat = allt
             }   
           }
     ans
}

