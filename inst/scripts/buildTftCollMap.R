
library(TFutils)
data(tftColl)
nmsig = names(tftColl)
alltok = strsplit(nmsig, "_") 
allhgnc = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype="SYMBOL")
allal = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype="ALIAS")
alsym = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
  keys=allal, column="SYMBOL", keytype="ALIAS")
names(alsym) = allal
names(allhgnc) = allhgnc
alltargs = c(allhgnc, alsym)
chk2 = sapply(alltok, function(x) sapply(x, function(z) 
  { 
  ans = try(alltargs[which(names(alltargs)==z)])
  if (inherits(ans, "try-error")) ans = NA_character_
  ans
  }))
newmap = cbind(tftname=names(tftColl), hgnc.heur=sapply(chk2, "[", 1), hgnc.fixed="no") 
# then edited by hand
# then mapped entries 556-615 using alis
