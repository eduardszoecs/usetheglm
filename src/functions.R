#' pairwise wilcox.test
#' @param y numeric; vector of data values
#' @param g factor; grouping vector
#' @param dunnett logical; if TRUE dunnett contrast, otherwise Tukey-contrasts
#' @param padj character; method for p-adjustment, see ?p.adjust.
pairwise_wilcox <- function(y, g, dunnett = TRUE, padj = 'holm'){
  tc <- t(combn(nlevels(g), 2))
  if(dunnett){
    tc <- tc[tc[ ,1] == levels(g)[1], ]
  }
  pval <- numeric(nrow(tc))
  for(i in seq_len(nrow(tc))){
    pval[i] <- wilcox.test(y[as.numeric(g) == tc[i, 1]], 
                           y[as.numeric(g) == tc[i, 2]])$p.value
  }
  pval <- p.adjust(pval, padj)
  names(pval) = paste(tc[,1], tc[,2], sep = '-')
  return(pval)
}
