### ----------------------------------------------------------------------------
### Motivating Example - Binomial data
### Written by Eduard Szöcs

### -------- Packages ----------------------------------------------------------
require(reshape2)
require(multcomp)

### --------- Custom Functions -------------------------------------------------
#' pairwise wilcox.test
#' @param y numeric; vector of data values
#' @param g factor; grouping vector
#' @param dunnett logical; if TRUE dunnett contrast, otherwise Tukey-contrasts
#' @param padj character; method for p-adjustment, see ?p.adjust.
pairwise_wilcox <- function(y, g, dunnett = TRUE, padj = 'holm'){
  tc <- t(combn(nlevels(g), 2))
  if(dunnett){
    tc <- tc[tc[ ,1] == 1, ]
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


### -------- Load + Clean ------------------------------------------------------
df <- read.table(header = TRUE, text = 'conc A B C D
0 1 1 0.9 0.9
32 0.8 0.8 1 0.8
64 0.9 1 1 1 
128 0.9 0.9 0.8 1
256 0.7 0.9 1 0.5
512 0.4 0.3 0.4 0.2')
df

dfm <- melt(df, id.vars = 'conc', value.name = 'y', variable.name = 'tank')
# conc as factor
dfm$conc <- factor(dfm$conc)
head(dfm)

boxplot(y ~ conc, data = dfm, 
        xlab = 'conc', ylab = 'Proportion surv.')


### -------- Methods -----------------------------------------------------------
### -------- Normal + Transformation
dfm$y_asin <- ifelse(dfm$y == 1, 
                     asin(sqrt(dfm$y)) - asin(sqrt(1/40)), 
                     asin(sqrt(dfm$y)) 
)

modlm <- aov(y_asin ~ conc, data = dfm)
summary(modlm)
# F-test
drop1(modlm, test = 'F')
# LOEC
summary(glht(mod, linfct = mcp(conc = 'Dunnett')), test = adjusted('holm'))


### -------- Binomial GLM
modglm <- glm(y ~ conc , data = dfm, family = binomial, weights = rep(10, nrow(dfm)))
# LR-test
drop1(modlm, test = 'Chisq')
summary(glht(modglm, linfct = mcp(conc = 'Dunnett')), test = adjusted('holm'))


