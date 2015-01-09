### ----------------------------------------------------------------------------
### Motivating Example - Count data
### Written by Eduard Szöcs

### -------- Setup ------------------------------------------------------------
prj <- "/home/edisz/Documents/Uni/Projects/PHD/6USETHEGLM/"
datadir <- file.path(prj, "data")   # data

### -------- Packages
require(MASS)
require(ggplot2)
require(scales)
require(reshape2)
require(plyr)
require(lmtest)
require(multcomp)

### --------- Custom Functions
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

#' Parametric bootstrap (PB)
#' @param m1 full model
#' @param m0 reduced model
#' @param data data object used to fit the models
#' @param npb number of boostraps

myPBmodcomp <- function(m1, m0, data, npb){
  ## create reference distribution of LR and coefs
  myPBrefdist <- function(m1, m0, data){
    # simulate from null
    x0 <- simulate(m0)
    # refit
    newdata0 <- data
    newdata0[ , as.character(formula(m0)[[2]])] <- x0
    m1r <-  try(update(m1, .~., data = newdata0))
    m0r <- try(update(m0, .~., data = newdata0))
    # simulate from model
    x1 <- simulate(m1)
    # refit
    newdata1 <- data
    newdata1[ , as.character(formula(m1)[[2]])] <- x1
    m1r1 <-  try(update(m1, .~., data = newdata1))
    
    # check convergence (otherwise return NA for LR)
    if(!is.null(m0r[['th.warn']]) | !is.null(m1r[['th.warn']]) | 
         inherits(m0r, "try-error") | inherits(m1r, "try-error")){
      LR <- 'convergence error'
    } else {
      LR <- -2 * (logLik(m0r) - logLik(m1r))
    }
    # convergence for coefs
    if(!is.null(m1r1[['th.warn']]) | inherits(m1r1, "try-error")){
      coefs <- 'convergence error'
    } else {
      coefs <- coef(m1r1)
    }
    # return LR and coefs.
    out <- list(LR = LR, coefs = coefs)
    return(out)
  }
  
  ## calculate reference distribution
  ref <- replicate(npb, myPBrefdist(m1 = m1, m0 = m0, data = data), 
                   simplify = FALSE)
  LR <- sapply(ref, function(x) x[['LR']])
  coefs <- lapply(ref, function(x) x[['coefs']])
  
  # check convergence
  nconv_LR <-  LR == 'convergence error'
  # rm those
  LR <- as.numeric(LR[!nconv_LR])
  nconv_LR <- sum(!nconv_LR)
  
  nconv_coefs <- sapply(coefs, function(x) any(x == 'convergence error'))
  # rm those
  coefs[nconv_coefs] <- NULL
  nconv_coefs <- sum(!nconv_coefs)
  coefs <- do.call(rbind, coefs)
  
  ## original stats
  LRo <- c(-2 * (logLik(m0) - logLik(m1)))
  COEFo <- coef(m1)
  
  ## LR p-value 
  p.pb <- mean(c(LR, LRo) >= LRo, na.rm = TRUE)
  
  ## p-value for coef
  # proportions of boostrap coefs > or < 0
  twosidep<-function(coef){
    p1 <- sum(coef > 0) / length(coef)
    p2 <- sum(coef < 0) / length(coef)
    p <- min(p1,p2) * 2
    return(p)
  }
  p.coef <- apply(coefs, 2, twosidep)
  return(list(nconv_LR = nconv_LR, nconv_coefs = nconv_coefs,
              p.pb = p.pb, p.coef = p.coef))
}



### -------- Load + Clean ------------------------------------------------------
df <- read.table(file.path(datadir, 'brock.csv'), header = TRUE, sep = ';')
df$Concentration <- factor(df$Concentration)
plot(Abundance ~ Concentration, data = df)



### -------- Methods -----------------------------------------------------------
### --------
### Normal + Transformation
## log(Ax+1) transformation
# Ax = 2 , for min(x) & x != 0
A <- 2 / min(df$Abundance[df$Abundance != 0])
df <- transform(df, log_Abundance = log(A*Abundance + 1))
plot(log_Abundance ~ Concentration, data = df)

## Fit Normal model
modlm <- lm(log_Abundance ~ Concentration, data = df)

## Inference, global test
# F-test
drop1(modlm, test = 'F')
# LR-test
drop1(modlm, test = 'Chisq')

## Inference, parametres
summary(modlm)
coef(summary(modlm))
# Wald-t
prettyNum(p.adjust(coef(summary(modlm))[ , 'Pr(>|t|)'], method = 'holm'))[2:6]
# same as
summary(glht(modlm, linfct = mcp(Concentration = 'Dunnett')),  
        test = adjusted('holm'))

## predicted mean values with back transformation
modlm_p <- predict(modlm, newdata = data.frame(Concentration = unique(df$Concentration)), 
                 se.fit = TRUE)
# mean (backtransfored)
modlm_fit <- (exp(modlm_p$fit) - 1) / A
# ci
modlm_lwr <- (exp(modlm_p$fit - 1.96*modlm_p$se.fit) - 1) / A
modlm_upr <- (exp(modlm_p$fit + 1.96*modlm_p$se.fit) - 1) / A



### --------
### Poisson GLM
modpois <- glm(Abundance ~ Concentration, data = df, family = poisson)
summary(modpois)
# residual deviance >> residual degrees of freedom
modpois$deviance / modpois$df.residual
# overdispersion!

## Inference, global test
# LR test
drop1(modpois, test = 'Chisq')
# # parametric bootstrap
# modpois0 <- glm(Abundance ~ 1, data = df, family = poisson)
# myPBmodcomp(modpois, modpois0, data = df, npb = 500)
# # p-value of boostraped LR is 0.001996008

## Inference, test of parameters
# Wald-z
summary(modpois)
prettyNum(p.adjust(coef(summary(modpois))[ , 'Pr(>|z|)'], method = 'holm'))[2:6]
summary(glht(modpois, linfct = mcp(Concentration = 'Dunnett')),  
        test = adjusted('holm'))
# parametric bootstrap
# see above, all p = 0

## predicted value on response scale
modpois_p <- predict(modpois, newdata = data.frame(Concentration = unique(df$Concentration)), 
                   se.fit = TRUE)
modpois_fit <- modpois$family$linkinv(modpois_p$fit)
modpois_lwr <- modpois$family$linkinv(modpois_p$fit - 1.96 * modpois_p$se.fit)
modpois_upr <- modpois$family$linkinv(modpois_p$fit + 1.96 * modpois_p$se.fit)



### --------
### quasi-Poisson GLM
modqpois <- glm(Abundance ~ Concentration, data = df, family = quasipoisson)

## Inference, global test
# F test
drop1(modqpois, test = 'F')

## Inference, test of parameters
# Wald-t
summary(modqpois)
prettyNum(p.adjust(coef(summary(modqpois))[ , 'Pr(>|t|)'], method = 'holm'))[2:6]
# multicomp performs a Wald-z
summary(glht(modqpois, linfct = mcp(Concentration = 'Dunnett')),  
        test = adjusted('holm'))

## predicted value on response scale
modqpois_p <- predict(modqpois, newdata = data.frame(Concentration = unique(df$Concentration)), 
                     se.fit = TRUE)
modqpois_fit <- modqpois$family$linkinv(modqpois_p$fit)
modqpois_lwr <- modqpois$family$linkinv(modqpois_p$fit - 1.96 * modqpois_p$se.fit)
modqpois_upr <- modqpois$family$linkinv(modqpois_p$fit + 1.96 * modqpois_p$se.fit)



### --------
### Negative binomial GLM
modnb <- glm.nb(Abundance ~ Concentration, data = df)

## Inference, global test
# LR test
# refit null model, to estimate theta of null
modnb0 <- glm.nb(Abundance ~ 1, data = df)
anova(modnb, modnb0, test = 'Chisq')
# parametric bootstrap
set.seed(1234)
myPBmodcomp(modnb, modnb0, data = df, npb = 2000)
# p-value of boostraped LR is 0.054
# 1 model did not converge

## Inference, test of parameters
# Wald-z
summary(modnb)
prettyNum(p.adjust(coef(summary(modnb))[ , 'Pr(>|z|)'], method = 'holm'))[2:6]
summary(glht(modnb, linfct = mcp(Concentration = 'Dunnett')),  
        test = adjusted('holm'))
# parametric bootstrap
# see above, 
# $p.coef
# (Intercept) Concentration0.1 Concentration0.3   Concentration1   Concentration3  Concentration10 
# 0.000            0.214            0.026            0.034            0.000            0.004 

## predicted value on response scale
modnb_p <- predict(modnb, newdata = data.frame(Concentration = unique(df$Concentration)), 
                     se.fit = TRUE)
modnb_fit <- modnb$family$linkinv(modnb_p$fit)
modnb_lwr <- modnb$family$linkinv(modnb_p$fit - 1.96 * modnb_p$se.fit)
modnb_upr <- modnb$family$linkinv(modnb_p$fit + 1.96 * modnb_p$se.fit)



### --------
### Non-parametric

# Kruskal-Wallis test
kruskal.test(Abundance ~ Concentration, data = df)
# pairwise Wilcox
pairwise_wilcox(df$Abundance, factor(df$Concentration), dunnett = TRUE, padj = 'holm')


### --------------- Summarize --------------------------------------------------
### Mean variance relationship
mv <- ddply(df, .(Concentration), summarise,
      m = mean(Abundance),
      var = var(Abundance))
plot(var ~ m, data = mv, xlab = 'Mean', ylab = 'Variance', pch = 16)
# Poisson
abline(0,1, col = 'darkred', lwd = 2)
# quasi-poisson
abline(0, summary(modqpois)$dispersion, col = 'steelblue', lwd = 2)
# negative binomial
m <- seq(min(mv$m), max(mv$m), 0.01)
lines(m, m + m^2 / modnb$theta, col = 'darkgreen', lwd = 3)
legend(80, 1000, 
       legend = c('Poisson','quasi-Poisson', 'negativ binomial'), 
       lwd = c(2,2,2), col = c('darkred', 'steelblue', 'darkgreen'))



### Plot data and model results
moddf <- data.frame(conc = rep(unique(df$Concentration), 4),
                    fit = c(modlm_fit, modpois_fit, modqpois_fit, modnb_fit),
                    upr = c(modlm_upr, modpois_upr, modqpois_upr, modnb_upr),
                    lwr = c(modlm_lwr, modpois_lwr, modqpois_lwr, modnb_lwr),
                    mod = rep(c('Normal', 'Poisson', 'quasi-Poisson', 'Negative binomial'), each = 6), 
                    strindsAsFactor = FALSE)
moddf$mod <- factor(moddf$mod, levels = c('Normal', 'Poisson', 'quasi-Poisson', 'Negative binomial'))

p <- ggplot() +
  # raw data
  geom_boxplot(data = df, aes(x = Concentration, y = Abundance), 
               width = 0.3, fill = 'grey85') +
  geom_point(data = df, aes(x = Concentration, y = Abundance)) +
  # estimated means + CI
  geom_point(data = moddf, 
             aes(x = 0.1 + rep(1:6,4) + rep(c(0.1, 0.2, 0.3, 0.4), each = 6), 
                 y = fit), size = 4, shape = 17) +
  geom_errorbar(data = moddf, 
                aes(x = 0.1 + rep(1:6,4) + rep(c(0.1, 0.2, 0.3, 0.4), each = 6), 
                    ymax = upr, ymin = lwr, linetype = mod), 
                width = 0.1, lwd = 1) + 
  # loec_bars
  # normal
  geom_path(aes(x = c(5-0.15, 5+0.55), y = c(220,220)), linetype = 'solid', size = 1) + 
  # pois
  geom_path(aes(x = c(2-0.15, 6+0.55), y = c(220-5,220-5)), linetype = '1111', size = 1) + 
  # neg bin
  geom_path(aes(x = c(3-0.15, 3+0.55), y = c(220-3*5,220-3*5)), linetype = 'dotdash', size = 1) + 
  geom_path(aes(x = c(5-0.15, 6+0.55), y = c(220-3*5,220-3*5)), linetype = 'dotdash', size = 1) + 
  scale_linetype_manual(values=c('solid', '1111', 'dashed', 'dotdash')) + 
  mytheme +
  guides(linetype=guide_legend(title='Model')) +
  xlab('Concetration [mg/L]') +
  theme(axis.title.x = element_text(face="bold",size=18),
        axis.text.x  = element_text(size=14),
        axis.text.y  = element_text(size=14),
        axis.title.y = element_text(face="bold",size=18),
        legend.position=c(0.8, 0.7),
        legend.key.width = unit(1, "cm"),
        legend.background = element_rect(color = 'black')
        )
p
ggsave(file.path(figdir, 'example.pdf'), p, height = 6, width = 8)
# ggsave(file.path(markdir, 'figure1.eps'), p)

