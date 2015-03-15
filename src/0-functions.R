## ----------------------------------------------------------------------------
### Custom function for simulation and plotting
### Written by Eduard Szöcs
### ----------------------------------------------------------------------------


### ----------------------------------------------------------------------------
#' pairwise wilcox.test
#' @param y numeric; vector of data values
#' @param g factor; grouping vector
#' @param dunnett logical; if TRUE dunnett contrast, otherwise Tukey-contrasts
#' @param padj character; method for p-adjustment, see ?p.adjust.
pairwise_wilcox <- function(y, g, dunnett = TRUE, padj = 'holm', alternative = 'less'){
  tc <- t(combn(nlevels(g), 2))
  # take on dunnett comparisons
  if(dunnett){
    tc <- tc[tc[ ,1] == 1, ]
  }
  pval <- numeric(nrow(tc))
  # use wilcox.exact (for tied data)
  for(i in seq_len(nrow(tc))){
    pval[i] <- wilcox.exact(y[as.numeric(g) == tc[i, 2]], 
                            y[as.numeric(g) == tc[i, 1]], exact = TRUE, 
                            alternative = alternative)$p.value
  }
  pval <- p.adjust(pval, padj)
  names(pval) = paste(tc[,1], tc[,2], sep = '-')
  return(pval)
}


### ----------------------------------------------------------------------------
### Parametric bootstrap (PB)
#' PB of LR statistic
#' @param m1 Full model
#' @param m0 reduced model
#' @param  data data used in the models
#' @return LR of boostrap
# generate reference distribution
myPBrefdist <- function(m1, m0, data){
  # simulate from null
  x0 <- simulate(m0)
  # refit with new data
  newdata0 <- data
  newdata0[ , as.character(formula(m0)[[2]])] <- x0
  m1r <-  try(update(m1, .~., data = newdata0), silent = TRUE)
  m0r <- try(update(m0, .~., data = newdata0), silent = TRUE)
  # check convergence (otherwise return NA for LR)
  if(inherits(m0r, "try-error") | inherits(m1r, "try-error")){
    LR <- 'convergence error'
  } else {
    if(!is.null(m0r[['th.warn']]) | !is.null(m1r[['th.warn']])){
      LR <- 'convergence error'
    } else {
      LR <- -2 * (logLik(m0r) - logLik(m1r))
    }
  }
  return(LR)
}

#' generate LR distribution and return p value
#' @param m1 Full model
#' @param m0 reduced model
#' @param data data used in m1 and m0
#' @param npb number of bootstrap samples
#' @return p-value of boostrapped LR values
myPBmodcomp <- function(m1, m0, data, npb){
  ## calculate reference distribution
  LR <- replicate(npb, myPBrefdist(m1 = m1, m0 = m0, data = data), 
                  simplify = TRUE)
  # rm those
  LR <- as.numeric(LR)
  nconv_LR <- sum(!is.na(LR))
  ## original stats
  LRo <- c(-2 * (logLik(m0) - logLik(m1)))
  ## p-value from parametric bootstrap
  p.pb <- mean(c(LR, LRo) >= LRo, na.rm = TRUE)
  return(list(nconv_LR = nconv_LR, p.pb = p.pb))
}



### ----------------------------------------------------------------------------
##### Simulation 1 -  Count data
### ----------------------
#' Function to create simulated data (6 groups)
#' @param N number of replicates per group
#' @param mu group means (vector of length 6)
#' @param theta overdispersion per group (vector of length 6) [variance = mu + mu^2/theta.]
#' @param nsims number of simulated datasets
#' @return list of two: x - factor of groups, y - simulated counts
dosim1 <- function(N, mu, theta, nsims = 100){
  Nj     <- rep(N, time = 6)                # 6 groups
  mus    <- rep(mu, times = Nj)             # vector of mus
  thetas <- rep(theta, times=Nj)            # vector of thetas
  x      <- factor(rep(1:6, times=Nj))      # factor
  y      <- replicate(nsims, rnegbin(sum(Nj), mus, thetas)) # draw from negative binomial distribution
  return(list(x = x, y = y))
}

### -----------------------------
#' Function to analyse simulated datasets
#' @param z a list as returned by dosim1
#' @param verbose logical; print output during run?
#' @param npb number of boostrap samples
#' @param nmax maximum number of simulations to analyes 
#' (should be always NULL = all simulations)
#' @return list of 9, with p-values and loecs from:
#' 1) F test for lm on ln(Ax+1) transformed variables [p_lm_f]
#' 2) LR test for neg.bin. glm [p_glm_lr]
#' 3) F test for quasi-poisson glm [p_qglm_f]
#' 4) parametric boostrap of neg.bin glm [p_glm_lrpb]
#' 5) kruskal-test
#' 6) LOEC from LM using one-sided Dunnett test [loeclm]
#' 7) LOEC from neg.bin. glm using one-sided Dunnett test [loecglm]
#' 8) LOEC from quasi-poisson glm using one-sided Dunnett test [loecqglm]
#' 9) LOEC from one-sided pairwise Wilcoxon with Holm-correction
resfoo1 <- function(z, verbose = TRUE, npb = 400, nmax = NULL){
  if(verbose){
    message('n: ', length(z$x) / 6, '; muc = ', mean(z$y[,1][z$x == 1]))
  }
  ana <- function(y, x, npb){
    # -------------
    # Transformations
    # ln(ax + 1) transformation
    A <- 1/min(y[y!=0])         
    yt <- log(A*y + 1)
    df <- data.frame(x, y, yt)
    
    # -------------
    # Models
    # gaussian
    modlm <- lm(yt ~ x, data = df)
    modlm.null <- lm(yt ~ 1, data = df)
    # negative binomial 
    modglm <- try(glm.nb(y ~ x, data = df), silent = TRUE)
    modglm.null <- try(glm.nb(y ~ 1, data = df), silent = TRUE)
    # quasipoisson (to tackle down convergence problems)
    modqglm <- glm(y ~ x, data = df, family = 'quasipoisson')
    modqglm.null <-  glm(y ~ 1, data = df, family = 'quasipoisson')
    # poisson
    modpglm <- glm(y ~ x, data = df, family = 'poisson')
    modpglm.null <-  glm(y ~ 1, data = df, family = 'poisson')
    
    # ------------- 
    # Test of effects
    # check convergence
    if(inherits(modglm, "try-error") | inherits(modglm.null, "try-error")){
      p_glm_lr <- 'convergence error'
      p_glm_lrpb <- 'convergence error'
    } else {
      if(!is.null(modglm[['th.warn']]) | !is.null(modglm.null[['th.warn']])){
        p_glm_lr <- 'convergence error'
        p_glm_lrpb <- 'convergence error'
      } else {
        p_glm_lr <- anova(modglm, modglm.null, test = 'Chisq')[2 , 'Pr(Chi)']
        # Parametric bootstrap for GLM LR
        glm_pb <- myPBmodcomp(modglm, modglm.null, data = df, npb = npb)
        p_glm_lrpb <- glm_pb$p.pb
      }
    }
    # F Tests
    p_lm_f <- anova(modlm, modlm.null, test = 'F')[2, 'Pr(>F)']
    p_qglm_f <- anova(modqglm, modqglm.null, test = 'F')[2, 'Pr(>F)']
    # LR Tests
    p_pglm_lr <- anova(modpglm, modpglm.null, test = 'Chisq')[2, 'Pr(>Chi)']
    # non-parametric test
    p_k <- kruskal.test(y ~ x, data = df)$p.value
    
    # ----------------
    # LOEC
    # lm
    mc_lm <- summary(glht(modlm, linfct = mcp(x = 'Dunnett'),  
                          alternative = 'less'), test = adjusted('holm'))$test$pvalues
    suppressWarnings( # intended warnings about no min -> no LOEC
      loeclm <- min(which(mc_lm < 0.05))
    )
    # negbin
    if(inherits(modglm, "try-error")){
      loecglm <- 'convergence error'
    } else {
      if(!is.null(modglm[['th.warn']])){
        loecglm <- 'convergence error'
      } else {
        mc_glm <- summary(glht(modglm, linfct = mcp(x = 'Dunnett'),  
                               alternative = 'less'), test = adjusted('holm'))$test$pvalues
        suppressWarnings(
          loecglm <- min(which(mc_glm  < 0.05))
        )
      }
    }
    # quasi
    mc_qglm <- summary(glht(modqglm, linfct = mcp(x = 'Dunnett'),  
                            alternative = 'less'), test = adjusted('holm'))$test$pvalues
    suppressWarnings( # intended warnings about no min -> no LOEC
      loecqglm <- min(which(mc_qglm < 0.05))
    ) 
    # pois
    mc_pglm <- summary(glht(modpglm, linfct = mcp(x = 'Dunnett'),  
                            alternative = 'less'), test = adjusted('holm'))$test$pvalues
    suppressWarnings( # intended warnings about no min -> no LOEC
      loecpglm <- min(which(mc_pglm < 0.05))
    ) 
    # pairwise wilcox
    suppressWarnings( # ties
      pw <- pairwise_wilcox(y, x, padj = 'holm', dunnett = TRUE)
    )
    suppressWarnings(
      loecpw <- min(which(pw < 0.05))
    )
    
    # ---------
    # return object
    return(list(p_lm_f = p_lm_f, p_glm_lr = p_glm_lr, p_qglm_f = p_qglm_f,
                p_glm_lrpb = p_glm_lrpb, p_pglm_lr = p_pglm_lr, p_k = p_k, 
                loeclm = loeclm, loecglm = loecglm, loecqglm = loecqglm, 
                loecpglm = loecpglm, loecpw = loecpw
    )
    )
  }
  # run on each simulated data
  if(!is.null(nmax)){
    res <- apply(z$y[ ,seq_len(nmax)], 2, ana, x = z$x, npb = npb)
  } else {
    res <- apply(z$y, 2, ana, x = z$x, npb = npb)
  }
  return(res)
}


### -----------------------------
#' Helper function to extract power for global tests from object returned by resfoo1
#' @param z an object as returned by resfoo1
#' @return a data.frame with power and convergence values
p_glob1 <- function(z){ 
  # extract p-values
  take <- c('p_lm_f', 'p_glm_lr','p_qglm_f', 'p_glm_lrpb', 'p_pglm_lr', 'p_k')
  ps <- ldply(z, function(w) as.numeric(unlist(w[take])))
  names(ps) <- take
  ps <- melt(ps)
  out <- ddply(ps, .(variable), summarize,
               power = sum(value < 0.05, na.rm = TRUE) / sum(!is.na(value)),
               conv = sum(!is.na(value)) / length(value))
  return(out)
}


### -----------------------------
#' Helper function to extract power for LEOC from object returned by resfoo1
#' @param u n object as returned by resfoo1
#' @param type what type of simulation was run ('t1' or 'power')
#' @return a data.frame with power values
#' @note for t1 error estimation LOEC should Inf (not there)
#' for power estimation loec should be at concentration 2
p_loec1 <- function(z, type = NULL){
  # extract p-values
  take <- c("loeclm", "loecglm", "loecqglm", "loecpglm", "loecpw")
  loecs <- ldply(z, function(w) as.numeric(unlist(w[take])))
  if(type == 't1'){
    # x should be Inf
    pow <- apply(loecs, 2, function(x) sum(x != Inf, na.rm = TRUE) / sum(!is.na(x)))
  } 
  if(type == 'power'){
    # x should be 2
    pow <- apply(loecs, 2, function(x) sum(x == 2, na.rm = TRUE) / sum(!is.na(x)))
  }
  names(pow) <- take
  return(pow)
}





### ----------------------------------------------------------------------------
##### Simulation 2 -  Proportions

### -----------------------------
#' Function to create simulated binomial data
#' @description Simulate data from binomial distribution.
#' @param N Number of replicates per group
#' @param pC probabilty in control groups
#' @param pE probabilty on effect groups
#' @param nsim number of simulated datasets
#' @param n_animals number of animals per replicate
#' 
#' @example 
#' # simulate 100 datasets. 
#' # Groups 1 & 2 with p = 0.9, Groups 3-6 with p = 0.3
#' sims <- dosim2(3)
#' # plot one realisation
#' plot(sims$x, sims$y[,1])
dosim2 <- function(N, pC = 0.95, pE = 0.3, nsim = 100, n_animals = 10){
  n_group <- 6        # number of groups
  p = c(rep(rep(pC, N), 2), rep(rep(pE, N), 4))    # expected proportions
  y <- replicate(nsim, rbinom(N * n_group, size = n_animals, prob = p))
  x      <- factor(rep(1:6, each = N))      
  return(list(x = x, y = y, n_animals = n_animals))
}


### -----------------------------
#' Function to analyse simulated datasets
#' @description Runs a normal model, a logistic model and a non-parametric 
#' tests on the simulated data
#' 
#' @param z simulated data, generated by dosim2()
#' @param verbose print status on the console?
#' @param asin Type of arcsine transformation. 'ecotox' or 'asin'. 
#' If 'ecotox' a special asin transformation is performed
#' @return list of 6, with p-values and loec from:
#' 1) F test for lm on arcsine transformed variables [lm_f]
#' 2) LR test for bin. glm [glm_lr]
#' 3) kruskal-test on untransformed values
#' 4) LOEC from LM using one-sided Dunnett test [loeclm]
#' 5) LOEC from bin. glm using one-sided Dunnett test [loecglm]]
#' 6) LOEC from one-sided pairwise Wilcoxon with Holm-correction
#' @example
#' sims <- dosim2(3)
#' resfoo2(sims)
resfoo2 <- function(z, verbose = TRUE, asin = 'ecotox'){
  if(verbose){
    message('n: ', length(z$x) / 6, '; muc = ', mean(z$y[,1][z$x == 1]) / 10)
  }
  
  ana <- function(y, x, n_animals, asin){
    # -------------
    # Transformations
    if(asin == 'ecotox'){
      y_asin <- ifelse(y  == 0, asin(sqrt(1 / (length(x) / 6 * n_animals))),
                       ifelse((y / n_animals) == 1, asin(1) - asin(sqrt(1 / (length(x) / 6 * n_animals))),
                              asin(sqrt(y / n_animals))))
    }
    if(asin == 'asin'){
      y_asin <- asin(sqrt(y / n_animals))
    }
    df <- data.frame(x, y, y_asin)
    
    # -------------
    # Models
    # Gaussian 
    modlm <- lm(y_asin ~ x, data = df)
    modlm.null <- lm(y_asin ~ 1, data = df)
    # binomial model
    modglm <- glm(cbind(y, n_animals - y) ~ x, data = df, 
                  family = binomial(link = 'logit'))
    modglm.null <- glm(cbind(y, n_animals - y) ~ 1, data = df, 
                       family = binomial(link = 'logit'))
    
    # ------------- 
    # Tests
    # LR Tests
    glm_lr <- anova(modglm, modglm.null, test = 'Chisq')[2 , 'Pr(>Chi)']
    # F Tests
    lm_f <- anova(modlm, modlm.null, test = 'F')[2, 'Pr(>F)']
    # kruskal test
    pk <- kruskal.test(y ~ x, data = df)$p.value
    
    # --------------
    # LOECs
    # multiple comparisons using one-sided Dunnett-contrasts
    mc_lm <- summary(glht(modlm, linfct = mcp(x = 'Dunnett'),  
                          alternative = 'less'), test = adjusted('holm'))$test$pvalues
    suppressWarnings( # intended warnings about no min -> no LOEC
      loeclm <- min(which(mc_lm < 0.05))
    )
    mc_glm <- summary(glht(modglm, linfct = mcp(x = 'Dunnett'),  
                           alternative = 'less'), test = adjusted('holm'))$test$pvalues
    suppressWarnings(
      loecglm <- min(which(mc_glm  < 0.05))
    )
    suppressWarnings(
      pw <- pairwise_wilcox(y, x, padj = 'holm', dunnett = TRUE)
    )
    suppressWarnings(
      loecpw <- min(which(pw < 0.05))
    )
    
    # --------------------------------------------
    # return object
    return(list(lm_f = lm_f, glm_lr = glm_lr, pk = pk, 
                loeclm = loeclm, loecglm = loecglm, loecpw = loecpw
    ))
    
  }
  # run models on simulated data
  res <- apply(z$y, 2, ana, x = z$x, n_animals = z$n_animals, asin = asin)
  res
}


### -----------------------------
#' Helper function to extract power for global tests from object returned by resfoo2
#' @param z an object as returned by resfoo2
#' @return a data.frame with power values
p_glob <- function(z){ 
  # extract p-values
  ps <- ldply(z, function(w) unlist(w)[1:3])
  pow <- apply(ps, 2, function(z) sum(z < 0.05, na.rm = TRUE)) / length(z)
  return(pow)
}

### -----------------------------
#' Helper function to extract power for LEOC from object returned by resfoo2
#' @param z n object as returned by resfoo2
#' @param type what type of simulation was run ('t1' or 'power')
#' @return a data.frame with power values
#' @note for t1 error estimation LOEC should Inf (not there)
#' for power estimation loec should be at concentration 2
p_loec <- function(z, type = NULL){
  # extract p-values
  loecs <- ldply(z, function(w) unlist(w)[4:6])
  if(type == 't1'){
    pow <- apply(loecs, 2, function(x) sum(x != Inf, na.rm = TRUE) / length(x))
  } 
  if(type == 'power'){
    pow <- apply(loecs, 2, function(x) sum(x == 2, na.rm = TRUE) / length(x))
  }
  return(pow)
}



### ----------------------------------------------------------------------------
### Plotting functions
# extract legend from ggplot
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# Function for strip labels
n_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="N") { 
    value <- paste0('n = ', value)
  }
  return(value)
}

# custom ggplot2 theme
mytheme <- theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size=14),
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=14,face="bold", vjust = 0),
        axis.title.y = element_text(size=14,face="bold", vjust = 1),
        legend.position="bottom",
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text= element_text(size=14, face = 'bold'))
