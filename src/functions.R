#! TODO:
#! Remove unused stuff from functions!

### ----------------------------------------------------------------------------
#' pairwise wilcox.test
#' @param y numeric; vector of data values
#' @param g factor; grouping vector
#' @param dunnett logical; if TRUE dunnett contrast, otherwise Tukey-contrasts
#' @param padj character; method for p-adjustment, see ?p.adjust.
pairwise_wilcox <- function(y, g, dunnett = TRUE, padj = 'holm'){
  tc <- t(combn(nlevels(g), 2))
  # take on dunnett comparisons
  if(dunnett){
    tc <- tc[tc[ ,1] == 1, ]
  }
  pval <- numeric(nrow(tc))
  # use wilcox.exact (for tied data)
  for(i in seq_len(nrow(tc))){
    pval[i] <- wilcox.exact(y[as.numeric(g) == tc[i, 1]], 
                           y[as.numeric(g) == tc[i, 2]], exact = TRUE)$p.value
  }
  pval <- p.adjust(pval, padj)
  names(pval) = paste(tc[,1], tc[,2], sep = '-')
  return(pval)
}


### ----------------------------------------------------------------------------
### Parametric bootstrap (PB)

### --------------------
### PB of LR statistic

# generate reference distribution
myPBrefdist <- function(m1, m0, data){
  # simulate from null
  x0 <- simulate(m0)
  # refit null
  newdata0 <- data
  newdata0[ , as.character(formula(m0)[[2]])] <- x0
  m1r <-  try(update(m1, .~., data = newdata0))
  m0r <- try(update(m0, .~., data = newdata0))
  # check convergence (otherwise return NA for LR)
  if(!is.null(m0r[['th.warn']]) | !is.null(m1r[['th.warn']]) | 
       inherits(m0r, "try-error") | inherits(m1r, "try-error")){
    LR <- 'convergence error'
  } else {
    LR <- -2 * (logLik(m0r) - logLik(m1r))
  }
  return(LR)
}

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
#' Function to create simulated data
dosim1 <- function(N, mu, theta, nsims = 100){
  Nj     <- rep(N, time = 6)                # 6 groups
  mus    <- rep(mu, times = Nj)             # vector of mus
  thetas <- rep(theta, times=Nj)            # vector of thetas
  x      <- factor(rep(1:6, times=Nj))      # factor
  y      <- replicate(nsims, rnegbin(sum(Nj), mus, thetas))
  return(list(x = x, y = y))
}


#' Function to analyse simulated datasets
#! TODO: Check convergence (th.warn.)
resfoo1 <- function(z, verbose = TRUE, npb = 400){
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
    modglm <- glm.nb(y ~ x, data = df)
    modglm.null <- glm.nb(y ~ 1, data = df)
    # quasipoisson (to tackle down convergence problems)
    modqglm <- glm(y ~ x, data = df, family = 'quasipoisson')
    modqglm.null <-  glm(y ~ 1, data = df, family = 'quasipoisson')
        
    # ------------- 
    # Test of effects
    # check convergence
    if(!is.null(modglm[['th.warn']]) | !is.null(modglm.null[['th.warn']])){
      p_glm_lr <- 'convergence error'
      p_glm_lrpb <- 'convergence error'
    } else {
      p_glm_lr <- anova(modglm, modglm.null, test = 'Chisq')[2 , 'Pr(Chi)']
      # Parametric bootstrap for GLM LR
      glm_pb <- myPBmodcomp(modglm, modglm.null, data = df, npb = npb)
      p_glm_lrpb <- glm_pb$p.pb
    }
    # F Tests
    p_lm_f <- anova(modlm, modlm.null, test = 'F')[2, 'Pr(>F)']
    p_qglm_f <- anova(modqglm, modqglm.null, test = 'F')[2, 'Pr(>F)']
    # non-parametric test
    p_k <- kruskal.test(y ~ x, data = df)$p.value
    
    # ----------------
    # LOEC
    # lm
    mc_lm <- summary(glht(modlm, linfct = mcp(x = 'Dunnett'),  
                          alternative = 'less'))$test$pvalues
    suppressWarnings( # intended warnings about no min -> no LOEC
      loeclm <- min(which(mc_lm < 0.05))
      )
    # negbin
    if(!is.null(modglm[['th.warn']])){
      loecglm <- 'convergence error'
    } else {
      mc_glm <- summary(glht(modglm, linfct = mcp(x = 'Dunnett'),  
                             alternative = 'less'))$test$pvalues
      suppressWarnings(
        loecglm <- min(which(mc_glm  < 0.05))
      )
    }
    # quasi
    mc_qglm <- summary(glht(modqglm, linfct = mcp(x = 'Dunnett'),  
                          alternative = 'less'))$test$pvalues
    suppressWarnings( # intended warnings about no min -> no LOEC
      loecqglm <- min(which(mc_qglm < 0.05))
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
                p_glm_lrpb = p_glm_lrpb, p_k = p_k, 
                loeclm = loeclm, loecglm = loecglm, loecqglm = loecqglm, 
                loecpw = loecpw
                )
           )
  }
  # run on simulated data
  res <- apply(z$y, 2, ana, x = z$x, npb = npb)
  return(res)
}

# Power
p_glob1 <- function(z){ 
  # extract p-values
  take <- c('p_lm_f', 'p_glm_lr','p_qglm_f', 'p_glm_lrpb', 'p_k')
  ps <- ldply(z, function(w) as.numeric(unlist(w[take])))
  names(ps) <- take
  ps <- melt(ps)
  out <- ddply(ps, .(variable), summarize,
        power = sum(value < 0.05, na.rm = TRUE) / sum(!is.na(value)),
        conv = sum(!is.na(value)) / length(value))
  return(out)
}


# loec
p_loec1 <- function(u, type = NULL){
  # extract p-values
  take <- c("loeclm", "loecglm", "loecqglm", "loecpw")
  loecs <- ldply(u, function(w) as.numeric(unlist(w[take])))
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
#' Function to create simulated binomial data
#' 
#' @description 
#' Simulate data from binomial distribution.
#' 
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

#' Function to analyse simulated datasets
#' 
#' @description Runs a normal model, a logistic model and a non-parametric tests on the simulated data
#' 
#' @param z simulated data, generated by dosim2()
#' @param verbose print status on the console?
#' @param asin Type of arcsine transformation. 'ecotox' or 'asin'. If 'ecotox' a special asin transformation is performed
#' 
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
      y_asin <- ifelse(y  == 0, asin(sqrt(1 / (4 * n_animals))),
                     ifelse((y / n_animals) == 1, asin(1) - asin(sqrt(1 / (4 * n_animals))),
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
#     # quasibinomial
#     modqglm <- glm(cbind(y, n_animals - y) ~ x, data = df, 
#                   family = quasibinomial(link = 'logit'))
#     modqglm.null <- glm(cbind(y, n_animals - y) ~ 1, data = df, 
#                   family = quasibinomial(link = 'logit'))

    # ------------- 
    # Tests
    # LR Tests
#     lm_lr <- lrtest(modlm, modlm.null)[2, 'Pr(>Chisq)']
    glm_lr <- lrtest(modglm, modglm.null)[2, 'Pr(>Chisq)']
    # F Tests
    lm_f <- anova(modlm, modlm.null, test = 'F')[2, 'Pr(>F)']
#     qglm_f <- anova(modqglm, modqglm.null, test = 'F')[2, 'Pr(>F)']
    # non-parametric test
    pk <- kruskal.test(y ~ x, data = df)$p.value
    
    # ------------------------------------------------------------
    # LOECs
    # multiple comparisons using Dunnett-contrasts
    pmclm <- p.adjust(coef(summary(modlm))[2:6 , 'Pr(>|t|)'], method = 'holm')
    pmcglm <- p.adjust(coef(summary(modglm))[2:6 , 'Pr(>|z|)'], method = 'holm')
#     pmcqglm <-  p.adjust(coef(summary(modqglm))[2:6 , 'Pr(>|t|)'], method = 'holm')
    # pairwise wilcox
    suppressWarnings(
      pw <- pairwise_wilcox(y, x, padj = 'holm', dunnett = TRUE))

    # extract LOEC (which level? 0 = Control)
    suppressWarnings(
      loeclm <- min(which(pmclm < 0.05)))
    suppressWarnings(
      loecglm <- min(which(pmcglm < 0.05)))
#     suppressWarnings(
#       loecqglm <- min(which(pmcqglm < 0.05)))
    suppressWarnings(
      loecpw <- min(which(pw < 0.05)))
    
    # --------------------------------------------
    # return object
    return(list(#lm_lr = lm_lr, 
                glm_lr = glm_lr, lm_f = lm_f, # qglm_f = qglm_f,
      pk = pk, 
      loeclm = loeclm, loecglm = loecglm, #loecqglm = loecqglm, 
      loecpw = loecpw
    ))
    
  }
  # run models on simulated data
  res <- apply(z$y, 2, ana, x = z$x, n_animals = z$n_animals, asin = asin)
  res
}


#### -----------------------------
### Extractor functions
p_glob <- function(z){ 
  # extract p-values
  ps <- ldply(z, function(w) unlist(w)[1:3])
  pow <- apply(ps, 2, function(z) sum(z < 0.05, na.rm = TRUE)) / length(z)
  return(pow)
}

# loec
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

# custom theme
mytheme <- theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size=14),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position="bottom",
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text= element_text(size=14, face = 'bold'))
