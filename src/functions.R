### ----------------------------------------------------------------------------
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
resfoo1 <- function(z, verbose = TRUE){
  if(verbose){
    message('n: ', length(z$x) / 6, '; muc = ', mean(z$y[,1][z$x == 1]))
  }
  ana <- function(y, x){
    # ln(ax + 1) transformation
    A <- 1/min(y[y!=0])         
    yt <- log(A*y + 1)
    df <- data.frame(x, y, yt)
#     dput(df)

    # ------------------------------
    # Maximum-Likelihood estimation
    # Gaussian model
    modlm <- lm(yt ~ x, data = df)
    modlm.null <- lm(yt ~ 1, data = df)
    
#     # Negative binomial model
#     modglm <- mle2(y ~ dnbinom(mu = exp(logmu), size = k),
#                       parameters = list(logmu ~ x),
#                       data = df,
#                       start = list(logmu = log(mean(df$y)), k = 10)
#                    )
#     modglm.null <- mle2(y ~ dnbinom(mu = exp(logmu), size = k),
#                         data = df,
#                         parameters = list(logmu ~ 1),
#                         start = list(logmu = log(mean(df$y)), k = 10)
#                         )
    modglm <- glm.nb(y ~ x, data = df)
    modglm.null <- glm.nb(y ~ 1, data = df)
        
    # ---------------------------------------------------------
    # F test
    plmf <- anova(modlm, modlm.null)[2, 'Pr(>F)']
    # = waldtest(modlm, modlm.null)
    # LR test
    plmlr <- lrtest(modlm, modlm.null)[2, 'Pr(>Chisq)']
    
    pglmwc <- waldtest(modglm, modglm.null, test = 'Chisq')[2, 'Pr(>Chisq)']
    pglmwf <- waldtest(modglm, modglm.null, test = 'F')[2, 'Pr(>F)']
    pglmlr <- lrtest(modglm, modglm.null)[2, 'Pr(>Chisq)']
    # non-parametric test
    pk <- kruskal.test(y ~ x, data = df)$p.value
    
    
    # ------------------------------------------------------------
    # LOECs
    # multiple comparisons using Dunnett-contrasts
    # no need for multcomp due to parametrisation
    # score tests
    pmclm <- p.adjust(coef(summary(modlm))[2:6 , 'Pr(>|t|)'], method = 'holm')
    pmcglm <- p.adjust(coef(summary(modglm))[2:6 , 'Pr(>|z|)'], method = 'holm')

    # pairwise wilcox
    suppressWarnings( # ties
      pw <- pairwise_wilcox(y, x, padj = 'holm', dunnett = TRUE)
    )
    
    # extract LOEC (which level? 0 = Control)
    suppressWarnings( # intended warnings about no min -> no LOEC
      loeclm <- min(which(pmclm < 0.05))
    )
    suppressWarnings(
      loecglm <- min(which(pmcglm < 0.05))
    )
    suppressWarnings(
      loecpw <- min(which(pw < 0.05))
    )
    
    # --------------------------------------------
    # return object
    return(list(
      plmf = plmf, plmlr = plmlr, 
      pglmwc = pglmwc, pglmwf = pglmwf, pglmlr = pglmlr, 
      pk = pk, 
      loeclm = loeclm, loecglm = loecglm, loecpw = loecpw
    ))
  }
  # run on simulated data
  res <- apply(z$y, 2, ana, x = z$x)
  res
}

p_glob1 <- function(z){ 
  # extract p-values
  ps <- ldply(z, function(w) unlist(w)[1:6])
  # calculate power
  pow <- apply(ps, 2, function(z) sum(z < 0.05, na.rm = TRUE)) / length(z)
  return(pow)
}


# loec
p_loec1 <- function(z, type = NULL){
  # extract p-values
  loecs <- ldply(z, function(w) unlist(w)[7:9]) 
  if(type == 't1'){
    pow <- apply(loecs, 2, function(x) sum(x != Inf, na.rm = TRUE) / length(x))
  } 
  if(type == 'power'){
    pow <- apply(loecs, 2, function(x) sum(x == 2, na.rm = TRUE) / length(x))
  }
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
    # switch for transformations
    if(asin == 'ecotox'){
      y_asin <- ifelse(y  == 0, asin(sqrt(1 / (4 * n_animals))),
                     ifelse((y / n_animals) == 1, asin(1) - asin(sqrt(1 / (4 * n_animals))),
                            asin(sqrt(y / n_animals))))
    }
    if(asin == 'asin'){
      y_asin <- asin(sqrt(y / n_animals))
    }
    # arcsin transformation of proportions
    df <- data.frame(x, y, y_asin)
    
    # ------------------------------
    # Maximum-Likelihood estimation
    # Gaussian model
    modlm <- lm(y_asin ~ x, data = df)
    modlm.null <- lm(y_asin ~ 1, data = df)
    # binomial model
    modglm <- glm(cbind(y, n_animals - y) ~ x, data = df, 
                  family = binomial(link = 'logit'))
    modglm.null <- glm(cbind(y, n_animals - y) ~ 1, data = df, 
                       family = binomial(link = 'logit'))

    
    # ---------------------------------------------------------
    # Global Test
    # F-test for LM
    plm <- anova(modlm, modlm.null)[2, 'Pr(>F)']
    # LRT-test for GLM
    pglm <- anova(modglm, modglm.null, test = 'Chisq')[2, 'Pr(>Chi)']
    # non-parametric test
    pk <- kruskal.test(y ~ x, data = df)$p.value
    
    
    # ------------------------------------------------------------
    # LOECs
    # multiple comparisons using Dunnett-contrasts
    pmclm <- p.adjust(coef(summary(modlm))[2:6 , 'Pr(>|t|)'], method = 'holm')
    pmcglm <- p.adjust(coef(summary(modglm))[2:6 , 'Pr(>|z|)'], method = 'holm')
    
    # pairwise wilcox
    suppressWarnings(
      pw <- pairwise_wilcox(y, x, padj = 'holm', dunnett = TRUE)
    )
    
    # extract LOEC (which level? 0 = Control)
    suppressWarnings(
      loeclm <- min(which(pmclm < 0.05))
    )
    suppressWarnings(
      loecglm <- min(which(pmcglm < 0.05))
    )
    suppressWarnings(
      loecpw <- min(which(pw < 0.05))
    )
    
    # --------------------------------------------
    # return object
    return(list(
      plm = plm, pglm = pglm, pk = pk,   
      loeclm = loeclm, loecglm = loecglm, loecpw = loecpw
    ))
    
  }
  # run models on simulated data
  res <- apply(z$y, 2, ana, x = z$x, n_animals = z$n_animals, asin = asin)
  res
}


#### -----------------------------
### Extractor functions
# global ps
p_glob <- function(z){ 
  # extract p-values
  ps <- ldply(z, function(w) c(lm = w$plm, glm = w$pglm, pk = w$pk))
  pow <- apply(ps, 2, function(z) sum(z < 0.05, na.rm = TRUE)) / length(z)
  return(pow)
}

# loec
p_loec <- function(z, type = NULL){
  # extract p-values
  loecs <- ldply(z, function(w) c(lm = w$loeclm, glm = w$loecglm, pw = w$loecpw
  ))
  if(type == 't1'){
    pow <- apply(loecs, 2, function(x) sum(x != Inf, na.rm = TRUE) / length(x))
  } 
  if(type == 'power'){
    pow <- apply(loecs, 2, function(x) sum(x == 2, na.rm = TRUE) / length(x))
  }
  return(pow)
}


#### -----------------------------
### Misc functions
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
