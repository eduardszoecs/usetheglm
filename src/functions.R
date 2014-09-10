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
#' Function to create simulated data
dosim <- function(N, mu, theta, nsims = 100){
  Nj     <- rep(N, time = 6)                # number of groups
  mus    <- rep(mu, times = Nj)             # vector of mus
  thetas <- rep(theta, times=Nj)            # vector of thetas
  x      <- factor(rep(1:6, times=Nj))      # factor
  y      <- replicate(nsims, rnegbin(sum(Nj), mus, thetas))
  return(list(x = x, y = y))
}


### ----------------------------------------------------------------------------
#' Function to analyse simulated datasets
resfoo <- function(z, verbose = TRUE){
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
    modlm <- mle2(yt ~ dnorm(mean = mu, sd = s),
                  parameters = list(mu ~ x, s ~ 1),
                  data = df,
                  start = list(mu = mean(yt), s = sd(yt)),
                  control = list(maxit = 500)
                  )
    modlm.null <- update(modlm, 
           parameters = list(mu ~ 1, sd ~ 1),
           start = list(mu = mean(yt), s = sd(yt)))
    
    # Negative binomial model
    modglm <- mle2(y ~ dnbinom(mu = exp(logmu), size = k),
                      parameters = list(logmu ~ x),
                      data = df,
                      start = list(logmu = log(mean(df$y)), k = 1)
                   )
    modglm.null <- update(modglm,
             parameters = list(logmu ~ 1),
             start = list(logmu = log(mean(df$y)), k = 1))
        
    # ---------------------------------------------------------
    # Likelihood-Ratio-Tests (global)
    plm <- anova(modlm, modlm.null)[2, 'Pr(>Chisq)']
    pglm <- anova(modglm, modglm.null)[2, 'Pr(>Chisq)']
    # non-parametric test
    pk <- kruskal.test(y ~ x, data = df)$p.value
    
    
    # ------------------------------------------------------------
    # LOECs
    # multiple comparisons using Dunnett-contrasts
    # no need for multcomp due to parametrisation
    # score tests
    pmclm <- p.adjust(coef(summary(modlm))[2:6 , 'Pr(z)'], method = 'holm')
    pmcglm <- p.adjust(coef(summary(modglm))[2:6 , 'Pr(z)'], method = 'holm')

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
      modlm = modlm, modlm.null = modlm.null, modglm = modglm, 
      modglm.null = modglm.null,
      plm = plm, pglm = pglm, pk = pk,   
      loeclm = loeclm, loecglm = loecglm, loecpw = loecpw
    ))
  }
  # run on simulated data
  res <- apply(z$y, 2, ana, x = z$x)
  res
}
