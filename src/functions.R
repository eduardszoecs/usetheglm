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


### ------------------------------------------------------------------------------
#' Function to analyse simulations
resfoo <- function(z){
  ana <- function(y, x){
    # ln(ax + 1) transformation
    A <- 1/min(y[y!=0])         
    yt <- log(A*y + 1)
    df <- data.frame(x, y, yt)
    
    # model using mle2
    modlm <- glm(yt ~ x, data = df)
    modlm.null <- glm(yt ~ 1, data = df)
    modglm <- glm.nb(y ~ x, data = df
                     # , control=glm.control(maxit=10,trace=3)
    )
    #     modglm.null <- glm.nb(y ~ 1, data = df)
    # using mle2
    # dput(y)
    modglm_bb <- mle2(y ~ dnbinom(mu = exp(logmu), size = k),
                      parameters = list(logmu ~ x),
                      data = df,
                      start = list(logmu = 1, k = 1),
                      lower = list(logmu = -Inf, k = 0),
                      upper = list(logmu = Inf, k = 1000),
                      method = 'L-BFGS-B')
    
    modglm_bb.null <- try(mle2(y ~ dnbinom(mu = exp(logmu), size = k),
                               parameters = list(logmu ~ 1),
                               data = df,
                               start = list(logmu = 1, k = 1),
                               lower = list(logmu = -Inf, k = 0),
                               upper = list(logmu = Inf, k = 1000), 
                               method = 'L-BFGS-B'),
                          silent = TRUE)
    if(is(modglm_bb.null, "try-error")){
      modglm_bb.null <- mle2(y ~ dnbinom(mu = exp(logmu), size = k),
                             parameters = list(logmu ~ 1),
                             data = df,
                             start = list(logmu = 0, k = 1),
                             lower = list(logmu = -Inf, k = 0),
                             upper = list(logmu = Inf, k = 1000), 
                             method = 'L-BFGS-B')
    }
    
    
    
    # global test
    plm <- lrtest(modlm, modlm.null)[2, 'Pr(>Chisq)']
    # same as
    # plm <- drop1(modlm, test = 'Chisq')["x", "Pr(>Chi)"]
    # pglm <- lrtest(modglm, modglm.null)[2, 'Pr(>Chisq)']
    # not the same as
    # drop1(modglm, test = 'Chisq')["x", "Pr(>Chi)"]
    # as no correct method of drop1 for glm.nb (theta is not refitted correctly!)
    pglm <- anova(modglm_bb, modglm_bb.null)[2, 'Pr(>Chisq)']
    pk <- kruskal.test(y ~ x, data = df)$p.value
    
    # multiple comparisons using Dunnett-contrasts
    smclm <- summary(glht(modlm, linfct = mcp(x = 'Dunnett')), test = adjusted('holm'))
    smcglm <- summary(glht(modglm, linfct = mcp(x = 'Dunnett')), test = adjusted('holm'))
    # nparcomp
    # np <- nparcomp(y ~ x, data = df, type = 'Dunnett', alternative = 'two.sided', asy.method = 'mult.t', info = FALSE)$Analysis[ , "p.Value"]
    # Holm needed?
    #     np <- p.adjust(np, method = 'holm')
    # pairwise wilcox
    suppressWarnings(
      pw <- pairwise_wilcox(y, x, padj = 'holm', dunnett = TRUE)
    )
    
    
    # LOEC (which level? 0 = Control)
    suppressWarnings(
      loeclm <- min(which(smclm$test$pvalues < 0.05))
    )
    suppressWarnings(
      loecglm <- min(which(smcglm$test$pvalues < 0.05))
    )
    suppressWarnings(
      loecpw <- min(which(pw < 0.05))
    )
    # loecnp <- min(which(np < 0.05))
    return(list(A = A, 
                # modlm = modlm, modglm=modglm,
                plm=plm, pglm=pglm, pk = pk,
                # smclm=smclm, smcglm=smcglm, 
                loeclm=loeclm, loecglm=loecglm, loecpw = loecpw
                # , loecnp = loecnp
    ))
  }
  # run on simulated data
  res <- apply(z$y, 2, ana, x = z$x)
  res
}
