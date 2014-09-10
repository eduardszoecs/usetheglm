#####--------------------------------------------------------------------------
### Simulation 2 - Proportions

nsims <- 100
N <- 3:6
pEs <- seq(0.6, 0.9, 0.01)

todo <- expand.grid(N = N, pE = pEs)

# function to create simulated data
dosim <- function(N, pC = 0.9, pE = 0.3, nsim){
  nanimals <-  10
  n_group <- 6
  p = c(rep(rep(pC, N), 2), rep(rep(pE, N), 4))
  surv <- replicate(nsim, rbinom(N * n_group, size = nanimals, prob = p))
  x      <- factor(rep(1:6, each = N))      
  y      <- surv/nanimals
  return(list(x = x, y = y))
}

# simulate data
sims <- NULL
set.seed(1234)
for(i in seq_len(nrow(todo))){
  sims[[i]] <- dosim(N = todo[i, 'N'], pE = todo[i, 'pE'], nsim = nsims)
}

# plot one realisation of simulated data
df <- data.frame(x = sims[[34]]$x, y = sims[[34]]$y[ , 2])
ggplot(df, aes(x = x, y = y)) +
  geom_boxplot(fill = 'grey80') +
  scale_x_discrete(labels = c('C', 'T1', 'T2', 'T3', 'T4', 'T5')) +
  labs(x = 'Treatment', y = 'Prop. surv.') +
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        text = element_text(size=14),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))



# analyse
resfoo <- function(z){
  ana <- function(y, x){
    # asin transformation
    y_asin <- ifelse(y == 0, asin(sqrt(1 / (4*10))),
                     ifelse(y == 1, asin(1) - asin(sqrt(1 / (4*10))),
                     asin(sqrt(y)))
    )
    
    df <- data.frame(x, y, y_asin)
    
    # models
    modlm <- glm(y ~ x, data = df)
    modlmasin <- glm(y_asin ~ x, data = df)
    modglm <- glm(y ~ x, data = df, family = binomial, weights = rep(10, nrow(df)))
    
    # global test
    plm <- drop1(modlm, test = 'Chisq')["x", "Pr(>Chi)"]
    plmasin <- drop1(modlmasin, test = 'Chisq')["x", "Pr(>Chi)"]
    pglm <- drop1(modglm, test = 'Chisq')["x", "Pr(>Chi)"]
    pnp <- kruskal.test(y ~ x, data = df)$p.value
    # same as:
#     pkasin <- kruskal.test(y_asin ~ x, data = df)$p.value
    
    # multiple comparisons using Dunnett-contrasts
    smclm <- summary(glht(modlm, linfct = mcp(x = 'Dunnett')), test = adjusted('holm'))
    smclmasin <- summary(glht(modlmasin, linfct = mcp(x = 'Dunnett')), test = adjusted('holm'))
    smcglm <- summary(glht(modglm, linfct = mcp(x = 'Dunnett')), test = adjusted('holm'))
    # pairwise wilcox
    smcnp <- pairwise_wilcox(y_asin, x, padj = 'holm', dunnett = TRUE)
    
    # LOEC (which level? 0 = Control)
    loeclm <- min(which(smclm$test$pvalues < 0.05))
    loeclmasin <- min(which(smclmasin$test$pvalues < 0.05))
    loecglm <- min(which(smcglm$test$pvalues < 0.05))
    loecnp <- min(which(smcnp < 0.05))
    return(list(plm=plm, plmasin = plmasin, pglm=pglm, pnp = pnp,
                loeclm=loeclm, loeclmasin = loeclmasin, loecglm=loecglm, loecnp = loecnp
    ))
  }
  # run models on simulated data
  res <- apply(z$y, 2, ana, x = z$x)
  res
}
res <- llply(sims, resfoo, .progress = 'text')


# global power
pow <- function(z){
  ps <- ldply(z, function(w) c(lm = w$plm, lmasin = w$plmasin, glm = w$pglm, np = w$pnp))
  apply(ps, 2, function(z) sum(z < 0.05)) / length(z)
}
pows <- ldply(res, pow)
pows$pE <- todo$pE
pows$N <- todo$N
powsm <- melt(pows, id.vars = c('pE', 'N'))

ggplot(powsm) +
  geom_line(aes(y = value, x = pE, group = variable)) +
  geom_point(aes(y = value, x = pE, fill = variable), size = 4, pch = 21, color = 'black') +
  coord_trans(xtrans = 'log10') +
  # use here
  scale_x_continuous(breaks = pEs) +
  facet_wrap(~N) + 
  # axes
  labs(x = expression('pE'), 
       y = expression(paste('Power (global test , ', alpha, ' = 0.05)'))) +
  # appearance
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size=14),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  # legend
  scale_fill_grey(name = '', 
                  breaks = c('lm', 'lmasin', 'glm', 'np'), 
                  labels = c('LM', 'LM + arcsin', 'GLM (bin.)', 'Kruskal'),
                  start = 0, end = 1) +
  theme(legend.position="bottom", legend.key = element_blank())
