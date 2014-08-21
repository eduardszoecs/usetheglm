if(!exists('ld')){
  source("/home/edisz/Downloads/donotlog/src/0-load.R")
}

#####--------------------------------------------------------------------------
#### Simulate data
nsims <- 100
# sample sizes
N <- c(3,6,9,12)
# mu control
ctrl <- exp(log(10) * seq(log10(1), log10(1000), by=0.25))
# both as grid
todo <- expand.grid(N = N, ctrl = ctrl)

# function to create simulated data
dosim <- function(N, mu, nsims = 100){
  Nj     <- rep(N, 6)                       # number of groups
  theta  <- rep(3.91, 6)                    # theta in groups
  mus    <- rep(mu, times = Nj)             # vector of mus
  thetas <- rep(theta, times=Nj)            # vector of thetas
  x      <- factor(rep(1:6, times=Nj))      # factor
  y      <- replicate(nsims, rnegbin(sum(Nj), mus, thetas))
  return(list(x = x, y = y))
}

# simulate data
sims <- NULL
set.seed(1234)
for(i in seq_len(nrow(todo))){
  N <- todo[i, 'N']
  takectrl <- todo[i, 'ctrl']
  taketrt <- takectrl * 0.5
  mu <- c(rep(takectrl, each = 2), rep(taketrt, each = 4))
  sims[[i]] <- dosim(N = N, mu = mu, nsims = nsims)
}

# plot one realisation of simulated data
todo[26, ]
df <- data.frame(x = sims[[26]]$x, y = sims[[26]]$y[ , 25])
df$yt <- log(1/min(df$y[df$y!=0]) * df$y + 1)
dfm <- melt(df)
levels(dfm$variable) <- c('y', 'ln(Ay + 1)')
p1 <- ggplot(dfm, aes(x = x, y = value)) +
  geom_boxplot(fill = 'grey80') +
  facet_wrap(~variable, scales = 'free_y') +
  scale_x_discrete(labels = c('C', 'T1', 'T2', 'T3', 'T4', 'T5')) +
  labs(x = 'Treatment', y = 'Abundance') +
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        text = element_text(size=14),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
p1
# ggsave(file.path(figdir, 'p1.pdf'), p1, width = 10, height = 6)


#####--------------------------------------------------------------------------
##### analyse simulations
res <- llply(sims, function(z){
  ana <- function(y, x){
#     print(y)
    A <- 1/min(y[y!=0])         # ln(ax + 1) transformation
    yt <- log(A*y + 1)
    
    df <- data.frame(x, y, yt)
    
    # models
    modlm <- glm(yt ~ x, data = df)
    modglm <- glm.nb(y ~ x, data = df)
    
    # global test
    plm <- anova(modlm, test = 'Chisq')["x", "Pr(>Chi)"]
    pglm <- anova(modglm, test = 'Chisq')["x", "Pr(>Chi)"]
    pk <- kruskal.test(y ~ x, data = df)$p.value
    
    # multiple comparisons using Dunnett-contrasts
    smclm <- summary(glht(modlm, linfct = mcp(x = 'Dunnett')), test = adjusted('holm'))
    smcglm <- summary(glht(modglm, linfct = mcp(x = 'Dunnett')), test = adjusted('holm'))
#     # coin
#     require(coin)
#     owt <- oneway_test(y ~ x, data = df, 
#                 xtrafo = function(data) {
#                   trafo(data, 
#                         factor_trafo = function(x){
#                           model.matrix(~x - 1) %*% t(contrMat(table(x), "Dunnett"))
#                         }
#                   )},
#                 ytrafo = function(data) trafo(data, numeric_trafo = rank)
#                 )
#     p.adjust(pvalue(owt,  method = "unadjusted"), 'holm')
#     # nparcomp
#     require(nparcomp)
#     np <- nparcomp(y ~ x, data = df, type = 'Dunnett', alternative = 'two.sided', asy.method = 'mult.t', info = FALSE)
#     np$Analysis[ , "p.Value"]
    # pairwise wilcox
    pw <- pairwise_wilcox(y, x, padj = 'holm', dunnett = TRUE)

    
    # LOEC (which level? 0 = Control)
    loeclm <- min(which(smclm$test$pvalues < 0.05))
    loecglm <- min(which(smcglm$test$pvalues < 0.05))
    loecpw <- min(which(pw < 0.05))
    return(list(A = A, 
                # modlm = modlm, modglm=modglm,
                plm=plm, pglm=pglm, pk = pk,
                # smclm=smclm, smcglm=smcglm, 
                loeclm=loeclm, loecglm=loecglm, loecpw = loecpw
                ))
  }
  # run models on simulated data
  res <- apply(z$y, 2, ana, x = z$x)
  res
}, .progress = 'text')

saveRDS(res, file.path(cachedir, 'res.rds'))
# res <- readRDS(file.path(cachedir, 'res.rds'))


#####--------------------------------------------------------------------------
##### compare methods
# global power
pow <- function(z){
  ps <- ldply(z, function(w) c(lm = w$plm, glm = w$pglm, pk = w$pk))
  apply(ps, 2, function(z) sum(z < 0.05)) / length(z)
}
pows <- ldply(res, pow)
pows$muc <- todo$ctrl
pows$N <- todo$N
powsm <- melt(pows, id.vars = c('muc', 'N'))

p2 <- ggplot(powsm) +
  geom_line(aes(y = value, x = muc, group = variable)) +
  geom_point(aes(y = value, x = muc, fill = variable), size = 4, pch = 21, color = 'black') +
  coord_trans(xtrans = 'log10') +
  scale_x_continuous(breaks = round(ctrl, 0)) +
  facet_wrap(~N) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Rejection rate of ', H[0], ' (global test , ', alpha, ' = 0.05)'))) +
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
                  breaks = c('lm', 'glm', 'pk'), 
                  labels = c('LM + log(Ay+1)', 'GLM (neg. bin.)', 'Kruskal'),
                  start = 0.1, end = 0.9) +
  theme(legend.position="bottom", legend.key = element_blank())
p2
ggsave(file.path(figdir, 'p2.pdf'), p2, width = 11, height = 11)

# loec
loec <- function(z){
  loecs <- ldply(z, function(w) c(lm = w$loeclm, glm = w$loecglm))
  out <- apply(loecs, 2, function(x) sum(x == 2))
  return(out)
}
loecs <- ldply(res, loec)
loecs$muc <- todo$ctrl
loecs$N <- todo$N
loecsm <- melt(loecs, id.vars = c('muc', 'N'))
loecsm$value <- loecsm$value / nsims

ggplot(loecsm) +
  geom_line(aes(y = value, x = muc, col = variable)) +
  coord_trans(xtrans = 'log10') +
  scale_x_continuous(breaks = round(unique(todo$ctrl), 0)) +
  facet_wrap(~N) + 
  theme_bw() +
  labs(x = 'mean Abundance in Control', y = 'Proportion correct LOEC')


