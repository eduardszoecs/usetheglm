if(!exists('ld')){
  source("/home/edisz/Downloads/donotlog/src/0-load.R")
}

#####--------------------------------------------------------------------------
### 1 - Counts
# Example data from Brock et al. 2014
df <- read.table(header = TRUE, text = '
                   C  0_1	0_3	1	3	10
175	29	27	36	26	20
65	114	78	11	13	37
154	72	27	105	33	NA
83 NA NA NA NA NA					
')

dfm <- melt(df)
names(dfm) <- c('x', 'y')
dfm$x <- factor(as.numeric(as.character(mapvalues(dfm$x, c('C', 'X0_1', 'X0_3', 'X1', 'X3', 'X10'), c(0, 0.1, 0.3, 1,3,10)))))
A <- min(dfm$y[dfm$y > 0], na.rm = TRUE)
dfm$yt <- log(A*dfm$y + 1)
dfm <- melt(dfm, value.name = 'y')
levels(dfm$variable) <- c('y', 'ln(Ay+1)')

p4 <- ggplot() + 
  geom_boxplot(data = dfm, aes(x = x, y = y), fill = 'grey80') +
  facet_wrap(~variable, scales = 'free_y') +
  scale_x_discrete(labels = c('C', 'T1', 'T2', 'T3', 'T4', 'T5')) +
  labs(x = 'Treatment', y = 'Abundance') +
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        text = element_text(size=14),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
p4
# ggsave(file.path(figdir, 'p4.pdf'), p4, width = 10, height = 6)

# combine p1 and p4
p4a <- p4 + ggtitle('Real data') + labs(x = '') 
p1a <- p1 + ggtitle('Simulated data')
p4_1 <- arrangeGrob(p4a, p1a)
p4_1
ggsave(file.path(figdir, 'p4_1.pdf'), p4_1, width = 12, height = 12)


### Split data
dfmt <- dfm[dfm$variable == "ln(Ay+1)", ]
dfmo <- dfm[dfm$variable == 'y', ]


### Fit different Models to the data
## lm + transformation
modlm <- glm(y ~ x, data = dfmt)
modlm2 <- lm(y ~ x, data = dfmo)
# global p-value
drop1(modlm, test = 'Chisq')
# Dunnett
summary(glht(modlm, linfct = mcp(x = 'Dunnett')), test = adjusted('holm'))


modnb <- glm.nb(y ~ x, data = dfmo)
summary(modnb)

# Mean-VarPlot
mv <- ddply(dfmo, .(x), summarise, 
            mean = mean(y, na.rm = TRUE), 
            var = var(y, na.rm = TRUE), 
            n = sum(!is.na(y)))
plot(var ~ mean, data = mv)
# Poission
abline(0,1, col = 'red', lwd = 2)
# NB
lines(30:110, 30:110 + (30:110)^2 / modnb$theta, col = 'steelblue', lwd = 3)

# plot( fitted(modnb), residuals(modnb, type="pearson"), pch = 16)
# var(residuals(modnb, type="pearson"))

coefplot2(list(modlm, modnb))
legend('topright', legend = c('lm', 'glm'), cex = 0.7, pch = 16, col = 1:2)



# #### Simulate base on properties of data
# prop <- ddply(dfmo[!is.na(dfmo$y), ], .(x), summarise, 
#       mu = fitdistr(y, densfun = 'negative binomial')$estimate[2],
#       theta = fitdistr(y, densfun = 'negative binomial')$estimate[1],
#       var = var(y),
#       n = length(y),
#       mulog = mean(log(y)),
#       varlog = var(log(y))
# )
# 
# # gamma for lognormal
# # modlmgamma <- glm(log(y) ~ x, data = dfmo, family = Gamma)
# # ss <- summary(modlmgamma, dispersion = 1)
# # sqrt(ss$deviance / ss$df.residual)
# 
# nsims  <- 100                          # number of simulations
# 
# dosim <- function() {                   # function to simulate data
#   Nj     <- prop$n                          # group sizes for 6 groups
#   mu     <- prop$mu                         # expected values in groups
#   theta  <- prop$theta                       # theta in groups
#   var    <- prop$var
#   mus    <- rep(mu, times = Nj)             # vector of mus
#   thetas <- rep(theta, times=Nj)            # vector of thetas
#   x      <- factor(rep(1:6, times=Nj))       # factor
#   # from negbin
#   y      <- replicate(nsims, rnegbin(sum(Nj), mus, thetas))
#   yg      <- replicate(nsims, rnegbin(sum(Nj), mus, theta = 3.91))
# #   # from lognormal
# #   yn      <- replicate(nsims, ceiling(rlnorm(sum(Nj), prop$mulog, sd = sqrt(prop$varlog))))
# #   yng     <- replicate(nsims, ceiling(rlnorm(sum(Nj), prop$mulog, sd = 0.6605833)))
# #   yout <- data.frame(rbind(y, yg, yn, yng))
# #   yout$methods <- rep(c('y', 'yg', 'yn', 'yng'), each = sum(Nj))
# 
#   yout <- data.frame(rbind(y, yg))
#   yout$methods <- rep(c('y', 'yg'), each = sum(Nj))
#   yout <- melt(yout, id.vars = 'methods')
#   return(yout)
# }
# sims <- dosim()
# 
# plot(x, sims$value[sims$methods == 'y' & sims$variable == 'X1'])
# 
# x      <- factor(rep(1:6, times=Nj))       # factor
# ana <- function(df){
#   y <- df$value
#   #     print(y)
#   A <- 1/min(y[y!=0])         # ln(ax + 1) transformation
#   yt <- log(A*y + 1)
# #   yt <- log(y)
#   
#   df <- data.frame(x, y, yt)
#   
#   # models
#   modlm <- glm(yt ~ x, data = df)
#   modglm <- glm.nb(y ~ x, data = df)
#   
#   # global test
#   plm <- drop1(modlm, test = 'Chisq')["x", "Pr(>Chi)"]
#   pglm <- drop1(modglm, test = 'Chisq')["x", "Pr(>Chi)"]
#   pk <- kruskal.test(y ~ x, data = df)$p.value
#   
#   # multiple comparisons using Dunnett-contrasts
#   smclm <- summary(glht(modlm, linfct = mcp(x = 'Dunnett')), test = adjusted('holm'))
#   smcglm <- summary(glht(modglm, linfct = mcp(x = 'Dunnett')), test = adjusted('holm'))
#   #     # coin
#   #     require(coin)
#   #     owt <- oneway_test(y ~ x, data = df, 
#   #                 xtrafo = function(data) {
#   #                   trafo(data, 
#   #                         factor_trafo = function(x){
#   #                           model.matrix(~x - 1) %*% t(contrMat(table(x), "Dunnett"))
#   #                         }
#   #                   )},
#   #                 ytrafo = function(data) trafo(data, numeric_trafo = rank)
#   #                 )
#   #     p.adjust(pvalue(owt,  method = "unadjusted"), 'holm')
#   # nparcomp
#   np <- nparcomp(y ~ x, data = df, type = 'Dunnett', alternative = 'two.sided', asy.method = 'mult.t', info = FALSE)$Analysis[ , "p.Value"]
#   # Holm needed?
#   #     np <- p.adjust(np, method = 'holm')
#   # pairwise wilcox
#   pw <- pairwise_wilcox(y, x, padj = 'holm', dunnett = TRUE)
#   
#   
#   # LOEC (which level? 0 = Control)
#   loeclm <- min(which(smclm$test$pvalues < 0.05))
#   loecglm <- min(which(smcglm$test$pvalues < 0.05))
#   loecpw <- min(which(pw < 0.05))
#   loecnp <- min(which(np < 0.05))
#   return(list(A = A, 
#               # modlm = modlm, modglm=modglm,
#               plm=plm, pglm=pglm, pk = pk,
#               # smclm=smclm, smcglm=smcglm, 
#               loeclm=loeclm, loecglm=loecglm, loecpw = loecpw, loecnp = loecnp,
#               method = unique(df$methods)
#   ))
# }
# 
# res <- dlply(sims, .(methods, variable), ana, .progress = 'text')
# 
# 
# # power
# pow <- function(z){
#   ps <- ldply(z, function(w) c(lm = w$plm, glm = w$pglm, pk = w$pk))
#   apply(ps[ , -1], 2, function(z) sum(z < 0.05)) / length(z)
# }
# 
# pow(res[1:100])
# pow(res[101:200])
# 
# # load
# loec <- function(z){
#   loecs <- ldply(z, function(w) c(lm = w$loeclm, glm = w$loecglm, pw = w$loecpw, np = w$loecnp))
#   out <- apply(loecs[ ,-1], 2, function(x) sum(x == 2))
#   return(out)
# }
# loec(res[1:100])
# loec(res[101:200])



#####--------------------------------------------------------------------------
# Proportions
# Data from Newman
df <- read.table(header = TRUE, text = 'conc A B C D
0 1 1 0.9 0.9
32 0.8 0.8 1 0.8
64 0.9 1 1 1 
128 0.9 0.9 0.8 1
256 0.7 0.9 1 0.5
512 0.4 0.3 0.4 0.2')
df

require(reshape2)
# to long
dfm <- melt(df, id.vars = 'conc', value.name = 'y', variable.name = 'tank')
# conc as factor
dfm$conc <- factor(dfm$conc)
head(dfm)
