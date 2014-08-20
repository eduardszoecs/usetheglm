#####--------------------------------------------------------------------------
### Setup project structure
rm(list =ls()) 

## Project Path
## You have to change this!
prj <- "/home/edisz/Downloads/donotlog/"

## Subfolder paths
srcdir <- file.path(prj, "src")     # source code
datadir <- file.path(prj, "data")   # data
figdir <- file.path(prj, "report/fig") # figures for latex


#####--------------------------------------------------------------------------
### install missing if needed!
#####

require(reshape2)
require(ggplot2)
require(plyr)
require(coefplot2)
require(multcomp)
require(MASS)

### To transform or to GLM?
# Example data from Brock et al. 2014
df <- read.table(header = TRUE, text = '
                   C	0_1	0_3	1	3	10
175	29	27	36	26	20
65	114	78	11	13	37
154	72	27	105	33	NA
83 NA NA NA NA NA					
')

df_m <- melt(df)
df_m$variable <- as.numeric(as.character(mapvalues(df_m$variable, c('C', 'X0_1', 'X0_3', 'X1', 'X3', 'X10'), c(0, 0.1, 0.3, 1,3,10))))
df_m$value_t <- log(2*df_m$value + 1)
df_m$variable <- factor(df_m$variable)

ggplot() + 
  geom_boxplot(data = df_m, aes(x = variable, y = value))

ggplot() + 
  geom_boxplot(data = df_m, aes(x = variable, y = value_t))


# LM
mod1 <- lm(value_t ~ variable, data = df_m)

coefplot2(mod1)
summary(mod1)
anova(mod1)
par(mfrow = c(2,2))
plot(mod1)

summary(glht(mod1, linfct = mcp(variable = 'Dunnett'), alternative = 'less'))
summary(glht(mod1, linfct = mcp(variable = 'Dunnett')))


# GLM
mod1 <- glm(value_t ~ variable, data = df_m, family = 'gaussian')
mod2a <- glm(value ~ variable, data = df_m, family = 'poisson')
summary(mod2a) # indication of overdisp
par(mfrow = c(2,2))
plot(mod2a)
mod2 <- glm.nb(value ~ variable, data = df_m)
plot(mod2)
coefplot2(mod2)
anova(mod2)
summary(mod2)
summary(glht(mod2, linfct = mcp(variable = 'Dunnett'), alternative = 'less'))
summary(glht(mod2, linfct = mcp(variable = 'Dunnett')))

par(mfrow = c(1,1))
coefplot2(list(mod1, mod2, mod2a))
AIC(mod1, mod2, mod2a)
logLik(mod1); logLik(mod2); logLik(mod2a)


### Mean-VarPlot
mv <- ddply(df_m, .(variable), summarise, mean = mean(value, na.rm = TRUE), var = var(value, na.rm = TRUE), n = sum(!is.na(value)))
plot(mv$mean, mv$var)
abline(0,1)
lines(30:120, 30:120 + (30:120)^2 / mod2$theta, pch = 16)

# These residuals should have constant variance = 1 if the mean variance relationship is correctly specified.
plot( fitted(mod2a), residuals(mod2a, type="pearson"))
var(residuals(mod2a, type="pearson"))
var(residuals(mod2, type="pearson"))

# properties per group
dlply(df_m, .(variable), function(x) fitdistr(x$value[!is.na(x$value)], densfun = 'negative binomial'))
ddply(df_m, .(variable), summarise, 
      mu = fitdistr(value[!is.na(value)], densfun = 'negative binomial')$estimate[2],
      theta = fitdistr(value[!is.na(value)], densfun = 'negative binomial')$estimate[1],
      var = var(value[!is.na(value)])
      )


### simulate data with properties from example above
Nj     <- c(4, 3, 3, 3, 3, 2)                            # group sizes for 6 groups
mu     <- c(119.25, 71.6667, 44, 50.6667, 24, 28.5)   # expected values in groups
theta  <- c(6.57, 3.82, 4.15, 1.50, 10.91, 17.8)                           # theta in groups
mus    <- rep(mu, times = Nj)             # vector of mus
thetas <- rep(theta, times=Nj)            # vector of thetas
x      <- factor(rep(1:6, times=Nj))       # factor
nsims  <- 100                          # number of simulations

dosim <- function() {                   # function to simulate data
  rnegbin(sum(Nj), mus, thetas)     
}
ysim <- replicate(nsims, dosim()) 
plot(x, ysim[,1])


ana <- function(y){
  A <- 1/min(y[y!=0])         # ln(ax + 1) transformation
  yt <- log(A*y + 1)
  # models
  modlm <- glm(yt ~ x)
  modglm <- glm.nb(y ~ x)
  # global test
  plm <- anova(modlm, test = 'Chisq')["x", "Pr(>Chi)"]
  pglm <- anova(modglm, test = 'Chisq')["x", "Pr(>Chi)"]
  # multiple comparisons
  smclm <- summary(glht(modlm, linfct = mcp(x = 'Williams')))
  smcglm <- summary(glht(modglm, linfct = mcp(x = 'Williams')))
  # NOEC (which level? 1 = Control)
  noeclm <- sum(!smclm$test$pvalues < 0.05) + 1
  noecglm <- sum(!smcglm$test$pvalues < 0.05) + 1
  return(list(A = A, modlm = modlm, modglm=modglm, plm=plm, pglm=pglm, 
              smclm=smclm, smcglm=smcglm, noeclm=noeclm, noecglm=noecglm))
}
# run models on simulated data
res <- apply(ysim, 2, ana)

# compare power
ps <- ldply(res, function(z) c(lm = z$plm, glm = z$pglm))
apply(ps, 2, function(z) sum(z < 0.05))

# tabulate NOECS
noecs <- ldply(res, function(z) c(lm = z$noeclm, glm = z$noecglm))
apply(noecs, 2, function(z) table(z) / 100)



#### Simulate data, change properties
# change sample sizes
N <- 2:12
nsims  <- 100                          # number of simulations
dosim <- function(N, mu = c(120, 120, 60, 60, 60, 60)){
  Nj     <- rep(N, 6)
  mu     <- mu   # expected values in groups
  theta  <- rep(3.91, 6)                # theta in groups = theta from mod2
  mus    <- rep(mu, times = Nj)             # vector of mus
  thetas <- rep(theta, times=Nj)            # vector of thetas
  x      <- factor(rep(1:6, times=Nj))       # factor
  nsims   <- 100
  y      <- replicate(nsims, rnegbin(sum(Nj), mus, thetas))
  return(list(x = x, y = y, N = N))
}
sims <- lapply(N, dosim)
plot(sims[[1]]$x, sims[[1]]$y[,4])

# analyse simulations
res <- lapply(sims, function(z){
  ana <- function(y, x){
    A <- 1/min(y[y!=0])         # ln(ax + 1) transformation
    yt <- log(A*y + 1)
    # models
    modlm <- glm(yt ~ x)
    modglm <- glm.nb(y ~ x)
    # global test
    plm <- anova(modlm, test = 'Chisq')["x", "Pr(>Chi)"]
    pglm <- anova(modglm, test = 'Chisq')["x", "Pr(>Chi)"]
    # multiple comparisons
    smclm <- summary(glht(modlm, linfct = mcp(x = 'Williams')))
    smcglm <- summary(glht(modglm, linfct = mcp(x = 'Williams')))
    # NOEC (which level? 1 = Control)
    noeclm <- sum(!smclm$test$pvalues < 0.05) + 1
    noecglm <- sum(!smcglm$test$pvalues < 0.05) + 1
    return(list(A = A, modlm = modlm, modglm=modglm, plm=plm, pglm=pglm, 
                smclm=smclm, smcglm=smcglm, noeclm=noeclm, noecglm=noecglm))
  }
  # run models on simulated data
  res <- apply(z$y, 2, ana, x = z$x)
  return(list(res = res, N = z$N))
})

# extract information
# Power
pow <- function(z){
  ps <- ldply(z$res, function(w) c(lm = w$plm, glm = w$pglm))
  c(N = z$N, apply(ps, 2, function(z) sum(z < 0.05)))
}
pows <- ldply(res, pow)
ggplot(melt(pows, id.vars = 'N')) +
  geom_line(aes(x = N, y = value, col = variable))


# NOEC
noec <- function(z){
  noecs <- ldply(z$res, function(w) c(lm = w$noeclm, glm = w$noecglm))
  out <- melt(dcast(melt(noecs), variable ~ value))
  names(out)[2] <- 'noec'
  out$N <- z$N
  return(out)
}
noecs <- ldply(res, noec)
ggplot(noecs[noecs$noec %in% c(1,2), ]) +
  geom_line(aes(x = N, y = value, col = variable, linetype = noec)) 


# change means
# change control mean
ctrl <- exp(log(10) * seq(log10(1), log10(1000), by=0.2))
ctrl
# reduction of 50%
trt <- ctrl * 0.5

mum <- cbind(matrix(rep(ctrl, each = 2), ncol = 2, byrow = TRUE),
             matrix(rep(trt, each = 4), ncol = 4, byrow = TRUE))


dosim <- function(N = 3, mu = mu){
  Nj     <- rep(N, 6)
  mu     <- mu    # expected values in groups
  theta  <- rep(3.91, 6)                            # theta in groups
  mus    <- rep(mu, times = Nj)             # vector of mus
  thetas <- rep(theta, times=Nj)            # vector of thetas
  x      <- factor(rep(1:6, times=Nj))       # factor
  nsims   <- 200
  y      <- replicate(nsims, rnegbin(sum(Nj), mus, thetas))
  return(list(x = x, y = y))
}
sims <- NULL
for(i in seq_len(nrow(mum))){
  sims[[i]] <- dosim(N = 3, mu = mum[i, ])
}
plot(sims[[8]]$x, sims[[8]]$y[,4])

# analyse simulations
res <- lapply(sims, function(z){
  ana <- function(y, x){
    A <- 1/min(y[y!=0])         # ln(ax + 1) transformation
    yt <- log(A*y + 1)
    # models
    modlm <- glm(yt ~ x)
    modglm <- glm.nb(y ~ x)
    # global test
    plm <- anova(modlm, test = 'Chisq')["x", "Pr(>Chi)"]
    pglm <- anova(modglm, test = 'Chisq')["x", "Pr(>Chi)"]
    # multiple comparisons
    smclm <- summary(glht(modlm, linfct = mcp(x = 'Williams')))
    smcglm <- summary(glht(modglm, linfct = mcp(x = 'Williams')))
    # NOEC (which level? 1 = Control)
    noeclm <- sum(!smclm$test$pvalues < 0.05) + 1
    noecglm <- sum(!smcglm$test$pvalues < 0.05) + 1
    return(list(A = A, 
                # modlm = modlm, modglm=modglm, 
                plm=plm, pglm=pglm, 
                smclm=smclm, smcglm=smcglm, 
                noeclm=noeclm, noecglm=noecglm
                ))
  }
  # run models on simulated data
  res <- apply(z$y, 2, ana, x = z$x)
  res
})



# extract information
# global power
pow <- function(z){
  ps <- ldply(z, function(w) c(lm = w$plm, glm = w$pglm))
  apply(ps, 2, function(z) sum(z < 0.05)) / nsims
}
pows <- ldply(res, pow)
pows$muc <- ctrl
pows
ggplot(melt(pows, id.vars = 'muc')) +
  geom_line(aes(x = muc, y = value, col = variable)) +
  coord_trans(xtrans = 'log10') +
  scale_x_continuous(breaks = round(ctrl, 0))


# NOEC
noec <- function(z){
  noecs <- ldply(z, function(w) c(lm = w$noeclm, glm = w$noecglm))
  melt(noecs)
}
noecs <- ldply(res, noec)
noecs$muc <- rep(ctrl, each = 400)

ggplot(noecs) +
  geom_bar(aes(x = value, fill = variable), position = 'dodge') +
  facet_wrap(~muc)




#### Change both
N <- 2:12
ctrl <- exp(log(10) * seq(log10(1), log10(1000), by=0.2))
todo <- expand.grid(N = N, ctrl = ctrl)

dosim <- function(N = N, mu = mu){
  Nj     <- rep(N, 6)
  mu     <- mu    # expected values in groups
  theta  <- rep(3.91, 6)                            # theta in groups
  mus    <- rep(mu, times = Nj)             # vector of mus
  thetas <- rep(theta, times=Nj)            # vector of thetas
  x      <- factor(rep(1:6, times=Nj))       # factor
  nsims   <- 200
  y      <- replicate(nsims, rnegbin(sum(Nj), mus, thetas))
  return(list(x = x, y = y))
}
sims <- NULL
for(i in seq_len(nrow(todo))){
  N <- todo[i, 'N']
  ctrl <- todo[i, 'ctrl']
  trt <- ctrl * 0.5
  mu <- c(rep(ctrl, each = 2), rep(trt, each = 4))
  sims[[i]] <- dosim(N = N, mu = mu)
}


# analyse simulations
res <- llply(sims, function(z){
  ana <- function(y, x){
    A <- 1/min(y[y!=0])         # ln(ax + 1) transformation
    yt <- log(A*y + 1)
    # models
    modlm <- glm(yt ~ x)
    modglm <- glm.nb(y ~ x)
    # global test
    plm <- anova(modlm, test = 'Chisq')["x", "Pr(>Chi)"]
    pglm <- anova(modglm, test = 'Chisq')["x", "Pr(>Chi)"]
    # multiple comparisons
    smclm <- summary(glht(modlm, linfct = mcp(x = 'Williams')))
    smcglm <- summary(glht(modglm, linfct = mcp(x = 'Williams')))
    # NOEC (which level? 1 = Control)
    noeclm <- sum(!smclm$test$pvalues < 0.05) + 1
    noecglm <- sum(!smcglm$test$pvalues < 0.05) + 1
    return(list(A = A, modlm = modlm, modglm=modglm, plm=plm, pglm=pglm, 
                smclm=smclm, smcglm=smcglm, noeclm=noeclm, noecglm=noecglm))
  }
  # run models on simulated data
  res <- apply(z$y, 2, ana, x = z$x)
  res
}, .progress = 'text')

saveRDS(res, '/home/edisz/res.rds')

# extract information
# global power
pow <- function(z){
  ps <- ldply(z, function(w) c(lm = w$plm, glm = w$pglm))
  apply(ps, 2, function(z) sum(z < 0.05)) / nsims
}
pows <- ldply(res, pow)
pows$muc <- ctrl
pows
ggplot(melt(pows, id.vars = 'muc')) +
  geom_line(aes(x = muc, y = value, col = variable)) +
  coord_trans(xtrans = 'log10') +
  scale_x_continuous(breaks = round(ctrl, 0))


# NOEC
noec <- function(z){
  noecs <- ldply(z, function(w) c(lm = w$noeclm, glm = w$noecglm))
  melt(noecs)
}
noecs <- ldply(res, noec)
noecs$muc <- rep(ctrl, each = 400)

ggplot(noecs) +
  geom_bar(aes(x = value, fill = variable), position = 'dodge') +
  facet_wrap(~muc)


