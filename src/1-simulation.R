if(!exists('ld')){
  source("/home/edisz/Documents/Uni/Projects/PHD/6USETHEGLM/src/0-load.R")
}

#####--------------------------------------------------------------------------
##### Simulation 1 -  Count data
### Power
#####------------------------------------
# Settings
# no of simulated datasets
nsims <- 100
# sample sizes
N <- c(3, 6 ,9 ,12)
ctrl <- 2^(c(1:5, 7, 9))
# both as grid
todo <- expand.grid(N = N, ctrl = ctrl)
# fixed theta
theta  <- rep(3.91, 6)  


#####------------------------------------
# simulate data
sims <- NULL
set.seed(seed)
for(i in seq_len(nrow(todo))){
  N <- todo[i, 'N']
  takectrl <- todo[i, 'ctrl']
  # reduce t2-t5 to 50%
  taketrt <- takectrl * 0.5
  mu <- c(rep(takectrl, each = 2), rep(taketrt, each = 4))
  sims[[i]] <- dosim(N = N, mu = mu, nsims = nsims, theta = theta)
}

# plot one realisation of simulated data
todo[15, ]
df <- data.frame(x = sims[[15]]$x, y = sims[[15]]$y[ , 2])
df$yt <- log(1 / min(df$y[df$y != 0]) * df$y + 1)
dfm <- melt(df)
levels(dfm$variable) <- c('y', 'ln(Ay + 1)')
p1 <- ggplot(dfm, aes(x = x, y = value)) +
  geom_boxplot(fill = 'grey80') +
  facet_wrap( ~variable, scales = 'free_y') +
  scale_x_discrete(labels = c('C', 'T1', 'T2', 'T3', 'T4', 'T5')) +
  labs(x = 'Treatment', y = 'Abundance') +
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face = "bold"))
p1
if(exp_plot){
  ggsave(file.path(figdir, 'p1.pdf'), p1, width = 10, height = 6)
}


#####------------------------------------
if(sim1){
# analyse simulations
#   res <- llply(sims[1:2], resfoo, .progress = 'text')
  res1 <- llply(sims, resfoo, .progress = 'text')
  saveRDS(res1, file.path(cachedir, 'res1.rds'))
} else {
  res1 <- readRDS(file.path(cachedir, 'res1.rds'))
}


#####------------------------------------
# Results
# global power
pow <- function(z){
  ps <- ldply(z, function(w) c(lm = w$plm, glm = w$pglm, pk = w$pk))
  apply(ps, 2, function(z) sum(z < 0.05)) / length(z)
}
pows <- ldply(res1, pow)
pows$muc <- todo$ctrl
pows$N <- todo$N
powsm <- melt(pows, id.vars = c('muc', 'N'))

p2 <- ggplot(powsm) +
  geom_line(aes(y = value, x = muc, group = variable)) +
  geom_point(aes(y = value, x = muc, fill = variable), size = 4, pch = 21, color = 'black') +
  coord_trans(xtrans = 'log2') +
  scale_x_continuous(breaks = round(ctrl, 0)) +
  facet_wrap(~N) + 
  # axes
  labs(x = expression(mu[C]), 
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
                  breaks = c('lm', 'glm', 'pk'), 
                  labels = c('LM + log(Ay+1)', 'GLM (neg. bin.)', 'Kruskal'),
                  start = 0, end = 1) +
  theme(legend.position = "bottom", legend.key = element_blank())
p2
if(exp_plot){
  ggsave(file.path(figdir, 'p2.pdf'), p2, width = 11, height = 11)
}



### loec
loec <- function(z){
  loecs <- ldply(z, function(w) c(lm = w$loeclm, glm = w$loecglm, pw = w$loecpw
                                  #, np = w$loecnp
                                  ))
  out <- apply(loecs, 2, function(x) sum(x == 2))
  return(out)
}
loecs <- ldply(res1, loec)
loecs$muc <- todo$ctrl
loecs$N <- todo$N
loecsm <- melt(loecs, id.vars = c('muc', 'N'))
loecsm$value <- loecsm$value / nsims

p3 <- ggplot(loecsm) +
  geom_line(aes(y = value, x = muc, group = variable)) +
  geom_point(aes(y = value, x = muc, fill = variable), size = 4, pch = 21, color = 'black') +
  coord_trans(xtrans = 'log2') +
  scale_x_continuous(breaks = round(unique(todo$ctrl), 0)) +
  facet_wrap( ~N) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Power (LOEC , ', alpha, ' = 0.05)'))) + 
  # appearance
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size = 14),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 14,face = "bold")) +
  # legend
  scale_fill_grey(name = '', 
                  breaks = c('lm', 'glm', 'pw', 'np'), 
                  labels = c('LM + log(Ay+1)', 'GLM (neg. bin.)', 'Wilcox', 'NRCE'),
                  start = 0, end = 1) +
  theme(legend.position="bottom", legend.key = element_blank())
p3
if(exp_plot){
  ggsave(file.path(figdir, 'p3.pdf'), p3, width = 11, height = 11)
}





### ----------------------------------------------------------------------------
# Type1 Error

#####------------------------------------
# settings
nsims <- 100
# sample sizes
N <- c(3, 6, 9, 12)
ctrl <- 2^(c(1:5, 7, 9))
# both as grid
todo <- expand.grid(N = N, ctrl = ctrl)
theta  <- rep(3.91, 6)  


#####------------------------------------
# simulate data
sims <- NULL
set.seed(seed)
for(i in seq_len(nrow(todo))){
  N <- todo[i, 'N']
  takectrl <- todo[i, 'ctrl']
  # all treatments with same mean
  taketrt <- takectrl * 1
  mu <- c(rep(takectrl, each = 2), rep(taketrt, each = 4))
  sims[[i]] <- dosim(N = N, mu = mu, nsims = nsims, theta = theta)
}


#####------------------------------------
# analyse data
res2 <- llply(sims, resfoo, .progress = 'text')



#####------------------------------------
# Results
# Global test (how often wrongly assigned an effect)
t1 <- function(z){
#   res2[[1]][[1]]$modglm@details$convergence == 0
  ps <- ldply(z, function(w) c(lm = w$plm, glm = w$pglm, pk = w$pk))
  apply(ps, 2, function(z) sum(z < 0.05)) / length(z)
}
t1s <- ldply(res2, t1)
t1s$muc <- todo$ctrl
t1s$N <- todo$N
t1sm <- melt(t1s, id.vars = c('muc', 'N'))

pt1 <- ggplot(t1sm) +
  geom_line(aes(y = value, x = muc, group = variable)) +
  geom_point(aes(y = value, x = muc, fill = variable), size = 4, pch = 21, color = 'black') +
  geom_segment(aes(x = 2, xend = 1024, y = 0.05, yend = 0.05), linetype = 'dashed') + 
  # maby try log2 transformation?
  coord_trans(xtrans = 'log2') +
  # use here
  scale_x_continuous(breaks = round(ctrl, 0)) +
  facet_wrap(~N) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Type 1 error (global test , ', alpha, ' = 0.05)'))) +
  # appearance
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size=14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face="bold")) +
  # legend
  scale_fill_grey(name = '', 
                  breaks = c('lm', 'glm', 'pk'), 
                  labels = c('LM + log(Ay+1)', 'GLM (neg. bin.)', 'Kruskal'),
                  start = 0, end = 1) +
  theme(legend.position = "bottom", legend.key = element_blank())
pt1
if(exp_plot){
  ggsave(file.path(figdir, 'pt1.pdf'), pt1, width = 11, height = 11)
}


# MC (how often assigned wrongly a LOEC)
loec <- function(z){
  loecs <- ldply(z, function(w) c(lm = w$loeclm, glm = w$loecglm, pw = w$loecpw
                                  #, np = w$loecnp
  ))
  out <- apply(loecs, 2, function(x) sum(x != Inf))
  return(out)
}
loecs <- ldply(res2, loec)
loecs$muc <- todo$ctrl
loecs$N <- todo$N
loecsm <- melt(loecs, id.vars = c('muc', 'N'))
loecsm$value <- loecsm$value / nsims



pt2 <- ggplot(loecsm) +
  geom_line(aes(y = value, x = muc, group = variable)) +
  geom_point(aes(y = value, x = muc, fill = variable), size = 4, pch = 21, color = 'black') +
  geom_segment(aes(x = 2, xend = 1024, y = 0.05, yend = 0.05), linetype = 'dashed') + 
  coord_trans(xtrans = 'log2') +
  scale_x_continuous(breaks = round(unique(todo$ctrl), 0)) +
  facet_wrap(~N) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Type 1 error (LOEC , ', alpha, ' = 0.05)'))) + 
  # appearance
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size = 14),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 14, face = "bold")) +
  # legend
  scale_fill_grey(name = '', 
                  breaks = c('lm', 'glm', 'pw'), 
                  labels = c('LM + log(Ay+1)', 'GLM (neg. bin.)', 'Wilcox'),
                  start = 0, end = 1) +
  theme(legend.position="bottom", legend.key = element_blank())
pt2
if(exp_plot){
  ggsave(file.path(figdir, 'pt2.pdf'), pt2, width = 11, height = 11)
}
