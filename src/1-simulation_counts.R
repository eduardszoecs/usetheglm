if(!exists("prj")){
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project'!")
} else {
  source(file.path(prj, "src", "0-load.R"))
}

### ----------------------------------------------------------------------------
### Simulation 1 -  Count data
### Written by Eduard Szöcs
### ----------------------------------------------------------------------------

#####--------------------------------------------------------------------------
### Power simulations

#####------------------------------------
### simulate data
# Settings
# no of simulated datasets
nsims <- 100
# sample sizes
N <- c(3, 6 ,9)
ctrl <- 2^(c(1:7))
# for testing
# nsims <- 20
# N <- 3
# ctrl <- 2^(c(1:3))
# both as grid
todo1_c <- expand.grid(N = N, ctrl = ctrl)
# fixed theta
theta  <- rep(4, 6)  

# create datasets
sims1_c <- NULL
# set.seed(seed)
for(i in seq_len(nrow(todo1_c))){
  N <- todo1_c[i, 'N']
  takectrl <- todo1_c[i, 'ctrl']
  # reduce t2-t5 to 50%
  taketrt <- takectrl * 0.5
  mu <- c(rep(takectrl, each = 2), rep(taketrt, each = 4))
  sims1_c[[i]] <- dosim1(N = N, mu = mu, nsims = nsims, theta = theta)
}

# plot one realisation of simulated data
todo1_c[15, ]
df <- data.frame(x = sims1_c[[15]]$x, y = sims1_c[[15]]$y[ , 2])
df$yt <- log(1 / min(df$y[df$y != 0]) * df$y + 1)
dfm <- melt(df)
levels(dfm$variable) <- c('y', 'ln(Ay + 1)')
ggplot(dfm, aes(x = x, y = value)) +
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



#####------------------------------------
# analyse simulations
if(sim1){
  res1_c <- llply(sims1_c, resfoo1, .progress = 'text', npb = 500)
  saveRDS(res1_c, file.path(cachedir, 'res1_c.rds'))
} else {
  res1_c <- readRDS(file.path(cachedir, 'res1_c.rds'))
}


#####------------------------------------
# Results
# global power
pow_glob_c <- ldply(res1_c, p_glob1)
pow_glob_c$muc <- rep(todo1_c$ctrl, each  = 5)
pow_glob_c$N <- rep(todo1_c$N, each = 5)
pow_glob_c$variable <-  factor(pow_glob_c$variable, unique(pow_glob_c$variable)[1:5], 
       labels = c('lm', 'glm_nb', 'glm_qp', 'glm_pb', 'np'))

plot_pow_glob_c <- ggplot(pow_glob_c) +
  geom_line(aes(y = power, x = log2(muc), group = variable, linetype = variable)) +
  geom_point(aes(y = power, x = log2(muc), shape = variable), color = 'black', size = 4) +
  facet_grid( ~N, labeller = n_labeller) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Power (global test , ', alpha, ' = 0.05)'))) +
  scale_x_continuous(breaks = log2(round(unique(todo1_c$ctrl), 0)), 
                     labels = round(unique(todo1_c$ctrl), 0)) +
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(16,2,4,0,17), 
                     labels = c('LM', expression(GLM[nb]), expression(GLM[qp]), 
                                expression(GLM[npb]), 'KW')) +
  scale_linetype_discrete('Method', labels = c('LM', expression(GLM[nb]), expression(GLM[qp]), 
                                               expression(GLM[npb]), 'KW')) +
  ylim(c(0,1))
plot_pow_glob_c

# # # plot for presentation
# pow_glob_c$glm <- ifelse(pow_glob_c$variable %in% c('lm', 'np'), 'Standard', 'GLM')
# library(extrafont)
# require(xkcd)
# loadfonts() # fonts appear only in pdf
# # for the dataman
# ratioxy <- 8
# mapping <- aes(x=x,
#                y=y,
#                scale=scale,
#                ratioxy=ratioxy,
#                angleofspine = angleofspine,
#                anglerighthumerus = anglerighthumerus,
#                anglelefthumerus = anglelefthumerus,
#                anglerightradius = anglerightradius,
#                angleleftradius = angleleftradius,
#                anglerightleg =  anglerightleg,
#                angleleftleg = angleleftleg,
#                angleofneck = angleofneck)
# dataman <- data.frame( x=4.5, y=0.4, N = 9,
#                        scale = 0.15 ,
#                        ratioxy = ratioxy,
#                        angleofspine = -pi/2 ,
#                        anglerighthumerus = pi/2 + pi/2,
#                        anglelefthumerus = pi/2 + pi/3, 
#                        anglerightradius = pi/2 +  pi/4,
#                        angleleftradius = pi/2 +  pi/4, 
#                        angleleftleg = 3*pi/2 + pi / 12 ,
#                        anglerightleg = 3*pi/2 - pi / 12,
#                        angleofneck = runif(1, 3*pi/2-pi/10, 3*pi/2+pi/10)
#                        )
# myp <- ggplot(pow_glob_c) +
#   geom_smooth(aes(y = power, x = log2(muc), group = variable, color = glm), 
#               size= 1, se = FALSE, span = 0.8, shape = 'X',
#               position = position_jitter(width = 0.025)) +
#   # geom_line(aes(y = power, x = log2(muc), group = variable, color = glm), linewidth = 4) +
#   # geom_point(aes(y = power, x = log2(muc), color = glm), size = 4, shape = 4) +
#   facet_grid( ~N, labeller = n_labeller) + 
#   # axes
#   labs(x = expression('Effect size'), 
#        y = 'Power') +
#   # xkcd appearance
#   theme_xkcd() +
#   xkcdaxis(c(1,7), c(0, 1)) + 
#   scale_y_continuous(breaks = c(0, 0.8, 1), labels = c('0', '80', '100')) + 
#   scale_x_continuous(breaks = c(1.5, 6.5), labels = c('small', 'big')) +
#   # colors
#   scale_color_manual('', values =  c('#FA6C00', 'gray25'), guide = guide_legend(override.aes = list(size = 2)))   +
#   theme(legend.justification=c(0,1), 
#         legend.position=c(0,1), 
#         axis.text=element_text(size=22),
#         axis.title=element_text(size=25,face="bold"),
#         strip.text.x = element_text(size = 25),
#         legend.text=element_text(size = 22)) +
#   # decoration
#   xkcdman(mapping, dataman) +
#   geom_text(data=data.frame(x=5.7, y = 0.55, N = 9, label = "Lo and behold!"),
#                             aes(x=x,y=y,label=label), size=6,
#             show_guide=F, family="xkcd") +
#   xkcdline(mapping=aes(xbegin=xb, ybegin=yb, xend =xe, yend= ye), 
#            data = data.frame(xb = 5.2, xe = 6, yb = 0.45, ye = 0.5, N = 9), 
#            xjitteramount = 0.4)
# myp 
# ggsave('/home/edisz/Documents/Uni/Projects/PHD/MISC/AG_presentations/fig/myp2.pdf', myp, width = 10, height = 6)


# max power for N = 3, excluding glm_nb
max(pow_glob_c[pow_glob_c$N == 3 & !pow_glob_c$variable %in% 'glm_nb', 'power'])
diffdf <- dcast(pow_glob_c, muc + N ~ variable, value.var = 'power')
diffdf$diff_lm_qp <- diffdf$lm - diffdf$glm_qp
diffdf

### loec
pow_loec_c <- ldply(res1_c, p_loec1, type = 'power')
pow_loec_c$muc <- todo1_c$ctrl
pow_loec_c$N <- todo1_c$N
pow_loec_c  <- melt(pow_loec_c, id.vars = c('muc', 'N'), value.name = 'power')
pow_loec_c$variable  <-  factor(pow_loec_c$variable, unique(pow_loec_c$variable)[1:4], 
                                labels = c('lm', 'glm_nb', 'glm_qp', 'np'))

plot_pow_loec_c <- ggplot(pow_loec_c) +
  geom_line(aes(y = power, x = log2(muc), group = variable, linetype = variable)) +
  geom_point(aes(y = power, x = log2(muc), shape = variable), color = 'black', size = 4) +
  facet_grid( ~N, labeller = n_labeller) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Power (LOEC , ', alpha, ' = 0.05)'))) + 
  scale_x_continuous(breaks = log2(round(unique(todo1_c$ctrl), 0)), 
                     labels = round(unique(todo1_c$ctrl), 0)) +
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(16,2,4,17), 
                     labels = c('LM', expression(GLM[nb]), expression(GLM[qp]), 
                                                               'WT')) +
  scale_linetype_discrete('Method',  labels = c('LM', expression(GLM[nb]), expression(GLM[qp]), 
                                                'WT')) +
  ylim(c(0,1))
plot_pow_loec_c

max(pow_loec_c[pow_loec_c$N == 3 & !pow_loec_c$variable %in% c('glm_nb'), 'power'])


# compare power
merged <- merge(pow_glob_c[ , c(1,2,4, 5)], pow_loec_c, by = c('N', 'muc', 'variable'), suffixes = c('glob', 'loec'))
merged$diff <- merged$powerloec - merged$powerglob
plot(merged$diff)
max(abs(merged$diff[!merged$variable %in% c('glm_nb', 'np')]))
max(abs(merged$diff[merged$variable %in% c('lm')]))
merged[merged$variable == 'lm', ]
diffdf <- dcast(pow_loec_c, muc + N ~ variable, value.var = 'power')
diffdf$diff_lm_qp <- diffdf$lm - diffdf$glm_qp
diffdf



### ----------------------------------------------------------------------------
# Type1 Error simulations

#####------------------------------------
### simulate data
# settings
nsims <- 100
# sample sizes
N <- c(3, 6, 9)
ctrl <- 2^(c(1:7))
# nsims <- 20
# N <- 3
# ctrl <- 2^(c(1, 3, 5, 7, 9))
# both as grid
todo2_c <- expand.grid(N = N, ctrl = ctrl)
theta  <- rep(3.91, 6)  

# create simulate data
sims2_c <- NULL
# set.seed(seed)
for(i in seq_len(nrow(todo2_c))){
  N <- todo2_c[i, 'N']
  takectrl <- todo2_c[i, 'ctrl']
  # all treatments with same mean
  taketrt <- takectrl * 1
  mu <- c(rep(takectrl, each = 2), rep(taketrt, each = 4))
  sims2_c[[i]] <- dosim1(N = N, mu = mu, nsims = nsims, theta = theta)
}


#####------------------------------------
# analyse data
if(sim1){
  res2_c <- llply(sims2_c, resfoo1, .progress = 'text', npb = 400)
  saveRDS(res2_c, file.path(cachedir, 'res2_c.rds'))
} else {
  res2_c <- readRDS(file.path(cachedir, 'res2_c.rds'))
}


#####------------------------------------
# Results
# Global test (how often wrongly assigned an effect)
t1_glob_c <- ldply(res2_c, p_glob1)
t1_glob_c$muc <- rep(todo2_c$ctrl, each  = 5)
t1_glob_c$N <- rep(todo2_c$N, each = 5)
t1_glob_c$variable  <-  factor(t1_glob_c$variable, unique(t1_glob_c$variable)[c(1, 2, 3, 4, 5)], 
                                labels = c('lm', 'glm_nb', 'glm_pb', 'glm_qp', 'np'))

plot_t1_glob_c  <- ggplot(t1_glob_c) +
  geom_line(aes(y = power, x = log2(muc), group = variable, linetype = variable)) +
  geom_point(aes(y = power, x = log2(muc), shape = variable), color = 'black', size = 4) +
  geom_segment(aes(x = -Inf, xend = Inf, y = 0.05, yend = 0.05), linetype = 'dashed') + 
  facet_grid( ~N, labeller = n_labeller) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Type 1 error (global test , ', alpha, ' = 0.05)'))) +
  scale_x_continuous(breaks = log2(round(unique(todo2_c$ctrl), 0)), 
                     labels = round(unique(todo2_c$ctrl), 0)) +
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(16,2,4,0,17)) +
  scale_linetype_discrete('Method')
plot_t1_glob_c


# LOEC
t1_loec_c <- ldply(res2_c, p_loec1, type = 't1')
t1_loec_c$muc <- todo2_c$ctrl
t1_loec_c$N <- todo2_c$N
t1_loec_c <- melt(t1_loec_c, id.vars = c('muc', 'N'), value.name = 't1')
t1_loec_c$variable  <-  factor(t1_loec_c$variable, unique(t1_loec_c$variable)[c(1, 2, 3, 4)], 
                                labels = c('lm', 'glm_nb', 'glm_qp', 'np'))

plot_t1_loec_c <- ggplot(t1_loec_c) +
  geom_line(aes(y = t1, x = log2(muc), group = variable, linetype = variable)) +
  geom_point(aes(y = t1, x = log2(muc), shape = variable), color = 'black', size = 4) +
  geom_segment(aes(x = -Inf, xend = Inf, y = 0.05, yend = 0.05), 
               linetype = 'dashed') + 
  facet_grid( ~N, labeller = n_labeller) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Type 1 error (LOEC , ', alpha, ' = 0.05)'))) + 
  scale_x_continuous(breaks = log2(round(unique(todo2_c$ctrl), 0)), 
                     labels = round(unique(todo2_c$ctrl), 0)) +
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(16,2,4,17)) +
  scale_linetype_discrete('Method')
plot_t1_loec_c 

