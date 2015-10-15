if(!exists("prj")){
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/usetheglm'!")
} else {
  source(file.path(prj, "src", "0-load.R"))
}

### ----------------------------------------------------------------------------
### Results -  Count data
### Written by Eduard Szöcs
### ----------------------------------------------------------------------------

## Load results
source(file.path(prj, "src", "1-simulations.R"))
# Power 
res1_c <- readRDS(file.path(cachedir, 'res1_c.rds'))
# Type I error
res2_c <- readRDS(file.path(cachedir, 'res2_c.rds'))


### ----------------------------------------------------------------------------
### Global test 
### ------------------------
# Plot Power
pow_glob_c <- ldply(res1_c, p_glob1)
pow_glob_c$muc <- rep(todo_c$ctrl, each  = 6)
pow_glob_c$N <- rep(todo_c$N, each = 6)
pow_glob_c$variable <-  factor(pow_glob_c$variable, unique(pow_glob_c$variable)[1:6], 
                               labels = c('lm', 'glm_nb', 'glm_qp', 'glm_pb', 'glm_p', 'np'))

plot_pow_glob_c <- ggplot(pow_glob_c) +
  geom_line(aes(y = power, x = log2(muc), group = variable, linetype = variable)) +
  geom_point(aes(y = power, x = log2(muc), shape = variable), color = 'black', size = 4) +
  facet_grid( ~N, labeller = n_labeller) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Power (global test , ', alpha, ' = 0.05)'))) +
  scale_x_continuous(breaks = log2(round(unique(todo_c$ctrl), 0)), 
                     labels = round(unique(todo_c$ctrl), 0)) +
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(16,2,4,0, 1, 17), 
                     labels = c('LM', expression(GLM[nb]), expression(GLM[qp]), 
                                expression(GLM[npb]), expression(GLM[p]), 'KW')) +
  scale_linetype_discrete('Method', labels = c('LM', expression(GLM[nb]), expression(GLM[qp]), 
                                               expression(GLM[npb]), expression(GLM[p]), 'KW')) +
  ylim(c(0,1))
plot_pow_glob_c


### ------------------------
### Tabular output
ldf <- dcast(pow_glob_c, N + muc ~ variable, value.var = "power")
colnames(ldf) <- c('N', '$\\mu_C$', 'LM', '$GLM_{nb}$', '$GLM_{qp}$', '$GLM_{p}$', 'np')
print(xtable(ldf, 
             caption = 'Count data simulations - Power to detect a treatment effect. N = sample sizes, 
             $\\mu_C$ = mean abundance in control, LM = Linear model after transformation, 
             $GLM_{nb}$ = negative binomial model, $GLM_{qp}$ = quasi-Poisson model, $GLM_{qp}$ = Poisson model,
            np = pairwise Wilcoxon test.', label = 'tab:pow_glob_c'),
      file = file.path(suppdir, "supp1", "tab", "pow_glob_c.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})




### ------------------------
## Plot T1-error
t1_glob_c <- ldply(res2_c, p_glob1)
t1_glob_c$muc <- rep(todo_c$ctrl, each  = 6)
t1_glob_c$N <- rep(todo_c$N, each = 6)
t1_glob_c$variable  <-  factor(t1_glob_c$variable, unique(t1_glob_c$variable)[c(1, 2, 3, 4, 5, 6)], 
                               labels = c('lm', 'glm_nb', 'glm_pb', 'glm_qp', 'glm_p', 'np'))

plot_t1_glob_c  <- ggplot(t1_glob_c) +
  geom_line(aes(y = power, x = log2(muc), group = variable, linetype = variable)) +
  geom_point(aes(y = power, x = log2(muc), shape = variable), color = 'black', size = 4) +
  geom_segment(aes(x = -Inf, xend = Inf, y = 0.05, yend = 0.05), linetype = 'dashed') + 
  facet_grid( ~N, labeller = n_labeller) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Type 1 error (global test , ', alpha, ' = 0.05)'))) +
  scale_x_continuous(breaks = log2(round(unique(todo_c$ctrl), 0)), 
                     labels = round(unique(todo_c$ctrl), 0)) +
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(16,2,4,0,1,17)) +
  theme(legend.position="none") + labs(x = NULL) +
  scale_y_log10(breaks = c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1))
  # coord_cartesian(ylim = c(-0.005, 0.25))
plot_t1_glob_c


# poster plot
# pp_t1  <- ggplot(t1_glob_c) +
#   geom_line(aes(y = power, x = log2(muc), group = variable, linetype = variable), color =  '#2A6491', size = 1) +
#   geom_point(aes(y = power, x = log2(muc), shape = variable), color =  '#2A6491', size = 5) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = 0.05, yend = 0.05), linetype = 'dashed') + 
#   facet_grid( ~N, labeller = n_labeller) + 
#   # axes
#   labs(x = expression(mu[C]), 
#        y = expression(paste('Type 1 error (global test , ', alpha, ' = 0.05)'))) +
#   scale_x_continuous(breaks = log2(round(unique(todo_c$ctrl), 0))[c(1, 3, 5, 7)], 
#                      labels = round(unique(todo_c$ctrl), 0)[c(1, 3, 5, 7)]
#                      )+
#   # appearance
#   mytheme + 
#   # legend title
#   scale_y_log10(breaks = c(0.01, 0.05, 0.25, 1)) +
#   scale_shape_manual('Method', values=c(16,2,4,0, 1, 17), 
#                      labels = c('LM', expression(GLM[nb]), expression(GLM[qp]), 
#                                 expression(GLM[npb]), expression(GLM[p]), 'KW')) +
#   scale_linetype_discrete('Method', labels = c('LM', expression(GLM[nb]), expression(GLM[qp]), 
#                                                expression(GLM[npb]), expression(GLM[p]), 'KW')) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         text = element_text(size=25),
#         axis.text = element_text(size=25),
#         axis.title.x = element_text(size=25,face="bold", vjust = 0),
#         axis.title.y = element_text(size=25,face="bold", vjust = 1),
#         legend.position="bottom",
#         legend.key = element_blank(),
#         strip.background = element_blank(),
#         strip.text= element_text(size=25, face = 'bold'))
# pp_t1




### ------------------------
### Tabular output
ldf <- dcast(t1_glob_c, N + muc ~ variable, value.var = "power")
ldf <- ldf[ , c(1,2,3,4,6,5,7, 8)]
colnames(ldf) <- c('N', '$\\mu_C$', 'LM', '$GLM_{nb}$', '$GLM_{qp}$', '$GLM_{pb}$', '$GLM_{p}$', 'np')
print(xtable(ldf, 
             caption = 'Count data simulations - Type 1 error to detect a global treatment effect. N = sample sizes, 
             $\\mu_C$ = mean abundance in control, LM = Linear model after transformation, 
             $GLM_{nb}$ = negative binomial model, $GLM_{qp}$ = quasi-Poisson model, 
             $GLM_{pb}$ = negative binomial model with parametric boostrap, 
             $GLM_{p}$ = Poisson model, np = Kruskal-Wallis test.', 
             label = 'tab:t1_glob_c'),
      file = file.path(suppdir, "supp1", "tab", "t1_glob_c.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})


### ------------------------
## Figure 2 
leg <- g_legend(plot_pow_glob_c + 
                  theme(legend.key = element_blank(),
                        legend.text = element_text(size = 16),
                        legend.title = element_text(size = 18, face = "bold"),
                        legend.key.size = unit(1.5, "cm"))
)
# combine plots and legend
p_glob_c <- arrangeGrob(
  plot_t1_glob_c + theme(legend.position="none") + labs(x = NULL),
  plot_pow_glob_c + theme(legend.position="none"),
  leg,
  nrow = 3,
  heights = c(10, 10, 2))
p_glob_c
if(exp_plot)
  ggsave(file.path(figdir, 'p_glob_c.pdf'), p_glob_c, width = 10, height = 8)


### ------------------------
# max power for N = 3, excluding glm_nb
max(pow_glob_c[pow_glob_c$N == 3 & !pow_glob_c$variable %in% c('glm_nb', 'glm_p'), 'power'])
diffdf <- dcast(pow_glob_c, muc + N ~ variable, value.var = 'power')
diffdf$diff_lm_qp <- diffdf$lm - diffdf$glm_qp
diffdf




### ---------------------------------------------------------------------------
### LOEC
### ------------------------
## Plot Power
pow_loec_c <- ldply(res1_c, p_loec1, type = 'power')
pow_loec_c$muc <- todo_c$ctrl
pow_loec_c$N <- todo_c$N
pow_loec_c  <- melt(pow_loec_c, id.vars = c('muc', 'N'), value.name = 'power')
pow_loec_c$variable  <-  factor(pow_loec_c$variable, unique(pow_loec_c$variable)[1:5], 
                                labels = c('lm', 'glm_nb', 'glm_qp', 'glm_p', 'np'))

plot_pow_loec_c <- ggplot(pow_loec_c) +
  geom_line(aes(y = power, x = log2(muc), group = variable, linetype = variable)) +
  geom_point(aes(y = power, x = log2(muc), shape = variable), color = 'black', size = 4) +
  facet_grid( ~N, labeller = n_labeller) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Power (LOEC , ', alpha, ' = 0.05)'))) + 
  scale_x_continuous(breaks = log2(round(unique(todo_c$ctrl), 0)), 
                     labels = round(unique(todo_c$ctrl), 0)) +
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(16,2,4,1,17), 
                     labels = c('LM', expression(GLM[nb]), expression(GLM[qp]), 
                                expression(GLM[p]), 'WT')) +
  scale_linetype_discrete('Method',  labels = c('LM', expression(GLM[nb]), 
                                                expression(GLM[qp]), expression(GLM[p]), 'WT')) +
  ylim(c(0,1))
# +
#   scale_y_log10(breaks = c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1), limit = c(0.01, 1))
plot_pow_loec_c

# ### Poster plot
# pppdf <- pow_loec_c[pow_loec_c$variable %in% c('lm', 'glm_qp', 'np'), ]
# pppdf <- droplevels(pppdf)
# ppp <- ggplot(pppdf) +
#   geom_line(aes(y = power, x = log2(muc), group = variable, linetype = variable), color = '#2A6491', size = 1) +
#   geom_point(aes(y = power, x = log2(muc), shape = variable), color = '#2A6491', size = 5) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = 0.8, yend = 0.8), linetype = 'dashed') + 
#   facet_grid( ~N, labeller = n_labeller) + 
#   # axes
#   labs(x = expression(mu[C]), 
#        y = expression(paste('Power (LOEC , ', alpha, ' = 0.05)'))) + 
#   scale_x_continuous(breaks = log2(round(unique(todo_c$ctrl), 0))[c(1, 3, 5, 7)], 
#                      labels = round(unique(todo_c$ctrl), 0)[c(1, 3, 5, 7)]
#   )+
#   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.8, 1), limits = c(0,1)) +
#   scale_shape_manual('Method', values=c(16,4,17), 
#                      labels = c('LM', expression(GLM[qp]), 'WT')) +
#   scale_linetype_discrete('Method',  labels = c('LM', 
#                                                 expression(GLM[qp]), 'WT')) +
#   mytheme + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         text = element_text(size=25),
#         axis.text = element_text(size=25),
#         axis.title.x = element_text(size=25,face="bold", vjust = 0),
#         axis.title.y = element_text(size=25,face="bold", vjust = 1),
#         legend.position="bottom",
#         legend.key = element_blank(),
#         strip.background = element_blank(),
#         strip.text= element_text(size=25, face = 'bold'))
# ppp

### ------------------------
### Tabular output
ldf <- dcast(pow_loec_c, N + muc ~ variable, value.var = "power")
colnames(ldf) <- c('N', '$\\mu_C$', 'LM', '$GLM_{nb}$', '$GLM_{qp}$', '$GLM_{p}$', 'np')
print(xtable(ldf, 
             caption = 'Count data simulations - Power to detect LOEC. N = sample sizes, 
             $\\mu_C$ = mean abundance in control, LM = Linear model after transformation, 
             $GLM_{nb}$ = negative binomial model, $GLM_{qp}$ = quasi-Poisson model, 
            $GLM_{p}$ = Poisson model, np = pairwise Wilcoxon test.', 
             label = 'tab:pow_loec_c'),
      file = file.path(suppdir, "supp1", "tab", "pow_loec_c.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})


### ------------------------
## Plot T1-error
t1_loec_c <- ldply(res2_c, p_loec1, type = 't1')
t1_loec_c$muc <- todo_c$ctrl
t1_loec_c$N <- todo_c$N
t1_loec_c <- melt(t1_loec_c, id.vars = c('muc', 'N'), value.name = 't1')
t1_loec_c$variable  <-  factor(t1_loec_c$variable, unique(t1_loec_c$variable)[c(1, 2, 3, 4, 5)], 
                               labels = c('lm', 'glm_nb', 'glm_qp', 'glm_p', 'np'))

plot_t1_loec_c <- ggplot(t1_loec_c) +
  geom_line(aes(y = t1, x = log2(muc), group = variable, linetype = variable)) +
  geom_point(aes(y = t1, x = log2(muc), shape = variable), color = 'black', size = 4) +
  geom_segment(aes(x = -Inf, xend = Inf, y = 0.05, yend = 0.05), 
               linetype = 'dashed') + 
  facet_grid( ~N, labeller = n_labeller) + 
  # axes
  labs(x = expression(mu[C]), 
       y = expression(paste('Type 1 error (LOEC , ', alpha, ' = 0.05)'))) + 
  scale_x_continuous(breaks = log2(round(unique(todo_c$ctrl), 0)), 
                     labels = round(unique(todo_c$ctrl), 0)) +
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(16,2,4,1, 17)) +
  theme(legend.position="none") + labs(x = NULL) +
  scale_y_log10(breaks = c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1), limit = c(0.01, 1))
plot_t1_loec_c 


### ------------------------
### Tabular output
ldf <- dcast(t1_loec_c, N + muc ~ variable, value.var = "t1")
colnames(ldf) <- c('N', '$\\mu_C$', 'LM', '$GLM_{nb}$', '$GLM_{qp}$', '$GLM_{p}$', 'np')
print(xtable(ldf, 
             caption = 'Count data simulations - Type 1 error to detect LOEC. N = sample sizes, 
             $\\mu_C$ = mean abundance in control, LM = Linear model after transformation, 
             $GLM_{nb}$ = negative binomial model, $GLM_{qp}$ = quasi-Poisson model, 
            $GLM_{p}$ = Poisson model, np = pairwise Wilcoxon.', 
             label = 'tab:t1_loec_c'),
      file = file.path(suppdir, "supp1", "tab", "t1_loec_c.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})

### ------------------------
## Figure 3
leg <- g_legend(plot_pow_loec_c + 
                  theme(legend.key = element_blank(),
                        legend.text = element_text(size = 16),
                        legend.title = element_text(size = 18, face = "bold"),
                        legend.key.size = unit(1.5, "cm")))
p_loec_c <- arrangeGrob(
  plot_t1_loec_c + theme(legend.position="none") + labs(x = NULL),
  plot_pow_loec_c + theme(legend.position="none"),
  leg,
  nrow = 3,
  heights = c(10, 10, 2))
p_loec_c
if(exp_plot)
  ggsave(file.path(figdir, 'p_loec_c.pdf'), p_loec_c, width = 10, height = 8)


### ------------------------
### Misc
max(pow_loec_c[pow_loec_c$N == 3 & !pow_loec_c$variable %in% c('glm_nb', 'glm_p'), 'power'])

# compare power
merged <- merge(pow_glob_c[ , c(1,2,4, 5)], 
                pow_loec_c, by = c('N', 'muc', 'variable'), 
                suffixes = c('glob', 'loec'))
merged$diff <- merged$powerloec - merged$powerglob
plot(merged$diff)
merged
ggplot(merged, aes(x = variable, y = diff)) +
  geom_boxplot()
a <- merged[!merged$variable %in% c('glm_nb', 'np', 'glm_p'), ]
# mean reduction
ddply(a, .(variable), summarise, mean(diff))

max(abs(merged$diff[!merged$variable %in% c('glm_nb', 'np')]))
max(abs(merged$diff[merged$variable %in% c('lm')]))
merged[merged$variable == 'lm', ]
diffdf <- dcast(pow_loec_c, muc + N ~ variable, value.var = 'power')
diffdf$diff_lm_qp <- diffdf$lm - diffdf$glm_qp
diffdf


# Convergence table
ldf <- dcast(pow_glob_c, N + muc ~ variable, value.var = "conv")
ldf <- ldf[ , -c(6, 8)]
colnames(ldf) <- c('N', '$\\mu_C$', 'LM', '$GLM_{nb}$', '$GLM_{qp}$', '$GLM_{p}$')
print(xtable(ldf, 
             caption = 'Count data simulations - Proportion of models converged. N = sample sizes, 
             $\\mu_C$ = mean abundance in control, LM = Linear model after transformation, 
             $GLM_{nb}$ = negative binomial model, $GLM_{qp}$ = quasi-Poisson model, 
             $GLM_{p}$ = Poisson model', 
             label = 'tab:conv'),
      file = file.path(suppdir, "supp1", "tab", "conv.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})

### ----------------------------------------------------------------------------
### Misc plots
# # XKCD-Plot
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
#   scale_color_manual('', values =  c('#FA6C00', 'gray25'), 
      # guide = guide_legend(override.aes = list(size = 2)))   +
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
# ggsave('/home/edisz/Documents/Uni/Projects/PHD/MISC/AG_presentations/fig/myp2.pdf', 
# myp, width = 10, height = 6)




### ----------------------------------------------------------------------------
### Results -  Binomial data
### Written by Eduard Szöcs
### ----------------------------------------------------------------------------
## load results
res1_p <- readRDS(file.path(cachedir, 'res1_p.rds'))
res2_p <- readRDS(file.path(cachedir, 'res2_p.rds'))


### ----------------------------------------------------------------------------
### Global test 
### ------------------------
# Plot Power
pow_glob_p <- ldply(res1_p, p_glob)
pow_glob_p$pE <- todo_p$pE
pow_glob_p$N <- todo_p$N
pow_glob_p <- melt(pow_glob_p, id.vars = c('pE', 'N'))
levels(pow_glob_p$variable) <- c('lm', 'glm', 'np')

plot_pow_glob_p <- ggplot(pow_glob_p) +
  geom_line(aes(y = value, x = pE, group = variable, linetype = variable)) +
  geom_point(aes(y = value, x = pE, shape = variable), size = 4, color = 'black') +
  facet_grid( ~N, labeller = n_labeller) +  
  # axes
  labs(x = expression(pi[C]), 
       y = expression(paste('Power (global test , ', alpha, ' = 0.05)'))) +
  # appearance
  mytheme +
  scale_shape_manual('Method', values=c(16,2,17), labels = c('LM', expression(GLM[bin]), 'KW')) +
  scale_linetype_discrete('Method', labels = c('LM', expression(GLM[bin]), 'KW')) +
  ylim(c(0,1))
plot_pow_glob_p


### ------------------------
### Tabular output
ldf <- dcast(pow_glob_p, N + pE ~ variable, value.var = "value")
colnames(ldf) <- c('N', '$p_E$', 'LM', '$GLM$', 'np')
print(xtable(ldf, 
             caption = 'Binomial data simulations - Power to detect a global treatment effect. N = sample sizes, 
             $p_E$ = probability in effect treatments, LM = Linear model after transformation, 
             $GLM$ = binomial model, np = Kruskal-Wallis test.', label = 'tab:pow_glob_p'),
      file = file.path(suppdir, "supp1", "tab", "pow_glob_p.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})


### ------------------------
## Plot T1-error
t1_glob_p <- ldply(res2_p, p_glob)
t1_glob_p$ps <- todo_p$pE
t1_glob_p$N <- todo_p$N
t1_glob_p <- melt(t1_glob_p, id.vars = c('ps', 'N'))
levels(t1_glob_p$variable) <- c('lm', 'glm', 'np')

plot_t1_glob_p <- ggplot(t1_glob_p) +
  geom_line(aes(y = value, x = ps, group = variable, linetype = variable)) +
  geom_point(aes(y = value, x = ps, shape = variable), 
             size = 4, color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = 0.05, yend = 0.05), 
               linetype = 'dashed') + 
  facet_grid( ~N, labeller = n_labeller) +  
  # axes
  labs(x =  expression(pi[C]), 
       y = expression(paste('Type 1 error (global test , ', alpha, ' = 0.05)'))) +
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(16,2,17)) +
  scale_linetype_discrete('Method') +
  coord_cartesian(ylim = c(-0.005, 0.1))
plot_t1_glob_p

### ------------------------
### Tabular output
ldf <- dcast(t1_glob_p, N + ps ~ variable, value.var = "value")
colnames(ldf) <- c('N', '$p$', 'LM', '$GLM$', 'np')
print(xtable(ldf, 
             caption = 'Binomial data simulations - Type 1 error to detect a global treatment effect. N = sample sizes, 
             $p$ = probability, LM = Linear model after transformation, 
             $GLM$ = binomial model, np = Kruskal-Wallis test.', label = 'tab:t1_glob_p'),
      file = file.path(suppdir, "supp1", "tab", "t1_glob_p.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})


### ------------------------
### Figure 4
leg <- g_legend(plot_pow_glob_p + 
                  theme(legend.key = element_blank(),
                        legend.text = element_text(size = 16),
                        legend.title = element_text(size = 18, face = "bold"),
                        legend.key.size = unit(1.5, "cm"))
)
p_glob_p <- arrangeGrob(
  plot_t1_glob_p + theme(legend.position="none") + labs(x = NULL),
  plot_pow_glob_p + theme(legend.position="none"),
  leg,
  nrow = 3,
  heights = c(10, 10, 2))
p_glob_p
if(exp_plot)
  ggsave(file.path(figdir, 'p_glob_p.pdf'), p_glob_p, width = 10, height = 8)


### ---------------------------------------------------------------------------
### LOEC
### ------------------------
## Plot Power
pow_loec_p <- ldply(res1_p, p_loec, type = 'power')
pow_loec_p$pE <- todo_p$pE
pow_loec_p$N <- todo_p$N
pow_loec_p  <- melt(pow_loec_p, id.vars = c('pE', 'N'))
levels(pow_loec_p$variable) <- c('lm', 'glm', 'np')

plot_pow_loec_p <- ggplot(pow_loec_p) +
  geom_line(aes(y = value, x = pE, group = variable, linetype = variable)) +
  geom_point(aes(y = value, x = pE, shape = variable), 
             size = 4,  color = 'black') +
  facet_grid( ~N, labeller = n_labeller) + 
  labs(x =  expression(pi[C]), 
       y = expression(paste('Power (LOEC , ', alpha, ' = 0.05)'))) + 
  mytheme +
  scale_shape_manual('Method', values=c(16,2,17), labels = c('LM', expression(GLM[bin]), 'WT')) +
  scale_linetype_discrete('Method', labels = c('LM', expression(GLM[bin]), 'WT')) +
  ylim(c(0,1))
plot_pow_loec_p

### ------------------------
### Tabular output
ldf <- dcast(pow_loec_p, N + pE ~ variable, value.var = "value")
colnames(ldf) <- c('N', '$p_E$', 'LM', '$GLM$', 'np')
print(xtable(ldf, 
             caption = 'Count data simulations - Power to detect LOEC. N = sample sizes, 
             $p_E$ = probability in effect treatments, LM = Linear model after transformation, 
             $GLM$ = binomial model, np = pairwise Wilcoxon.', label = 'tab:pow_loec_p'),
      file = file.path(suppdir, "supp1", "tab", "pow_loec_p.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})

### ------------------------
## Plot T1-error
t1_loec_p <- ldply(res2_p, p_loec, type = 't1')
t1_loec_p$ps <- todo_p$pE
t1_loec_p$N <- todo_p$N
t1_loec_p <- melt(t1_loec_p, id.vars = c('ps', 'N'))
levels(t1_loec_p$variable) <- c('lm', 'glm', 'np')


plot_t1_loec_p <- ggplot(t1_loec_p) +
  geom_line(aes(y = value, x = ps, group = variable, linetype = variable)) +
  geom_point(aes(y = value, x = ps, shape = variable), 
             size = 4, color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = 0.05, yend = 0.05), 
               linetype = 'dashed') + 
  facet_grid( ~N, labeller = n_labeller) + 
  # axes
  labs(x =  expression(pi[C]), 
       y = expression(paste('Type 1 error (LOEC , ', alpha, ' = 0.05)'))) + 
  # appearance
  mytheme + 
  # legend title
  scale_shape_manual('Method', values=c(16,2,17)) +
  scale_linetype_discrete('Method') +
  coord_cartesian(ylim = c(-0.005, 0.1))
plot_t1_loec_p

### ------------------------
### Tabular output
ldf <- dcast(t1_loec_p, N + ps ~ variable, value.var = "value")
colnames(ldf) <- c('N', '$p_E$', 'LM', '$GLM$', 'np')
print(xtable(ldf, 
             caption = 'Binomial data simulations - Type 1 error to detect LOEC. N = sample sizes, 
             $p$ = probability, LM = Linear model after transformation, 
             $GLM$ = binomial model, np = pairwise Wilcoxon.', label = 'tab:t1_loec_p'),
      file = file.path(suppdir, "supp1", "tab", "t1_loec_p.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})


### ------------------------
### Figure 5
leg <- g_legend(plot_pow_loec_p + 
                  theme(legend.key = element_blank(),
                        legend.text = element_text(size = 16),
                        legend.title = element_text(size = 18, face = "bold"),
                        legend.key.size = unit(1.5, "cm"))
)
p_loec_p <- arrangeGrob(
  plot_t1_loec_p + theme(legend.position="none") + labs(x = NULL),
  plot_pow_loec_p + theme(legend.position="none"),
  leg,
  nrow = 3,
  heights = c(10, 10, 2))
p_loec_p
if(exp_plot)
  ggsave(file.path(figdir, 'p_loec_p.pdf'), p_loec_p, width = 10, height = 8)

### Misc
compdf <- dcast(pow_glob_p, pE + N ~ variable)
compdf$diff <- compdf$lm - compdf$glm
compdf [compdf$N == 9, ]
compdf[which.max(abs(compdf$diff)), ]

