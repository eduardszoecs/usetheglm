if(!exists("prj")){
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project'!")
} else {
  source(file.path(prj, "src", "0-load.R"))
}

### ----------------------------------------------------------------------------
### Compile results
### Written by Eduard Szöcs
### ----------------------------------------------------------------------------

#####--------------------------------------------------------------------------
# Count simulations

# run simulations
source(file.path(srcdir, '1-simulation_counts.R'))
# global test
plot_pow_glob_c
plot_t1_glob_c
leg <- g_legend(plot_pow_glob_c + 
                  theme(legend.key = element_blank(),
                        legend.text = element_text(size = 16),
                        legend.title = element_text(size = 18, face = "bold"),
                        legend.key.size = unit(1.5, "cm"))
                )
plot(leg)
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

# loec
plot_pow_loec_c
plot_t1_loec_c 
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


## --------- Tabular output ----------------------------------------------------

# Convergence table

ldf <- dcast(pow_glob_c, N + muc ~ variable, value.var = "conv")
ldf <- ldf[ , -c(6, 7)]
colnames(ldf) <- c('N', '$\\mu_C$', 'LM', '$GLM_{nb}$', '$GLM_{qp}$')
print(xtable(ldf, 
             caption = 'Count data simulations - Proportion of models converged. N = sample sizes, 
             $\\mu_C$ = mean abundance in control, LM = Linear model after transformation, 
             $GLM_{nb}$ = negative binomial model, $GLM_{qp}$ = quasi-Poisson model.', label = 'tab:conv'),
      file = file.path(prj, "supplement", "supp1", "tab", "conv.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})

# power table - global test
ldf <- dcast(pow_glob_c, N + muc ~ variable, value.var = "power")
colnames(ldf) <- c('N', '$\\mu_C$', 'LM', '$GLM_{nb}$', '$GLM_{qp}$', '$GLM_{pb}$', 'np')
print(xtable(ldf, 
             caption = 'Count data simulations - Power to detect a global treatment effect. N = sample sizes, 
             $\\mu_C$ = mean abundance in control, LM = Linear model after transformation, 
             $GLM_{nb}$ = negative binomial model, $GLM_{qp}$ = quasi-Poisson model, 
             $GLM_{pb}$ = negative binomial model with parametric boostrap, np = Kruskal-Wallis test.', label = 'tab:pow_glob_c'),
      file = file.path(prj, "supplement", "supp1", "tab", "pow_glob_c.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})

# power table - loec
ldf <- dcast(pow_loec_c, N + muc ~ variable, value.var = "power")
colnames(ldf) <- c('N', '$\\mu_C$', 'LM', '$GLM_{nb}$', '$GLM_{qp}$', 'np')
print(xtable(ldf, 
             caption = 'Count data simulations - Power to detect LOEC. N = sample sizes, 
             $\\mu_C$ = mean abundance in control, LM = Linear model after transformation, 
             $GLM_{nb}$ = negative binomial model, $GLM_{qp}$ = quasi-Poisson model, 
            np = pairwise Wilcoxon test.', label = 'tab:pow_loec_c'),
      file = file.path(prj, "supplement", "supp1", "tab", "pow_loec_c.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})

# type 1 error - global test
ldf <- dcast(t1_glob_c, N + muc ~ variable, value.var = "power")
ldf <- ldf[ , c(1,2,3,4,6,5,7)]
colnames(ldf) <- c('N', '$\\mu_C$', 'LM', '$GLM_{nb}$', '$GLM_{qp}$', '$GLM_{pb}$', 'np')
print(xtable(ldf, 
             caption = 'Count data simulations - Type 1 error to detect a global treatment effect. N = sample sizes, 
             $\\mu_C$ = mean abundance in control, LM = Linear model after transformation, 
             $GLM_{nb}$ = negative binomial model, $GLM_{qp}$ = quasi-Poisson model, 
             $GLM_{pb}$ = negative binomial model with parametric boostrap, np = Kruskal-Wallis test.', label = 'tab:t1_glob_c'),
      file = file.path(prj, "supplement", "supp1", "tab", "t1_glob_c.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})


# type 1 error - loec
ldf <- dcast(t1_loec_c, N + muc ~ variable, value.var = "t1")
colnames(ldf) <- c('N', '$\\mu_C$', 'LM', '$GLM_{nb}$', '$GLM_{qp}$', 'np')
print(xtable(ldf, 
             caption = 'Count data simulations - Type 1 error to detect LOEC. N = sample sizes, 
             $\\mu_C$ = mean abundance in control, LM = Linear model after transformation, 
             $GLM_{nb}$ = negative binomial model, $GLM_{qp}$ = quasi-Poisson model, 
            np = pairwise Wilcoxon.', label = 'tab:t1_loec_c'),
      file = file.path(prj, "supplement", "supp1", "tab", "t1_loec_c.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})

#####--------------------------------------------------------------------------
# Binomial simulations
source(file.path(srcdir, '1-simulation_survival.R'))
# global
plot_pow_glob_p
plot_t1_glob_p
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


# loec
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



## --------- Tabular output ----------------------------------------------------
# power table - global test
ldf <- dcast(pow_glob_p, N + pE ~ variable, value.var = "value")
colnames(ldf) <- c('N', '$p_E$', 'LM', '$GLM$', 'np')
print(xtable(ldf, 
             caption = 'Binomial data simulations - Power to detect a global treatment effect. N = sample sizes, 
             $p_E$ = probability in effect treatments, LM = Linear model after transformation, 
             $GLM$ = binomial model, np = Kruskal-Wallis test.', label = 'tab:pow_glob_p'),
      file = file.path(prj, "supplement", "supp1", "tab", "pow_glob_p.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})

# power table - loec
ldf <- dcast(pow_loec_p, N + pE ~ variable, value.var = "value")
colnames(ldf) <- c('N', '$p_E$', 'LM', '$GLM$', 'np')
print(xtable(ldf, 
             caption = 'Count data simulations - Power to detect LOEC. N = sample sizes, 
             $p_E$ = probability in effect treatments, LM = Linear model after transformation, 
             $GLM$ = binomial model, np = pairwise Wilcoxon.', label = 'tab:pow_loec_p'),
      file = file.path(prj, "supplement", "supp1", "tab", "pow_loec_p.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})

# type 1 error - global test
ldf <- dcast(t1_glob_p, N + ps ~ variable, value.var = "value")
colnames(ldf) <- c('N', '$p$', 'LM', '$GLM$', 'np')
print(xtable(ldf, 
             caption = 'Binomial data simulations - Type 1 error to detect a global treatment effect. N = sample sizes, 
             $p$ = probability, LM = Linear model after transformation, 
             $GLM$ = binomial model, np = Kruskal-Wallis test.', label = 'tab:t1_glob_p'),
      file = file.path(prj, "supplement", "supp1", "tab", "t1_glob_p.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})


# type 1 error - loec
ldf <- dcast(t1_loec_p, N + ps ~ variable, value.var = "value")
colnames(ldf) <- c('N', '$p_E$', 'LM', '$GLM$', 'np')
print(xtable(ldf, 
             caption = 'Binomial data simulations - Type 1 error to detect LOEC. N = sample sizes, 
             $p$ = probability, LM = Linear model after transformation, 
             $GLM$ = binomial model, np = pairwise Wilcoxon.', label = 'tab:t1_loec_p'),
      file = file.path(prj, "supplement", "supp1", "tab", "t1_loec_p.tex"), 
      table.placement = 'H', caption.placement = 'top', size = 'footnotesize',
      include.rownames = FALSE, sanitize.text.function=function(x){x})
