if(!exists('ld')){
  source("/home/edisz/Documents/Uni/Projects/PHD/6USETHEGLM/src/0-load.R")
}
#####--------------------------------------------------------------------------
# run simulation scripts
source(file.path(srcdir, '1-simulation_counts.R'))
source(file.path(srcdir, '1-simulation_survival.R'))

# Counts
# global
plot_pow_glob_c
plot_t1_glob_c
leg <- g_legend(plot_pow_glob_c + 
                  theme(legend.key = element_blank()) +
                  guides(fill = guide_legend(override.aes = list(size=5)))
                )

p_glob_c <- arrangeGrob(
  arrangeGrob(plot_pow_glob_c + theme(legend.position="none"),
                  plot_t1_glob_c + theme(legend.position="none"),
                  nrow = 2),
  leg, ncol = 2, widths = c(10,2))
p_glob_c

# loec
leg <- g_legend(plot_pow_loec_c + 
                  theme(legend.key = element_blank()) +
                  guides(fill = guide_legend(override.aes = list(size=5)))
)

p_loec_c <- arrangeGrob(
  arrangeGrob(plot_pow_loec_c + theme(legend.position="none"),
              plot_t1_loec_c + theme(legend.position="none"),
              nrow = 2),
  leg, ncol = 2, widths = c(10,2))
p_loec_c
