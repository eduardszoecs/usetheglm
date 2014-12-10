if(!exists('ld')){
  source("/home/edisz/Documents/Uni/Projects/PHD/6USETHEGLM/src/0-load.R")
}


### --------------------------------------------
# data to nice graphs
# Counts
source(file.path(srcdir, '1-simulation_counts.R'))
# global
plot_pow_glob_c
plot_t1_glob_c
leg <- g_legend(plot_pow_glob_c + 
                  theme(legend.key = element_blank())
                )
plot(leg)

p_glob_c <- arrangeGrob(
  plot_pow_glob_c + theme(legend.position="none") + labs(x = NULL),
  plot_t1_glob_c + theme(legend.position="none"),
  leg,
  nrow = 3, 
  heights = c(10, 10, 2))
p_glob_c
# ggsave(file.path(figdir, 'p_glob_c.pdf'), p_glob_c, width = 14, height = 8)
# ggsave(file.path(figdir, 'p_glob_c.jpeg'), p_glob_c, width = 14, height = 8)

# loec
leg <- g_legend(plot_pow_loec_c + 
                  theme(legend.key = element_blank()))
p_loec_c <- arrangeGrob(
  plot_pow_loec_c + theme(legend.position="none") + labs(x = NULL),
  plot_t1_loec_c + theme(legend.position="none"),
  leg,
  nrow = 3,
  heights = c(10, 10, 2))
p_loec_c
# ggsave(file.path(figdir, 'p_loec_c.pdf'), p_loec_c, width = 14, height = 8)


### --------------------------------------------
# Binomial data
source(file.path(srcdir, '1-simulation_survival.R'))
# global
plot_pow_glob_p
plot_t1_glob_p
leg <- g_legend(plot_pow_glob_p + 
                  theme(legend.key = element_blank()) +
                  guides(fill = guide_legend(override.aes = list(size=5)))
)
p_glob_p <- arrangeGrob(
  plot_pow_glob_p + theme(legend.position="none") + labs(x = NULL),
  plot_t1_glob_p + theme(legend.position="none"),
  leg,
  nrow = 3,
  heights = c(10, 10, 2))
p_glob_p
ggsave(file.path(figdir, 'p_glob_p.pdf'), p_glob_p, width = 14, height = 8)


# loec
leg <- g_legend(plot_pow_loec_p + 
                  theme(legend.key = element_blank()) +
                  guides(fill = guide_legend(override.aes = list(size=5)))
)
p_loec_p <- arrangeGrob(
  plot_pow_loec_p + theme(legend.position="none") + labs(x = NULL),
  plot_t1_loec_p + theme(legend.position="none"),
  leg,
  nrow = 3,
  heights = c(10, 10, 2))
p_loec_p
ggsave(file.path(figdir, 'p_loec_p.pdf'), p_loec_p, width = 14, height = 8)

