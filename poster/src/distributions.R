require(ggplot2)
require(tikzDevice)

prj <- '/home/edisz/Documents/Uni/Projects/PHD/CONFERENCES/SETAC_Barcelona/presentation/'
figdir <- file.path(prj, 'fig')

# set default theme
mytheme <- theme_set(theme_minimal())
# modify default theme
mytheme <- theme_update(
  axis.text.x = element_text(size=25), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  axis.ticks = element_blank()
)

# annotate size
size <- 8

pth <- '/home/edisz/Documents/Uni/Projects/PHD/CONFERENCES/SETAC2014_Barcelona/poster/setac_poster/fig/'
# Normal distribution
pn <- ggplot(data.frame(x = seq(-1, 5, 0.1), y = dnorm(seq(-1,5,0.1), 2, 1)),
  aes(x=x,y=y)) + 
  geom_area(fill="grey65") +
  geom_line() +
  scale_y_continuous("", breaks=c()) + 
  scale_x_continuous(breaks = seq(-1, 5)) + 
  geom_hline(aes(yintercept = 0)) +
  labs(x = '')
pn
ggsave(file.path(pth, 'pn.svg'), pn)

# Poisson distribution
pp <- ggplot(data.frame(x = 0:5, y = dpois(0:5, 2)),
       aes(x=x,y=y)) + 
  geom_bar(stat = "identity", width = 0.2, fill="grey65") + 
  scale_y_continuous("", breaks=c()) + 
  geom_hline(aes(yintercept = 0)) +
  labs(x = '')
ggsave(file.path(pth, 'pp.svg'), pp)

# Binomial distribution
pb <- ggplot(data.frame(x = 0:8, y = dbinom(0:8, size = 8, 0.8)),
       aes(x=x,y=y)) + 
  geom_bar(stat = "identity", width = 0.3, fill="grey65") + 
  scale_y_continuous("", breaks=c()) + 
  scale_x_continuous(breaks = 0:8)+ 
  geom_hline(aes(yintercept = 0)) +
  labs(x = '')
pb
ggsave(file.path(pth, 'pb.svg'), pb)

