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

### ------------------------
### Power simulations

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

# # plot one realisation of simulated data
# todo1_c[15, ]
# df <- data.frame(x = sims1_c[[15]]$x, y = sims1_c[[15]]$y[ , 2])
# df$yt <- log(1 / min(df$y[df$y != 0]) * df$y + 1)
# dfm <- melt(df)
# levels(dfm$variable) <- c('y', 'ln(Ay + 1)')
# ggplot(dfm, aes(x = x, y = value)) +
#   geom_boxplot(fill = 'grey80') +
#   facet_wrap( ~variable, scales = 'free_y') +
#   scale_x_discrete(labels = c('C', 'T1', 'T2', 'T3', 'T4', 'T5')) +
#   labs(x = 'Treatment', y = 'Abundance') +
#   theme_bw(base_size = 12, 
#            base_family = "Helvetica") +
#   theme(panel.grid.major = element_blank(),
#         text = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14,face = "bold"))

# run methods
if(sim1){
  res1_c <- llply(sims1_c, resfoo1, .progress = 'text', npb = 500)
  saveRDS(res1_c, file.path(cachedir, 'res1_c.rds'))
} 


### ------------------------
# Type1 Error simulations

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

# run methods
if(sim1){
  res2_c <- llply(sims2_c, resfoo1, .progress = 'text', npb = 400)
  saveRDS(res2_c, file.path(cachedir, 'res2_c.rds'))
} 



### ----------------------------------------------------------------------------
### Simulation 2 -  Binomial data
### Written by Eduard Szöcs
### ----------------------------------------------------------------------------

### ------------------------
### Power Simulations

# Settings
# no of simulated datasets
nsims <- 250
# sample sizes
N <- c(3, 6 ,9)
n_animals <- 10
# proportions in effect groups
pEs <- seq(0.6, 0.95, 0.05)
# both as grid
todo1_p <- expand.grid(N = N, pE = pEs)

# create simulate data
sims1_p <- NULL
set.seed(1234)
for(i in seq_len(nrow(todo1_p))){
  sims1_p[[i]] <- dosim2(N = todo1_p[i, 'N'], 
                         pC = 0.95, pE = todo1_p[i, 'pE'], 
                         nsim = nsims, n_animals = n_animals)
}

# # plot one realisation of simulated data
# df <- data.frame(x = sims1_p[[22]]$x, y = sims1_p[[22]]$y[ , 7] / 10)
# ggplot(df, aes(x = x, y = y)) +
#   geom_boxplot(fill = 'grey80') +
#   scale_x_discrete(labels = c('C', 'T1', 'T2', 'T3', 'T4', 'T5')) +
#   labs(x = 'Treatment', y = 'Prop. surv.') +
#   theme_bw(base_size = 12, 
#            base_family = "Helvetica") +
#   theme(panel.grid.major = element_blank(),
#         text = element_text(size=14),
#         axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"))

# run methods
if(sim2){
  res1_p <- llply(sims1_p, resfoo2, .progress = 'text')
  saveRDS(res1_p, file.path(cachedir, 'res1_p.rds'))
}

### ------------------------
# Type1 Error simulations

# Settings
# sample sizes
N <- c(3, 6 ,9)
n_animals <- 10
# proportions 
ps <- seq(0.6, 0.95, 0.05)
# both as grid
todo2_p <- expand.grid(N = N, ps = ps)

# create simulate data
sims2_p <- NULL
for(i in seq_len(nrow(todo2_p))){
  sims2_p[[i]] <- dosim2(N = todo2_p[i, 'N'], 
                         pC = todo2_p[i, 'ps'], pE = todo2_p[i, 'ps'],    # pC = CE
                         nsim = nsims, n_animals = n_animals)
}

# run methods
if(sim2){
  res2_p <- llply(sims2_p, resfoo2, asin = 'ecotox', .progress = 'text')
  saveRDS(res2_p, file.path(cachedir, 'res2_p.rds'))
}
  
