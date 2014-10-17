if(!exists("prj")){
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project'!")
} else {
  source(file.path(prj, "src", "0-load.R"))
}

#####--------------------------------------------------------------------------
### Simulation 3 - Shapiro-Wilk Normality Test (count data)
nsims <- 100
N <- c(2, 3, 4, 6, 8, 10, 20)
mu <-  2^(c(1:5, 7, 9))
todo <- expand.grid(N = N, mu = mu)


#####------------------------------------
# simulate data
sims3_c <- NULL
set.seed(seed)
for(i in seq_len(nrow(todo))){
  sims3_c[[i]] <- dosim3(N = todo[i, 'N'], mu = todo[i, 'mu'], nsim = nsims)
}

sims3_c[[1]]


#####------------------------------------
# analyse simulations
if(sim3){
  res3_c <- ldply(sims3_c, resfoo3, .progress = 'text')
  saveRDS(res3_c, file.path(cachedir, 'res3_c.rds'))
} else {
  res3_c <- readRDS(file.path(cachedir, 'res3_c.rds'))
}


#####------------------------------------
# Results
res <- cbind(todo, res3_c)
res <- melt(res, id.vars = names(todo))

ggplot(res, aes(x = mu, y = value, col = variable, group = variable)) +
  geom_point() +
  geom_line() +
  facet_wrap(~N) +
  coord_trans(xtrans = 'log2') +
  scale_x_continuous(breaks = round(unique(todo$mu), 0)) +
  theme_bw(base_size = 12, 
           base_family = "Helvetica") +
  geom_segment(aes(x = 2, xend = 1024, y = 0.05, yend = 0.05), 
               linetype = 'dashed', col = 'black') 


