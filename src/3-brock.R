if(!exists("prj")){
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project'!")
} else {
  source(file.path(prj, "src", "0-load.R"))
}
### ----------------------------------------------------------------------------
### Motivating Example - Count data
### Written by Eduard Szöcs
### ----------------------------------------------------------------------------

### -------- Load + Clean ------------------------------------------------------
df <- read.table(file.path(datadir, 'brock.csv'), header = TRUE, sep = ';')
df$Concentration <- factor(df$Concentration)
plot(Abundance ~ Concentration, data = df)

ddply(df, .(Concentration), summarize, mean(Abundance))



### -------- Methods -----------------------------------------------------------
### -----------------------------
### Normal + Transformation
## log(Ax+1) transformation
# Ax = 2 , for min(x) & x != 0
min(df$Abundance[df$Abundance != 0])
A <- 2 / min(df$Abundance[df$Abundance != 0])
min(df$Abundance[df$Abundance != 0]) * A
df <- transform(df, log_Abundance = log(A*Abundance + 1))
plot(log_Abundance ~ Concentration, data = df)

## Fit Normal model
modlm <- lm(log_Abundance ~ Concentration, data = df)

## Inference, global test
# F-test
drop1(modlm, test = 'F')

## Inference, LOEC
# one-sided Dunnett test
summary(glht(modlm, linfct = mcp(Concentration = 'Dunnett'),  
             alternative = 'less'), test = adjusted('holm'))

## predicted mean values with back transformation
modlm_p <- predict(modlm, newdata = data.frame(Concentration = unique(df$Concentration)), 
                 se.fit = TRUE)
# mean (backtransfored)
modlm_fit <- (exp(modlm_p$fit) - 1) / A
# ci
modlm_lwr <- (exp(modlm_p$fit - 1.96*modlm_p$se.fit) - 1) / A
modlm_upr <- (exp(modlm_p$fit + 1.96*modlm_p$se.fit) - 1) / A



### -----------------------------
### Poisson GLM
modpois <- glm(Abundance ~ Concentration, data = df, family = poisson)

# residual deviance >> residual degrees of freedom
modpois$deviance / modpois$df.residual
sum(resid(modpois, type = 'pearson')^2) / modpois$df.residual
# overdispersion!

## Inference, global test
# LR test
drop1(modpois, test = 'Chisq')

## Inference, LOEC
## one-sided Dunnett test
summary(glht(modpois, linfct = mcp(Concentration = 'Dunnett'), 
             alternative = 'less'), test = adjusted('holm'))

## predicted value on response scale
modpois_p <- predict(modpois, newdata = data.frame(Concentration = unique(df$Concentration)), 
                   se.fit = TRUE)
modpois_fit <- modpois$family$linkinv(modpois_p$fit)
modpois_lwr <- modpois$family$linkinv(modpois_p$fit - 1.96 * modpois_p$se.fit)
modpois_upr <- modpois$family$linkinv(modpois_p$fit + 1.96 * modpois_p$se.fit)


### -----------------------------
### quasi-Poisson GLM
modqpois <- glm(Abundance ~ Concentration, data = df, family = quasipoisson)
summary(modqpois)
# # disp same as
# summary(modpois)
# disp <- sum(resid(modpois, type = 'pearson')^2) / modpois$df.residual
# sqrt(diag(vcov(modpois))) * sqrt(disp)

## Inference, global test
# F test
drop1(modqpois, test = 'F')

## Inference, LOEC
# # one-sided Dunnett test
summary(glht(modqpois, linfct = mcp(Concentration = 'Dunnett'),  
             alternative = 'less'), test = adjusted('holm'))

## predicted value on response scale
modqpois_p <- predict(modqpois, newdata = data.frame(Concentration = unique(df$Concentration)), 
                     se.fit = TRUE)
modqpois_fit <- modqpois$family$linkinv(modqpois_p$fit)
modqpois_lwr <- modqpois$family$linkinv(modqpois_p$fit - 1.96 * modqpois_p$se.fit)
modqpois_upr <- modqpois$family$linkinv(modqpois_p$fit + 1.96 * modqpois_p$se.fit)



### -----------------------------
### Negative binomial GLM
modnb <- glm.nb(Abundance ~ Concentration, data = df)
summary(modnb)
# alternatively:
# require(COUNT)
# cc <- nbinomial(Abundance ~ Concentration, data = df)
# summary(cc)

## Inference, global test
# LR test
# refit null model, to estimate theta of null
modnb0 <- glm.nb(Abundance ~ 1, data = df)
anova(modnb, modnb0, test = 'Chisq')

### parametric bootstrap
set.seed(1234)
myPBmodcomp(modnb, modnb0, data = df, npb = 500)



## Inference, LOEC
# # one-sided Dunnett test
summary(glht(modnb, linfct = mcp(Concentration = 'Dunnett'),  
             alternative = 'less'), test = adjusted('holm'))

## predicted value on response scale
modnb_p <- predict(modnb, newdata = data.frame(Concentration = unique(df$Concentration)), 
                     se.fit = TRUE)
modnb_fit <- modnb$family$linkinv(modnb_p$fit)
modnb_lwr <- modnb$family$linkinv(modnb_p$fit - 1.96 * modnb_p$se.fit)
modnb_upr <- modnb$family$linkinv(modnb_p$fit + 1.96 * modnb_p$se.fit)


### -----------------------------
### non-parametric
kruskal.test(df$Abundance, df$Concentration)
pairwise_wilcox(df$Abundance, df$Concentration)



### --------------- Summarize --------------------------------------------------
### -----------------------------
### Mean variance relationships
mv <- ddply(df, .(Concentration), summarise,
      m = mean(Abundance),
      var = var(Abundance))
plot(var ~ m, data = mv, xlab = 'Mean', ylab = 'Variance', pch = 16)
# Poisson
abline(0,1, col = 'darkred', lwd = 2)
# quasi-poisson
abline(0, summary(modqpois)$dispersion, col = 'steelblue', lwd = 2)
# negative binomial
m <- seq(min(mv$m), max(mv$m), 0.01)
lines(m, m + m^2 / modnb$theta, col = 'darkgreen', lwd = 3)
legend(80, 1000, 
       legend = c('Poisson','quasi-Poisson', 'negativ binomial'), 
       lwd = c(2,2,2), col = c('darkred', 'steelblue', 'darkgreen'))



### -----------------------------
### Plot data and model results
moddf <- data.frame(conc = rep(unique(df$Concentration), 4),
                    fit = c(modlm_fit, modpois_fit, modqpois_fit, modnb_fit),
                    upr = c(modlm_upr, modpois_upr, modqpois_upr, modnb_upr),
                    lwr = c(modlm_lwr, modpois_lwr, modqpois_lwr, modnb_lwr),
                    mod = rep(c('Normal', 'Poisson', 'quasi-Poisson', 
                                'Negative binomial'), each = 6), 
                    strindsAsFactor = FALSE)
moddf$mod <- factor(moddf$mod, levels = c('Normal', 'Poisson', 'quasi-Poisson', 
                                          'Negative binomial'))

p <- ggplot() +
  # raw data
#   geom_boxplot(data = df, aes(x = Concentration, y = Abundance), 
#                width = 0.3, fill = 'grey85') +
  geom_point(data = df, aes(x = Concentration, y = Abundance), size = 4) +
  # estimated means + CI
  geom_point(data = moddf, 
             aes(x = 0.1 + rep(1:6,4) + rep(c(0.1, 0.2, 0.3, 0.4), each = 6), 
                 y = fit), size = 4, shape = 17) +
  geom_errorbar(data = moddf, 
                aes(x = 0.1 + rep(1:6,4) + rep(c(0.1, 0.2, 0.3, 0.4), each = 6), 
                    ymax = upr, ymin = lwr, linetype = mod), 
                width = 0.1, lwd = 1) + 
  # loec_bars
  # normal
  geom_path(aes(x = c(5-0.15, 5+0.55), y = c(220,220)), 
            linetype = 'solid', size = 2) + 
  # pois
  geom_path(aes(x = c(2-0.15, 6+0.55), y = c(220-5,220-5)), 
            linetype = '1111', size = 2) + 
  # qpois 
  geom_path(aes(x = c(5-0.15, 5+0.55), y = c(220-2*5,220-2*5)), 
            linetype = 'dashed', size = 2) +
  # neg bin
  geom_path(aes(x = c(3-0.15, 3+0.55), y = c(220-3*5,220-3*5)), 
            linetype = 'dotdash', size = 2) + 
  geom_path(aes(x = c(5-0.15, 6+0.55), y = c(220-3*5,220-3*5)), 
            linetype = 'dotdash', size = 2) + 
  scale_linetype_manual(values=c('solid', '1111', 'dashed', 'dotdash')) + 
  mytheme +
  guides(linetype=guide_legend(title='Model')) +
  xlab('Concentration [mg/L]') +
  theme(axis.title.x = element_text(face="bold",size=18),
        axis.text.x  = element_text(size=14),
        axis.text.y  = element_text(size=14),
        axis.title.y = element_text(face="bold",size=18),
        legend.position=c(0.8, 0.7),
        legend.key.width = unit(1, "cm"),
        legend.background = element_rect(color = 'black')
        )
p
if(exp_plot)
  ggsave(file.path(figdir, 'example.pdf'), p, height = 6, width = 8)

