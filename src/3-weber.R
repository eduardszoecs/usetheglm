### ----------------------------------------------------------------------------
### Motivating Example - Binomial data
### Written by Eduard Szöcs

### -------- Packages ----------------------------------------------------------
require(reshape2)
require(multcomp)


### -------- Load + Clean ------------------------------------------------------
df <- read.table(header = TRUE, text = 'conc A B C D
0 1 1 0.9 0.9
32 0.8 0.8 1 0.8
64 0.9 1 1 1 
128 0.9 0.9 0.8 1
256 0.7 0.9 1 0.5
512 0.4 0.3 0.4 0.2')
df

dfm <- melt(df, id.vars = 'conc', value.name = 'y', variable.name = 'tank')
# conc as factor
dfm$conc <- factor(dfm$conc)
head(dfm)

boxplot(y ~ conc, data = dfm, 
        xlab = 'conc', ylab = 'Proportion surv.')



### -------- Methods -----------------------------------------------------------
### -------- Normal + Transformation
dfm$y_asin <- ifelse(dfm$y == 1, 
                     asin(sqrt(dfm$y)) - asin(sqrt(1/10)), 
                     asin(sqrt(dfm$y)) 
)

boxplot(y_asin ~ conc, data = dfm, 
        xlab = 'conc', ylab = 'Proportion surv.')

modlm <- lm(y_asin ~ conc, data = dfm)
summary(modlm)
plot(fitted(modlm), residuals(modlm)); abline(h = 0)

# F-test
drop1(modlm, test = 'F')
# LOEC
# one-sided dunnett
summary(glht(modlm, linfct = mcp(conc = 'Dunnett'), alternative = 'less'))


### -------- Binomial GLM
modglm <- glm(y ~ conc , data = dfm, family = binomial(link = 'logit'), weights = rep(10, nrow(dfm)))
modglm.null <- glm(y ~ 1 , data = dfm, family = binomial, weights = rep(10, nrow(dfm)))
summary(modglm)
# LR-test
drop1(modglm, test = 'Chisq')
# one-sided Dunnett contrasts
summary(glht(modglm, linfct = mcp(conc = 'Dunnett'), alternative = 'less'))


