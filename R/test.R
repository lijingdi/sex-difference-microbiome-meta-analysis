# Shinichi's code to show the models

# pacakges

library(tidyverse)
library(metafor)
library(orchaRd)
library(here)

# load the data

alpha_es <- read.csv(here("data", "alpha_es_new.csv"))


# adding observation level effect

alpha_es$obs <- 1:nrow(alpha_es)

# model


# the original code

###homoscedastic model
fit_ho <- rma.mv(yi, vi, slab = Authors,random = list(~ 1 | Protocol/Authors/subgroup), mod = ~metric -1,
                 data = alpha_es[alpha_es$stats == "lnCVR",], method = "REML")

summary(fit_ho)

# Shinichi's code to show the models

fit_ho1 <- rma.mv(yi, vi, 
                  random = list(~ 1 | Authors/subgroup, # you cannot put both Authors and Protocls - they are not identifiable (20 and 22 levels)
                                # also need to ask how subgroups is coded - is it correct?
                                ~1 | obs), # is this nested - if not code like this
                  mod = ~metric -1,
                  data = alpha_es[alpha_es$stats == "lnCVR",], 
                  test = "t", # please add this
                  method = "REML",
                  sparse = TRUE) # please add this - make it faster

summary(fit_ho1)

orchard_plot(fit_ho, mod = "metric", group = "Authors", xlab = "lnCVR")
orchard_plot(fit_ho1, mod = "metric", group = "Authors", xlab = "lnCVR")


###heteroscedastic model
fit_he <- rma.mv(yi, vi, slab = Authors,random = list(~1|metric, ~ 1 | Protocol/Authors/subgroup), mod = ~metric -1,
                 data = alpha_es[alpha_es$stats == "lnCVR",], method = "REML")


# Shinichi's code to show the models - hetero model - allow to residual varinace to vary

fit_he1 <- rma.mv(yi, vi, 
                  random = list(~ 1 | Authors/subgroup, # you cannot put both Authors and Protocls - they are not identifiable (20 and 22 levels)
                                # also need to ask how subgroups is coded - is it correct?
                                ~ metric | obs), # is this nested - if not code like this
                  struct = "DIAG",
                  mod = ~metric -1,
                  data = alpha_es[alpha_es$stats == "lnCVR",], 
                  test = "t", # please add this
                  method = "REML",
                  sparse = TRUE) # please add this - make it faster

summary(fit_he1)



orchard_plot(fit_he, mod = "metric", group = "Authors", xlab = "lnCVR")
orchard_plot(fit_he1, mod = "metric", group = "Authors", xlab = "lnCVR")

anova(fit_ho, fit_he) ##no significant difference between the models?
# Shinichi's
anova(fit_ho1, fit_he1) # Hetero is better model

