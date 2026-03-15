# Model comparison for alpha-diversity effect sizes

library(tidyverse)
library(metafor)
library(orchaRd)
library(here)

alpha_es <- read.csv(here("data", "alpha_es_new.csv"))
alpha_es$obs <- seq_len(nrow(alpha_es))

alpha_es_lncvr <- alpha_es %>% filter(stats == "lnCVR")

# Homoscedastic model
fit_ho <- rma.mv(
  yi,
  vi,
  random = list(
    ~1 | Authors/subgroup,
    ~1 | obs
  ),
  mod = ~metric - 1,
  data = alpha_es_lncvr,
  test = "t",
  method = "REML",
  sparse = TRUE
)

# Heteroscedastic model
fit_he <- rma.mv(
  yi,
  vi,
  random = list(
    ~1 | Authors/subgroup,
    ~metric | obs
  ),
  struct = "DIAG",
  mod = ~metric - 1,
  data = alpha_es_lncvr,
  test = "t",
  method = "REML",
  sparse = TRUE
)

summary(fit_ho)
summary(fit_he)

orchaRd::orchard_plot(fit_ho, mod = "metric", group = "Authors", xlab = "lnCVR")
orchaRd::orchard_plot(fit_he, mod = "metric", group = "Authors", xlab = "lnCVR")

anova(fit_ho, fit_he)
