suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(stringr)
})

fits <- readRDS("results/meta_model_fits.rds")

extract_meta_model_results <- function(fit_obj, domain_label = NULL, stats_label = NULL) {
  coef_mat <- as.matrix(fit_obj$b)
  term_names <- rownames(coef_mat)
  if (is.null(term_names)) {
    term_names <- paste0("term_", seq_len(nrow(coef_mat)))
  }

  out <- tibble(
    term = term_names,
    metric = stringr::str_remove(term_names, "^metric"),
    estimate = as.numeric(coef_mat[, 1]),
    se = as.numeric(fit_obj$se),
    ci_lb = as.numeric(fit_obj$ci.lb),
    ci_ub = as.numeric(fit_obj$ci.ub),
    pval = as.numeric(fit_obj$pval)
  )

  stat_vec <- suppressWarnings(as.numeric(out$estimate / out$se))
  stat_name <- if (!is.null(fit_obj$test) && fit_obj$test %in% c("t", "knha", "hksj")) "tval" else "zval"
  out[[stat_name]] <- stat_vec

  if (!is.null(domain_label)) out$domain <- domain_label
  if (!is.null(stats_label)) out$stats <- stats_label

  out %>% dplyr::select(any_of(c("domain", "stats")), dplyr::everything())
}

alpha_model_results <- bind_rows(lapply(
  names(fits$alpha_fits),
  function(stat_name) extract_meta_model_results(fits$alpha_fits[[stat_name]], "alpha", stat_name)
))

beta_model_results <- extract_meta_model_results(fits$beta_fit, "beta", "lnRR")

write.csv(alpha_model_results, "results/alpha_meta_model_results.csv", row.names = FALSE)
write.csv(beta_model_results, "results/beta_meta_model_results.csv", row.names = FALSE)

cat("rows alpha:", nrow(alpha_model_results), "\n")
cat("rows beta:", nrow(beta_model_results), "\n")
