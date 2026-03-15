##### Functions for sex-difference meta-analysis #####

lnCVR2 <- function(mod.list, CMean, CSD, CN, EMean, ESD, EN) {
  yi <- log((ESD / EMean) / (CSD / CMean)) + (1 / (2 * (EN - 1))) - (1 / (2 * (CN - 1)))
  yi <- yi + 0.5 * (((CSD^2) / (CN * CMean^2)) - ((ESD^2) / (EN * EMean^2)))

  vi <- CSD^2 / (CN * CMean^2) + CSD^4 / (2 * CN^2 * CMean^4) +
    1 / (2 * (CN - 1)) + 1 / (2 * (CN - 1)^2) +
    ESD^2 / (EN * EMean^2) + ESD^4 / (2 * EN^2 * EMean^4) +
    1 / (2 * (EN - 1)) + 1 / (2 * (EN - 1)^2)

  data.frame(mod.list, CN, EN, yi, vi)
}

lnVR2 <- function(mod.list, CMean, CSD, CN, EMean, ESD, EN) {
  yi <- log(ESD / CSD) + (1 / (2 * (EN - 1))) - (1 / (2 * (CN - 1)))
  vi <- 1 / (2 * (CN - 1)) + 1 / (2 * (CN - 1)^2) +
    1 / (2 * (EN - 1)) + 1 / (2 * (EN - 1)^2)

  data.frame(mod.list, CN, EN, yi, vi)
}

lnRR2 <- function(mod.list, CMean, CSD, CN, EMean, ESD, EN) {
  yi <- log(EMean / CMean) +
    0.5 * (((ESD^2) / (EN * EMean^2)) - ((CSD^2) / (CN * CMean^2)))

  vi <- CSD^2 / (CN * CMean^2) + CSD^4 / (2 * CN^2 * CMean^4) +
    ESD^2 / (EN * EMean^2) + ESD^4 / (2 * EN^2 * EMean^4)

  data.frame(mod.list, CN, EN, yi, vi)
}

calculate_ES <- function(dat, mod) {
  stopifnot(all(c("metric", "mean_male", "sd_male", "n_male", "mean_female", "sd_female", "n_female") %in% names(dat)))

  out <- lapply(unique(dat$metric), function(metric_name) {
    dat_sub <- dat[dat$metric == metric_name, , drop = FALSE]

    cvr <- lnCVR2(dat_sub[mod], dat_sub$mean_male, dat_sub$sd_male, dat_sub$n_male, dat_sub$mean_female, dat_sub$sd_female, dat_sub$n_female)
    vr <- lnVR2(dat_sub[mod], dat_sub$mean_male, dat_sub$sd_male, dat_sub$n_male, dat_sub$mean_female, dat_sub$sd_female, dat_sub$n_female)
    rr <- lnRR2(dat_sub[mod], dat_sub$mean_male, dat_sub$sd_male, dat_sub$n_male, dat_sub$mean_female, dat_sub$sd_female, dat_sub$n_female)

    dplyr::bind_rows(
      dplyr::mutate(cvr, stats = "lnCVR"),
      dplyr::mutate(vr, stats = "lnVR"),
      dplyr::mutate(rr, stats = "lnRR")
    )
  })

  dplyr::bind_rows(out)
}

calculate_ES_treatment <- function(dat, mod) {
  stopifnot(all(c("metric", "mean_control", "sd_control", "n_control", "mean_treatment", "sd_treatment", "n_treatment") %in% names(dat)))

  out <- lapply(unique(dat$metric), function(metric_name) {
    dat_sub <- dat[dat$metric == metric_name, , drop = FALSE]

    cvr <- lnCVR2(dat_sub[mod], dat_sub$mean_control, dat_sub$sd_control, dat_sub$n_control, dat_sub$mean_treatment, dat_sub$sd_treatment, dat_sub$n_treatment)
    vr <- lnVR2(dat_sub[mod], dat_sub$mean_control, dat_sub$sd_control, dat_sub$n_control, dat_sub$mean_treatment, dat_sub$sd_treatment, dat_sub$n_treatment)
    rr <- lnRR2(dat_sub[mod], dat_sub$mean_control, dat_sub$sd_control, dat_sub$n_control, dat_sub$mean_treatment, dat_sub$sd_treatment, dat_sub$n_treatment)

    dplyr::bind_rows(
      dplyr::mutate(cvr, stats = "lnCVR"),
      dplyr::mutate(vr, stats = "lnVR"),
      dplyr::mutate(rr, stats = "lnRR")
    )
  })

  dplyr::bind_rows(out)
}

get_lower_tri <- function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  cormat
}

get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  cormat
}

theme_MicrobeR <- function() {
  ggplot2::theme_classic(base_size = 20, base_family = "Helvetica") +
    ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", linewidth = 1, fill = NA)) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), strip.background = ggplot2::element_blank())
}
