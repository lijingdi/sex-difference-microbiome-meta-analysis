#####functions for sex-difference meta-analysis#####
# Create lnCVR2 function, higher lnCRV2 means higher variability among females
#log-transformed ratio of two coefficients of variation 
lnCVR2 <- function(mod.list, CMean, CSD, CN, EMean, ESD, EN)
{
  yi<-log((ESD / EMean) / (CSD / CMean)) + (1 / (2 * (EN - 1))) - (1 / (2 * (CN - 1)))
  yi<-yi + 0.5 * (((CSD^2) / (CN * CMean^2)) - ((ESD^2) / (EN * EMean^2)))	
  vi<-CSD^2 / (CN * CMean^2) + CSD^4 / (2 * CN^2 * CMean^4) +
    1 / (2 * (CN - 1)) + 1 / (2 * (CN - 1)^2) +
    ESD^2 / (EN * EMean^2) + ESD^4 / (2 * EN^2 * EMean^4) +
    1 / (2 * (EN - 1)) + 1 / (2 * (EN - 1)^2)
  result_df <- data.frame(mod.list, CN, EN, yi,vi)
#  result_df <- cbind(result_df, mod.list)
  return(result_df)
}

#lnVR2 (better estimator)
#log-tranformed ratio of variability, likely a more useful measure where a mean-variance relationship is likely to exist
lnVR2 <- function(mod.list, CMean, CSD, CN, EMean, ESD, EN){
  yi<-log(ESD / CSD) + (1 / (2 * (EN - 1))) - (1 / (2 * (CN - 1))) 
  vi<- 1 / (2 * (CN - 1)) + 1 / (2 * (CN - 1)^2) + 
    1 / (2 * (EN - 1)) + 1 / (2 * (EN - 1)^2)
  result_df <- data.frame(mod.list, CN, EN, yi,vi)
#  result_df <- cbind(result_df, mod.list)
  
  return(result_df)
}

# Create lnRR2 function 
lnRR2 <- function(mod.list, CMean, CSD, CN, EMean, ESD, EN){
  yi<-log( EMean / CMean) + 
    0.5 * (((ESD^2) / (EN * EMean^2)) - ((CSD^2) / (CN * CMean^2)))	
  vi <- CSD^2 / (CN * CMean^2) + CSD^4 / (2 * CN^2 * CMean^4)  + 
    ESD^2 / (EN * EMean^2) + ESD^4 / (2 * EN^2 * EMean^4) 
  result_df <- data.frame(mod.list,CN, EN, yi,vi)
#  result_df <- cbind(result_df, mod.list)
  
  return(result_df)
}

# this function computes corrected lnCVR, lnVR and lnRR to generate effect size for alpha-diversity
# lnRR2 (better estimator)
calculate_ES <- function(dat, mod, metric=c("Chao1.Chao1", "Evenness", "Observed", "PD",
                                                                     "Shannon", "Simpson")) {
  mydata_CVR <- data.frame(CN=numeric(0), EN=numeric(0),
                       yi=numeric(0), vi=numeric(0), Authors=character(0), metric=character(0))
  mydata_VR <- data.frame(CN=numeric(0), EN=numeric(0),
                          yi=numeric(0), vi=numeric(0), Authors=character(0), metric=character(0))
  mydata_RR <- data.frame(CN=numeric(0), EN=numeric(0),
                          yi=numeric(0), vi=numeric(0), Authors=character(0), metric=character(0))
  for (i in unique(dat$metric)) {
      dat_sub <- dat[dat$metric == i,]
      add_CVR <- lnCVR2(dat_sub[mod], dat_sub$mean_male, dat_sub$sd_male, dat_sub$n_male, dat_sub$mean_female, dat_sub$sd_female, dat_sub$n_female)
      add_VR <- lnVR2(dat_sub[mod], dat_sub$mean_male, dat_sub$sd_male, dat_sub$n_male, dat_sub$mean_female, dat_sub$sd_female, dat_sub$n_female)
      add_RR <- lnRR2(dat_sub[mod], dat_sub$mean_male, dat_sub$sd_male, dat_sub$n_male, dat_sub$mean_female, dat_sub$sd_female, dat_sub$n_female)
      mydata_CVR <- rbind(mydata_CVR, add_CVR)
      mydata_VR <- rbind(mydata_VR, add_VR)
      mydata_RR <- rbind(mydata_RR, add_RR)
  }
  mydata_CVR <- mydata_CVR %>% mutate(stats = "lnCVR")
  mydata_VR <- mydata_VR %>% mutate(stats = "lnVR")
  mydata_RR <- mydata_RR %>% mutate(stats = "lnRR")
  mydata_all <- do.call("bind_rows", list(mydata_CVR, mydata_VR, mydata_RR))
  return(mydata_all)
}


# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

theme_MicrobeR <- function () { 
  theme_classic(base_size=20, base_family="Helvetica") +
    theme(panel.border = element_rect(color="black", size=1, fill=NA)) +
    theme(axis.line = element_blank(), strip.background = element_blank())
}


