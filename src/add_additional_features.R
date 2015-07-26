################################################################################################
################################################################################################
####            Code for deriving additional features utilized in the analysis             #####
####                                     Team JayHawks                                      ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
## add_additional_features: Derive additional features by summing and collapsing original     ##
##                          features from clean_core_data.                                    ##
##                                                                                            ##
################################################################################################


add_additional_features <- function(fileName, data){
  input_data <- read.csv(fileName)
  out_data <- data
  
  out_data$AGE <- as.numeric(out_data$AGEGRP2)
  out_data$ECOG <- as.numeric(out_data$ECOG_C)
  out_data$DISEASE_SITE_NEW = ifelse(out_data$DISEASE_SITE %in% c(1,2), "G12", "G3")
  out_data$RACE_NEW = ifelse(out_data$RACE_C %in% c("White", "Other"), "White/Other", "NonWhite")
  
  sub_dat <- out_data
  ##sub_dat$metastases <- as.matrix(make_metastases_data( input_data )[,c((1:8)[-7],18)])
  ##sum_weights=c(3.0665884,2.0926818,1.3665459,0.4740491,4.3143560,1.1118549,3.0536851,1.2068681)
  sub_dat$metastases <- as.matrix(make_metastases_data( input_data )[, c("BONE" , "LYMPH_NODES", "LUNGS", "PROSTATE", "LIVER", "ADRENAL", "OTHER")])
  sum_weights=c(2.8720038, 2.1705066 , 0.8278407 ,-0.3270099  ,4.4245667 ,2.9434344 ,0.7459380 )
  #names(sum_weights)=c("metastases.BONE","metastases.LYMPH_NODES","metastases.LUNGS","metastases.PROSTATE","metastases.LIVER","metastases.OTHER","metastases.ADRENAL","metastases.DISEASE_SITE")
  names(sum_weights)=c("metastases.BONE","metastases.LYMPH_NODES","metastases.LUNGS","metastases.PROSTATE","metastases.LIVER", "metastases.ADRENAL","metastases.OTHER")
  sub_dat$metastases_sum <- sub_dat$metastases %*% sum_weights
  out_data$metastases_sum = sub_dat$metastases_sum
  
  sum_weights2 <- c(0.1678374, -0.4059502)
  names(sum_weights2) <- c("sg_harmful","sg_protective")
  out_data$harm_pro <- as.matrix(out_data[names(sum_weights2)]) %*% sum_weights2
  
  sum_weights3 <- c(2.883853, -2.140534)
  names(sum_weights3) <- c("sg_harmful","sg_protective")
  out_data$harm_pro2 <- as.matrix(out_data[names(sum_weights3)]) %*% sum_weights3
  
  # hardcoded weights of the first principle component from the principle component analysis
  pc1_weights <- c(0.4660384, -0.2602531, -0.3944890, 0.5011690, 0.3759285, 0.4086103)
  names(pc1_weights) <- c("rm_LDH", "rm_ALB", "rm_HB", "rm_ALP", "rm_PSA", "rm_AST")
  # add pc1 as a new feature
  out_data$pc1 <- as.matrix(out_data[names(pc1_weights)]) %*% pc1_weights
  
  # add disease history derived features
  out_data$MetabolicSyndrome = ifelse((out_data$rr_DIAB == "YES" & (out_data$sg_BETA_BLOCKING == 1 | out_data$sg_ACE_INHIBITORS == 1)), 1, 0)
  out_data$Obesity = ifelse(out_data$BMI >=25, 1, 0)
  out_data$StatinObese = ifelse(out_data$sg_HMG_COA_REDUCT == 1 & out_data$Obesity == 1, 1, 0)
  out_data$StainMetabolicSyndrome = ifelse(out_data$sg_HMG_COA_REDUCT == 1 & out_data$MetabolicSyndrome == 1, 1, 0)
  out_data$MetforminObese = ifelse(out_data$rr_metforminInfo == 1 & out_data$Obesity == 1, 1, 0)
  out_data$MetforminMetabolicSyndrome = ifelse(out_data$rr_metforminInfo == 1 & out_data$MetabolicSyndrome == 1, 1, 0)
  out_data$StatinMetforminObese = ifelse(out_data$sg_HMG_COA_REDUCT == 1 & out_data$rr_metforminInfo == 1 & out_data$Obesity == 1, 1, 0)
  out_data$StatinMetforminMetabolicSyndrome = ifelse(out_data$sg_HMG_COA_REDUCT == 1 & out_data$rr_metforminInfo == 1 & out_data$MetabolicSyndrome == 1, 1, 0)
  
  return(out_data)
}