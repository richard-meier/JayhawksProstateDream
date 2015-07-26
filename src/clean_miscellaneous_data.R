################################################################################################
################################################################################################
####       Code to process Response and Baseline/Demographic data from Core Table           ####
####                                     Joseph Usset                                       ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
## Input: Training, Test, or Validation Data 												  ##
## Output: Response and Numberic/imputed Baseline/Demographic variables 					  ##
##                                                                                            ##
## New variables:                                    										  ##
## 				Total: sum of all metastases locations per patient                            ##
##              Visceral: sum of visceral metastases                    					  ##
##              Other: sum of rare metastases (not bone, lymphnode, lungs prostate liver)     ##
##				Disease Site variable: Same as Halabi 										  ##	
################################################################################################

### FUNCTION TO CLEAN THE DATA!!!! ###
clean_miscellaneous_data <- function( fileName ){
  options(stringsAsFactors=TRUE)
  
  input_data <- read.csv( fileName )
  
  ### Extract outcomes for sub-challenge 1 and 2, as well as baseline/demographic or categorical variables considered by Halabi ###
  base_ids=input_data[,c("STUDYID", "RPT", "LKADT_REF", "LKADT_PER", "DEATH", "DISCONT", "AGEGRP2", "ECOG_C", "ANALGESICS", "RACE_C", "PRIOR_RADIOTHERAPY")]
  
  ### Clean the baseline/non-laboratory data ###
  base_non_lab_data <- input_data[,c("LKADT_P", "BMI")]
  base_non_lab_data$LKADT_P <- as.character(as.matrix(base_non_lab_data$LKADT_P))
  base_non_lab_data$BMI <- as.character(as.matrix(base_non_lab_data$BMI))
  
  #### Get rid of dots ##3
  base_non_lab_data <- reformatTableColumnsToProperNumericVectors(base_non_lab_data)
  
  #### Add back IDs
  base_non_lab_data <- data.frame(base_ids, base_non_lab_data)
  
  ### Make all factors numeric ###
  base_non_lab_data$AGEGRP2 <- as.numeric(base_non_lab_data$AGEGRP2) - 1
  base_non_lab_data$ANALGESICS <- as.numeric(base_non_lab_data$ANALGESICS) - 1
  suppressWarnings(	# this introduces NAs for "." charactersm which will be imputed later
	base_non_lab_data$ECOG_C <- as.numeric(as.character(base_non_lab_data$ECOG_C))
  )
  base_non_lab_data$ECOG_C[which(base_non_lab_data$ECOG_C == 3)] <- 2
  base_non_lab_data$DEATH <- as.numeric(as.factor(base_non_lab_data$DEATH)) - 1
  base_non_lab_data$PRIOR_RADIOTHERAPY <- as.numeric(base_non_lab_data$PRIOR_RADIOTHERAPY)-1
  
  ### Impute ECOG and BMI with Median ###
  base_non_lab_data$ECOG_C <- with(base_non_lab_data, impute(ECOG_C, median))
  base_non_lab_data$BMI <- with(base_non_lab_data, impute(BMI, median))
  
  ### Scale BMI ###
  base_non_lab_data$z_BMI <- scale(base_non_lab_data$BMI)
  
  ###### Create and Combine with Metastases Data ########
  metastases_data <- make_metastases_data(input_data)
  out <- data.frame( base_non_lab_data , metastases_data )
  
  return( out )
}