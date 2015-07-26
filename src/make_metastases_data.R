################################################################################################
################################################################################################
####            		Code to process Metastases Data from Core Table                         ####
####                                     Joseph Usset                                       ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
## Input: Training, Test, or Validation Data 											                        	  ##
## Output: Numeric variables for metastases location + some new derived variables 			      ##
##                                                                                            ##
## New variables:                                    										                      ##
## 				Total: sum of all metastases locations per patient                                  ##
##              Visceral: sum of visceral metastases                    					            ##
##              Other: sum of rare metastases (not bone, lymphnode, lungs prostate liver)     ##
##				Disease Site variable: Same as Halabi 										                          ##	
################################################################################################



make_metastases_data <- function(input_data){
  
  ### Convert metastases to Numeric ### 
  input_data$BONE <- as.numeric(input_data$BONE) - 1
  if( sum(input_data$BONE) == 0 ) input_data$BONE <- rep(1, nrow(input_data))
  
  input_data$LYMPH_NODES <- as.numeric(input_data$LYMPH_NODES) - 1
  input_data$LUNGS <- as.numeric(input_data$LUNGS) - 1
  input_data$PROSTATE <- as.numeric(input_data$PROSTATE) - 1
  input_data$LIVER <- as.numeric(input_data$LIVER) - 1
  input_data$OTHER <- as.numeric(input_data$OTHER) - 1
  input_data$ADRENAL <- as.numeric(input_data$ADRENAL) - 1
  input_data$BLADDER <- as.numeric(input_data$BLADDER) - 1
  input_data$SOFT_TISSUE <- as.numeric(input_data$SOFT_TISSUE) - 1
  input_data$PLEURA <- as.numeric(input_data$PLEURA) - 1
  input_data$RECTAL <- as.numeric(input_data$RECTAL) - 1
  input_data$KIDNEYS <- as.numeric(input_data$KIDNEYS) - 1
  input_data$PERITONEUM <- as.numeric(input_data$PERITONEUM) - 1
  input_data$ABDOMINAL <- as.numeric(input_data$ABDOMINAL) - 1
  input_data$COLON <- as.numeric(input_data$COLON) - 1
  
  metastases_data <- input_data[, c("BONE", "LYMPH_NODES", "LUNGS", "PROSTATE", "LIVER", "OTHER", 
                                    "SOFT_TISSUE", "ADRENAL", "BLADDER", "PLEURA", "ABDOMINAL", 
                                    "KIDNEYS", "PERITONEUM", "RECTAL", "COLON")]
  
  ### Create derived variables ###
  total <- apply( metastases_data , 1, sum, na.rm = TRUE)
  visceral <- apply( metastases_data [,-c(1,2)], 1, sum, na.rm = TRUE)
  other <- apply( metastases_data [,8:ncol(metastases_data)], 1, sum, na.rm = TRUE)
  
  metastases_data$TOTAL <- total
  metastases_data$VISCERAL <- visceral
  metastases_data$OTHER <- other
  
  ### Make variable Used by Lahabi ###
  disease_site <- rep(NA, nrow(input_data))
  disease_site[which(metastases_data$LYMPH_NODES == 1)] <- 3
  disease_site[intersect(which(metastases_data$BONE == 1), which(metastases_data$VISCERAL < 1))] <- 2
  disease_site[which(metastases_data$VISCERAL >= 1)] <- 1
  metastases_data$DISEASE_SITE <- disease_site
  metastases_data$DISEASE_SITE[which(is.na(metastases_data$DISEASE_SITE))] <- 2
  
  return(data.frame(metastases_data))
  
}