################################################################################################
################################################################################################
####                Code to load and derive all features used for prediction.               ####
####                 Joseph Usset, Rama Ragavan, Richard Meier, Stefan Graw                 ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
## clean_core_data: loads the core data table and derives additional features                 ##
##                                                                                            ##
## clean_all_data: loads all utilized tables and stores their features in a data frame        ##
##                                                                                            ##
################################################################################################

library(survival)
library(Hmisc)

setwd(".../GitHub/JayhawksProstateDream")
source("src/make_metastases_data.R")
source("src/clean_miscellaneous_data.R")
source("src/parsing_utils.R")
source("src/clean_lab_value_data.R")
source("src/add_additional_features.R")
source("src/clean_medical_history_data.R")
source("src/cleanPriorTreatment.R")
source("src/lesion_volume.R")

clean_core_data <- function(fileName, priorMedFile,imputation="median"){
  data.misc = clean_miscellaneous_data(fileName)
  data.lab = clean_lab_value_data(fileName,imputationApproach=imputation)
  data.priorTreatment = cleanTreatment_AddNewFeatures(fileName)
  data.medHist = clean_medical_history_data(coreFile=fileName, priorMedFile=priorMedFile)
  data.complete = cbind(data.misc,data.lab,data.priorTreatment,data.medHist)
  data.complete = add_additional_features(fileName=fileName, data.complete)
  return(data.complete)
}

clean_all_data <- function(coreFileName, lesionFileName, priorMedFile, isTraining, imputation="median"){
  data.out = clean_core_data(fileName = coreFileName, priorMedFile=priorMedFile, imputation=imputation)
  data.lesion = read.csv(lesionFileName)
  data.out = create_target_lesion_volume(lesion_data=data.lesion, core_data=data.out, training=isTraining)
  return(data.out)
}
