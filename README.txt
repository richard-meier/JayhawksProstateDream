######################################################################################
##################              Jayhawks - Read Me file             ##################
####     This file provides:                                                      ####
####        (1) INSTRUCTIONS to run our R code/prediction algorithm  #################
####        (2) descriptions of all the helper files used                         ####
######################################################################################


#### (1) INSTRUCTIONS ####

In order to make the main scripts work you should first make sure that all dependencies are installed. You can do this by running the script:
installRequirements.R
... which is located in the project directory. 

Next you have to make sure that the paths to files and scripts are setup correctly. First open the script:
main.R
... and change setwd(".../GitHub/JayhawksProstateDream") so that it points to where the software project folder is located on your local hard-drive.
Also change setwd(".../Prostate_DREAM/data") so that it points to the directory where the data provided for the challenge is located. The program assumes that the training data is inside this directory and that the data for the final scoring round is located in a sub-folder with the name "finalScoringSet".

Next open the script:
dataCleaningMain.R
... and change setwd(".../GitHub/JayhawksProstateDream") the same way you did in the previous files.

Next open the script:
deriveHardcodedWeights.R
... and change setwd(".../GitHub/JayhawksProstateDream") the same way you did in the previous files.

After you have made and saved these changes, you can use... 
source(".../main.R") to run the main program
source(".../deriveHardcodedWeights.R") to run the supplementary program showing how we obtained hardcoded weights for collapsed features.


#### (2) DESCRIPTION OF HELPER FILES ####

## main project files
dataCleaningMain.R - main routines for loading data, cleaning data and collapsing features
deriveHardcodedWeights.R - shows how to derive hard coded weights used to generate novel features
installRequirements.R - script to install required packages
main.R - main program used for submission

## utility and helper files
add_additional_features.R - creates many additional derived variables (described in write-up), uses hard-coded weights for variables (derived in deriveHardcodeWeights.R)
clean_lab_value_data.R - script to clean lab value data
clean_medical_history.R - script to clean medical history data
clean_miscellaneous_data.R - reads in/cleans main response variables + ECOG, AGE, RACE, and BMI (in core data)
cleanPriorTreatment.R - script to clean prior treatment data (in core Data)
lesion_volume.R - derived target lesion volume variable (from supplementary lesion data sets)
make_metastases_data.R - cleans/imputes/and creates new metastastses variables (e.g. total, disease site..) from core data 
ModelCombination.R - script used to create coxph model ensembles and perform predictions
modelTuning.R - script containing functions for the greedy, curated forward/backward procedure, as well as cross-validation functions
parsing_utils.R - script containing utility functions for parsing the data
