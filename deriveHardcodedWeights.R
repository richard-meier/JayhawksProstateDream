################################################################################################
################################################################################################
####      Code to demonstrate how hardcoded weights for collapsed features are derived      ####
####                                     Team JayHawks                                      ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
##                                                                                            ##
################################################################################################

setwd(".../GitHub/JayhawksProstateDream")
source("dataCleaningMain.R")

### load data
train_data <- clean_all_data(coreFileName="R:/BIO RSCH/genomics/bfridley/Prostate_DREAM/data/CoreTable_training.csv", lesionFileName="R:/BIO RSCH/genomics/bfridley/Prostate_DREAM/data/LesionMeasure_training.csv", priorMedFile="R:/BIO RSCH/genomics/bfridley/Prostate_DREAM/data/PriorMed_training.csv", isTraining=TRUE)
test_data <- clean_all_data(coreFileName="R:/BIO RSCH/genomics/bfridley/Prostate_DREAM/data/leaderboard data/prostate_cancer_challenge_data_leaderboard/CoreTable_leaderboard.csv", lesionFileName="R:/BIO RSCH/genomics/bfridley/Prostate_DREAM/data/leaderboard data/prostate_cancer_challenge_data_leaderboard/LesionMeasure_leaderboard.csv", priorMedFile="R:/BIO RSCH/genomics/bfridley/Prostate_DREAM/data/leaderboard data/prostate_cancer_challenge_data_leaderboard/PriorMed_leaderboard.csv", isTraining=FALSE)
validation_data <- clean_all_data(coreFileName="R:/BIO RSCH/genomics/bfridley/Prostate_DREAM/data/finalScoringSet/CoreTable_validation.csv", lesionFileName="R:/BIO RSCH/genomics/bfridley/Prostate_DREAM/data/finalScoringSet/LesionMeasure_validation.csv", priorMedFile="R:/BIO RSCH/genomics/bfridley/Prostate_DREAM/data/finalScoringSet/PriorMed_validation.csv", isTraining=FALSE)



### filter medical history data 
medicalHistoryDataColumns = c(
	"rr_ARTTHROM" , "rr_CEREBACC" , "rr_CHF" , "rr_DVT" , "rr_DIAB" , "rr_GASTREFL" , "rr_GIBLEED" , 
	"rr_MI" , "rr_PUD" , "rr_PULMEMB" , "rr_PATHFRAC" , "rr_SPINCOMP" , "rr_COPD" , "rr_MHBLOOD" , 
	"rr_MHCARD" , "rr_MHCONGEN" , "rr_MHEAR" , "rr_MHENDO" , "rr_MHEYE" , "rr_MHGASTRO" , "rr_MHGEN" , 
	"rr_MHHEPATO" , "rr_MHIMMUNE" , "rr_MHINFECT" , "rr_MHINJURY" , "rr_MHINVEST" , "rr_MHMETAB" , 
	"rr_MHMUSCLE" , "rr_MHNEOPLA" , "rr_MHNERV" , "rr_MHPSYCH" , "rr_MHRENAL" , "rr_MHRESP" , 
	"rr_MHSKIN" , "rr_MHSOCIAL" , "rr_MHSURG" , "rr_MHVASC"
)
medHistoryDataBinary <- as.matrix(train_data[names(train_data) %in% medicalHistoryDataColumns])
medHistoryDataBinary[medHistoryDataBinary == "YES"] = 1
medHistoryDataBinary[medHistoryDataBinary == "NO"] = 0
medHistoryDataBinary = apply(apply(medHistoryDataBinary, FUN=as.character,MARGIN=2), FUN=as.numeric,MARGIN=2)
medHistoryDataBinaryWithResponse = data.frame(train_data[c("LKADT_P" , "DEATH")], medHistoryDataBinary)

#### beneficial vs non-beneficial & z-score weighted medical hostiry data
zScoresForMedHistory=c()
rr_harmfulVariables=c()
rr_protectiveVariables=c()
for(i in 3:length(medHistoryDataBinaryWithResponse))
{
	model = coxph(Surv(time = LKADT_P, event = DEATH) ~ medHistoryDataBinaryWithResponse[,i], data = medHistoryDataBinaryWithResponse, na.action = na.omit)
	zScoresForMedHistory = c(zScoresForMedHistory, summary(model)$coefficient[,4])
	if(summary(model)$coefficients[5]<0.1)
	{
		if(exp(model$coefficients)>=1)
		{
			rr_harmfulVariables=c(rr_harmfulVariables,colnames(medHistoryDataBinaryWithResponse)[i])
		} 
		else
		{
			rr_protectiveVariables=c(rr_protectiveVariables,colnames(medHistoryDataBinaryWithResponse)[i])
		}
	}
}



### Combine to all common variables from train, test, and validation set ###
train_test_common <- intersect(names(train_data), names(test_data))
common_variables <- intersect( train_test_common, names(validation_data))
all_challenge_data <- rbind(train_data[,common_variables], test_data[,common_variables], validation_data[,common_variables])

### Derive first principal component for lab values ###
X_main_lab <- as.matrix(all_challenge_data[,c("rm_LDH", "rm_ALB", "rm_HB", "rm_ALP", "rm_PSA", "rm_AST")])
svd_X_main_lab <- svd(X_main_lab)
pc1_weights <- svd_X_main_lab$v[,1]



### beneficial vs non-beneficial & z-score weighted prior treatment data
priorColumns = c(
	"sg_ORCHIDECTOMY", "sg_PROSTATECTOMY", "sg_TURP", "sg_LYMPHADENECTOMY",  "sg_BILATERAL_ORCHIDECTOMY", "sg_PRIOR_RADIOTHERAPY",
	"sg_ANALGESICS", "sg_ANTI_ANDROGENS", "sg_GLUCOCORTICOID", "sg_GONADOTROPIN", "sg_BISPHOSPHONATE", "sg_CORTICOSTEROID", 
	"sg_IMIDAZOLE", "sg_ACE_INHIBITORS", "sg_BETA_BLOCKING", "sg_HMG_COA_REDUCT", "sg_ESTROGENS"
)
responTreat = train_data[c("LKADT_P" , "DEATH",priorColumns)]
zScoresForPriorTreatment=c()
sg_harmfulVariables=c()
sg_protectiveVariables=c()
for(i in 3:length(responTreat)){
	model = coxph(Surv(time = LKADT_P, event = DEATH) ~ responTreat[,i], data = responTreat, na.action = na.omit)
	zScoresForPriorTreatment=c(zScoresForPriorTreatment, summary(model)$coefficient[,4])
	if(summary(model)$coefficients[5]<0.1){
		if(exp(model$coefficients)>=1){
			sg_protectiveVariables=c(sg_protectiveVariables,colnames(responTreat)[i])
		} else{
			sg_harmfulVariables=c(sg_harmfulVariables,colnames(responTreat)[i])
		}
	}
}



### create weights for the metastases overall score variable
train_data$metastases_main <- as.matrix( train_data [, c("BONE" , "LYMPH_NODES", "LUNGS", "PROSTATE", "LIVER", "ADRENAL", "OTHER")])
cox_metastases = coxph(Surv(time = LKADT_P, event = DEATH) ~  metastases_main, data = train_data, na.action = na.omit)
metastases_z_score_weights <- summary(cox_metastases)$coefficients[,"z"] 
names(metastases_z_score_weights) <- NULL



### create weights for the harmpro feature
model = coxph(Surv(time = LKADT_P, event = DEATH) ~ sg_harmful + sg_protective, data = train_data, na.action = na.omit)
harmpro_weights = summary(model)$coefficient[,1]



cat("Weights for pca1 of lab values:\n","\t",pc1_weights,"\n",sep=" ")
cat("Weights for sum of z-score weighted medical history features:\n","\t",zScoresForMedHistory,"\n",sep=" ")
cat("Protective medical history features:\n","\t",rr_protectiveVariables,"\n",sep=" ")
cat("Harmful medical history features:\n","\t",rr_harmfulVariables,"\n",sep=" ")
cat("Weights for sum of z-score weighted prior treatment features:\n","\t",zScoresForPriorTreatment,"\n",sep=" ")
cat("Protective prior treatment features:\n","\t",sg_protectiveVariables,"\n",sep=" ")
cat("Harmful prior treatment features:\n","\t",sg_harmfulVariables,"\n",sep=" ")
cat("Weights for sum of z-score weighted metastases features:\n","\t",metastases_z_score_weights,"\n",sep=" ")
cat("Weights for feature harmpro:\n","\t",harmpro_weights,"\n",sep=" ")



