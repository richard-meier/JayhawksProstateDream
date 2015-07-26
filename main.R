################################################################################################
################################################################################################
####              Main program for fitting the final model and predicting risk              ####
####                                     Team JayHawks                                      ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
##                                                                                            ##
################################################################################################

setwd(".../GitHub/JayhawksProstateDream")
source("src/modelCombination.R")
source("dataCleaningMain.R")
library(mgcv)

setwd(".../Prostate_DREAM/data")

#### load data ####
train_data <- clean_all_data(
	coreFileName="CoreTable_training.csv", 
	lesionFileName="LesionMeasure_training.csv", 
	priorMedFile="PriorMed_training.csv", 
	isTraining=TRUE
)
validation_data <- clean_all_data(
	coreFileName="finalScoringSet/CoreTable_validation.csv", 
	lesionFileName="finalScoringSet/LesionMeasure_validation.csv", 
	priorMedFile="finalScoringSet/PriorMed_validation.csv", 
	isTraining=FALSE
)



#### Challenge 1A - Predicting Risk Scores ####

## ensemble model fit
emod = fitModelCombination(
	models=c(
		"Surv(time = LKADT_P, event = DEATH) ~ pspline(pc1,df=2) + AGE + rr_MI + AGE*rr_MI + ECOG + metastases_sum + RACE_NEW + metastases_sum*RACE_NEW + rr_protective + rm_NA. + metastases_sum*rm_NA. + sg_ESTROGENS + rm_HB + metastases_sum*rm_HB",
		"Surv(time = LKADT_P, event = DEATH) ~ pspline(pc1, df = 2) + rm_NEU + rr_z_score_weighted_medHistory + AGE + rr_MI + AGE*rr_MI + ECOG + metastases_sum + RACE_NEW + metastases_sum*RACE_NEW + rr_ormedHistory + tlv + sg_ESTROGENS + sg_BISPHOSPHONATE + rm_NEU*sg_BISPHOSPHONATE",
		"Surv(time = LKADT_P, event = DEATH) ~ tlv + harm_pro2 + pc1 + rr_z_score_weighted_medHistory + rr_ormedHistory + PROSTATE + ECOG_C + sg_ESTROGENS + rm_PHOS + rr_ormedHistory:rm_PHOS + rm_ALP:sg_GLUCOCORTICOID",
		"Surv(time = LKADT_P, event = DEATH) ~ rm_AST + rm_ALP + tlv + rm_HB + rr_z_score_weighted_medHistory + rm_LDH + rr_ormedHistory + rm_NEU + rm_PHOS + rr_ormedHistory*rm_PHOS + ECOG + rm_ALP*ECOG + harm_pro + metastases_sum",
		"Surv(time = LKADT_P, event = DEATH) ~ pc1 + rr_z_score_weighted_medHistory + sg_z_score_weighted_priorMed + tlv + rm_NEU:tlv + rr_ormedHistory + rm_PHOS + rr_ormedHistory*rm_PHOS + rm_ALP + pc1*rm_ALP + rm_HB + sg_z_score_weighted_priorMed*rm_HB"
	),
	mWeights=c(0.2,0.2,0.2,0.2,0.2),
	train_data = train_data
)

## make predictions
v_scores=predictEnsembleRisk(emod, t_data=validation_data)
v_submission_out <- data.frame( validation_data$RPT, v_scores)
colnames(v_submission_out) = c("RPT", "RISK")
# write.csv(v_submission_out,"SUBMISSIONS/final/q1a_risk_final.csv",row.names=FALSE)



#### Challenge 1B - Predicting Time To event ####
require(mgcv)

## append risk scores for training and validation ##
train_data$risk <- predictEnsembleRisk(emod, t_data=train_data)
validation_data$risk <- v_scores

## Restrict partipants in training data who died, fit non-parametric regression ##
idx_death <- which(train_data$DEATH == 1)
mod_tte <- gam( LKADT_P ~ s(risk, k = 5) , data = train_data[idx_death,])

## Make the TTE model predictions for validation!!
tte_validation <- predict( mod_tte, newdata = validation_data)

## One prediction extrapolates with estimated TTE less than 0, we will truncate at 100 days
tte_validation[which(tte_validation < 100)] <- 100

## make predictions for question q1b ###
TTE <- data.frame(validation_data$RPT, tte_validation)
colnames(TTE) = c("RPT", "TIMETOEVENT")
# write.csv( TTE, file = "SUBMISSIONS/final/q1b_tte_final.csv" , row.names = FALSE )
