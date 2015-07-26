################################################################################################
################################################################################################
####                     Code for predicting and tuning model performance.                  ####
####         Prostate Cancer Dream Challenge Team (scoring methods), Richard Meier          ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
## score_q1a: provided by the Prostate Cancer Dream Challenge Team                            ##
##                                                                                            ##
## score_q1b: provided by the Prostate Cancer Dream Challenge Team                            ##
##                                                                                            ##
## score_q2: provided by the Prostate Cancer Dream Challenge Team                             ##
##                                                                                            ##
## clean_all_data: loads all utilized tables and stores their features in a data frame        ##
##                                                                                            ##
## determineFeatureToAdd: evaluates performance after trying to individually add each         ##
##                        available feature to a target prior model                           ##
##                                                                                            ##
## crossValidate: cross-validate performance of a target coxph model                          ##
##                                                                                            ##
## determineFeatureToRemove: evaluates performance after individually excluding components    ##
##                           of a target model                                                ##
##                                                                                            ##
################################################################################################


library(Bolstad2); library(ROCR); library(survival); library(timeROC); library(devtools); library(dismo);

# question 1A
# riskScoreGlobal is the global risk score
# riskScore12, riskScore18, riskScore24 are the risk scores at 12, 18 and 24 months
# time is called LKADT_P in the CoreTable meaning the last known follow up time in days
# death is last known follow up status (F=survival, T=death)
# all input parameters are vectors
# returned value is a vector containing:
#  * concordance index
#  * AUC of ROC at 12, 18, and 24 months
#  * integrated AUC (integrated over all time points)
## from Justin: https://mail.google.com/mail/u/1/#inbox/14c8f82a26d23b80
score_q1a<-function(time, death, riskScoreGlobal, riskScore12, riskScore18, riskScore24)
{

  if (missing(riskScore12)) {
    riskScore12 <- riskScoreGlobal
  }
  if (missing(riskScore18)) {
    riskScore18 <- riskScoreGlobal
  }
  if (missing(riskScore24)) {
    riskScore24 <- riskScoreGlobal
  }

  auc12 <- timeROC(T=time,
                  delta=death,
                  marker=riskScore12,
                  cause=1,
                  weighting="marginal",
                  times=12 * 30.5,
                  iid=FALSE)$AUC[2]

  auc18 <- timeROC(T=time,
                   delta=death,
                   marker=riskScore18,
                   cause=1,
                   weighting="marginal",
                   times=18 * 30.5,
                   iid=FALSE)$AUC[2]

  auc24 <- timeROC(T=time,
                   delta=death,
                   marker=riskScore24,
                   cause=1,
                   weighting="marginal",
                   times=24 * 30.5,
                   iid=FALSE)$AUC[2]

  # compute global concordance index
  surv <- Surv(time,death)
  cIndex <- survConcordance(surv ~ riskScoreGlobal)$concordance

  # compute iAUC from 6 to 30 months
  times <- seq(6,30,by=1) * 30.5
  aucs <- timeROC(T=time,
                  delta=death,
                  marker=riskScoreGlobal,
                  cause=1,
                  weighting="marginal",
                  times=times,
                  iid=FALSE)$AUC

  # Simpsons rules for integrating under curve
  iAUC <- sintegral(times, aucs)$int / (max(times) - min(times))

  return (list(cIndex=cIndex, auc12=auc12, auc18=auc18, auc24=auc24, iAUC=iAUC))
}

# question 1B
# predTime is the predicted exact survival time for all patients in days
# LKADT_P is last known follow up time in days
# DEATH is last known follow up status (F=survival, T=death)
# all input parameters are vectors
# returned value is RMSE
score_q1b<-function(predTime,LKADT_P,DEATH)
{
  x=LKADT_P[DEATH==T]
  y=predTime[DEATH==T]
  sqrt(sum((x-y)^2)/length(x))
}

# question 2
# pred is the predicted risk with a higher value corresponding to a higher
# probility of discontinue due to AE before 3 months
# y is the censored true label of discontinue due to AE before 3 month
# y=1 if it happens, 0, otherwise
# both input parameters are vectors
# returned value is AUC of PR
score_q2<-function(pred,y)
{
  ## remove subjects whose discontinuation status is NA
  pred <- pred[!is.na(y)]
  y <- y[!is.na(y)]

  prf=ROCR::performance(ROCR::prediction(pred,y),"prec","rec")
  x=prf@x.values[[1]]
  y=prf@y.values[[1]]
  auc=sum((x[-1]-x[-length(x)])*y[-1])
  auc
}

determineFeatureToAdd <- function(data, initial_model, initial_exclusion, vFold=500,c_range=534){
	exclude=c(initial_exclusion)
	feature_set=setdiff(names(data),exclude)
	
	maxModel=""
	maxAUC=0
	
	model_str=c()
	model_iAUC_mean=c()
	model_iAUC_sd=c()
	
	for(addedFeature in feature_set){
		cat("PROCESSING FEATURE:",addedFeature,"... ")
		modelString=paste(initial_model[1],sep="")
		for(i in 2:length(initial_model)){
			modelString=paste(modelString," + ",initial_model[i],sep="")
		}
		modelString=paste(modelString," + ",addedFeature,sep="")
		
		subM0=paste("Surv(time = LKADT_P, event = DEATH) ~ ",modelString,sep="")
		Submodels=c(subM0)
		frequencyTable = sort(table(train_data[addedFeature]))
		abundanceRate = frequencyTable[length(frequencyTable)] / sum(frequencyTable)
		cat("abundance rate =",abundanceRate,"\n")
		if(abundanceRate<0.77){
			for(i in 1:length(initial_model)){
				submodel=paste(modelString," + ",initial_model[i],"*",addedFeature,sep="")
				submodel=paste("Surv(time = LKADT_P, event = DEATH) ~ ",submodel,sep="")
				Submodels=c(Submodels, submodel)
			}
		}
		else{
			cat("--> HIGH ABUNDANCE: SKIPPING INTERACTIONS!\n",sep="")
		}
		for(submodel in Submodels){
			aucs=c()
			for( ii in 1:vFold){
				index <- sample(1:1600)[1:c_range]
				sub_train <- train_data[-index,]
				sub_test <- train_data[index,]
				cox_sub = coxph(
					eval(parse(text=submodel)),
					data = sub_train, na.action = na.omit
				)
				risk_scores <- predict( cox_sub, type = "risk", newdata = sub_test )
				out_s = score_q1a( time = sub_test$LKADT_P, death = sub_test$DEATH, riskScoreGlobal = risk_scores)
				aucs=c(aucs, out_s$iAUC)
			}
			AUC = mean(aucs)
			AUCSD = sd(aucs)
			
			model_str=c(model_str, submodel)
			model_iAUC_mean=c(model_iAUC_mean, AUC)
			model_iAUC_sd=c(model_iAUC_sd, AUCSD)
			if(maxAUC<AUC){
				maxAUC=AUC
				maxModel=submodel
				cat("New max model!!!\n\t",maxModel,"\n\tiAUC=",maxAUC,"\n",sep="")
			}
		}
	}
	return( cbind(model_str, model_iAUC_mean, model_iAUC_sd) )
}

crossValidate = function(model, train_data, vFold=10000,c_range=534){
	aucs1=c()
	aucs2=c()
	for( ii in 1:vFold){
		index <- sample(1:1600)[1:c_range]
		sub_train <- train_data[-index,]
		sub_test <- train_data[index,]
		cox_sub = coxph(
			eval(parse(text=model)),
			data = sub_train, na.action = na.omit
		)
		risk_scores <- predict( cox_sub, type = "risk", newdata = sub_test )
		
		out_s = score_q1a( time = sub_test$LKADT_P, death = sub_test$DEATH, riskScoreGlobal = risk_scores)
		aucs1=c(aucs1, out_s$iAUC)
		
		rocs = timeROC( sub_test$LKADT_P, delta = sub_test$DEATH, marker = risk_scores, weighting = "cox", cause = 1, times = c(365, 548, 730 ))
		aucs2=c(aucs2, mean(rocs$AUC))
	}
	cat(model, "\tiAUC=", mean(aucs1), "\taAUC=", mean(aucs2), "\n",sep="")
	return(aucs1)
}

determineFeatureToRemove <- function(data, initial_model, vFold=10000,c_range=534){
	exclusion=c()
	models=c()
	iAUCs=c()
	aAUCs=c()
	for(removedFeature in initial_model){
		modelString="0"
		for(i in 1:length(initial_model)){
			if(initial_model[i] != removedFeature){
				modelString=paste(modelString," + ",initial_model[i],sep="")
			}
		}
		modelString=paste("Surv(time = LKADT_P, event = DEATH) ~ ", modelString, sep="")
		aucs1=c()
		aucs2=c()
		for( ii in 1:vFold ){
			index <- sample(1:1600)[1:c_range]
			sub_train <- train_data[-index,]
			sub_test <- train_data[index,]
			cox_sub = coxph(
				eval(parse(text=modelString)),
				data = sub_train, na.action = na.omit
			)
			risk_scores <- predict( cox_sub, type = "risk", newdata = sub_test )
			
			out_s = score_q1a( time = sub_test$LKADT_P, death = sub_test$DEATH, riskScoreGlobal = risk_scores)
			aucs1=c(aucs1, out_s$iAUC)
			
			rocs = timeROC( sub_test$LKADT_P, delta = sub_test$DEATH, marker = risk_scores, weighting = "cox", cause = 1, times = c(365, 548, 730 ))
			aucs2=c(aucs2, mean(rocs$AUC))
		}
		iAUC = mean(aucs1)
		AUC = mean(aucs2)
		cat(modelString,"\tiAUC=",iAUC,"\taAUC=",AUC,"\n",sep="")
		iAUCs=c(iAUCs,iAUC)
		aAUCs=c(aAUCs,AUC)
		models=c(models,modelString)
		exclusion=c(exclusion,removedFeature)
	}
	return(cbind(exclusion,models,iAUCs,aAUCs))
}

retrieveCoefficientMatrix = function(coxModel){
	smry=summary(cox_test)
	significanceTest=data.frame(smry$coefficients)
	significanceTest$significant=significanceTest$p<0.05
	return(significanceTest)
}



### example: feature addition analysis
#init = c("pspline(pc1, df = 2)","rm_NEU","rr_z_score_weighted_medHistory","sg_totalPriorMedication","AGE","rr_MI","AGE*rr_MI","ECOG","metastases_sum","RACE_NEW","metastases_sum*RACE_NEW","rm_PSA","rr_ormedHistory")
#exclude=setdiff(names(train_data), names(test_data))
#exclude=c(exclude,"STUDYID","RPT","LKADT_P","LKADT_REF","LKADT_PER","DEATH","PER_REF","DISCONT","ENDTRS_C","ENTRT_PC")
#exclude=c(exclude,"LDH","ALB","HB","ALP","PSA","AST")
#exclude=c(exclude,"pc1","BONE","metastases")

#ADD_RESULT_MATRIX = determineFeatureToAdd(data=train_data, initial_model=init, initial_exclusion=exclude)
#save(ADD_RESULT_MATRIX, file=".../TRMT.RData")

#test=data.frame(ADD_RESULT_MATRIX,stringsAsFactors=FALSE)
#test$model_iAUC_mean=as.numeric(test$model_iAUC_mean)
#test$model_iAUC_sd=as.numeric(test$model_iAUC_sd)
#test=test[with(test,order(model_iAUC_mean, decreasing = TRUE)),]

#for(rwnmb in 1:30){
#	empty=crossValidate(model=test[rwnmb,1], train_data)
#	cat("\n")
#}


### example: crossvalidation of a model average
#averagingTest = crossValidateAverage(
#	model1="Surv(time = LKADT_P, event = DEATH) ~ pspline(pc1,df=2) + AGE + rr_MI + AGE*rr_MI + ECOG + metastases_sum + RACE_NEW + metastases_sum*RACE_NEW + z_PSA + rr_protective + z_PSA*rr_protective + rm_NA. + metastases_sum*rm_NA. + rm_NEU + rm_NA.*rm_NEU",
#	model2="Surv(time = LKADT_P, event = DEATH) ~ pspline(pc1, df = 2) + rm_NEU + rr_z_score_weighted_medHistory + sg_totalPriorMedication + AGE + rr_MI + AGE:rr_MI + ECOG + metastases_sum + RACE_NEW + metastases_sum:RACE_NEW + z_PSA + rr_ormedHistory + rm_PHOS + rr_ormedHistory*rm_PHOS",
#	train_data
#)


### example crossvalidation of a specific model
#currentModel="Surv(time = LKADT_P, event = DEATH) ~ pspline(pc1, df = 2) + rm_NEU + rr_z_score_weighted_medHistory + sg_totalPriorMedication + AGE + rr_MI + AGE:rr_MI + ECOG + metastases_sum + RACE_NEW + metastases_sum:RACE_NEW + rm_PSA + rr_ormedHistory"
#mah_test=crossValidate(model=currentModel, train_data)


### example searching for best model combination
#cvmc = searchForBestModelCombination(
#	models=c(
#		"Surv(time = LKADT_P, event = DEATH) ~ pspline(pc1,df=2) + AGE + rr_MI + AGE*rr_MI + ECOG + metastases_sum + RACE_NEW + metastases_sum*RACE_NEW + rr_protective + rm_NA. + metastases_sum*rm_NA. + sg_ESTROGENS + rr_ormedHistory",
#		"Surv(time = LKADT_P, event = DEATH) ~ pspline(pc1, df = 2) + rm_NEU + rr_z_score_weighted_medHistory + AGE + rr_MI + AGE*rr_MI + ECOG + metastases_sum + RACE_NEW + metastases_sum*RACE_NEW + rr_ormedHistory + tlv + sg_ESTROGENS + sg_BISPHOSPHONATE + rm_NEU*sg_BISPHOSPHONATE",
#		"Surv(time = LKADT_P, event = DEATH) ~ rm_ALP + tlv + harm_pro2 + pc1 + rr_z_score_weighted_medHistory + rr_ormedHistory + PROSTATE + ECOG_C + sg_ESTROGENS + rm_PHOS + rr_ormedHistory*rm_PHOS",
#		"Surv(time = LKADT_P, event = DEATH) ~ rm_AST + rm_ALP + tlv + harm_pro2 + rm_HB + rr_z_score_weighted_medHistory + rm_LDH + rr_ormedHistory + rm_NEU + rm_PHOS + rr_ormedHistory*rm_PHOS + ECOG + rm_ALP*ECOG"
#	),
#	train_data=train_data, rn=10, vFold=15000, modelBias=c(1,1,1,1)
#)
#test=cvmc[with(cvmc,order(iaucs, decreasing = TRUE)),]


### example: feature removal analysis
#init=c("pspline(pc1, df = 2)","rm_NEU","rr_z_score_weighted_medHistory","sg_totalPriorMedication","AGE","rr_MI","AGE:rr_MI","ECOG","metastases_sum","RACE_NEW","metastases_sum:RACE_NEW","z_PSA","rm_WBC")
#rmv_run_test = determineFeatureToRemove(data=train_data, initial_model=init)


### example: coefficient matrix retrieval
#cox_test = coxph(
#	Surv(time = LKADT_P, event = DEATH) ~ pspline(pc1, df = 2) + rm_NEU + rr_z_score_weighted_medHistory + sg_totalPriorMedication + AGE + rr_MI + AGE:rr_MI + ECOG + metastases_sum + RACE_NEW + metastases_sum:RACE_NEW + z_PSA + rr_ormedHistory,
#	data = train_data, na.action = na.omit
#)
#retrieveCoefficientMatrix(cox_test)

