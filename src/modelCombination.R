################################################################################################
################################################################################################
####                  Code to create and utilize a coxph model ensemble                     ####
####                                     Richard Meier                                      ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
## fitModelCombination: fit a target coxph model ensemble and store the fit in a custom       ##
##                      CoxphModelEnsemble object.                                            ##
##                                                                                            ##
## predictEnsembleRisk: predict the risk for a new data set utilizing a fitted                ##
##                      CoxphModelEnsemble object.                                            ##
##                                                                                            ##
## crossValidateModelCombination: cross-validate performance of a target coxph model          ##
##                      ensemble utilizing a target dataset                                   ##
##                                                                                            ##
################################################################################################

source("src/modelTuning.R")

fitModelCombination <- function(models, mWeights, train_data){
	fMods <- list()
	for(i in 1:length(models)){
		model=models[i]
		cox_sub = coxph(
			eval(parse(text=model)), data = train_data, na.action = na.omit
		)
		assignmentString = paste("fMods$M",i," <- cox_sub",sep="")
		eval(parse(text=assignmentString))
	}
	ensemble <- list(
		modelStrings = models,
		mWeights = mWeights,
		fittedModels = c(fMods)
	)
	class(ensemble) <- 'CoxphModelEnsemble'
	return(ensemble)
}

predictEnsembleRisk <- function(coxEnsemble, t_data){
	risk_scores <- rep(0,length(t_data[,1]))
	for(i in 1:length(coxEnsemble$modelStrings)){
		retrievalString = paste("coxEnsemble$fittedModels$M",i,sep="")
		cox_sub <- eval(parse(text=retrievalString))
		risk_sub <- predict( cox_sub, type = "risk", newdata = t_data )
		risk_scores <- risk_scores + coxEnsemble$mWeights[i]*risk_sub
	}
	return(risk_scores)
}

crossValidateModelCombination = function(models,mWeights,train_data,vFold=10000,c_range=534){
	aucs1=c()
	aucs2=c()
	for( ii in 1:vFold){
		index <- sample(1:1600)[1:c_range]
		sub_train <- train_data[-index,]
		sub_test <- train_data[index,]
		emod=fitModelCombination(models=models, mWeights=mWeights, train_data=sub_train)
		risk_scores <- predictEnsembleRisk(emod, t_data=sub_test)
		
		out_s = score_q1a( time = sub_test$LKADT_P, death = sub_test$DEATH, riskScoreGlobal = risk_scores)
		aucs1=c(aucs1, out_s$iAUC)
		
		rocs = timeROC( sub_test$LKADT_P, delta = sub_test$DEATH, marker = risk_scores, weighting = "cox", cause = 1, times = c(365, 548, 730 ))
		aucs2=c(aucs2, mean(rocs$AUC))
	}
	cat("iAUC=", mean(aucs1), "\taAUC=", mean(aucs2), "\n",sep="")
	return(aucs1)
}
