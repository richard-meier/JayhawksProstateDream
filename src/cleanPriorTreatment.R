################################################################################################
################################################################################################
####                          Code to process Prior Treament data                           ####
####                                     Stefan Graw                                        ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
## cleanPriorTreatment: changes Y/YES into 1 and "", which stands for NO, into 0              ##
##                      also removes vectors, that only contain NA's                          ##
##                                                                                            ##
## prior treatments: 6 procedures, 1 therapy, 12 medications                                  ##
## cleanTreatment_AddNewFeatures: create the following new variables:                         ##
##              totalPrior: sum of all prior treatments for each patient                      ##
##              totalPriorProcedures: sum of all prior procedures for each patient            ##
##              totalPriorTherapy: sum of all prior therapies (just PRIOR_RADIOTHERAPY)       ##
##              totalPriorMedication: sum of all prior medication for each patient            ##
##              protective: sum of all marginally significant parameters with a negative HR   ##
##              harmful: sum of all marginally significant parameters with a positive HR      ##
##              pca_weighted_priorMed: weighted principal component analysis                  ##
##              z_score_weighted_priorMed: sum of z-score weighted priorMed features          ##
################################################################################################

cleanPriorTreatment <- function(coreDataTable){
## extracting prior data
  priorColumns = c(
	"ORCHIDECTOMY", "PROSTATECTOMY", "TURP", "LYMPHADENECTOMY", "SPINAL_CORD_SURGERY", "BILATERAL_ORCHIDECTOMY", 
	"PRIOR_RADIOTHERAPY", "ANALGESICS", "ANTI_ANDROGENS", "GLUCOCORTICOID", "GONADOTROPIN", "BISPHOSPHONATE", "CORTICOSTEROID", 
	"IMIDAZOLE", "ACE_INHIBITORS", "BETA_BLOCKING", "HMG_COA_REDUCT", "ESTROGENS", "ANTI_ESTROGENS"
  )

  ###subCore_priorData = coreDataTable[,76:94]
  subCore_priorData = coreDataTable[,priorColumns]
  for(i in 1:length(subCore_priorData)){
# remove vectors with only NA's 
    if(all(is.na(subCore_priorData[,i]))){
      subCore_priorData[,i]=NULL
      next
    }
    
##changing "Y" and "YES" into "1" 
## and "" into "0"
    levels(subCore_priorData[,i]) <- c(levels(subCore_priorData[,i]),"1", "0")
    subCore_priorData[subCore_priorData[,i] == 'Y',i] <- '1'
    subCore_priorData[subCore_priorData[,i] == 'YES',i] <- '1'
    subCore_priorData[subCore_priorData[,i] == '',i] <- '0'
    subCore_priorData[,i] <- factor(subCore_priorData[,i])
  }
  return(subCore_priorData)
}


cleanTreatment_AddNewFeatures <- function(filePath){
### Read in the core data###
  options(stringsAsFactors=TRUE)
  core_data <- read.csv(filePath)
  
### prior treatment ###
# 6 procedures, 1 therapy, 12 medications
  subCore_prior = cleanPriorTreatment(core_data)
  
# remove two factors (almost all patients didn't have it)
  subCore_prior_reduced = subCore_prior
  subCore_prior_reduced$SPINAL_CORD_SURGERY = subCore_prior_reduced$ANTI_ESTROGENS = NULL
  
  #                             YES   NO
  # 1   LKADT_P
  # 2   DEATH                   937  663
  # 
  # 3   ORCHIDECTOMY            269 1331
  # 4   PROSTATECTOMY           448 1152
  # 5   TURP                    142 1458
  # 6   LYMPHADENECTOMY         144 1456
  # 7   SPINAL_CORD_SURGERY       4 1596
  # 8   BILATERAL_ORCHIDECTOMY  190 1410
  # 9   PRIOR_RADIOTHERAPY      926  674
  # 10  ANALGESICS              496 1104
  # 11  ANTI_ANDROGENS         1441  159
  # 12  GLUCOCORTICOID          723  877
  # 13  GONADOTROPIN           1407  193
  # 14  BISPHOSPHONATE          642  958
  # 15  CORTICOSTEROID          200 1400
  # 16  IMIDAZOLE               168 1432
  # 17  ACE_INHIBITORS          336 1264
  # 18  BETA_BLOCKING           342 1258
  # 19  HMG_COA_REDUCT          406 1194
  # 20  ESTROGENS               155 1445
  # 21  ANTI_ESTROGENS           10 1590
  
###################   creating different sums   #################################################
  
### getting total of all treatments and procedures, therapy, medications
  subCore_prior_reduced_binary = apply(apply(subCore_prior_reduced, FUN=as.character,MARGIN=2), FUN=as.numeric,MARGIN=2)
  
  totalPrior = rowSums(subCore_prior_reduced_binary)
  totalPriorProcedures = rowSums(subCore_prior_reduced_binary[,which(colnames(subCore_prior_reduced_binary)=="ORCHIDECTOMY") : which(colnames(subCore_prior_reduced_binary)=="BILATERAL_ORCHIDECTOMY")])
  totalPriorTherapy = subCore_prior_reduced_binary[,which(colnames(subCore_prior_reduced_binary)=="PRIOR_RADIOTHERAPY")]
  totalPriorMedication = rowSums(subCore_prior_reduced_binary[,which(colnames(subCore_prior_reduced_binary)=="ANALGESICS") : which(colnames(subCore_prior_reduced_binary)=="ESTROGENS")])
  
  protective = apply(cbind(subCore_prior_reduced$TURP, subCore_prior_reduced$ANALGESICS),FUN=sum,MARGIN = 1)
  harmful = apply(cbind(subCore_prior_reduced$PROSTATECTOMY, subCore_prior_reduced$GONADOTROPIN, subCore_prior_reduced$ACE_INHIBITORS),FUN=sum,MARGIN = 1)
  
  subCorePriorReducedBinary = as.matrix(subCore_prior_reduced)
  subCorePriorReducedBinary = apply(apply(subCorePriorReducedBinary, FUN=as.character,MARGIN=2), FUN=as.numeric,MARGIN=2)

##-1.021091  2.02266  -2.489221  1.091006  -0.8004286  0.208643  -5.290268  -0.5304798  0.1181942  1.790251  0.5947606  -0.9124802  0.9156636  2.220504  0.05435133  0.5340498  1.38837
  
  z_score_weighted_priorMed = rep(0,nrow(subCorePriorReducedBinary))
  
  z_scores_vector<-c(-1.021091,2.02266,-2.489221,1.091006,-0.8004286,0.208643,-5.290268,-0.5304798,0.1181942,1.790251,0.5947606,-0.9124802,0.9156636,2.220504,0.05435133,0.5340498,1.38837)
  
  for(i in 1:ncol(subCorePriorReducedBinary))
  {
    z_score_weighted_priorMed = z_score_weighted_priorMed + as.numeric(as.character(subCorePriorReducedBinary[,i]*z_scores_vector[i]))
  }
    
  out=data.frame(subCore_prior_reduced,totalPrior, totalPriorProcedures, totalPriorTherapy, totalPriorMedication, harmful,protective, z_score_weighted_priorMed)
  colnames(out)=getNewColumnNames(out,prefix="sg_")
  return(out)				
}