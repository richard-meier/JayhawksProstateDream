################################################################################################
################################################################################################
####                          Code to process Medical History data                          ####
####                                     Rama Raghavan                                      ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
## DATA:    Total - 37 medical history data variables split as following:                     ## 
##            13 diagnosis history variables                                                  ##
##            21 disease history in body system variables                                     ##  
##            1 number of surgery procedures variable                                         ##
##            1 social issues variable                                                        ##
##            1 lab test condition variable                                                   ##           
##                                                                                            ##
## clean_medical_history_data: create the following new variables:                            ##
##              totalmedHistory: sum of all total medical history for each patient            ##
##              ormedHistory: if patient has any medical history then 1 else 0                ##
##              diseaseBodySystemHistory: sum of all disease history in body system           ##
##              labTestResult: lab test condition variable in binary                          ##
##              socialIssues: social issues variable in binary                                ##
##              surgeryProcedures: number of surgery procedures variable in binary            ##
##              sum_CV_Diseases: sum of cardio vascular disease conditions                    ##
##                            (CV Diseases = MI, MHCARD, ARTTHROM, CEREBACC, DVT)             ##
##              or_CV_Diseases: if patient has any cardio vascular disease then 1 else 0      ##
##              harmful: sum of all marginally significant parameters with a positive HR      ##
##              protective: sum of all marginally significant parameters with a negative HR   ##
##              z_score_weighted_medHistory: sum of z-score weighted medical history features ##
##              pca_weighted_medHistory: weighted principal component analysis for medHistory ##
##              metforminInfo: if metformin information for patient is present in CMDECOD     ## 
##                             column of PriorMed table then 1 else 0                         ##
################################################################################################

## needs path to the data file and priorMed file
clean_medical_history_data <- function(coreFile, priorMedFile)
{
    library(survival)
    options(stringsAsFactors=TRUE)
    
    ## read data cvs file
    dataToBeCleaned <- read.csv( file=coreFile)

    medicalHistoryDataColumns = c("ARTTHROM" , "CEREBACC" , "CHF" , "DVT" , "DIAB" , "GASTREFL" , "GIBLEED" , "MI" , "PUD" , "PULMEMB" , "PATHFRAC" , "SPINCOMP" , "COPD" , "MHBLOOD" , "MHCARD" , "MHCONGEN" , "MHEAR" , "MHENDO" , "MHEYE" , "MHGASTRO" , "MHGEN" , "MHHEPATO" , "MHIMMUNE" , "MHINFECT" , "MHINJURY" , "MHINVEST" , "MHMETAB" , "MHMUSCLE" , "MHNEOPLA" , "MHNERV" , "MHPSYCH" , "MHRENAL" , "MHRESP" , "MHSKIN" , "MHSOCIAL" , "MHSURG" , "MHVASC")
    
    ## filter medical history data 
    medHistoryData<-dataToBeCleaned[names(dataToBeCleaned) %in% medicalHistoryDataColumns]
    
    ## make all yes/no columns standard with 'YES' and 'NO' values
    i=1
    lastItem=ncol(medHistoryData)
    while(i<=ncol(medHistoryData))
    {
        if(!all(is.na(medHistoryData[,i])))
        {
            levels(medHistoryData[,i]) <- c(levels(medHistoryData[,i]),"YES", "NO")
            medHistoryData[medHistoryData[,i] == 'Y',i] <- 'YES'
            medHistoryData[medHistoryData[,i] == '',i] <- 'NO'
            medHistoryData[,i] <- factor(medHistoryData[,i])
            i=i+1
        }
        else
        {
            medHistoryData[,i] = NULL
            medicalHistoryDataColumns <- medicalHistoryDataColumns [! medicalHistoryDataColumns %in% colnames(medHistoryData)[i]]
        }
    }

    ## code YES and NO to be 1 and 0 (binary) for further calculations
    medHistoryDataBinary = as.matrix(medHistoryData)
    medHistoryDataBinary[medHistoryDataBinary == "YES"] = 1
    medHistoryDataBinary[medHistoryDataBinary == "NO"] = 0
    medHistoryDataBinary = apply(apply(medHistoryDataBinary, FUN=as.character,MARGIN=2), FUN=as.numeric,MARGIN=2)
    
    ## is total number of medical history complaints reported
    totalmedHistory = rowSums(medHistoryDataBinary) 
    
    ## 1 if any one medical history complaint reported, 0 otherwise
    ormedHistory = totalmedHistory 
    ormedHistory[ormedHistory>0]<-1
    
    ## history of diagnosis sub-section (sum of all history of diagnosis)
    diagnosisHistory = rowSums(medHistoryDataBinary[,colnames(medHistoryDataBinary) %in% intersect(c("ARTTHROM" , "CEREBACC" , "CHF" , "DVT" , "DIAB" , "GASTREFL" , "GIBLEED" , "MI" , "PUD" , "PULMEMB" , "PATHFRAC" , "SPINCOMP" , "COPD"), medicalHistoryDataColumns)]) 
    
    ## history of disease in this body system sub-section, contains sum of number of history of disease in body systems
    diseaseBodySystemHistory1 = rowSums(medHistoryDataBinary[,colnames(medHistoryDataBinary) %in% intersect(c("MHBLOOD" , "MHCARD" , "MHCONGEN" , "MHEAR" , "MHENDO" , "MHEYE" , "MHGASTRO" , "MHGEN" , "MHHEPATO" , "MHIMMUNE" , "MHINFECT" , "MHINJURY"), medicalHistoryDataColumns)]) 
    diseaseBodySystemHistory2 = rowSums(medHistoryDataBinary[,colnames(medHistoryDataBinary) %in% intersect(c("MHMETAB" , "MHMUSCLE" , "MHNEOPLA" , "MHNERV" , "MHPSYCH" , "MHRENAL" , "MHRESP" , "MHSKIN"), medicalHistoryDataColumns)]) 
    diseaseBodySystemHistory3 = medHistoryDataBinary[,colnames(medHistoryDataBinary) %in% intersect(c("MHVASC"), medicalHistoryDataColumns)]
    diseaseBodySystemHistory = diseaseBodySystemHistory1+diseaseBodySystemHistory2+diseaseBodySystemHistory3 
    
    ## condition reported as a result of lab test
    labTestResult = medHistoryDataBinary[,colnames(medHistoryDataBinary) %in% intersect(c("MHINVEST"), medicalHistoryDataColumns)]
     
    ## personal issues that could impact condition reported
    socialIssues = medHistoryDataBinary[,colnames(medHistoryDataBinary) %in% intersect(c("MHSOCIAL"), medicalHistoryDataColumns)]
     
    ## having history of surgicel/medical procedure (1/0)
    surgeryProcedures = medHistoryDataBinary[,colnames(medHistoryDataBinary) %in% intersect(c("MHSURG"), medicalHistoryDataColumns)]
    
    ## variable to consider grouping - CV Diseases = MI OR MHCARD OR ARTTHROM OR CEREBACC OR DVT 
    sum_CV_Diseases = rowSums(medHistoryDataBinary[,colnames(medHistoryDataBinary) %in% intersect(c("ARTTHROM","CEREBACC","MI","MHCARD","DVT"), medicalHistoryDataColumns)]) 
    
    ## CV_Diseases - applied 'or' operation, so if any one complaint present
    or_CV_Diseases = sum_CV_Diseases
    or_CV_Diseases[or_CV_Diseases>0]<-1 
    
    ## create protective, harmful and z-score weighted variables based on training dataset 
    ## harmful = sum of all marginally significant parameters with a positive HR
    ## protective: sum of all marginally significant parameters with a negative HR
    ## z_score_weighted_medHistory: sum of z-score weighted medical history features
    
    medHistoryDataBinary <- data.frame(medHistoryDataBinary)
        
    protective = rep(0,length(medHistoryDataBinary[,1]))
    harmful = rep(0,length(medHistoryDataBinary[,1]))
    z_score_weighted_medHistory = rep(0,length(medHistoryDataBinary[,1]))
    
    # COXPH hazards ratio based classification on training data:
    #harmful:  CHF 
    #harmful:  MI 
    #harmful:  MHCARD 
    #harmful:  MHGEN 
    #harmful:  MHPSYCH 
    #protective:  MHVASC
    
    protective = medHistoryDataBinary$MHVASC
    harmful = apply(cbind(medHistoryDataBinary$CHF, medHistoryDataBinary$MI, medHistoryDataBinary$MHCARD, medHistoryDataBinary$MHGEN, medHistoryDataBinary$MHPSYCH),FUN=sum,MARGIN = 1)
    
    # z-score weights obtained from training set:
    #-0.1039634  0.3999269  1.686504  0.1523668  0.3880475  0.3379828  1.189789  2.546904  -0.7108761  1.275867  1.043461  -0.6048098  0.6788655  0.9067699  2.471922  -0.4939589  -0.3376898  -0.7965857  -0.8095489  0.4019387  2.619782  0.09382656  -0.095713  0.7533073  0.0313078  -0.2604924  0.7935342  -0.9298297  -1.346463  -0.5689742  2.060213  0.7715366  1.348532  0.2417617  0.9718068  -0.2436993  -1.798341  
        
    z_scores_vector <- c(-0.1039634,0.3999269,1.686504,0.1523668,0.3880475,0.3379828,1.189789,2.546904,-0.7108761,1.275867,1.043461,-0.6048098,0.6788655,0.9067699,2.471922,-0.4939589,-0.3376898,-0.7965857,-0.8095489,0.4019387,2.619782,0.09382656,-0.095713,0.7533073,0.0313078,-0.2604924,0.7935342,-0.9298297,-1.346463,-0.5689742,2.060213,0.7715366,1.348532,0.2417617,0.9718068,-0.2436993,-1.798341)
        
    for(i in 1:length(medHistoryDataBinary))
    {
        z_score_weighted_medHistory = z_score_weighted_medHistory + as.numeric(as.character(medHistoryDataBinary[,i]*z_scores_vector[i]))
    }
    
    ## create new metformin variable using PriorMed table 
    ## if metformin information is present in CMDECOD column of PriorMed table then get RPT (patient ID). 
    ## Find those RPTs in medHistoryData table and add a new column metforminInfo (value = 1 if metformin info present else value = 0)
    
    medHistoryDataRPT<-dataToBeCleaned[names(dataToBeCleaned) %in% c("RPT")]
    
    priorMedicationData <- read.csv( priorMedFile, header=T)
    index<-with(priorMedicationData,grepl("METFORMIN",CMDECOD))
    metforminRows<-priorMedicationData[index,]
    metforminInfo = ifelse(medHistoryDataRPT$RPT %in% metforminRows$RPT, 1, 0)
    
    ## combine cleaned columns and summarized new columns and return new dataset
    medHistoryDataCleaned<-cbind(medHistoryData,totalmedHistory,ormedHistory,diseaseBodySystemHistory,labTestResult,socialIssues,surgeryProcedures,sum_CV_Diseases,or_CV_Diseases, harmful, protective, z_score_weighted_medHistory, metforminInfo)
    
    colnames(medHistoryDataCleaned)=getNewColumnNames(medHistoryDataCleaned,prefix="rr_")
	return(medHistoryDataCleaned)
}
