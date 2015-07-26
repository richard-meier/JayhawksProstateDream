################################################################################################
################################################################################################
####                          Code to process baseline lab value data                       ####
####                                     Richard Meier                                      ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
## clean_lab_value_data: extract lab values from core table, reformat data,                   ##
##                       filter out features with more than 2/3 missing values,               ##
##                       impute missing values, transform features to normal                  ##
################################################################################################

clean_lab_value_data <- function(fileName,imputationApproach){
  # parse data
  core_table=read.csv(file=fileName, stringsAsFactors = FALSE)
  baseline_columns=c(
    "ALP","ALT","AST","CA","CREAT","HB","LDH","NEU","PLT","PSA","TBILI","TESTO","WBC",
    "CREACL","NA.","MG","PHOS","ALB","TPRO","RBC","LYM","BUN","CCRC","GLU","CREACLCA"
  )
  baseline_table=core_table[baseline_columns]
  baseline_table=reformatTableColumnsToProperNumericVectors(baseline_table)
  
  # filter high NA rates
  filtered_columns=findColumnsWithLessMissingValuesThanThreshold(baseline_table,threshold=2/3)
  baseline_table=baseline_table[filtered_columns]
  
  baseline.raw = baseline_table
  baseline.transformed = baseline_table
  baseline.TI.median = baseline_table
  baseline.TI.min = baseline_table
  baseline.TI.LC = baseline_table
  baseline.TI.MI = baseline_table
  
  # process median imputation approach #
  for(i in 1:length(baseline.TI.median[1,])){
    baseline.TI.median[,i]=impute(baseline.TI.median[,i], median)
  }
  
  # process minimum imputation approach #
  for(i in 1:length(baseline.TI.min[1,])){
    baseline.TI.min[,i]=impute(baseline.TI.min[,i], min)
  }
  
  # process linear combination imputation approach #
  if(imputationApproach == "linear"){
    lcImputationObject = transcan(
      # ~    ALT+AST   +CREAT+HB+LDH+NEU+PLT    +TBILI+TESTO+WBC    +MG     +ALB+TPRO+BUN+CCRC+GLU,
      ~ALP+ALT+AST+CA+CREAT+HB+LDH+NEU+PLT+PSA+TBILI+TESTO+WBC+NA.+MG+PHOS+ALB+TPRO+BUN+CCRC+GLU,
      iter=100,
      data=baseline.TI.LC, transformed=TRUE, imputed=TRUE, n.impute=5
    )
    linear_combination.imputations = impute.transcan(lcImputationObject, imputation=5, data=baseline.TI.LC , list.out=TRUE, pr=FALSE, check=FALSE)
    for(feature in names(linear_combination.imputations)){
      FEAT=unlist(linear_combination.imputations[feature])
      baseline.TI.LC[feature]=FEAT
    }
    baseline.TFTS.LC=lcImputationObject$transformed
  }
  
  # process multiple imputation approach #
  if(imputationApproach == "multiple"){
    miImputationObject=aregImpute(
      ~ALP+ALT+AST+CA+CREAT+HB+LDH+NEU+PLT+PSA+TBILI+TESTO+WBC+NA.+MG+PHOS+ALB+TPRO+BUN+CCRC+GLU,
      data=baseline.TI.MI, n.impute=5, nk=4, match='closest'
    )
    multiple.imputations=impute.transcan(miImputationObject, imputation=5, data=baseline.TI.MI, list.out=TRUE, pr=FALSE, check=FALSE)
    for(feature in names(multiple.imputations)){
      FEAT=unlist(multiple.imputations[feature])
      baseline.TI.MI[feature]=FEAT
    }
  }
  
  # transform tables (to z_scores) which apply #
  baseline.transformed=normalizeFeatures(baseline.transformed,filtered_columns)
  baseline.TI.median=normalizeFeatures(baseline.TI.median,filtered_columns)
  baseline.TI.min=normalizeFeatures(baseline.TI.min,filtered_columns)
  baseline.TI.LC=normalizeFeatures(baseline.TI.LC,filtered_columns)
  baseline.TI.MI=normalizeFeatures(baseline.TI.MI,filtered_columns)
  
  # change column names
  colnames(baseline.TI.median)=getNewColumnNames(baseline.TI.median,prefix="rm_")
  colnames(baseline.TI.min)=getNewColumnNames(baseline.TI.min,prefix="rm_")
  colnames(baseline.TI.MI)=getNewColumnNames(baseline.TI.MI,prefix="rm_")
  colnames(baseline.TI.LC)=getNewColumnNames(baseline.TI.LC,prefix="rm_")
  colnames(baseline.raw)=getNewColumnNames(baseline.raw,prefix="rm_")
  colnames(baseline.transformed)=getNewColumnNames(baseline.transformed,prefix="rm_")
  
  # return output
  if(imputationApproach == "median"){
    return(baseline.TI.median)
  }
  else if(imputationApproach == "minimum"){
    return(baseline.TI.min)
  }
  else if(imputationApproach == "multiple"){
    return(baseline.TI.MI)
  }
  else if(imputationApproach == "linear"){
    return(baseline.TI.LC)
  }
  else if(imputationApproach == "non_raw"){
    return(baseline.raw)
  }
  else if(imputationApproach == "non_transformed"){
    return(baseline.transformed)
  }
  else{
    cat("WARNING: Invalid imputation approach! Defaulting to median imputation!\n")
    return(baseline.TI.median)
  }
}