################################################################################################
################################################################################################
####                          Code to process baseline lab value data                       ####
####                                     Richard Meier                                      ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
## replaceDotsWithNaSigns: if input is "." turn it into NA, otherwise do nothing              ##
##                                                                                            ##
## reformatTableColumnsToProperNumericVectors: for each column in a data frame replace all    ##
##                    "." with NA and type cast values to numeric                             ##
##                                                                                            ##
## normalizeFeatures: Transform each column in the lab value subset to make its distribution  ##
##                    approximately standard normal.                                          ##
##                                                                                            ##
## findColumnsWithLessMissingValuesThanThreshold: Return all columns of a data frame with     ##
##                    a lower rate of missing values than a target threshold.                 ##
##                                                                                            ##
## getNewColumnNames: Return column names of a target data frame after adding a target prefix ##
##                                                                                            ##
################################################################################################

replaceDotsWithNaSigns=function(variable){
  if(variable==".") {
    return(NA)
  }
  else{
    return(variable)
  }
}

reformatTableColumnsToProperNumericVectors=function(data){
  for(i in 1:length(data[1,])){
    columnWithReplacedDots=sapply(data[,i], FUN=replaceDotsWithNaSigns)
    data[,i]=as.numeric(columnWithReplacedDots)
  }
  return(data)
}

normalizeFeatures=function(data, filtered_columns){
  columns_already_normal=c("CA", "HB", "NA.", "PHOS", "TPRO", "MG")
  columns_left_skewed=c("ALB")
  baseline.transformed=data
  for(i in 1:length(filtered_columns)){
    current_subset=unlist(na.omit(data[filtered_columns[i]]))
    if(filtered_columns[i] %in% columns_already_normal){
      baseline.transformed[filtered_columns[i]]=scale(unlist(data[filtered_columns[i]]))
    }
    else if(filtered_columns[i] %in% columns_left_skewed){
      baseline.transformed[filtered_columns[i]]=scale(unlist(data[filtered_columns[i]])^3)
    }
    else{
      offset=abs(min(current_subset))+1;
      baseline.transformed[filtered_columns[i]]=scale(log(unlist(data[filtered_columns[i]])+offset))
    }
  }
  return(baseline.transformed)
}

findColumnsWithLessMissingValuesThanThreshold=function(data,threshold){
  na_rate=c()
  baseline_columns=colnames(data)
  for(i in 1:length(baseline_columns)){
    na_rate=c(na_rate,sum(is.na(data[baseline_columns[i]]))/length(data[,1]))
  }
  na_rates_for_columns=data.frame(baseline_columns,na_rate)
  names=na_rates_for_columns[na_rates_for_columns$na_rate<threshold,]$baseline_columns
  return(as.character(as.matrix(names)))
}

getNewColumnNames=function(data,prefix){
  newNames=c()
  for(i in 1:length(colnames(data))){
    featureName=paste(prefix,colnames(data)[i],sep="")
    newNames=c(newNames,featureName)
  }
  return(newNames)
}