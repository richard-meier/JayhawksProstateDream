################################################################################################
################################################################################################
####            		Code to process Create Target Lesion Volume variable 				####
####                                     Joseph Usset                                       ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
## Input: Supplementary Lesion Measures Data Set for Training, Test, Validation               ##
## Output: Target Lesion Volume Variable, Non-Target Lesion Total Variable					  ##
################################################################################################



lesion_volume <- function( subject, data){ 
	
  #### Helper function that get lesion volume data for one subject

  ### INPUT:
	  ### subject ID and RAW lesion data
  ### OUTPUT:
  	### Total number of non-target, target lesions, and target lesion volume

  idx_subject <- which(c(data$RPT == subject)) #### index subject
  lesion_subj <- data[idx_subject,]            #### extract lesion data for that subject
  
  subj_non_target <- lesion_subj[which(c(lesion_subj$LSCAT != "TARGET") ),] ### separate into target and non-target
  subj_target <- lesion_subj[which(c(lesion_subj$LSCAT == "TARGET")),]     ### 
  
  ## Sum volume of target lesions 
  total_lesion_volume <- sum(as.numeric(as.character(subj_target[, "LSSTRESC"] )), na.rm = TRUE)
  
  ### output is # lesions, #non-target, # target, total target lesion volumes
  c( nrow(lesion_subj), nrow(subj_non_target), nrow(subj_target), total_lesion_volume  )
  
}  

create_target_lesion_volume <- function( lesion_data, core_data, training = TRUE){  
	
	## Function extracts target lesion volume across all subjects 
		## Input: Supplementary Lesion Data, and core data (to help input for Ascent trial)
		## Output: Target Lesion Volume, and number of non-target lesions per subject 
  
  #### Extract main variables from SCREENING ###
  if(training == TRUE){
    idx_screen <- which(lesion_data$VISIT == "SCREENING")
  }else{
    idx_screen <- which(lesion_data$VISIT == 1)      #### validation data uses different notation     
  }

  ### Extract lesion volume and number non-target lesions for subjects with non-missing data! ###
  lesion_base <- lesion_data[ idx_screen , c("STUDYID", "RPT", "VISIT", "LSCAT", "LSTEST", "LSLOC", "LSSTRESC", "LSSTRESU")]
  subjects <- unique(lesion_base$RPT) ## Find unique subject 
  lesion_volumes <- t(sapply( subjects, lesion_volume, data = lesion_base)  ) ### get lesion volumes for subjects 
  colnames(lesion_volumes) <- c("Total Lesions", "Total Non-Target Lesions", "Total Target-Lesions","Target Lesion Volume")

  
  if(training == TRUE){ 
  	
  	#### Then we need to impute data for the ASCENT TRIAL !!!! ####
    #### Target lesion volumes also need to be converted from MM to CM in training data ###
    
    ############# Extract all the lesion data in training ###
    idx_les <- which(core_data$RPT %in% subjects) 				 ### Correspond to all subjects in EFC and celgene trials ### 
    dat_les <- core_data[idx_les,]
    
    ### Exract important variables 
    dat_les$lesion_total <- lesion_volumes[,"Total Lesions"]
    dat_les$non_target_lesion_total <- lesion_volumes[,"Total Non-Target Lesions"]
    dat_les$target_lesion_total <- lesion_volumes[,"Total Target-Lesions"]
    dat_les$target_lesion_volume <- lesion_volumes[,"Target Lesion Volume"]
    #############################################
    
    ### Index for ascent 
    idx_asc <- which(core_data$STUDYID == "ASCENT2")
    
    ### Impute total lesions from Core table for the ASCENT data ###
    core_data$non_target_lesion_total <- rep(0,nrow(core_data))
    core_data$non_target_lesion_total[idx_les] <- lesion_volumes[,"Total Non-Target Lesions"]
    core_data$non_target_lesion_total[idx_asc] <- core_data$TOTAL[idx_asc]
    
    ### IMPUTATION OF VOLUME VARIABLES for ASCENT TRIAL in TRAINING DATA  ###
    core_data$target_lesion_volume <- rep(0, nrow(core_data))
    core_data$target_lesion_volume[idx_les] <- lesion_volumes[,"Target Lesion Volume"]
    
    ### get average size of target lesion within survived/died ###
    idx_surv <- which(dat_les$DEATH == 0); ### indices for survival   
    idx_death <- which(dat_les$DEATH == 1)
    avg_lesion_volumes_surv <- mean(lesion_volumes[idx_surv,"Target Lesion Volume"]/lesion_volumes[idx_surv,"Total Target-Lesions"], na.rm = TRUE)
    avg_lesion_volumes_die <- mean(lesion_volumes[idx_death,"Target Lesion Volume"]/lesion_volumes[idx_death,"Total Target-Lesions"], na.rm = TRUE)
    
    ### Get indices for survived/died in ASCENT TRIAL ###
    idx_asc_d <- intersect(idx_asc, which(core_data$DEATH == 1))
    idx_asc_s <- intersect(idx_asc, which(core_data$DEATH == 0))
    
    ### Calculate Number of target non-bone lesions from core data for ASCENT ###
    ascent_data <- core_data[idx_asc,]
    number_non_bone_lesion_ascent <- ascent_data$TOTAL - ascent_data$BONE
    core_data[idx_asc_s,c("target_lesion_volume")] <- number_non_bone_lesion_ascent[idx_asc_s]*avg_lesion_volumes_surv*2   ### 2 chosen to make 3rd quartile (survivors) match celgene trial 
    core_data[idx_asc_d,c("target_lesion_volume")] <- number_non_bone_lesion_ascent[idx_asc_d]*avg_lesion_volumes_die*3.1  ### 3.1 chosen to make 3rd quartile (deceased) match celgene trial
    
    ########## Convert training data from MM to CM!!! ###########################
    core_data$tlv <- core_data$target_lesion_volume/10
    
    ### Remove target_lesion_volume variable, and only keep the scaled tlv variable ###
    core_data <- core_data[ , -which(names(core_data) %in% c("target_lesion_volume"))]
    
  } else { ### if this is the test or validation data  
    core_data$non_target_lesion_total <- lesion_volumes[,"Total Non-Target Lesions"]
    core_data$tlv <- lesion_volumes[,"Target Lesion Volume"] ### no need to scale in validaiton!!
  }
  
  core_data
}


