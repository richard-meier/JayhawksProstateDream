################################################################################################
################################################################################################
####                           Code to install required packages.                           ####
####                                     Team JayHawks                                      ####
####                            Prostate Cancer DREAM Challenge                             ####
################################################################################################
##                                                                                            ##
################################################################################################

if(!require(survival)){
	install.packages("survival")
}
if(!require(Bolstad2)){
	install.packages("Bolstad2")
}
if(!require(ROCR)){
	install.packages("ROCR")
}
if(!require(timeROC)){
	install.packages("timeROC")
}
if(!require(devtools)){
	install.packages("devtools")
}
if(!require(dismo)){
	install.packages("dismo")
}
if(!require(mgcv)){
	install.packages("mgcv")
}
