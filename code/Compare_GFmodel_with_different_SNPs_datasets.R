require(raster)
require(geosphere)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(gradientForest)
require(fields)
require(sf)
#####SNPs both located in selective regions and associated with environmental variables#####
SNPsENVdata <- read.csv("CMX_snp_env.csv",header = TRUE)
ENVdata <- SNPsENVdata[,4:22]
SNPsdata <- SNPsENVdata[,23:1659]
maxLevel <- log2(0.368*nrow(ENVdata)/2)
CMXmod_bio <- gradientForest(cbind(ENVdata, SNPsdata), predictor.vars=colnames(ENVdata),
							response.vars=colnames(SNPsdata), ntree=2000, compact=T, nbin =1001,
							maxLevel=maxLevel, trace=T, corr.threshold=0.5)
							
#####random SNPs located in selective regions#####
#Replicate 1			
SNPsENVdata <- read.csv("CMX_snp_env_select_ran1.csv",header = TRUE)
ENVdata <- SNPsENVdata[,4:22]
SNPsdata <- SNPsENVdata[,23:1547]
maxLevel <- log2(0.368*nrow(ENVdata)/2)
CMXmod_selR1 <- gradientForest(cbind(ENVdata, SNPsdata), predictor.vars=colnames(ENVdata),
								response.vars=colnames(SNPsdata), ntree=2000, compact=T, nbin =1001,
								maxLevel=maxLevel, trace=T, corr.threshold=0.5)
#Replicate 2
SNPsENVdata <- read.csv("CMX_snp_env_select_ran2.csv",header = TRUE)
ENVdata <- SNPsENVdata[,4:22]
SNPsdata <- SNPsENVdata[,23:1547]
maxLevel <- log2(0.368*nrow(ENVdata)/2)
CMXmod_selR2 <- gradientForest(cbind(ENVdata, SNPsdata), predictor.vars=colnames(ENVdata),
								response.vars=colnames(SNPsdata), ntree=2000, compact=T, nbin =1001,
								maxLevel=maxLevel, trace=T, corr.threshold=0.5)
#Replicate 3
SNPsENVdata <- read.csv("CMX_snp_env_select_ran3.csv",header = TRUE)
ENVdata <- SNPsENVdata[,4:22]
SNPsdata <- SNPsENVdata[,23:1547]
maxLevel <- log2(0.368*nrow(ENVdata)/2)
CMXmod_selR3 <- gradientForest(cbind(ENVdata, SNPsdata), predictor.vars=colnames(ENVdata),
								response.vars=colnames(SNPsdata), ntree=2000, compact=T, nbin =1001,
								maxLevel=maxLevel, trace=T, corr.threshold=0.5)

#####random SNPs associated with environmental variables#####
#Replicate 1
SNPsENVdata <- read.csv("CMX_snp_env_climate_ran1.csv",header = TRUE)
ENVdata <- SNPsENVdata[,4:22]
SNPsdata <- SNPsENVdata[,23:1547]
maxLevel <- log2(0.368*nrow(ENVdata)/2)
CMXmod_cliR1 <- gradientForest(cbind(ENVdata, SNPsdata), predictor.vars=colnames(ENVdata),
								response.vars=colnames(SNPsdata), ntree=2000, compact=T, nbin =1001,
								maxLevel=maxLevel, trace=T, corr.threshold=0.5)
#Replicate 2
SNPsENVdata <- read.csv("CMX_snp_env_climate_ran2.csv",header = TRUE)
ENVdata <- SNPsENVdata[,4:22]
SNPsdata <- SNPsENVdata[,23:1547]
maxLevel <- log2(0.368*nrow(ENVdata)/2)
CMXmod_cliR2 <- gradientForest(cbind(ENVdata, SNPsdata), predictor.vars=colnames(ENVdata),
								response.vars=colnames(SNPsdata), ntree=2000, compact=T, nbin =1001,
								maxLevel=maxLevel, trace=T, corr.threshold=0.5)
#Replicate 3
SNPsENVdata <- read.csv("CMX/CMX_snp_env_climate_ran3.csv",header = TRUE)
ENVdata <- SNPsENVdata[,4:22]
SNPsdata <- SNPsENVdata[,23:1547]
maxLevel <- log2(0.368*nrow(ENVdata)/2)
CMXmod_cliR3 <- gradientForest(cbind(ENVdata, SNPsdata), predictor.vars=colnames(ENVdata),
								response.vars=colnames(SNPsdata), ntree=2000, compact=T, nbin =1001,
								maxLevel=maxLevel, trace=T, corr.threshold=0.5)


#####random all SNPs#####
#Replicate 1
SNPsENVdata <- read.csv("CMX_snp_env_all_ran1.csv",header = TRUE)
ENVdata <- SNPsENVdata[,4:22]
SNPsdata <- SNPsENVdata[,23:1547]
maxLevel <- log2(0.368*nrow(ENVdata)/2)
CMXmod_allR1 <- gradientForest(cbind(ENVdata, SNPsdata), predictor.vars=colnames(ENVdata),
								response.vars=colnames(SNPsdata), ntree=2000, compact=T, nbin =1001,
								maxLevel=maxLevel, trace=T, corr.threshold=0.5)
#Replicate 2
SNPsENVdata <- read.csv("CMX_snp_env_all_ran2.csv",header = TRUE)
ENVdata <- SNPsENVdata[,4:22]
SNPsdata <- SNPsENVdata[,23:1547]
maxLevel <- log2(0.368*nrow(ENVdata)/2)
CMXmod_allR2 <- gradientForest(cbind(ENVdata, SNPsdata), predictor.vars=colnames(ENVdata),
								response.vars=colnames(SNPsdata), ntree=2000, compact=T, nbin =1001,
								maxLevel=maxLevel, trace=T, corr.threshold=0.5)
#Replicate 3
SNPsENVdata <- read.csv("CMX_snp_env_all_ran3.csv",header = TRUE)
ENVdata <- SNPsENVdata[,4:22]
SNPsdata <- SNPsENVdata[,23:1547]
maxLevel <- log2(0.368*nrow(ENVdata)/2)
CMXmod_allR3 <- gradientForest(cbind(ENVdata, SNPsdata), predictor.vars=colnames(ENVdata),
								response.vars=colnames(SNPsdata), ntree=2000, compact=T, nbin =1001,
								maxLevel=maxLevel, trace=T, corr.threshold=0.5)

###count model performance 
write.csv(CMXmod_selR1$overall.imp, file="CMXmod_selR1_imp.csv")
write.csv(CMXmod_selR1$X, file="CMXmod_selR1_Y.csv")
write.csv(CMXmod_selR1$Y, file="CMXmod_selR1_X.csv")
write.csv(CMXmod_selR1$imp.rsq, file="CMXmod_selR1_impRsq.csv")
write.csv(CMXmod_selR1$result, file="CMXmod_selR1_result.csv")
write.csv(CMXmod_selR1$res.u, file="CMXmod_selR1_res_u.csv")
write.csv(CMXmod_selR1$res, file="CMXmod_selR1_res.csv")

write.csv(CMXmod_selR2$overall.imp, file="CMXmod_selR2_imp.csv")
write.csv(CMXmod_selR2$X, file="CMXmod_selR2_Y.csv")
write.csv(CMXmod_selR2$Y, file="CMXmod_selR2_X.csv")
write.csv(CMXmod_selR2$imp.rsq, file="CMXmod_selR2_impRsq.csv")
write.csv(CMXmod_selR2$result, file="CMXmod_selR2_result.csv")
write.csv(CMXmod_selR2$res.u, file="CMXmod_selR2_res_u.csv")
write.csv(CMXmod_selR2$res, file="CMXmod_selR2_res.csv")

write.csv(CMXmod_selR3$overall.imp, file="CMXmod_selR3_imp.csv")
write.csv(CMXmod_selR3$X, file="CMXmod_selR3_Y.csv")
write.csv(CMXmod_selR3$Y, file="CMXmod_selR3_X.csv")
write.csv(CMXmod_selR3$imp.rsq, file="CMXmod_selR3_impRsq.csv")
write.csv(CMXmod_selR3$result, file="CMXmod_selR3_result.csv")
write.csv(CMXmod_selR3$res.u, file="CMXmod_selR3_res_u.csv")
write.csv(CMXmod_selR3$res, file="CMXmod_selR3_res.csv")

write.csv(CMXmod_allR1$overall.imp, file="CMXmod_allR1_imp.csv")
write.csv(CMXmod_allR1$X, file="CMXmod_allR1_Y.csv")
write.csv(CMXmod_allR1$Y, file="CMXmod_allR1_X.csv")
write.csv(CMXmod_allR1$imp.rsq, file="CMXmod_allR1_impRsq.csv")
write.csv(CMXmod_allR1$result, file="CMXmod_allR1_result.csv")
write.csv(CMXmod_allR1$res.u, file="CMXmod_allR1_res_u.csv")
write.csv(CMXmod_allR1$res, file="CMXmod_allR1_res.csv")

write.csv(CMXmod_allR2$overall.imp, file="CMXmod_allR2_imp.csv")
write.csv(CMXmod_allR2$X, file="CMXmod_allR2_Y.csv")
write.csv(CMXmod_allR2$Y, file="CMXmod_allR2_X.csv")
write.csv(CMXmod_allR2$imp.rsq, file="CMXmod_allR2_impRsq.csv")
write.csv(CMXmod_allR2$result, file="CMXmod_allR2_result.csv")
write.csv(CMXmod_allR2$res.u, file="CMXmod_allR2_res_u.csv")
write.csv(CMXmod_allR2$res, file="CMXmod_allR2_res.csv")

write.csv(CMXmod_allR3$overall.imp, file="CMXmod_allR3_imp.csv")
write.csv(CMXmod_allR3$X, file="CMXmod_allR3_Y.csv")
write.csv(CMXmod_allR3$Y, file="CMXmod_allR3_X.csv")
write.csv(CMXmod_allR3$imp.rsq, file="CMXmod_allR3_impRsq.csv")
write.csv(CMXmod_allR3$result, file="CMXmod_allR3_result.csv")
write.csv(CMXmod_allR3$res.u, file="CMXmod_allR3_res_u.csv")
write.csv(CMXmod_allR3$res, file="CMXmod_allR3_res.csv")

write.csv(CMXmod_cliR1$overall.imp, file="CMXmod_cliR1_imp.csv")
write.csv(CMXmod_cliR1$X, file="CMXmod_cliR1_Y.csv")
write.csv(CMXmod_cliR1$Y, file="CMXmod_cliR1_X.csv")
write.csv(CMXmod_cliR1$imp.rsq, file="CMXmod_cliR1_impRsq.csv")
write.csv(CMXmod_cliR1$result, file="CMXmod_cliR1_result.csv")
write.csv(CMXmod_cliR1$res.u, file="CMXmod_cliR1_res_u.csv")
write.csv(CMXmod_cliR1$res, file="CMXmod_cliR1_res.csv")

write.csv(CMXmod_cliR2$overall.imp, file="CMXmod_cliR2_imp.csv")
write.csv(CMXmod_cliR2$X, file="CMXmod_cliR2_Y.csv")
write.csv(CMXmod_cliR2$Y, file="CMXmod_cliR2_X.csv")
write.csv(CMXmod_cliR2$imp.rsq, file="CMXmod_cliR2_impRsq.csv")
write.csv(CMXmod_cliR2$result, file="CMXmod_cliR2_result.csv")
write.csv(CMXmod_cliR2$res.u, file="CMXmod_cliR2_res_u.csv")
write.csv(CMXmod_cliR2$res, file="CMXmod_cliR2_res.csv")

write.csv(CMXmod_cliR3$overall.imp, file="CMXmod_cliR3_imp.csv")
write.csv(CMXmod_cliR3$X, file="CMXmod_cliR3_Y.csv")
write.csv(CMXmod_cliR3$Y, file="CMXmod_cliR3_X.csv")
write.csv(CMXmod_cliR3$imp.rsq, file="CMXmod_cliR3_impRsq.csv")
write.csv(CMXmod_cliR3$result, file="CMXmod_cliR3_result.csv")
write.csv(CMXmod_cliR3$res.u, file="CMXmod_cliR3_res_u.csv")
write.csv(CMXmod_cliR3$res, file="CMXmod_cliR3_res.csv")