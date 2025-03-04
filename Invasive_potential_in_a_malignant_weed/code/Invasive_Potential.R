require(raster)
require(geosphere)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(gradientForest)
require(fields)
require(sf)
#Load adjusted_cosine_similarit.R
source("adjusted_cosine_similarity.R")
#Load climate data and allele frequencies 
SNPsENVdata <- read.csv("CMX_snp_env.csv",header = TRUE)
ENVdata <- SNPsENVdata[,4:22]
SNPsdata <- SNPsENVdata[,23:1547]
maxLevel <- log2(0.368*nrow(ENVdata)/2)
#GF model
CMXmod <- gradientForest(cbind(ENVdata, SNPsdata), predictor.vars=colnames(ENVdata),
							response.vars=colnames(SNPsdata), ntree=2000, compact=T, nbin =1001,
							maxLevel=maxLevel, trace=T, corr.threshold=0.5)
#Load climate data in current invasive ranges
rpENVdata <- read.csv("CMX_current_invasive_ranges_random_points.csv",header = TRUE)
top_variables<-c("bio15","bio19","bio3","bio17","bio14","bio8")
CMX_present_dis_tgrid=cbind(rpENVdata[,c("lon","lat")], predict(CMXmod,rpENVdata[,top_variables]))

#####Present climate in study areas##### 
#transform climate data with GF model
present_CMX<-read.csv("CMX_study_areas_present_climate_randon_points.csv",header = TRUE)
CMX_present=cbind(present_CMX[,c("lon","lat")], predict(CMXmod,present_CMX[,top_variables]))
CMX_present_spilt <- split(CMX_present, seq(nrow(CMX_present)))
#Calculation invasive potential, Multithreading implementation refer to "Maladaptation, migration and extirpation fuel climate change risk in a forest tree species"
cl <- makeCluster(20)
registerDoParallel(cl)
present_InP <- foreach(i = 1:length(CMX_present_spilt), .packages=c("fields","gdm","geosphere")) %dopar%{
present_oneInP <- CMX_present_spilt[[i]]
present_combinedInP <- CMX_present_dis_tgrid[,c("lon","lat")]
#Calculate adjusted cosine similarity between each random point within study areas and all random points within current invasive range
present_combinedInP["InP"] <- c(adjusted_cosine_similarity(present_oneInP[,top_variables], CMX_present_dis_tgrid[,top_variables]))
present_coordInP <- present_oneInP[,c("lon","lat")]
#Select maxium adjusted cosine similarity as invasive potential
present_maxCoordsInP <- present_combinedInP[which(present_combinedInT$InT == max(present_combinedInT$InT)),]
present_maxValInP <- present_maxCoordsInP$InP
present_maxPtInP <- present_maxCoordsInP[,c("lon","lat")]
present_outInP <- c(x1=present_coordInP[[1]], y1=present_coordInP[[2]], InvasivePotential=present_maxValInP, x2=present_maxPtInP[[1]],y2=present_maxPtInP[[2]])
}
stopCluster(cl)
present_InP_rbind <- do.call(rbind, present_InP)
write.csv(present_InP_rbind, paste0("CMX_present_InP.csv"), row.names=FALSE)

#####Future climate ACCESS-CM2 model(SSP126) for 2021-2040 in study areas#####
#transform climate data with GF model
futureA126_CMX<-read.csv("CMX_study_areas_future_climateA126_random_points.csv",header = TRUE)
CMX_futureA126=cbind(futureA126_CMX[,c("lon","lat")], predict(CMXmod,futureA126_CMX[,top_variables]))
CMX_futureA126_spilt <- split(CMX_futureA126, seq(nrow(CMX_futureA126)))
#Calculation invasive potential, Multithreading implementation refer to "Maladaptation, migration and extirpation fuel climate change risk in a forest tree species"
cl <- makeCluster(20)
registerDoParallel(cl)
futureA126_InP <- foreach(i = 1:length(CMX_futureA126_spilt), .packages=c("fields","gdm","geosphere")) %dopar%{
futureA126_oneInP <- CMX_futureA126_spilt[[i]]
futureA126_combinedInP <- CMX_futureA126_dis_tgrid[,c("lon","lat")]
#Calculate adjusted cosine similarity between each random point within study areas and all random points within current invasive range
futureA126_combinedInP["InP"] <- c(adjusted_cosine_similarity(futureA126_oneInP[,top_variables], CMX_futureA126_dis_tgrid[,top_variables]))
futureA126_coordInP <- futureA126_oneInP[,c("lon","lat")]
#Select maxium adjusted cosine similarity as invasive potential
futureA126_maxCoordsInP <- futureA126_combinedInP[which(futureA126_combinedInT$InT == max(futureA126_combinedInT$InT)),]
futureA126_maxValInP <- futureA126_maxCoordsInP$InP
futureA126_maxPtInP <- futureA126_maxCoordsInP[,c("lon","lat")]
futureA126_outInP <- c(x1=futureA126_coordInP[[1]], y1=futureA126_coordInP[[2]], InvasivePotential=futureA126_maxValInP, x2=futureA126_maxPtInP[[1]],y2=futureA126_maxPtInP[[2]])
}
stopCluster(cl)
futureA126_InP_rbind <- do.call(rbind, futureA126_InP)
write.csv(futureA126_InP_rbind, paste0("CMX_futureA126_InP.csv"), row.names=FALSE)

#####Future climate ACCESS-CM2 model(SSP585) for 2021-2040 in study areas#####
#transform climate data with GF model
futureA585_CMX<-read.csv("CMX_study_areas_future_climateA585_random_points.csv",header = TRUE)
CMX_futureA585=cbind(futureA585_CMX[,c("lon","lat")], predict(CMXmod,futureA585_CMX[,top_variables]))
CMX_futureA585_spilt <- split(CMX_futureA585, seq(nrow(CMX_futureA585)))
#Calculation invasive potential, Multithreading implementation refer to "Maladaptation, migration and extirpation fuel climate change risk in a forest tree species"
cl <- makeCluster(20)
registerDoParallel(cl)
futureA585_InP <- foreach(i = 1:length(CMX_futureA585_spilt), .packages=c("fields","gdm","geosphere")) %dopar%{
futureA585_oneInP <- CMX_futureA585_spilt[[i]]
futureA585_combinedInP <- CMX_futureA585_dis_tgrid[,c("lon","lat")]
#Calculate adjusted cosine similarity between each random point within study areas and all random points within current invasive range
futureA585_combinedInP["InP"] <- c(adjusted_cosine_similarity(futureA585_oneInP[,top_variables], CMX_futureA585_dis_tgrid[,top_variables]))
futureA585_coordInP <- futureA585_oneInP[,c("lon","lat")]
#Select maxium adjusted cosine similarity as invasive potential
futureA585_maxCoordsInP <- futureA585_combinedInP[which(futureA585_combinedInT$InT == max(futureA585_combinedInT$InT)),]
futureA585_maxValInP <- futureA585_maxCoordsInP$InP
futureA585_maxPtInP <- futureA585_maxCoordsInP[,c("lon","lat")]
futureA585_outInP <- c(x1=futureA585_coordInP[[1]], y1=futureA585_coordInP[[2]], InvasivePotential=futureA585_maxValInP, x2=futureA585_maxPtInP[[1]],y2=futureA585_maxPtInP[[2]])
}
stopCluster(cl)
futureA585_InP_rbind <- do.call(rbind, futureA585_InP)
write.csv(futureA585_InP_rbind, paste0("CMX_futureA585_InP.csv"), row.names=FALSE)

#####Future climate BCC-CSM2-MR model(SSP126) for 2021-2040 in study areas#####
#transform climate data with GF model
futureB126_CMX<-read.csv("CMX_study_areas_future_climateB126_random_points.csv",header = TRUE)
CMX_futureB126=cbind(futureB126_CMX[,c("lon","lat")], predict(CMXmod,futureB126_CMX[,top_variables]))
CMX_futureB126_spilt <- split(CMX_futureB126, seq(nrow(CMX_futureB126)))
#Calculation invasive potential, Multithreading implementation refer to "Maladaptation, migration and extirpation fuel climate change risk in a forest tree species"
cl <- makeCluster(20)
registerDoParallel(cl)
futureB126_InP <- foreach(i = 1:length(CMX_futureB126_spilt), .packages=c("fields","gdm","geosphere")) %dopar%{
futureB126_oneInP <- CMX_futureB126_spilt[[i]]
futureB126_combinedInP <- CMX_futureB126_dis_tgrid[,c("lon","lat")]
#Calculate adjusted cosine similarity between each random point within study areas and all random points within current invasive range
futureB126_combinedInP["InP"] <- c(adjusted_cosine_similarity(futureB126_oneInP[,top_variables], CMX_futureB126_dis_tgrid[,top_variables]))
futureB126_coordInP <- futureB126_oneInP[,c("lon","lat")]
#Select maxium adjusted cosine similarity as invasive potential
futureB126_maxCoordsInP <- futureB126_combinedInP[which(futureB126_combinedInT$InT == max(futureB126_combinedInT$InT)),]
futureB126_maxValInP <- futureB126_maxCoordsInP$InP
futureB126_maxPtInP <- futureB126_maxCoordsInP[,c("lon","lat")]
futureB126_outInP <- c(x1=futureB126_coordInP[[1]], y1=futureB126_coordInP[[2]], InvasivePotential=futureB126_maxValInP, x2=futureB126_maxPtInP[[1]],y2=futureB126_maxPtInP[[2]])
}
stopCluster(cl)
futureB126_InP_rbind <- do.call(rbind, futureB126_InP)
write.csv(futureB126_InP_rbind, paste0("CMX_futureB126_InP.csv"), row.names=FALSE)

#####Future climate BCC-CSM2-MR model(SSP585) for 2021-2040 in study areas#####
#transform climate data with GF model
futureB585_CMX<-read.csv("CMX_study_areas_future_climateB585_random_points.csv",header = TRUE)
CMX_futureB585=cbind(futureB585_CMX[,c("lon","lat")], predict(CMXmod,futureB585_CMX[,top_variables]))
CMX_futureB585_spilt <- split(CMX_futureB585, seq(nrow(CMX_futureB585)))
#Calculation invasive potential, Multithreading implementation refer to "Maladaptation, migration and extirpation fuel climate change risk in a forest tree species"
cl <- makeCluster(20)
registerDoParallel(cl)
futureB585_InP <- foreach(i = 1:length(CMX_futureB585_spilt), .packages=c("fields","gdm","geosphere")) %dopar%{
futureB585_oneInP <- CMX_futureB585_spilt[[i]]
futureB585_combinedInP <- CMX_futureB585_dis_tgrid[,c("lon","lat")]
#Calculate adjusted cosine similarity between each random point within study areas and all random points within current invasive range
futureB585_combinedInP["InP"] <- c(adjusted_cosine_similarity(futureB585_oneInP[,top_variables], CMX_futureB585_dis_tgrid[,top_variables]))
futureB585_coordInP <- futureB585_oneInP[,c("lon","lat")]
#Select maxium adjusted cosine similarity as invasive potential
futureB585_maxCoordsInP <- futureB585_combinedInP[which(futureB585_combinedInT$InT == max(futureB585_combinedInT$InT)),]
futureB585_maxValInP <- futureB585_maxCoordsInP$InP
futureB585_maxPtInP <- futureB585_maxCoordsInP[,c("lon","lat")]
futureB585_outInP <- c(x1=futureB585_coordInP[[1]], y1=futureB585_coordInP[[2]], InvasivePotential=futureB585_maxValInP, x2=futureB585_maxPtInP[[1]],y2=futureB585_maxPtInP[[2]])
}
stopCluster(cl)
futureB585_InP_rbind <- do.call(rbind, futureB585_InP)
write.csv(futureB585_InP_rbind, paste0("CMX_futureB585_InP.csv"), row.names=FALSE)