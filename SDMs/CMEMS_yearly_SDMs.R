########################## R code for running Lottia yearly SDMs ############################

# These models use SST data from Copernicus (https://marine.copernicus.eu/)
# Specifically the "SST_GLO_SST_L4_REP_OBSERVATIONS_010_011" dataset, with mean daily SSTs
# Models were built on average SST from 1981-2011, and projected yearly from 2012-2019

####################### Getting layers #######################
library(ncdf4)
library(rgdal)
library(raster)

## To crop/mask SST layers to coast
setwd("~/Desktop/PD_stuffies/SDMS/shpfiles")
myshp <- readOGR("coast.buff.shp")
setwd("~/Desktop/PD_stuffies/SDMS/CMEM_data")
dat.fname <-"~/Desktop/PD_stuffies/SDMS/CMEM_data/cmem.sst.2012.nc"
sst.12 <- stack(dat.fname,varname="analysed_sst")
sst.12.m <- mean(sst.12)
sst.12.m.c <- crop(sst.12.m, extent(myshp))
sst.12.m.c.m <- mask(sst.12.m.c, myshp)
plot(sst.12.m.c.m)
writeRaster(sst.12.m.c.m, "sst.12.CMEMs.tif", format="GTiff", overwrite=TRUE)


dat.fname <-"~/Desktop/PD_stuffies/SDMS/CMEM_data/cmem.sst.2013.nc"
sst.13 <- stack(dat.fname,varname="analysed_sst")
sst.13.m <- mean(sst.13)
sst.13.m.c <- crop( sst.13.m, extent(myshp))
sst.13.m.c.m <- mask(sst.13.m.c, myshp)
plot(sst.13.m.c.m)
writeRaster(sst.13.m.c.m, "sst.13.CMEMs.tif", format="GTiff", overwrite=TRUE)


dat.fname <-"~/Desktop/PD_stuffies/SDMS/CMEM_data/CMEMs.sst.2014.nc"
sst.14 <- stack(dat.fname,varname="analysed_sst")
sst.14.m <- mean(sst.14)
sst.14.m.c <- crop( sst.14.m, extent(myshp))
sst.14.m.c.m <- mask(sst.14.m.c, myshp)
plot(sst.14.m.c.m)
writeRaster(sst.14.m.c.m, "sst.14.CMEMs.tif", format="GTiff", overwrite=TRUE)


dat.fname <-"~/Desktop/PD_stuffies/SDMS/CMEM_data/CMEMs.sst.2015.nc"
sst.15 <- stack(dat.fname,varname="analysed_sst")
sst.15.m <- mean(sst.15)
sst.15.m.c <- crop( sst.15.m, extent(myshp))
sst.15.m.c.m <- mask(sst.15.m.c, myshp)
plot(sst.15.m.c.m)
writeRaster(sst.15.m.c.m, "sst.15.CMEMs.tif", format="GTiff", overwrite=TRUE)


dat.fname <-"~/Desktop/PD_stuffies/SDMS/CMEM_data/CMEMs.sst.2016.nc"
sst.16 <- stack(dat.fname,varname="analysed_sst")
sst.16.m <- mean(sst.16)
sst.16.m.c <- crop( sst.16.m, extent(myshp))
sst.16.m.c.m <- mask(sst.16.m.c, myshp)
plot(sst.16.m.c.m)
writeRaster(sst.16.m.c.m, "sst.16.CMEMs.tif", format="GTiff", overwrite=TRUE)


dat.fname <-"~/Desktop/PD_stuffies/SDMS/CMEM_data/CMEMs.sst.2017.nc"
sst.17 <- stack(dat.fname,varname="analysed_sst")
sst.17.m <- mean(sst.17)
sst.17.m.c <- crop( sst.17.m, extent(myshp))
sst.17.m.c.m <- mask(sst.17.m.c, myshp)
plot(sst.17.m.c.m)
writeRaster(sst.17.m.c.m, "sst.17.CMEMs.tif", format="GTiff", overwrite=TRUE)


dat.fname <-"~/Desktop/PD_stuffies/SDMS/CMEM_data/CMEMs.sst.2018.nc"
sst.18 <- stack(dat.fname,varname="analysed_sst")
sst.18.m <- mean(sst.18)
sst.18.m.c <- crop( sst.18.m, extent(myshp))
sst.18.m.c.m <- mask(sst.18.m.c, myshp)
plot(sst.18.m.c.m)
writeRaster(sst.18.m.c.m, "sst.18.CMEMs.tif", format="GTiff", overwrite=TRUE)


dat.fname <-"~/Desktop/PD_stuffies/SDMS/CMEM_data/CMEMs.sst.2019.nc"
sst.19 <- stack(dat.fname,varname="analysed_sst")
sst.19.m <- mean(sst.19)
sst.19.m.c <- crop( sst.19.m, extent(myshp))
sst.19.m.c.m <- mask(sst.19.m.c, myshp)
plot(sst.19.m.c.m)
writeRaster(sst.19.m.c.m, "sst.19.CMEMs.tif", format="GTiff", overwrite=TRUE)


####################### Running models #######################
library(biomod2)

env.curr.c <- 'sst.81.11.CMEMs.tif'
env.curr.c=stack(env.curr.c)
env.curr.c.2 <- aggregate(env.curr.c, fact=2)

LG.xy <- read.csv("LG.pres.data.all.csv", sep=";")
myResp<-as.numeric(LG.xy$Lottia)
myRespName <- 'Lottia'

SPC_PresAbs <- BIOMOD_FormatingData(resp.var = myResp,
                                    expl.var = env.curr.c.2,
                                    resp.xy = LG.xy[,c('x', 'y')],
                                    resp.name = myRespName,
                                    PA.nb.rep = 3,
                                    PA.nb.absences = 500,
                                    PA.strategy = 'random')

SPC_PresAbs
plot(SPC_PresAbs)

# Set modeling options
MySpc_options <- BIOMOD_ModelingOptions(
  GLM = list( type = 'quadratic', interaction.level = 0 ),
  GAM = list( algo = 'GAM_gam' ),
  MARS = list( type = 'quadratic' ))

MySpc_models <- BIOMOD_Modeling( data = SPC_PresAbs,
                                 models = c("GLM","GAM", "RF","MARS", "FDA"),
                                 models.options = MySpc_options,
                                 NbRunEval = 10,
                                 DataSplit = 70,
                                 VarImport = 3,
                                 models.eval.meth=c('TSS','ROC'),
                                 do.full.models = F )

# get models evaluation scores
MyModels_scores <- get_evaluations(MySpc_models)

## MyModels_scores is a 5 dimension array containing the scores of the models
dim(MyModels_scores)
dimnames(MyModels_scores)

#graphically see model scores
models_scores_graph(MySpc_models, by = "models" , metrics = c("ROC","TSS"), xlim = c(0.5,1), ylim = c(0.5,1))
models_scores_graph(MySpc_models, by = "cv_run" , metrics = c("ROC","TSS"), xlim = c(0.5,1), ylim = c(0.5,1))
models_scores_graph(MySpc_models, by = "data_set" , metrics = c("ROC","TSS"), xlim = c(0.5,1), ylim = c(0.5,1))

ROC<- MyModels_scores["ROC","Testing.data",,,]
TSS<- MyModels_scores["TSS","Testing.data",,,]

write.table(ROC, TSS, file="LG.curr.model.evals.txt", sep="\t", quote = FALSE)

MySpc_models_proj_current <- BIOMOD_Projection( modeling.output = MySpc_models,
                                                new.env = env.curr.c,
                                                proj.name = "current",
                                                binary.meth = "ROC",
                                                output.format = ".img",
                                                do.stack = FALSE )

plot(MySpc_models_proj_current,  str.grep="PA1_RUN1")

MySpc_ensemble_models <- BIOMOD_EnsembleModeling( modeling.output = MySpc_models,
                                                  chosen.models ='all',
                                                  em.by = 'all',  #combine all models
                                                  eval.metric = 'all',
                                                  eval.metric.quality.threshold = c(0.55,0.8),
                                                  models.eval.meth = c('TSS','ROC'),
                                                  prob.mean = FALSE,
                                                  prob.cv = TRUE, #coefficient of variation across predictions
                                                  committee.averaging = TRUE,
                                                  prob.mean.weight = TRUE,
                                                  VarImport = 0 )

MySpc_ensemble_models_scores <- get_evaluations(MySpc_ensemble_models)
write.table(MySpc_ensemble_models_scores, file="LG.curr.model.ensmbl.evals.txt", sep="\t", quote = FALSE)


MySpc_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_current,
  binary.meth = "ROC",  #make binary predictions (pres/abs) based on ROC score
  output.format = ".img",
  do.stack = FALSE )

#### Doing for all other time points

env.12 <- 'sst.12.CMEMs.tif'
env.12=stack(env.12)
names(env.1) <- c('env.0.1')

MySpc_models_proj_12 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                          new.env = env.12,
                                          proj.name = "2012",
                                          binary.meth = c("ROC"),
                                          output.format = ".img",
                                          do.stack = FALSE)

MySpc_ensemble_models_proj_12 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_12,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)


env.13 <- 'sst.13.CMEMs.tif'
env.13=stack(env.13)
names(env.1) <- c('env.0.1')

MySpc_models_proj_13 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                           new.env = env.13,
                                           proj.name = "2013",
                                           binary.meth = c("ROC"),
                                           output.format = ".img",
                                           do.stack = FALSE)

MySpc_ensemble_models_proj_13 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_13,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)

env.14 <- 'sst.14.CMEMs.tif'
env.14=stack(env.14)
names(env.1) <- c('env.0.1')

MySpc_models_proj_14 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                           new.env = env.14,
                                           proj.name = "2014",
                                           binary.meth = c("ROC"),
                                           output.format = ".img",
                                           do.stack = FALSE)

MySpc_ensemble_models_proj_14 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_14,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)


env.15 <- 'sst.15.CMEMs.tif'
env.15=stack(env.15)
names(env.1) <- c('env.0.1')

MySpc_models_proj_15 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                           new.env = env.15,
                                           proj.name = "2015",
                                           binary.meth = c("ROC"),
                                           output.format = ".img",
                                           do.stack = FALSE)

MySpc_ensemble_models_proj_15 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_15,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)


env.16 <- 'sst.16.CMEMs.tif'
env.16=stack(env.16)
names(env.1) <- c('env.0.1')

MySpc_models_proj_16 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                           new.env = env.16,
                                           proj.name = "2016",
                                           binary.meth = c("ROC"),
                                           output.format = ".img",
                                           do.stack = FALSE)

MySpc_ensemble_models_proj_16 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_16,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)


env.17 <- 'sst.17.CMEMs.tif'
env.17=stack(env.17)
names(env.1) <- c('env.0.1')

MySpc_models_proj_17 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                           new.env = env.17,
                                           proj.name = "2017",
                                           binary.meth = c("ROC"),
                                           output.format = ".img",
                                           do.stack = FALSE)

MySpc_ensemble_models_proj_17 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_17,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)


env.18 <- 'sst.18.CMEMs.tif'
env.18=stack(env.18)
names(env.1) <- c('env.0.1')

MySpc_models_proj_18 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                           new.env = env.18,
                                           proj.name = "2018",
                                           binary.meth = c("ROC"),
                                           output.format = ".img",
                                           do.stack = FALSE)

MySpc_ensemble_models_proj_18 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_18,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)


env.19 <- 'sst.19.CMEMs.tif'
env.19=stack(env.19)
names(env.1) <- c('env.0.1')

MySpc_models_proj_19 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                           new.env = env.19,
                                           proj.name = "2019",
                                           binary.meth = c("ROC"),
                                           output.format = ".img",
                                           do.stack = FALSE)

MySpc_ensemble_models_proj_19 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_19,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)


#### Plotting SDMs #####
library(sf)
library(ggplot2)
library(rasterVis)
library(gridExtra)
library(rgdal)
library(tmap)

setwd("~/Desktop/PD_stuffies/SDMS")
MyBinLG_81.11 <- raster::stack("Lottia/proj_current_CMEM/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_2012 <- raster::stack("Lottia/proj_2012/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_2013 <- raster::stack("Lottia/proj_2013/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_2014<- raster::stack("Lottia/proj_2014/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_2015<- raster::stack("Lottia/proj_2015/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_2016<- raster::stack("Lottia/proj_2016/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_2017<- raster::stack("Lottia/proj_2017/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_2018<- raster::stack("Lottia/proj_2018/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_2019<- raster::stack("Lottia/proj_2019/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")

LG.sst.yearly <- stack(MyBinLG_81.11, MyBinLG_2012, MyBinLG_2013, MyBinLG_2014, MyBinLG_2015, MyBinLG_2016, MyBinLG_2017, MyBinLG_2018, MyBinLG_2019)

LG.sst.yearly.2 <- aggregate(LG.sst.yearly, fact=2)

pdf(file=paste("LG.sst.yearly.pdf", sep=""))
levelplot(LG.sst.yearly.2, par.settings = RdBuTheme, row=1, main="Yearly SDMs", names.attr=paste0(c('1981-2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019')))
dev.off()
