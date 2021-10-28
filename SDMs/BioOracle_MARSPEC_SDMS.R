################################################################################
#Scripts written to run SDMs on Lottia gigantea

# Adapted and/or written by ES Nielsen, UC Davis, USA
#
# Code is provided as is, without support 
################################################################################



######################## DOWNLOAD PREDICTOR VARIABLES ############################
##################################################################################

LIB <- c("rgbif", "biomod2", "ggplot2", "gridExtra", "knitr", "raster", 
         "ade4", "rworldmap", "cleangeo", "maptools", "rasterVis", "rgdal","rgeos", "sdmpredictors")
#for(i in LIB) { install.packages(i, repos="http://ftp.sun.ac.za/ftp/pub/mirrors/cran.za.r/") ; library(i, character.only=T) }
for(i in LIB) { library(i, character.only=T) }

library(sdmpredictors)

#download landscape features
a.contp <- load_layers( layercodes = c("WC_bio5", "WC_bio6", "WC_bio1") , equalarea=FALSE, rasterstack=TRUE)

#download marine features
o.contp <- load_layers( layercodes = c("MS_biogeo08_sss_mean_5m", "MS_biogeo09_sss_min_5m", "MS_biogeo10_sss_max_5m", "MS_biogeo13_sst_mean_5m", "MS_biogeo14_sst_min_5m", "MS_biogeo15_sst_max_5m") , equalarea=FALSE, rasterstack=TRUE)

#resample because air and ocean vars have different resolution/extent
a<-extent(-180, 180, -90, 90) #input layer extent
o<-extent(-180, 180, -90, 90)
extent_list<-list(a, o)
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]), min(matrix_extent[2,]), max(matrix_extent[4,]))
ranges<-apply(as.matrix(best_extent), 1, diff)
reso<-res(a.contp) #choose layer you want to keep resolution
nrow_ncol<-ranges/reso
s2<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=a.contp@crs) #choose layer crs you want to keep
a.c.r.2 <-resample(a.contp, s2, method="ngb") #resample by s2
o.c.r.2 <-resample(o.contp, s2, method="ngb") #resample by s2
contemp.r.2=stack(a.c.r.2, o.c.r.2) #stack resampled layers

# Load species occurrence data
LG.xy <- read.csv("~/Desktop/PD_stuffies/SDMS/LG.pres.data.all.csv", sep=";")

# Extract env values for our species locations:
envdata <- extract(x=contemp.r.2, y=cbind(LG.xy$x, LG.xy$y))
data<- cbind2 (x=LG.xy, y=envdata)

#Remove rows with NA
envdataf <- envdata[complete.cases(envdata), ]

#Crop layers to CA extent
CA.ext <- extent(-130, -107, 20, 50)
env.curr.c <- crop(contemp.r.2, CA.ext)

# Principal Component Analysis.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# The `dudi.pca` function allows to perform the PCA over the whole study area.
#this requires library(ade4)
# We decide to keep only 2 principal component axes to summarize the whole environmental niche.
pca1 <- dudi.pca(envdataf, scannf = F, nf = 2)
round(pca1$eig/sum(pca1$eig)*100, 2)

# A preliminary test is to look for potential outliers in the environmental data.
plot(pca1$li[,1:2]) # PCA scores on first two axes
summary(pca1$li)

#plot a correlation circle
s.corcircle(pca1$co)

# calculate Pearson correlations between pairs of variables 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

( cor_current <- cor(envdataf, method="spearman") )

# reformat correlation table for graphical analyse
#here you require the package reshape
cor_current[ upper.tri(cor_current, diag = T) ] <- NA
cor_current_resh <- na.omit( melt( cor_current) )
colnames(cor_current_resh) <- c("var1", "var2", "correlation")

# only consider absolute value of correlations
cor_current_resh$correlation <- abs(cor_current_resh$correlation)
cor_current_resh[order(cor_current_resh$correlation),]

# make a correlation plot
gg_cor <- ggplot(cor_current_resh, aes(x = var1, y = var2 , fill = correlation) )
gg_cor + geom_tile() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle=90, vjust=0.5))


# Subset env variables to choose those that are non-correlated
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
VARSEL <- c("WC_bio6", "MS_biogeo08_sss_mean_5m", "MS_biogeo13_sst_mean_5m")
contemp.3vars <- stack(subset(contemp.crop, VARSEL))
vif(contemp.3vars)
#Variables      VIF
#1                 WC_bio6 5.310549
#2 MS_biogeo08_sss_mean_5m 5.332209
#3 MS_biogeo13_sst_mean_5m 6.895266

writeRaster(env.curr.c, "LG.mintemp.sss.sst.curr.tif", format="GTiff", overwrite=TRUE)

######## Building SDMs with biomod2 ###############
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(biomod2)

env.curr.c <- 'LG.mintemp.sss.sst.curr.tif'
env.curr.c=stack(env.curr.c)

LG.xy <- read.csv("~/Desktop/PD_stuffies/SDMS/LG.pres.data.all.csv", sep=";")
myResp<-as.numeric(LG.xy$Lottia)
myRespName <- 'Lottia'

SPC_PresAbs <- BIOMOD_FormatingData(resp.var = myResp,
                                    expl.var = env.curr.c,
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

######## Runing SDMs into past/future with biomod2 ###########
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##### LOAD MID-HOLOCENE LAYERS #####

a.miroc.6 <- load_layers(layercodes = c("WC_bio6_mrmid"), equalarea=FALSE, rasterstack=TRUE, 
                         datadir=NULL)

a.ccsm.6 <- load_layers(layercodes = c("WC_bio6_ccmid"), equalarea=FALSE, rasterstack=TRUE, 
                        datadir=NULL)

o.ens.6 <- load_layers(layercodes = c("MS_biogeo08_sss_mean_6kya", "MS_biogeo13_sst_mean_6kya"), equalarea=FALSE, rasterstack=TRUE, 
                       datadir=NULL)

a.miroc.6.2 <-resample(a.miroc.6, s2, method="ngb")
a.ccsm.6.2 <-resample(a.ccsm.6, s2, method="ngb") 
o.ens.6.2 <-resample(o.ens.6, s2, method="ngb")

miroc.6 =stack(a.miroc.6.2, o.ens.6.2 ) 
ccsm.6 =stack(a.ccsm.6.2, o.ens.6.2 )

miroc.6.crop <- crop(miroc.6, CA.ext)
ccsm.6.crop <- crop(ccsm.6, CA.ext)
env.MH.c <- overlay(miroc.6.crop, ccsm.6.crop, fun=mean)

vif(env.MH.c)
#Variables      VIF
#1   layer.1 5.065931
#2   layer.2 5.106875
#3   layer.3 6.083624

##### LOAD LGM LAYERS #####

a.miroc.21 <- load_layers(layercodes = c("WC_bio6_cclgm"), equalarea=FALSE, rasterstack=TRUE, 
                          datadir=NULL)

a.ccsm.21 <- load_layers(layercodes = c("WC_bio6_mrlgm"), equalarea=FALSE, rasterstack=TRUE, 
                         datadir=NULL)

o.ens.adj.21 <- load_layers(layercodes = c("MS_biogeo08_sss_mean_21kya_adjCCSM", "MS_biogeo13_sst_mean_21kya_adjCCSM"), equalarea=FALSE, rasterstack=TRUE, 
                            datadir=NULL)

o.ens.nc.21 <- load_layers(layercodes = c("MS_biogeo08_sss_mean_21kya_noCCSM", "MS_biogeo13_sst_mean_21kya_noCCSM"), equalarea=FALSE, rasterstack=TRUE, 
                           datadir=NULL)

a.miroc.21.2 <-resample(a.miroc.21, s2, method="ngb")
a.ccsm.21.2 <-resample(a.ccsm.21, s2, method="ngb") 
o.ens.adj.21.2 <-resample(o.ens.adj.21, s2, method="ngb")
o.ens.nc.21.2 <-resample(o.ens.nc.21, s2, method="ngb")

miroc.nc.21 =stack(a.miroc.21.2, o.ens.nc.21.2 ) 
ccsm.21 =stack(a.ccsm.21.2, o.ens.adj.21.2 )

miroc.nc.21.crop <- crop(miroc.nc.21, CA.ext)
ccsm.21.crop <- crop(ccsm.21, CA.ext)

env.LGM.c <- overlay(miroc.nc.21.crop, ccsm.21.crop, fun=mean)

plot(env.LGM.c)
vif(env.LGM.c)
#Variables      VIF
#1   layer.1 3.922244
#2   layer.2 9.971863
#3   layer.3 6.830728


names(env.MH.c) <- c('WC_bio6','MS_biogeo08_sss_mean_5m','MS_biogeo13_sst_mean_5m')
names(env.LGM.c) <- c('WC_bio6','MS_biogeo08_sss_mean_5m','MS_biogeo13_sst_mean_5m')


env.MH.c <- stack(env.MH.c)
env.LGM.c <- stack(env.LGM.c)

##### LOAD BIO CURR LAYERS #####
## need to use these instead of MARSPEC present day layers, as they are at same resolution as future env layers
a.contp <- load_layers( layercodes = c("WC_bio6") , equalarea=FALSE, rasterstack=TRUE)
environment.contp <- load_layers( layercodes = c("BO2_salinitymean_bdmin", "BO2_tempmean_bdmin"))

a.c.r.2 <-resample(a.contp, s2, method="ngb")
o.c.r.2 <-resample(environment.contp, s2, method="ngb")
contemp.r.2=stack(a.c.r.2, o.c.r.2)
writeRaster(contemp.r.2, "3env.BIO.curr.tif", format="GTiff", overwrite=TRUE)

##### LOAD FUTURE LAYERS #####

##RCP45 2050
mr.50.45.t <- load_layers( layercodes = c("WC_bio6_mr45_2050") , equalarea=FALSE, rasterstack=TRUE)
cc.50.45.t <- load_layers( layercodes = c("WC_bio6_cc45_2050") , equalarea=FALSE, rasterstack=TRUE)

t.50.45 <- overlay(mr.50.45.t, cc.50.45.t, fun=mean)

o.50.45 <- load_layers( layercodes = c("BO2_RCP45_2050_salinitymean_bdmin", "BO2_RCP45_2050_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
a.mid.2 <-resample(t.50.45, s2, method="ngb")
o.mid.2 <-resample(o.50.45, s2, method="ngb")
env.50.45 <- stack(a.mid.2, o.mid.2)

##RCP85 2050
mr.50.85.t <- load_layers( layercodes = c("WC_bio6_mr85_2050"))
cc.50.85.t <- load_layers( layercodes = c("WC_bio6_cc85_2050"))

t.50.85 <- overlay(mr.50.85.t, cc.50.85.t, fun=mean)

o.50.85 <- load_layers( layercodes = c( "BO2_RCP85_2050_salinitymean_bdmin", "BO2_RCP85_2050_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)

a.mid.2 <-resample(t.50.85, s2, method="ngb")
o.mid.2 <-resample(o.50.85, s2, method="ngb")
env.50.85 <- stack(a.mid.2, o.mid.2)

## RCP 45 2100/2070
mr.90.45.t <- load_layers( layercodes = c("WC_bio6_mr45_2070"))
cc.90.45.t <- load_layers( layercodes = c("WC_bio6_cc45_2070"))

t.90.45 <- overlay(mr.90.45.t, cc.90.45.t, fun=mean)

o.90.45 <- load_layers( layercodes = c("BO2_RCP45_2100_salinitymean_bdmin", "BO2_RCP45_2100_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
a.mid.2 <-resample(t.90.45, s2, method="ngb")
o.mid.2 <-resample(o.90.45, s2, method="ngb")
env.90.45 <- stack(a.mid.2, o.mid.2)

## RCP 85 2100/2070
mr.90.85.t <- load_layers( layercodes = c("WC_bio6_mr85_2070"))
cc.90.85.t <- load_layers( layercodes = c("WC_bio6_cc85_2070"))

t.90.85 <- overlay(mr.90.85.t, cc.90.85.t, fun=mean)

a.mid.2 <-resample(t.90.85, s2, method="ngb")
o.mid.2 <-resample(o.90.85, s2, method="ngb")

env.90.85 <- stack(a.mid.2, o.mid.2)

writeRaster(env.50.45, "env.50.45.tif", format="GTiff", overwrite=TRUE)
writeRaster(env.50.85, "env.50.85.tif", format="GTiff", overwrite=TRUE)
writeRaster(env.90.45, "env.90.45.tif", format="GTiff", overwrite=TRUE)
writeRaster(env.90.85, "env.90.85.tif", format="GTiff", overwrite=TRUE)

curr.B.c <- crop(contemp.r.2, CA.ext)
curr.B.c <- stack(curr.B.c)

env.50.45.c <- crop(env.50.45, CA.ext)
env.50.45.c <- stack(env.50.45.c)

env.50.85.c <- crop(env.50.85, CA.ext)
env.50.85.c <- stack(env.50.85.c)

env.90.45.c <- crop(env.90.45, CA.ext)
env.90.45.c <- stack(env.90.45.c)

env.90.85.c <- crop(env.90.85, CA.ext)
env.90.85.c <- stack(env.90.85.c)

names(env.50.45.c) <- c('WC_bio6','BO2_2100_salinitymean_bdmin','BO2_2100_tempmean_bdmin')
env.50.45.c <- stack(env.50.45.c)

SPC_PresAbs <- BIOMOD_FormatingData(resp.var = myResp,
                                    expl.var = curr.B.c,
                                    resp.xy = LG.xy[,c('x', 'y')],
                                    resp.name = myRespName,
                                    PA.nb.rep = 3,
                                    PA.nb.absences = 5000,
                                    PA.strategy = 'random') 


######## Runing SDMs just using SST ###########
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sst.contp <- load_layers( layercodes = c( "MS_biogeo13_sst_mean_5m") , equalarea=FALSE, rasterstack=TRUE)
CA.ext <- extent(-130, -107, 20, 50)
sst.contp.c <- crop(sst.contp, CA.ext)

MyBinCA_1 <- raster::stack("Lottia/proj_BIO.current/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
sst.c.m <- mask(sst.contp.c, MyBinCA_1)

names(sst.6.c.m) <- c('MS_biogeo13_sst_mean_5m')
names(sst.21.c.m) <- c('MS_biogeo13_sst_mean_5m')

sst.21.c.m <-stack(sst.21.c.m)
sst.6.c.m <-stack(sst.6.c.m)

## Running future models
sst.B.curr <- 'sst.BIO.curr.tif'
sst.B.curr=stack(sst.B.curr)

sst.B.curr.c <- crop(sst.B.curr, CA.ext)
sst.B.curr.c=stack(sst.B.curr.c)

sst.50.45 <- 'sst.50.45.tif'
sst.50.45=stack(sst.50.45)

sst.50.45.c <- crop(sst.50.45, CA.ext)
sst.50.45.c=stack(sst.50.45.c)

sst.50.85 <- 'sst.50.85.tif'
sst.50.85=stack(sst.50.85)

sst.50.85.c <- crop(sst.50.85, CA.ext)
sst.50.85.c=stack(sst.50.85.c)

sst.90.85 <- 'sst.90.85.tif'
sst.90.85=stack(sst.90.85)

sst.90.85.c <- crop(sst.90.85, CA.ext)
sst.90.85.c=stack(sst.90.85.c)

sst.90.45 <- 'sst.90.45.tif'
sst.90.45=stack(sst.90.45)

sst.90.45.c <- crop(sst.90.45, CA.ext)
sst.90.45.c=stack(sst.90.45.c)

sst.B.curr.c.m <- mask(sst.B.curr.c, MyBinCA_1)
sst.50.45.c.m <- mask(sst.50.45.c, MyBinCA_1)
sst.50.85.c.m <- mask(sst.50.85.c, MyBinCA_1)
sst.90.45.c.m <- mask(sst.90.45.c, MyBinCA_1)
sst.90.85.c.m <- mask(sst.90.85.c, MyBinCA_1)

names(sst.50.45.c.m) <- c('sst.BIO.curr')
names(sst.50.85.c.m) <- c('sst.BIO.curr')
names(sst.90.45.c.m) <- c('sst.BIO.curr')
names(sst.90.85.c.m) <- c('sst.BIO.curr')

sst.50.45.c.m=stack(sst.50.45.c.m)
sst.50.85.c.m=stack(sst.50.85.c.m)
sst.90.45.c.m=stack(sst.90.45.c.m)
sst.90.85.c.m=stack(sst.90.85.c.m)
sst.B.curr.c.m=stack(sst.B.curr.c.m)


# Set modeling options
MySpc_options <- BIOMOD_ModelingOptions(
  GBM = list( n.trees = 1000 ),
  GAM = list( algo = 'GAM_mgcv' ) )


MySpc_models_proj_sst.50.85 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                          new.env = sst.50.85.c.m,
                                          proj.name = "sst.50.85",
                                          binary.meth = c("ROC"),
                                          output.format = ".img",
                                          do.stack = FALSE)

MySpc_ensemble_models_proj_sst.50.85 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_sst.50.85,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)


MySpc_models_proj_sst.90.85 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                                  new.env = sst.90.85.c.m,
                                                  proj.name = "sst.90.85",
                                                  binary.meth = c("ROC"),
                                                  output.format = ".img",
                                                  do.stack = FALSE)

MySpc_ensemble_models_proj_sst.90.85 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_sst.90.85,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)

MySpc_models_proj_sst.90.45 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                                  new.env = sst.90.45.c.m,
                                                  proj.name = "sst.90.45",
                                                  binary.meth = c("ROC"),
                                                  output.format = ".img",
                                                  do.stack = FALSE)

MySpc_ensemble_models_proj_sst.90.45 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_sst.90.45,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)


## PAST models

names(env.MH.c) <- c('LG.mintemp.sss.sst.curr.1','LG.mintemp.sss.sst.curr.2','LG.mintemp.sss.sst.curr.3')
names(env.LGM.c) <- c('LG.mintemp.sss.sst.curr.1','LG.mintemp.sss.sst.curr.2','LG.mintemp.sss.sst.curr.3')

MySpc_models_proj_MH <- BIOMOD_Projection( modeling.output = MySpc_models,
                                                  new.env = env.MH.c,
                                                  proj.name = "3env.2nd.MH",
                                                  binary.meth = c("ROC"),
                                                  output.format = ".img",
                                                  do.stack = FALSE)

MySpc_ensemble_models_proj_MH <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_MH,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)

MySpc_models_proj_LGM <- BIOMOD_Projection( modeling.output = MySpc_models,
                                           new.env = env.LGM.c,
                                           proj.name = "3env.2nd.LGM",
                                           binary.meth = c("ROC"),
                                           output.format = ".img",
                                           do.stack = FALSE)

MySpc_ensemble_models_proj_LGM <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_LGM,
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


MyBinLG_B_curr <- raster::stack("Lottia/proj_BIO.current/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_M_curr <- raster::stack("Lottia/proj_LG.3.current/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_3.MH <- raster::stack("Lottia/proj_LG.3env.MH/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_3.LGM <- raster::stack("Lottia/proj_LG.3env.LGM/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_3.50.45<- raster::stack("Lottia/proj_LG.3env.50.45/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_3.50.85<- raster::stack("Lottia/proj_LG.3env.50.85/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_3.90.45<- raster::stack("Lottia/proj_LG.3env.90.45/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_3.90.85<- raster::stack("Lottia/proj_LG.3env.90.85/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")


MyBinLG_B_sst_curr <- raster::stack("Lottia/proj_sst.curr.BIO/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_M_sst_curr <- raster::stack("Lottia/proj_SST.mar.curr/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_1.MH <- raster::stack("Lottia/proj_sst.only.6/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_1.LGM <- raster::stack("Lottia/proj_sst.only.21/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_1.50.45<- raster::stack("Lottia/proj_sst.50.45/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_1.50.85<- raster::stack("Lottia/proj_sst.50.85/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_1.90.45<- raster::stack("Lottia/proj_sst.90.45/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinLG_1.90.85<- raster::stack("Lottia/proj_sst.90.85/individual_projections/Lottia_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")


LG.3.past <- stack(MyBinLG_M_curr, MyBinLG_3.MH, MyBinLG_3.LGM)
LG.3.fut.45 <- stack(MyBinLG_B_curr, MyBinLG_3.50.45, MyBinLG_3.90.45)
LG.3.fut.85 <- stack(MyBinLG_B_curr, MyBinLG_3.50.85, MyBinLG_3.90.85)

LG.3.past.2 <- aggregate(LG.3.past, fact=2)
LG.3.fut.45.2 <- aggregate(LG.3.fut.45, fact=2)
LG.3.fut.85.2 <- aggregate(LG.3.fut.85, fact=2)

LG.1.past <- stack(MyBinLG_M_sst_curr, MyBinLG_1.MH, MyBinLG_1.LGM)
LG.1.fut.45 <- stack(MyBinLG_B_sst_curr, MyBinLG_1.50.45, MyBinLG_1.90.45)
LG.1.fut.85 <- stack(MyBinLG_B_sst_curr, MyBinLG_1.50.85, MyBinLG_1.90.85)

LG.1.past.2 <- aggregate(LG.1.past, fact=2)
LG.1.fut.45.2 <- aggregate(LG.1.fut.45, fact=2)
LG.1.fut.85.2 <- aggregate(LG.1.fut.85, fact=2)


#pdf(file=paste("LG.3env.past.SDMs.pdf", sep=""))
#levelplot(LG.3.past.2, par.settings = RdBuTheme, names.attr=paste0(c('Present', 'Mid-Holocene', 'LGM')))
#dev.off()

#pdf(file=paste("LG.3env.fut.45.SDMs.pdf", sep=""))
#levelplot(LG.3.fut.45.2, par.settings = RdBuTheme, names.attr=paste0(c('Present', '2050', '2070/2100')))
#dev.off()

#pdf(file=paste("LG.3env.fut.85.SDMs.pdf", sep=""))
#levelplot(LG.3.fut.85.2, par.settings = RdBuTheme, names.attr=paste0(c('Present', '2050', '2070/2100')))
#dev.off()

#pdf(file=paste("LG.sst.past.SDMs.pdf", sep=""))
#levelplot(LG.1.past.2, par.settings = RdBuTheme, names.attr=paste0(c('Present', 'Mid-Holocene', 'LGM')))
#dev.off()

#pdf(file=paste("LG.sst.fut.45.SDMs.pdf", sep=""))
#levelplot(LG.1.fut.45.2, par.settings = RdBuTheme, names.attr=paste0(c('Present', '2050', '2100')))
#dev.off()

#pdf(file=paste("LG.sst.fut.85.SDMs.pdf", sep=""))
#levelplot(LG.1.fut.85.2, par.settings = RdBuTheme, row=1, names.attr=paste0(c('Present', '2050', '2100')))
#dev.off()

###
### Plotting hab suit maps ####
###

LG.past <- stack(LG.3.past.2, LG.1.past.2)
LG.fut45 <- stack(LG.3.fut.45.2, LG.1.fut.45.2)
LG.fut85 <- stack(LG.3.fut.85.2, LG.1.fut.85.2)

pdf(file=paste("LG.past.pdf", sep=""))
levelplot(LG.past, par.settings = RdBuTheme, row=1, main="Past SDMs", names.attr=paste0(c('3env-Present', '3env-MH', '3env-LGM', 'sst-Present', 'sst-MH', 'sst-LGM')))
dev.off()

pdf(file=paste("LG.fut.45.pdf", sep=""))
levelplot(LG.fut45, par.settings = RdBuTheme, row=1, main="Future RCP 4.5", names.attr=paste0(c('3env-Present', '3env-2050', '3env-2070/2100', 'sst-Present', 'sst-2050', 'sst-2100')))
dev.off()
                
pdf(file=paste("LG.fut.85.pdf", sep=""))
levelplot(LG.fut85, par.settings = RdBuTheme, row=1, main="Future RCP 8.5", names.attr=paste0(c('3env-Present', '3env-2050', '3env-2070/2100', 'sst-Present', 'sst-2050', 'sst-2100')))
dev.off()   

###
### Plotting occurence pts over annual SSTs ####
###
library(rgdal)
library(rasterVis)
setwd("~/Desktop/PD_stuffies/SDMS/shpfiles")
map1 <- readOGR("coast.clip.shp")
setwd("~/Desktop/PD_stuffies/SDMS")
LG.xy <- read.csv("~/Desktop/PD_stuffies/SDMS/LG.pres.data.all.csv", sep=";")
coordinates(LG.xy) <- ~x+y

o.contp <- load_layers( layercodes = c("MS_biogeo13_sst_mean_5m") , equalarea=FALSE, rasterstack=TRUE)
CA.ext <- extent(-130, -107, 20, 50)
sstmap <- crop(o.contp, CA.ext)

pdf(file=paste("LG.xy.ssts.pdf", sep=""))
levelplot(sstmap, par.settings = BuRdTheme, margin=FALSE, main="Occurence Pts & Mean SSTs") + layer(sp.polygons(map1, fill='white', alpha=0.2)) + layer(sp.points(LG.xy, pch=16, cex=0.5, col="black"))
dev.off()


