#R project -- assessing if more bands are always better
#Lee Chengfa Benjamin
#2018-04-11

##Readme 
#Includes a hard coded removal of anomalous results (line 231, object = meta_md)
#To remove hard coding if codes are reused (ie. set meta_md=meta in code)

##libraries required
library(raster)
library(rgdal)
library(RStoolbox)
#also require "dplyr" and "graphics" packages, not loaded to avoid masking

##Import data
#setwd("C:/Users/leech/OneDrive/Documents/Uni/Winter 2017/MB2 - Programming and Geostatistics/Project")
ras <-brick("data/p224r63_2011.grd") 
train <- readOGR("data/training_2011.shp")
valid <- readOGR("data/validation_2011.shp")
#projection checking
  identical(crs(train),crs(ras))
  identical(crs(valid),crs(ras))
  
##Parameters for supervised classification
rep <- 10 #set number of reps
sc.model <- "mlc" #rf (random forest) or mlc (max likelihood)
#models from caret package also usable (eg. "kknn" for k-nearest neighbour)
seedno <- 30 #standardising execution of supervised classification

#############################################################
################  Execution script  #########################
#############################################################


##creating directories for output
if(!dir.exists("output")){dir.create("output",FALSE)}
if(!dir.exists("output/map")){dir.create("output/map",FALSE)}
if(!dir.exists("output/stats")){dir.create("output/stats",FALSE)}


##Doing a benchmark classification for the original raster stack
maxband <- length(names(ras))
for (x in 1:rep){
  cat(paste0("\n ","Running for x= ",x," of ",rep," for repetitive classification."))
  set.seed(seedno) #to ensure supervised classification is similar
  ptm <- proc.time() #start timer
  sc.valid <- superClass(ras,trainData=train,
                         valData=valid,
                         responseCol="class_name",
                         model=sc.model
  )
  timing <-  proc.time() - ptm #end timer
  if(x==1){
    meta <- as.list(c(timing[3],
                      sc.valid$validation$performance$overall[1:2], #overall accuracy & kappa
                      "benchmark","benchmark" #creation of column for no. of bands removed and name of bands removed
    ))
    names(meta) <- c(names(timing[3]),
                     names(sc.valid$validation$performance$overall[1:2]),
                     "n_bands","bands.removed")
  }else{
    compile <- as.list(c(timing[3],
                         sc.valid$validation$performance$overall[1:2], #overall accuracy & kappa
                         "benchmark","benchmark"))
    names(compile) <- c(names(timing[3]),
                        names(sc.valid$validation$performance$overall[1:2]),
                        "n_bands","bands.removed")
    meta <- dplyr::bind_rows(meta,compile)
  }
}
writeRaster(sc.valid$map,filename=paste0("output/map/","class_original_stack",".tif"),
            format="GTiff",geoinfo=geoinfo)


##Raster PCA of the original stack
##Feature reduction
pca.ras <- rasterPCA(ras,maskCheck=TRUE)
  summary(pca.ras$model) #summary of data
  loadings(pca.ras$model)

  
##Raster inflation (creating more bands)
##Band creation (spectral indices)
spec.ras <- spectralIndices(ras,
                            blue="B1_sre",green="B2_sre",red="B3_sre",nir="B4_sre",
                            swir2="B5_sre",swir3="B7_sre")
spec.ras <- spec.ras[[-c(2:4,15,19)]] #removing indices with infinity values or other issues
##Band creation (tassel cap)
tas.ras <- tasseledCap(ras[[-6]],sat="Landsat7ETM") #band 6 not to be included
##Band creation (band ratios)
bandno <- seq(1,maxband,1)
for (c in 1:maxband){
  for (d in bandno[-c]){
    new.band <- ras[[c]]/ras[[d]]
    new.band[is.infinite(new.band)] <- 0 #removes infinity values and convert to zero
    names(new.band) <- paste0(names(ras)[c],".",names(ras)[d])#special characters cannot be coded into names
    ras <- stack(ras,new.band)
  }
}
##Combining spectral indices and tasseled cap for inflating raster stack
ras <- stack(ras,spec.ras,tas.ras)
rm(spec.ras,tas.ras,new.band) #removing old objects to declutter environment


maxband <- length(names(ras)) #overwrite new maxband


##Supervised classification and timing
##Supervised classification of the inflated stack (without removal)
for (x in 1:rep){
  cat(paste0("\n ","Running for x= ",x," of ",rep," for repetitive classification."))
  set.seed(seedno) #to ensure supervised classification is similar
  ptm <- proc.time() #start timer
  sc.valid <- superClass(ras,trainData=train,
                         valData=valid,
                         responseCol="class_name",
                         model=sc.model
                         )
  timing <-  proc.time() - ptm #end timer
  compile <- as.list(c(timing[3],
                       sc.valid$validation$performance$overall[1:2], #overall accuracy & kappa
                       length(names(ras)),"null"))
  names(compile) <- c(names(timing[3]),
                      names(sc.valid$validation$performance$overall[1:2]),
                      "n_bands","bands.removed")
  meta <- dplyr::bind_rows(meta,compile)
}
writeRaster(sc.valid$map,filename=paste0("output/map/","class_",length(names(new.ras)),".tif"),
            format="GTiff",geoinfo=geoinfo)


##Supervised classification of the inflated stack (with progressive band removal)
for (n in 1:(maxband-2)){
  cat(paste0("\n ","Running for n = ",n," of ",maxband-2," for band removing (overall)."))
  for (m in 1:rep){
    print(paste0("Running for m= ",m," of ",rep," for repetitive classification."))
    set.seed(m*n) #ensure band removal is different/pseudorandom for all reiterations, preemptive measure against fixed seed setting later in loop
    x <- sort(sample(1:maxband,n)) 
    new.ras <- ras[[-x]] #reduced raster
    removed <- names(ras[[x]]) #tracking removed bands
    set.seed(seedno) #to ensure supervised classification is similar
    #superClass
    ptm <- proc.time() #start timer
    sc.valid <- superClass(new.ras,trainData=train,
                           valData=valid,
                           responseCol="class_name",
                           model=sc.model
                           )
    timing <-  proc.time() - ptm #end timer
    compile <- as.list(c(timing[3],
                         sc.valid$validation$performance$overall[1:2], #overall accuracy & kappa
                         length(names(new.ras)), #count number of layers
                         do.call("paste", as.list(removed))))
    names(compile) <- c(names(timing[3]),
                        names(sc.valid$validation$performance$overall[1:2]),
                        "n_bands","bands.removed")
    meta <- dplyr::bind_rows(meta,compile)
  }
writeRaster(sc.valid$map,filename=paste0("output/map/","class_",length(names(new.ras)),".tif"),
            format="GTiff",geoinfo=geoinfo)
}


rm(new.ras) #declutter RAM


#Supervised classification of PCA of original raster stack
for (x in 1:rep){
  cat(paste0("\n ","Running for x= ",x," of ",rep," for PCA7 repetitive classification."))
  set.seed(seedno)
  ptm <- proc.time() #start timer
  sc.valid <- superClass(pca.ras$map,trainData=train,
                         valData=valid,
                         responseCol="class_name",
                         model=sc.model
  )
  timing <-  proc.time() - ptm #end timer
  compile <- as.list(c(timing[3],
                       sc.valid$validation$performance$overall[1:2], #overall accuracy & kappa
                       "PCA7","null"))
  names(compile) <- c(names(timing[3]),
                      names(sc.valid$validation$performance$overall[1:2]),
                      "n_bands","bands.removed")
  meta <- dplyr::bind_rows(meta,compile)
}
writeRaster(sc.valid$map,filename=paste0("output/map/","class_","PCA_original",".tif"),
            format="GTiff",geoinfo=geoinfo)


rm(pca.ras) #declutter RAM


#Raster PCA of the inflated raster stack
pca.ras2 <- rasterPCA(ras,maskCheck=TRUE)
  summary(pca.ras2$model) #summary of data
  loadings(pca.ras2$model)
rm(ras) #declutter RAM
pca2.name <- paste0("PCA",maxband)


#Supervised classification for PCA of the inflated raster stack
for (x in 1:rep){
  cat(paste0("\n ","Running for x= ",x," of ",rep," for PCA49 repetitive classification."))
  set.seed(seedno)
  ptm <- proc.time() #start timer
  sc.valid <- superClass(pca.ras2$map,trainData=train,
                         valData=valid,
                         responseCol="class_name",
                         model=sc.model
  )
  timing <-  proc.time() - ptm #end timer
  compile <- as.list(c(timing[3],
                       sc.valid$validation$performance$overall[1:2], #overall accuracy & kappa
                       pca2.name,"null"))
  names(compile) <- c(names(timing[3]),
                      names(sc.valid$validation$performance$overall[1:2]),
                      "n_bands","bands.removed")
  meta <- dplyr::bind_rows(meta,compile)
}
writeRaster(sc.valid$map,filename=paste0("output/map/","class_","PCA_expanded",".tif"),
            format="GTiff",geoinfo=geoinfo)


#Pre-processing of compiled metadata
meta$elapsed <- as.numeric(meta$elapsed) #attribute change to numeric
meta$Accuracy<- as.numeric(meta$Accuracy) #attribute change to numeric
write.csv(meta,file=paste0("output/meta_data_",seedno,".csv"),col.names = TRUE)


#Doing boxplot descriptive statistics
meta_md <- meta[-which(meta$elapsed>250),] #hard coded removal anomalous result
meta_md$n_bands <- factor(meta_md$n_bands,
                         levels=c(benchmarks,
                                  rev(unique(meta_md$n_bands)[-match(benchmarks,unique(meta_md$n_bands))])))
benchmarks <- c("benchmark","PCA7",pca2.name)
graphics.par <- function(data.box,formula.box,title.box,ylab.box,#create function for boxplot
                         xlab.off,lab.off.b,lab.off.p1,lab.off.p2){ 
  box.output <- graphics::boxplot(formula=formula.box,data=data.box,las=2,
                                  main=title.box,xlab="No. of bands",ylab=ylab.box,
                                  cex.label=1.3,cex.main=1.8,axes=FALSE)
  axis(2,cex.axis=1.3)
  axis(1,cex.axis=1.3,at=c(4,11,21,31,41,51,61,70),
       labels=c("3","10","20","30","40","50","60","69"))
  rect(xleft=0.4,xright=3.5,ybottom=min((box.output$stats[1,])*0.999), 
       ytop=(max(box.output$stats[5,])*1.002),#note that $stats do not include anomalous points
       col="grey",density=20)
  abline(h=box.output$stats[3,match(benchmarks,box.output$names)],
         col=c("red","darkgreen","blue"),lty=5,lwd=2)
  text(xlab.off,y=(box.output$stats[3,match(benchmarks,box.output$names)]+
               c(lab.off.b,lab.off.p1,lab.off.p2)),
       labels=benchmarks,col=c("red","darkgreen","blue"),cex=1)
}
#Processing time
windows(width=12)
  graphics.par(data.box=meta_md,formula.box=elapsed~n_bands,
               title.box="Processing Time",ylab.box="Time (s)",
               xlab.off=68,lab.off.b = 4.5,lab.off.p1 = 4.5,lab.off.p2 = 4.5)
  savePlot(filename="output/stats/processingtime_boxplot.jpg",type="jpeg")
#Overall accuracy
windows(width=12)
  graphics.par(data.box=meta_md,formula.box=Accuracy~n_bands,
               title.box="Overall Accuracy",ylab.box="Accuracy (%)",
               xlab.off=1,lab.off.b = -0.0015,lab.off.p1 = 0.002,lab.off.p2 = 0.002)
  savePlot(filename="output/stats/accuracy_boxplot.jpg",type="jpeg")
#Ratio of overall accuracy over processing time
windows(width=12)
  graphics.par(data.box=meta_md,formula.box=Accuracy/elapsed~n_bands,
               title.box="Ratio of Overall Accuracy over Processing Time",
               ylab.box=expression(paste("Accuracy / Time (% ",s^-1,")")),
               xlab.off=68,lab.off.b = 0.00025,lab.off.p1 = 0.00025,lab.off.p2 = -0.00025)
  savePlot(filename="output/stats/ratioaccuracytime_boxplot.jpg",type="jpeg")  
  
  
#MANOVA for overall accuracy and processing time
fit.man <- manova(cbind(meta$elapsed,meta$Accuracy)~meta$n_bands)
summary.aov(fit.man) #classification difference with stack size
  sink("output/stats/MANOVA_summary.txt") #saving R console output
  summary.aov(fit.man) 
  sink()


##Linear regression
meta.position <- c(which(meta$n_bands=="benchmark"),which(meta$n_bands=="PCA7"),
                   which(meta$n_bands==pca2.name))
meta2 <- meta[-meta.position,]
meta2$n_bands <- as.numeric(meta2$n_bands)
lm.ras <- lm(cbind(meta2$elapsed,meta2$Accuracy)~
               meta2$n_bands,data=meta2)
summary(lm.ras) #interpretation based on gradient/coefficient of n_bands
  sink("output/stats/Regression_summary.txt") #saving R console output
  summary(lm.ras)
  sink()

#Percentage increase in processing time and accuracy per band added to stack
coeff.time <- lm.ras$coefficients[2,1] #time was the first column in rbind
coeff.acc <- lm.ras$coefficients[2,2]
prop.time <- coeff.time/mean(meta$elapsed[meta$n_bands==3])*100 #smallest stack size = 3
prop.acc <- coeff.acc/mean(meta$Accuracy[meta$n_bands==3])*100
