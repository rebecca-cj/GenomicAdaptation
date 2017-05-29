###################################################
#                                                 #
#  Allele freq changes associated with predicted  #
#           future climate changes                #
#    ----------------------------------------     #
#   Using abs climate predicitions data from      #   
#   Briscoe et al 2016 (Global Change Biol,       #
#   22:2425-2439;  baseline AWAP 1990-2009)       #
#                   &                             #
#    linear slope only to calculate allele        #
#            frequency change                     #
#                                                 #
#   Jan 2017  Rebecca Jordan                      #
#                                                 #
###################################################

# This analysis uses the absolute future climate projections from AWAP data as "Future" climate 
# estimate (data from Briscoe et al 2016, Global Change Biology, 22:2425-2439) and Atlas of Living
# Australia (www.ala.org.au) data as the "Current" baseline (data that was used in enviro association analysis)
# "Climate change" therefore = AWAP Future projection - ALA current prediction

load("CC_allele_frq.RData")

##--- LIBRARIES ---##
library(geosphere)
library(scales)

#-------------------#
#    IMPORT DATA    #
#-------------------#

#-- Environmental data for 10 enviro1-subset2 variables, by population

enviro.data <- t(read.table("~/Documents/GBSanalysis_2/Outlier05_Apr2016/envirodata/enviro1-subset2.txt",header=T, row.names=1))

#-- Alternate allele frequency, by population, for 81 outliers (2+/4 tests) with strong+ association (BF>20) with at least one of 10 enviro1-subset2 environmental variables

outlier.frq <- read.table("~/Documents/GBSanalysis_2/Outlier07_Oct2016/allelefrq/comout2plus.2ndfilt.strongenviro.frq", header=T, row.names = 1)
outlier.frq <- outlier.frq[,c(12,1:10,14:39)]
colnames(outlier.frq)[2:11] <- rownames(enviro.data)

#-- Site latitude and longitude 
tmp.latlongs<-read.csv("~/Documents/Fieldwork 2013/PointsForPlotting/Nat_Outlier_Testsites.csv",header=T,stringsAsFactors = F)
site.latlongs<-tmp.latlongs[2:3]
row.names(site.latlongs) <- tmp.latlongs[,"Site"]
rm(tmp.latlongs)

#-- Projected temperature and rainfall (2050 & 2070), based on calculations from 1990-2009 AWAP data

MAT.pred <- read.table("~/Documents/GBSanalysis_2/Outlier07_Oct2016/ClimateChangeData/AWAP_envirodata/MAT_predictions_MK12Nov2016.txt", header=T, row.names=1, stringsAsFactors = F)
SumRain.pred <- read.table("~/Documents/GBSanalysis_2/Outlier07_Oct2016/ClimateChangeData/AWAP_envirodata/SumRain_predictions_MK12Nov2016.txt", header=T, row.names=1, stringsAsFactors =F)
WinRain.pred <- read.table("~/Documents/GBSanalysis_2/Outlier07_Oct2016/ClimateChangeData/AWAP_envirodata/WinRain_predictions_MK12Nov2016.txt", header=T, row.names=1, stringsAsFactors =F)

  #-- get projected change based on 'future' AWAP climate projections and 'current' ALA climate data
  #     absolute change for Temperature, % change for rainfall
MAT.change <- MAT.pred - enviro.data["bio01",]
SumRain.change <- ((SumRain.pred / enviro.data["rain_sum",])*100 ) - 100
WinRain.change <- ((WinRain.pred / enviro.data["rain_win",])*100 ) - 100

summary.change <- cbind( "Tempmean" = apply(MAT.change, 2, mean), "Tempmin" = apply(MAT.change, 2, min), "Tempmax" = apply(MAT.change, 2, max),
       "SumRainmean" = apply(SumRain.change, 2, mean), "SumRainmin" = apply(SumRain.change, 2, min), "SumRainmax" = apply(SumRain.change, 2, max),
       "WinRainmean" = apply(WinRain.change, 2, mean), "WinRainmin" = apply(WinRain.change, 2, min), "WinRainmax" = apply(WinRain.change, 2, max))
write.table(summary.change, file = "summary_projected_future_change.txt", quote = F, sep = "\t")

summary.change.sd <- cbind( "TempMean" = apply(MAT.change, 2, mean), "TempSD" = apply(MAT.change, 2, sd), 
                         "SumRainMean" = apply(SumRain.change, 2, mean), "SumRainSD" = apply(SumRain.change, 2, sd),
                         "WinRainMean" = apply(WinRain.change, 2, mean), "WinRainSD" = apply(WinRain.change, 2, sd))
write.table(summary.change.sd, file = "summary_projected_future_change_sd.txt", quote = F, sep = "\t")

summary.boxplot.data <- cbind( "MAT" = MAT.change[,c(1:2,7:8,11:12)], "SumRain" = SumRain.change[,c(1:2,7:8,11:12)], "WinRain" = WinRain.change[,c(1:2,7:8,11:12)])
summary.boxplot.data <- summary.boxplot.data[,c(1,5,3,2,6,4,7,11,9,8,12,10,13,17,15,14,18,16)]

summary.boxplot.data.2070 <- summary.boxplot.data[,c(4:6,10:12,16:18)]
pdf("FutureClimateChange_Boxplots_2070.pdf", width = 9, height = 6)
tiff("FutureClimateChange_Boxplots_2070.tiff", width = 9, height = 6, units = "in", res = 600)
{
  par(mar = c(2,4,0.5,0.5), fig = c(0,0.35,0,1))
  boxplot(summary.boxplot.data.2070[,1:3], col = c("red", "blue", "goldenrod"), names = NA, axes = F, ylim = c(0,3.5), xlim = c(0,4), range = 0)
  axis(2, at = seq(0,3.5,0.5), las = 2, lwd = 2, cex.axis = 1.5, pos = 0.25)
  axis(1, at = c(0.5,3.5), labels = NA, pos = 0, lwd = 2)
  mtext("Mean Annual Temp.", 1, line = 0.5, cex = 1.5)
  mtext(bquote("Temperature Change ("*degree*"C)"), 2, cex = 1.5, line = 2.25)
  par(mar = c(2,4,0.5,0.5), fig = c(0.35,1,0,1), new = T)
  boxplot(summary.boxplot.data.2070[,4:9], col = c("red", "blue", "goldenrod"), names = NA, at = c(1:3,5:7), axes = F, ylim = c(-30,30), xlim = c(0,8), range = 0)
  axis(2, at = seq(-30,30,5), las = 2, lwd = 2, cex.axis = 1.5, pos = 0.25)
  axis(1, at = c(0.5,3.5), labels = NA, pos = -30, lwd = 2)
  axis(1, at = c(4.5,7.5), labels = NA, pos = -30, lwd = 2)
  mtext("Summer", 1, line = -0.5, cex = 1.5, at = 2)
  mtext("Winter", 1, line = -0.5, cex = 1.5, at = 6)
  mtext("Precipitation", 1, line = 0.8, cex = 1.5, at = 4)
  mtext("Precipitation Change (%)", 2, line = 2, cex = 1.5)
  legend(1,30, c("ACCESS 1.0", "HadGEM2-ES", "GDFL-CM3"), pch = rep(15,3), col = c("red","blue","goldenrod"), cex = 1, pt.cex = 1.2)
}
dev.off()


#-- Site plotting colours

plot.col <- read.table("../plot_colours.txt", header=T, row.names = 1, stringsAsFactors = F)

#---------------------------------#
#    GEOGRAPHIC DISTANCE MATRIX   #
#---------------------------------#

#-- create geographic distance matrix using 'geosphere'

natsite.dist <- distm(site.latlongs[,1:2], fun=distVincentyEllipsoid)/1000
geog.dist <- as.matrix(natsite.dist)
colnames(geog.dist) <- row.names(site.latlongs)
rownames(geog.dist) <- row.names(site.latlongs)
diag(geog.dist) <- NA


#--------------------------------#
#    EXPLORE PREDICTED CHANGE    #
#--------------------------------#

# plot predicted change for all 6 models per site to get idea of variation between models

#-- MAT
tmp <- cbind("MAT" = enviro.data["bio01",], MAT.change)
plot.data <- tmp[order(tmp$MAT, decreasing = T),]
maxY <- ceiling(max(plot.data[,2:13])/0.5)*0.5
minY <- floor(min(plot.data[,2:13])/0.5)*0.5
rm(tmp)

plot(1:26, plot.data[,2], pch = 1, col = "blue", lwd = 1.25, ylim = c(minY, maxY), axes = F, main = expression("Change in Mean Annual Temperature ("*degree*"C)"),
     xlab = "Sites (Left to Right = Warmest to Coolest)", ylab = expression("Change from 'current' (ALA) to 'future' (AWAP) climate ("*degree*")"))
axis(1, at=seq(1,26), labels = row.names(plot.data), las = 2)
axis(2, las = 2)
pchX = 2
for( i in seq(4,12,2)){
  points(1:26, plot.data[,i], pch = pchX, col = "blue", lwd = 1.25)
  pchX = pchX + 1
}
pchX = 1
for( i in seq(3,13,2)){
  points(1:26, plot.data[,i], pch = pchX, col = "red", lwd = 1.25)
  pchX = pchX + 1
}
legend(23, 3.5, c('Access1.0','Access1.3','CanESM2','GDFLCM3','HadGEM2-CC','HadGEM2-ES', '-------------', '2050', '2070'), pch = c(seq(1,6),rep(1,3)), 
       col = c(rep("black",6), "white", "blue", "red"), pt.lwd = 1.25, cex = 0.8)

#-- Summer Rainfall
tmp <- cbind("SumRain" = enviro.data["rain_sum",], SumRain.change)
plot.data <- tmp[order(tmp$SumRain, decreasing = T),]
maxY <- ceiling(max(plot.data[,2:13])/5)*5
minY <- floor(min(plot.data[,2:13])/5)*5
rm(tmp)

par(fig=c(0,1,0.5,1), mar = c(2,4,3.5,1))
plot(1:26, plot.data[,1], pch = 1, col = "blue", lwd = 1.25, ylim = c(minY, maxY), axes = F, 
     ylab = "% change from 'current' (ALA) to 'future' (AWAP) climate", xlab = NA, main = "Projected % change in Summer Rainfall")  
axis(1, at=seq(1,26), labels = NA)
axis(2, at=seq(minY, maxY, 5), las = 2)
abline(h = 0)
text(0.2, 30, "a) 2050", pos = 4)
pchX = 2
for( i in seq(3,11,2)){
  points(1:26, plot.data[,i], pch = pchX, col = "blue", lwd = 1.25)
  pchX = pchX + 1
}
par(fig=c(0,1,0,0.5), mar = c(5,4,0.5,1), new = T)
pchX = 1
plot(1:26, plot.data[,2], pch = 1, col = "red", lwd = 1.25, ylim = c(minY, maxY), axes = F,  
     xlab = "Sites (Left to Right = Wettest to Driest)", ylab = "% change from 'current' (ALA) to 'future' (AWAP) climate")
axis(1, at=seq(1,26), labels = row.names(plot.data), las = 2)
axis(2, at=seq(minY, maxY, 5), las = 2)
abline(h = 0)
text(0.2, 30, "b) 2070", pos = 4)
pchX = 2
for( i in seq(4,12,2)){
  points(1:26, plot.data[,i], pch = pchX, col = "red", lwd = 1.25)
  pchX = pchX + 1
}
legend(0.8, -15, c('Access1.0','Access1.3','CanESM2','GDFLCM3','HadGEM2-CC','HadGEM2-ES'), pch = seq(1,6), pt.lwd = 1.25, cex = 0.8)


#-- Winter Rainfall
tmp <- cbind("WinRain" = enviro.data["rain_win",], WinRain.change)
plot.data <- tmp[order(tmp$WinRain, decreasing = T),]
maxY <- ceiling(max(plot.data[,2:13])/5)*5
minY <- floor(min(plot.data[,2:13])/5)*5
rm(tmp)

par(fig=c(0,1,0.5,1), mar = c(2,4,3.5,1))
plot(1:26, plot.data[,2], pch = 1, col = "blue", lwd = 1.25, ylim = c(minY, maxY), axes = F, 
     ylab = "% change from 'current' (ALA) to 'future' (AWAP) climate", xlab = NA, main = "Projected % change in Winter Rainfall")  
axis(1, at=seq(1,26), labels = NA)
axis(2, at=seq(minY, maxY, 5), las = 2)
abline(h = 0)
text(0.2, maxY, "a) 2050", pos = 4)
pchX = 2
for( i in seq(4,12,2)){
  points(1:26, plot.data[,i], pch = pchX, col = "blue", lwd = 1.25)
  pchX = pchX + 1
}
legend(4, 40, c('Access1.0','Access1.3','CanESM2','GDFLCM3','HadGEM2-CC','HadGEM2-ES'), pch = seq(1,6), pt.lwd = 1.25, cex = 0.8)
par(fig=c(0,1,0,0.5), mar = c(5,4,0.5,1), new = T)
pchX = 1
plot(1:26, plot.data[,2], pch = 1, col = "red", lwd = 1.25, ylim = c(minY, maxY), axes = F,  
     xlab = "Sites (Left to Right = Wettest to Driest)", ylab = "% change from 'current' (ALA) to 'future' (AWAP) climate")
axis(1, at=seq(1,26), labels = row.names(plot.data), las = 2)
axis(2, at=seq(minY, maxY, 5), las = 2)
abline(h = 0)
text(0.2, maxY, "b) 2070", pos = 4)
pchX = 2
for( i in seq(3,13,2)){
  points(1:26, plot.data[,i], pch = pchX, col = "red", lwd = 1.25)
  pchX = pchX + 1
}


#-----------------------------#
#   MEAN ANNUAL TEMPERATURE   #
#-----------------------------#

##-- Minimum distance to closest population with future temperture --##

# using absolute future projected climate (from AWAP data)

min.dist.MAT <- function(CCmodel, plot.colours){
  # uses absolute future projected climate (from AWAP data)
  # input objects = site.latlongs, enviro.data, MAT.pred / MAT.change, geog.dist
  
  # create data.frame with current and predicted future MAT as well as
  # distance to nearest site with predicted future MAT
  output <- site.latlongs[,c("Latitude","Longitude")]
  row.names(output) <- rownames(site.latlongs)
  output[,"MAT"] <- enviro.data["bio01",]
  output[,CCmodel] <- MAT.pred[,CCmodel]  # uses absolute future climate prediction
  for(pops in rownames(output)){
    tmp <- cbind("MAT" = enviro.data["bio01",], "Dist" = geog.dist[pops,])
    min.geog.dist <- min(tmp[which(tmp[,"MAT"] >= output[pops,CCmodel]), 2])  # uses absolute future climate prediction
    if(is.finite(min.geog.dist)){
      output[pops,paste("MinDist_",CCmodel,sep="")] <- min.geog.dist
    } else {
      if(is.infinite(min.geog.dist)){
        output[pops,paste("MinDist_",CCmodel,sep="")] <- NA
    }}}    
  
  # plot
  pdf(paste("MinGeogDist_",CCmodel,".pdf",sep=""))
  plot.obj <- output[,c("MAT", paste("MinDist_",CCmodel,sep=""))]
  plot.obj[,"PlotCol"] <- plot.colours
  maxY <- ceiling(max(plot.obj[,paste("MinDist_",CCmodel,sep="")], na.rm=T)/50)*50
  minY <- floor(min(plot.obj[,paste("MinDist_",CCmodel,sep="")], na.rm=T)/50)*50
  plot(plot.obj$MAT, plot.obj[,paste("MinDist_",CCmodel,sep="")], xlab = "Mean Annual Temp (MAT; Bio01)", ylab = "Min. distance to site matching future MAT (km)", 
       pch = 16, cex = 1.5, col = plot.obj$PlotCol, axes = F, xlim = c(13.5,18), ylim = c(minY,(maxY+25)))
  axis(1, at=seq(13.5,18,0.5))
  axis(2, at=seq(minY,maxY,50))
  tmp <- plot.obj[which(is.na(plot.obj[,paste("MinDist_",CCmodel,sep="")])),]
  points(tmp$MAT, rep((maxY+25),nrow(tmp)), col = tmp[,"PlotCol"], pch = 3, cex = 1.5, lwd = 1.5)
  legend(13.5,(maxY+25), c("Central NSW","Southern NSW", "Central Vic", "Western Vic", "------------", "No site match"), 
         pch = c(rep(16,5), 3), col = c("red","goldenrod","green","blue","white","black"), cex = 0.8)
  dev.off()                        
  return(output)
}

min.dist.MAT.Access1.0.2050 <- min.dist.MAT("MAT_Access1.0_2050", plot.col[,2])
min.dist.MAT.Access1.0.2070 <- min.dist.MAT("MAT_Access1.0_2070", plot.col[,2])

min.dist.MAT.HadGEM2.ES.2050 <- min.dist.MAT("MAT_HadGEM2.ES_2050", plot.col[,2])
min.dist.MAT.HadGEM2.ES.2070 <- min.dist.MAT("MAT_HadGEM2.ES_2070", plot.col[,2])

min.dist.MAT.GDFLCM3.2050 <- min.dist.MAT("MAT_GDFLCM3_2050", plot.col[,2])
min.dist.MAT.GDFLCM3.2070 <- min.dist.MAT("MAT_GDFLCM3_2070", plot.col[,2])

##-- Allele frequency increase required to match future predicted climate --##

# using difference between AWAP 'future' and ALA 'current' climate data

delta.MAT.increase <- function(varX, CCmodel, CCprediction_df){
  # specificially designed for 'outlier.frq' data.frame structure
  # rows = loci
  # column 1 = alternate allele, column 2-11 = enviro variables
  # column 12-37 = alternate allele frequency for 26 populations in E. microcarpa outlier/EAA study
  
  # CCpredicition_df = change in temperature
  
  # for each SNP, determine direction of association between Mean Annual Temp and alt allele freq ( + / - )
  # and change to positive association (if negative) 
  # by converting alt allele freq to ref allele frequency
  
  #-- get loci with strong or very strong association for input environmental variable

  tmp <- outlier.frq[grep("strong", outlier.frq[,varX]),]
  
  #-- create new output dataframes (empty)
  output.current <- matrix(nrow=nrow(tmp), ncol=26)
  row.names(output.current) <- row.names(tmp)
  colnames(output.current) <- colnames(tmp)[12:37]
  
  output.change <- matrix(nrow=nrow(tmp), ncol=26)
  row.names(output.change) <- row.names(tmp)
  colnames(output.change) <- colnames(tmp)[12:37]
  
  #-- calculate allele freq change, using slope of linear model only.  Add data to data frames.
  for(snpX in rownames(tmp)){
    
    # calculate linear model for allele with postive relationship with environmental variable  
    tmp.frq <- tmp[snpX,12:37]
    tmp.lm <- lm(unlist(tmp.frq) ~ enviro.data[varX,])
    if(summary(tmp.lm)$coefficients[2,1] < 0){
      tmp.frq = 1-tmp.frq
      tmp.lm <- lm(unlist(tmp.frq) ~ enviro.data[varX,])
      }
      
    output.current[snpX,] <- as.matrix(tmp.frq)
    
    # get environmental coefficient (change in allele frq over 1 unit of environmental variable)
    enviro.coef <- summary(tmp.lm)$coefficients[2,1]
    
    # calculate change for different increases in environmental variable
    for(site in colnames(output.current)){
      output.change[snpX,site] <- enviro.coef * CCprediction_df[site,CCmodel]  # using slope only
      }
  }
  output.list <- list()
  output.list[["Current"]] <- output.current
  output.list[["FutureFrqChange"]] <- output.change
 
  #-- plot current frq vs change required, by population 
  pdf(paste("AlleleFrqChange_bypop_",CCmodel,".pdf",sep=""))
  x.max <- ceiling(max(output.list$FutureFrqChange)/0.02)*0.02
  x.min <- floor(min(output.list$FutureFrqChange)/0.02)*0.02
  for(i in colnames(output.list$Current)){
    temp_change <- CCprediction_df[i,CCmodel]  # used projected change
    plot(output.list$FutureFrqChange[,i], output.list$Current[,i], pch = 16, ylim = c(0,1), xlim = c(x.min,x.max), ylab = "Current allele frequency", xlab = bquote( Delta~"allele freq associated with projected MAT increase ["*.(round(temp_change,2))*degree~"C]"), main = paste(i, CCmodel, sep="   "))
    for (j in rownames(output.list$Current)){
      # "mark" allele not currently present in population
      if( output.list$Current[j,i] == 0 ){
        points(output.list$FutureFrqChange[j,i], output.list$Current[j,i], pch = 15, col = "red", cex = 1.2)
      }
      # "mark" allele where future increase = fixation
      if( output.list$Current[j,i] + output.list$FutureFrqChange[j,i] > 1){
        points(output.list$FutureFrqChange[j,i], output.list$Current[j,i], pch = 15, col = "goldenrod", cex = 1.2)
      }}
  }
  dev.off()
  
  return(output.list)
}

frq.change.MAT.Access1.0.2050 <- delta.MAT.increase("bio01","MAT_Access1.0_2050", MAT.change)
frq.change.MAT.Access1.0.2070 <- delta.MAT.increase("bio01","MAT_Access1.0_2070", MAT.change)

frq.change.MAT.HadGEM2.ES.2050 <- delta.MAT.increase("bio01","MAT_HadGEM2.ES_2050", MAT.change)
frq.change.MAT.HadGEM2.ES.2070 <- delta.MAT.increase("bio01","MAT_HadGEM2.ES_2070", MAT.change)

frq.change.MAT.GDFLCM3.2050 <- delta.MAT.increase("bio01","MAT_GDFLCM3_2050", MAT.change)
frq.change.MAT.GDFLCM3.2070 <- delta.MAT.increase("bio01","MAT_GDFLCM3_2070", MAT.change)

##-- Count alleles reaching fixation or frq not found in current (sampled) populations --##

# using allele frequence change associated with difference between AWAP 'future' and ALA 'current' climate data

mal.adapt.count <- function(dataIN, varX, plotcolours, CCmodel){
  # dataIN = list created by "delta.MAT.increase" function containing two lists
  # 1 = "Current" = current allele freq
  # 2 = "FutureFrqChange" = change in allele freq associated with climate change
  
  output <- matrix(ncol = 26, nrow = 6)
  row.names(output) <- c("fixation", "highFrq", "absent","TotalNoLoci", varX, "PlotCol")
  colnames(output) <- colnames(dataIN$Current)
  output["TotalNoLoci",] <- rep(nrow(dataIN$Current),26)
  output[varX,] <- enviro.data[varX,]
  output["PlotCol",] <- plotcolours
  
  for( site in colnames(dataIN$Current)){
    fix.allele = 0
    high.allele = 0
    absent.allele = 0
    for( locus in row.names(dataIN$Current) ){
      if( dataIN$Current[locus,site] == 0){
        absent.allele = absent.allele + 1
      }
      loc.increase <- dataIN$Current[locus,site] + dataIN$FutureFrqChange[locus,site]
      if( loc.increase > 1 ){
        fix.allele = fix.allele + 1
      } else {
        if (loc.increase > max(dataIN$Current[locus,])){
          high.allele = high.allele + 1
        }}}
    output["fixation",site] <- fix.allele
    output["highFrq",site] <- high.allele
    output["absent",site] <- absent.allele
  }

  # plot
  maxX = max(as.numeric(output[1:2,]))
  y.limit = ceiling(maxX/5)*5
  maxB = max(as.numeric(output["absent",]))
  
  pdf(paste("MalAdaptCount_",CCmodel,".pdf",sep=""))
  par(mar=c(5,4,1,1),fig=c(0,1,0,0.75))
  plot(output[varX,], output["fixation",], pch = 3, col = output["PlotCol",], lwd = 1.5,
       xlab = "Mean Annual Temperature (Bio01)", ylab = "No. of loci (Fixation or High Frequency)", ylim = c(0, y.limit), xlim = c(13,18), axes=F)
  points(output[varX,], output["highFrq",], pch = 16, col = output["PlotCol",])
  axis(1)
  axis(2, las = 2)
  text(13, y.limit, labels = "b)", cex = 1.25)
  legend(13, y.limit-3, c("Central NSW","Southern NSW", "Central Vic", "Western Vic"), pch = rep(16,4), col = c("red","goldenrod","green","blue"), cex = 0.7, pt.cex=0.9)
  legend(13, y.limit-10, c("Fixation", "High Frq"), pch = c(3,16), cex = 0.7, pt.cex = 0.9)
  par(mar=c(2,4,1,1), fig=c(0,1,0.7,1), new=T)
  plot(output[varX,], output["absent",], pch = 17, col = output["PlotCol",], xlab = NA, ylab = "No. of loci (Absent)",
       xlim = c(13,18), axes = F)
  axis(1, labels = NA)
  axis(2, las = 2)
  text(13, maxB, labels = "a)", cex = 1.25)
  dev.off()

  return( output )
}

mal.count.MAT.Access1.0.2050 <- mal.adapt.count(frq.change.MAT.Access1.0.2050, "bio01", plot.col[,2], "MAT_Access1.0_2050")
mal.count.MAT.Access1.0.2070 <- mal.adapt.count(frq.change.MAT.Access1.0.2070, "bio01", plot.col[,2], "MAT_Access1.0_2070")

mal.count.MAT.HadGEM2.ES.2050 <- mal.adapt.count(frq.change.MAT.HadGEM2.ES.2050, "bio01", plot.col[,2], "MAT_HadGEM2.ES_2050")
mal.count.MAT.HadGEM2.ES.2070 <- mal.adapt.count(frq.change.MAT.HadGEM2.ES.2070, "bio01", plot.col[,2], "MAT_HadGEM2.ES_2070")

mal.count.MAT.GDFLCM3.2050 <- mal.adapt.count(frq.change.MAT.GDFLCM3.2050, "bio01", plot.col[,2], "MAT_GDFLCM3_2050")
mal.count.MAT.GDFLCM3.2070 <- mal.adapt.count(frq.change.MAT.GDFLCM3.2070, "bio01", plot.col[,2], "MAT_GDFLCM3_2070")
  

#------------------------------------------------#
#    SEASONAL PRECIPITATION (Summer or Winter)   #
#------------------------------------------------#

##-- Minimum distance to closest population with future seasonal rainfall --##

# using absolute future projected climate (from AWAP data)

min.dist.rain <- function(varX, CCprediction_df, CCmodel, plot.colours){
  # uses absolute future projected climate (from AWAP data)
  # input objects = site.latlongs, enviro.data, xxxRain.change, geog.dist
  
  # create data.frame with current and predicted future rainfall,
  # distance to nearest site with predicted future rainfall, and
  # whether the predicted change is an increase or decrease in rainfall
  
  #varX = "rain_sum"
  #CCprediction_df = SumRain.pred
  #CCmodel = "SumRain_Access1.0_2050"
  #plot.colours = plot.col[,2]
  
  output <- site.latlongs[,c("Latitude","Longitude")]
  row.names(output) <- rownames(site.latlongs)
  output[,varX] <- enviro.data[varX,]
  output[,CCmodel] <- CCprediction_df[,CCmodel]
  
  for(pops in rownames(output)){
    #pops = "NAL"
    tmp <- cbind(enviro.data[varX,], "Dist" = geog.dist[pops,])
    colnames(tmp)[1] = varX
    if( output[pops,CCmodel] < output[pops,varX] ){
      # decrease in rainfall
      output[pops,"Direction"] <- "decrease"
      min.geog.dist <- min(tmp[which(tmp[,varX] <= output[pops,CCmodel]), 2]) * -1  # "negative" distance in indicate movement toward decreased rainfall
    } else {
      if( output[pops,CCmodel] > output[pops,varX]){ 
        # increase in rainfall
        output[pops,"Direction"] <- "increase"
        min.geog.dist <- min(tmp[which(tmp[,varX] >= output[pops,CCmodel]), 2])
    }}
    
    if(is.finite(min.geog.dist)){
      output[pops,paste("MinDist_",CCmodel,sep="")] <- min.geog.dist
    } else {
      if(is.infinite(min.geog.dist)){
        output[pops,paste("MinDist_",CCmodel,sep="")] <- NA
      }}}    

  # plot
  pdf(paste("MinGeogDist_",CCmodel,".pdf",sep=""))
  plot.obj <- output[,c(varX, paste("MinDist_",CCmodel,sep=""), "Direction")]
  plot.obj[,"Direction"] <- as.factor(plot.obj[,"Direction"])
  plot.obj[,"PlotCol"] <- plot.colours
  
  maxY <- ceiling(max(plot.obj[,paste("MinDist_",CCmodel,sep="")], na.rm=T)/25)*25
  minY <- floor(min(plot.obj[,paste("MinDist_",CCmodel,sep="")], na.rm=T)/25)*25
  maxX <- ceiling(max(plot.obj[,varX])/10)*10
  minX <- floor(min(plot.obj[,varX])/10)*10
  bufferX <- (maxX-minX)*0.3  # to create extra space on X axis for labels and legend
  if(varX == "rain_sum"){plot.var = "Summer"} else {if(varX == "rain_win"){plot.var = "Winter"}}
  
  par(mar = c(5,4,2.5,0.5))
  plot(plot.obj[,varX], plot.obj[,paste("MinDist_",CCmodel,sep="")], xlab = paste(plot.var, " Rain (mm)", sep=""), ylab = "Distance (km)",
       pch = 16, cex = 1.5, col = plot.obj$PlotCol, axes = F, xlim = c(minX,maxX+bufferX), ylim = c((minY-10),(maxY+10)), lwd = 1.5)
  axis(1, at=seq(minX,maxX,10))
  axis(2, at=seq(minY,maxY,25), las = 2, cex = 0.8)
  abline(h = 0)
  mtext(paste("Min. distance to site with predicted future ", plot.var, " rainfall",sep=""), side = 3, line = 1, cex = 1.25, font = 2)
  mtext(paste("[",CCmodel,"]",sep=""), cex = 1, line = 0)
  text(maxX, 10, "Increased rainfall", pos = 4, cex = 1)
  text(maxX, -10, "Decreased rainfall", pos = 4, cex = 1)
  
  tmp.increase.NA <- plot.obj[which(is.na(plot.obj[,paste("MinDist_",CCmodel,sep="")]) & plot.obj$Direction == "increase"),]
  tmp.decrease.NA <- plot.obj[which(is.na(plot.obj[,paste("MinDist_",CCmodel,sep="")]) & plot.obj$Direction == "decrease"),]
  points(tmp.increase.NA[,varX], rep((maxY+10),nrow(tmp.increase.NA)), col = tmp.increase.NA$PlotCol, pch = 3, cex = 1.25, lwd = 1.5)
  points(tmp.decrease.NA[,varX], rep((minY-10),nrow(tmp.decrease.NA)), col = tmp.decrease.NA$PlotCol, pch = 3, cex = 1.25, lwd = 1.5)
  legend(maxX+5,(maxY+10), c("Central NSW","Southern NSW", "Central Vic", "Western Vic", "------------", "No site match"), 
         pch = c(rep(16,5), 3), col = c("red","goldenrod","green","blue","white","black"), cex = 0.8, pt.cex = 1)
  dev.off()                        
  return(output)
}

min.dist.SumRain.Access1.0.2050 <- min.dist.rain("rain_sum", SumRain.pred, "SumRain_Access1.0_2050", plot.col[,2])
min.dist.SumRain.HadGEM2.ES.2050 <- min.dist.rain("rain_sum", SumRain.pred, "SumRain_HadGEM2.ES_2050", plot.col[,2])
min.dist.SumRain.GDFLCM3.2050 <- min.dist.rain("rain_sum", SumRain.pred, "SumRain_GDFLCM3_2050", plot.col[,2])

min.dist.SumRain.Access1.0.2070 <- min.dist.rain("rain_sum", SumRain.pred, "SumRain_Access1.0_2070", plot.col[,2])
min.dist.SumRain.HadGEM2.ES.2070 <- min.dist.rain("rain_sum", SumRain.pred, "SumRain_HadGEM2.ES_2070", plot.col[,2])
min.dist.SumRain.GDFLCM3.2070 <- min.dist.rain("rain_sum", SumRain.pred, "SumRain_GDFLCM3_2070", plot.col[,2])

min.dist.WinRain.Access1.0.2050 <- min.dist.rain("rain_win", WinRain.pred, "WinRain_Access1.0_2050", plot.col[,2])
min.dist.WinRain.HadGEM2.ES.2050 <- min.dist.rain("rain_win", WinRain.pred, "WinRain_HadGEM2.ES_2050", plot.col[,2])
min.dist.WinRain.GDFLCM3.2050 <- min.dist.rain("rain_win", WinRain.pred, "WinRain_GDFLCM3_2050", plot.col[,2])

min.dist.WinRain.Access1.0.2070 <- min.dist.rain("rain_win", WinRain.pred, "WinRain_Access1.0_2070", plot.col[,2])
min.dist.WinRain.HadGEM2.ES.2070 <- min.dist.rain("rain_win", WinRain.pred, "WinRain_HadGEM2.ES_2070", plot.col[,2])
min.dist.WinRain.GDFLCM3.2070 <- min.dist.rain("rain_win", WinRain.pred, "WinRain_GDFLCM3_2070", plot.col[,2])

##-- Allele frequency change required to match future predicted climate --##

# using difference between AWAP 'future' and ALA 'current' climate data

delta.Rain.increase <- function(varX, CCmodel, CCprediction_df){
  # specificially designed for 'outlier.frq' data.frame structure
  # rows = loci
  # column 1 = alternate allele, column 2-11 = enviro variables
  # column 12-37 = alternate allele frequency for 26 populations in E. microcarpa outlier/EAA study
  
  # CCpredicition_df = absolute future projected rainfall (absolute difference between this and 'current' calculated in script below)
  
  # for each SNP, determine direction of association between Precipitation and alt allele freq ( + / - )
  # and change to positive association (if negative) 
  # by converting alt allele freq to ref allele frequency
  
  #varX = "rain_sum"
  #CCmodel = "SumRain_Access1.0_2050"
  #CCpredicition_df = SumRain.pred
  
    #-- get loci with strong or very strong association for input environmental variable
  tmp <- outlier.frq[grep("strong", outlier.frq[,varX]),]
  
  #-- create new output dataframes (empty)
  output.current <- matrix(nrow=nrow(tmp), ncol=26)
  row.names(output.current) <- row.names(tmp)
  colnames(output.current) <- colnames(tmp)[12:37]
  
  output.change <- matrix(nrow=nrow(tmp), ncol=26)
  row.names(output.change) <- row.names(tmp)
  colnames(output.change) <- colnames(tmp)[12:37]
  
  #-- calculate allele freq change, using slope of linear model only.  Add data to data frames.
  for(snpX in rownames(tmp)){
    
    # calculate linear model for allele with positive relationship with rainfall variable
    tmp.frq <- tmp[snpX,12:37]
    tmp.lm <- lm(unlist(tmp.frq) ~ enviro.data[varX,])
    if(summary(tmp.lm)$coefficients[2,1] < 0){
      tmp.frq = 1-tmp.frq
      tmp.lm <- lm(unlist(tmp.frq) ~ enviro.data[varX,])
    }
    
    output.current[snpX,] <- as.matrix(tmp.frq)
    
    # get environmental coefficient (change in allele frq over 1 unit of environmental variable)
    enviro.coef <- summary(tmp.lm)$coefficients[2,1]
    
    # calculate change for different increases in environmental variable
    for(site in colnames(output.current)){
      #rainfall.change <- CCpred_df - enviro.data[varX,site]
      output.change[snpX,site] <- enviro.coef * (CCprediction_df[site,CCmodel] - enviro.data[varX,site])
    }
  }
  output.list <- list()
  output.list[["Current"]] <- output.current
  output.list[["FutureFrqChange"]] <- output.change
  
  # plot current frq vs change required, by population 
  pdf(paste("AlleleFrqChange_bypop_",CCmodel,"_increase.pdf",sep=""), width = 6, height = 4)
  for(i in colnames(output.list$Current)){
    rainfall_change <- ((CCprediction_df[i,CCmodel] / enviro.data[varX,i])*100)-100
    x.max <- ceiling(max(output.list$FutureFrqChange)*100)/100
    x.min <- floor(min(output.list$FutureFrqChange)*100)/100
    par(mar=c(3,2,2,1))
    plot(output.list$FutureFrqChange[,i], output.list$Current[,i], pch = 16, ylim = c(0,1.04), xlim = c(x.min,x.max), axes=F,
         ylab = NA, xlab = NA,main = paste(i, CCmodel, sep="   "))
    mtext("Current allele frequency (wet-associated allele)", side = 2, line = 0.5)
    mtext(bquote( Delta~"freq wet-associated allele for predicted rainfall change ["*.(round(rainfall_change,2))*"%]"), side = 1, line = 1.75)
    axis(1, at = seq(x.min,x.max,0.02), pos = 0)
    axis(2, pos = 0, at=seq(0,1,0.2), labels=c(NA,0.2,0.4,0.6,0.8,1), las = 2)
    if(x.min < 0){
      text(x.min, 1.04, "Increase dry-associated allele", pos = 4, cex = 0.8)
    }
    if(x.max > 0){
      text(x.max, 1.04, "Increase wet-associated allele", pos = 2, cex = 0.8)
    }
    for (j in rownames(output.list$Current)){
      # "mark" allele not currently present in population
      if( output.list$Current[j,i] == 0 || output.list$Current[j,i] == 1){
        points(output.list$FutureFrqChange[j,i], output.list$Current[j,i], pch = 15, col = "red", cex = 1.2)
      } else {
      # "mark" allele where future increase = fixation (excluding alleles already fixed)
      if( output.list$Current[j,i] + output.list$FutureFrqChange[j,i] > 1 || output.list$Current[j,i] - output.list$FutureFrqChange[j,i] < 0){
          points(output.list$FutureFrqChange[j,i], output.list$Current[j,i], pch = 15, col = "goldenrod", cex = 1.2)
      }}}} 
  dev.off()
  
  return(output.list)
}

frq.change.SumRain.Access1.0.2050 <- delta.Rain.increase("rain_sum","SumRain_Access1.0_2050", SumRain.pred)
frq.change.SumRain.HadGEM2.ES.2050 <- delta.Rain.increase("rain_sum","SumRain_HadGEM2.ES_2050", SumRain.pred)
frq.change.SumRain.GDFLCM3.2050 <- delta.Rain.increase("rain_sum", "SumRain_GDFLCM3_2050", SumRain.pred)

frq.change.SumRain.Access1.0.2070 <- delta.Rain.increase("rain_sum","SumRain_Access1.0_2070", SumRain.pred)
frq.change.SumRain.HadGEM2.ES.2070 <- delta.Rain.increase("rain_sum","SumRain_HadGEM2.ES_2070", SumRain.pred)
frq.change.SumRain.GDFLCM3.2070 <- delta.Rain.increase("rain_sum", "SumRain_GDFLCM3_2070", SumRain.pred)

frq.change.WinRain.Access1.0.2050 <- delta.Rain.increase("rain_win","WinRain_Access1.0_2050", WinRain.pred)
frq.change.WinRain.HadGEM2.ES.2050 <- delta.Rain.increase("rain_win","WinRain_HadGEM2.ES_2050", WinRain.pred)
frq.change.WinRain.GDFLCM3.2050 <- delta.Rain.increase("rain_win", "WinRain_GDFLCM3_2050", WinRain.pred)

frq.change.WinRain.Access1.0.2070 <- delta.Rain.increase("rain_win","WinRain_Access1.0_2070", WinRain.pred)
frq.change.WinRain.HadGEM2.ES.2070 <- delta.Rain.increase("rain_win","WinRain_HadGEM2.ES_2070", WinRain.pred)
frq.change.WinRain.GDFLCM3.2070 <- delta.Rain.increase("rain_win", "WinRain_GDFLCM3_2070", WinRain.pred)

##-- Count alleles reaching fixation or frq not found in current (sampled) populations --##

# using allele frequence change associated with difference between AWAP 'future' and ALA 'current' climate data

rain.mal.adapt.count <- function(dataIN, varX, plotcolours, CCmodel){
  # dataIN = list created by "delta.Rain.increase" function containing two lists
  # 1 = "Current" = current allele freq
  # 2 = "FutureFrqChange" = change in wet-associated allele freq with climate change (includes + and - changes depending on direction of rainfall change)
  
  output <- matrix(ncol = 26, nrow = 6)
  row.names(output) <- c("change_to_fixation", "highFrq", "currently_fixed","TotalNoLoci", varX, "PlotCol")
  colnames(output) <- colnames(dataIN$Current)
  output["TotalNoLoci",] <- rep(nrow(dataIN$Current),26)
  output[varX,] <- enviro.data[varX,]
  output["PlotCol",] <- plotcolours
  
  for( site in colnames(dataIN$Current)){
    become.fixed.allele = 0
    high.allele = 0
    fixed.allele = 0
    for( locus in row.names(dataIN$Current) ){
      if( dataIN$Current[locus,site] == 0 || dataIN$Current[locus,site] == 1 ){
        fixed.allele = fixed.allele + 1
      } else {
      loc.change <- dataIN$Current[locus,site] + dataIN$FutureFrqChange[locus,site]
      if( loc.change > 1 || loc.change < 0){
        become.fixed.allele = become.fixed.allele + 1
      } else {
        if (loc.change > max(dataIN$Current[locus,]) || loc.change < min(dataIN$Current[locus,])){
          high.allele = high.allele + 1
        }}}}
    output["change_to_fixation",site] <- become.fixed.allele
    output["highFrq",site] <- high.allele
    output["currently_fixed",site] <- fixed.allele
  }
  
  
  # plot
  maxX = max(as.numeric(output[1:2,]), na.rm=T)
  if(maxX > 0){ y.limit = ceiling(maxX/5)*5 } else { y.limit = 1}
  maxB = max(as.numeric(output["currently_fixed",]), na.rm=T)
  x.min.limit <- floor(min(enviro.data[varX,])/10)*10
  x.max.limit <- ceiling(max(enviro.data[varX,])/10)*10
  
  pdf(paste("MalAdaptCount_",CCmodel,".pdf",sep=""))
  par(mar=c(5,4,1,1),fig=c(0,1,0,0.75))
  plot(output[varX,], output["change_to_fixation",], pch = 3, col = output["PlotCol",], lwd = 1.5,
       xlab = varX, ylab = "No. of loci (Changing to Fixation or High Frequency)", ylim = c(0, y.limit), xlim = c(x.min.limit,x.max.limit), axes=F)
  points(output[varX,], output["highFrq",], pch = 16, col = output["PlotCol",])
  axis(1)
  axis(2, las = 2)
  text(x.min.limit, y.limit, labels = "b)", cex = 1.25)
  legend(x.min.limit, y.limit-1, c("Central NSW","Southern NSW", "Central Vic", "Western Vic"), pch = rep(16,4), col = c("red","goldenrod","green","blue"), cex = 0.7, pt.cex=0.9)
  legend(x.min.limit+20, y.limit-1, c("Fixation", "High Frq"), pch = c(3,16), cex = 0.7, pt.cex = 0.9)
  par(mar=c(2,4,1,1), fig=c(0,1,0.7,1), new=T)
  plot(output[varX,], output["currently_fixed",], pch = 17, col = output["PlotCol",], xlab = NA, ylab = "No. of loci (Currently Fixed)",
       xlim = c(x.min.limit,x.max.limit), axes = F)
  axis(1, labels = NA)
  axis(2, las = 2)
  text(x.min.limit, maxB, labels = "a)", cex = 1.25)
  dev.off()
  
  return( output )
}

mal.count.SumRain.Access1.0.2050 <- rain.mal.adapt.count(frq.change.SumRain.Access1.0.2050, "rain_sum", plot.col[,2], "SumRain_Access1.0_2050")
mal.count.SumRain.Access1.0.2070 <- rain.mal.adapt.count(frq.change.SumRain.Access1.0.2070, "rain_sum", plot.col[,2], "SumRain_Access1.0_2070")

mal.count.SumRain.HadGEM2.ES.2050 <- rain.mal.adapt.count(frq.change.SumRain.HadGEM2.ES.2050, "rain_sum", plot.col[,2], "SumRain_HadGEM2.ES_2050")
mal.count.SumRain.HadGEM2.ES.2070 <- rain.mal.adapt.count(frq.change.SumRain.HadGEM2.ES.2070, "rain_sum", plot.col[,2], "SumRain_HadGEM2.ES_2070")

mal.count.SumRain.GDFLCM3.2050 <- rain.mal.adapt.count(frq.change.SumRain.GDFLCM3.2050, "rain_sum", plot.col[,2], "SumRain_GDFLCM3_2050")
mal.count.SumRain.GDFLCM3.2070 <- rain.mal.adapt.count(frq.change.SumRain.GDFLCM3.2070, "rain_sum", plot.col[,2], "SumRain_GDFLCM3_2070")

mal.count.WinRain.Access1.0.2050 <- rain.mal.adapt.count(frq.change.WinRain.Access1.0.2050, "rain_win", plot.col[,2], "WinRain_Access1.0_2050")
mal.count.WinRain.Access1.0.2070 <- rain.mal.adapt.count(frq.change.WinRain.Access1.0.2070, "rain_win", plot.col[,2], "WinRain_Access1.0_2070")

mal.count.WinRain.HadGEM2.ES.2050 <- rain.mal.adapt.count(frq.change.WinRain.HadGEM2.ES.2050, "rain_win", plot.col[,2], "WinRain_HadGEM2.ES_2050")
mal.count.WinRain.HadGEM2.ES.2070 <- rain.mal.adapt.count(frq.change.WinRain.HadGEM2.ES.2070, "rain_win", plot.col[,2], "WinRain_HadGEM2.ES_2070")

mal.count.WinRain.GDFLCM3.2050 <- rain.mal.adapt.count(frq.change.WinRain.GDFLCM3.2050, "rain_win", plot.col[,2], "WinRain_GDFLCM3_2050")
mal.count.WinRain.GDFLCM3.2070 <- rain.mal.adapt.count(frq.change.WinRain.GDFLCM3.2070, "rain_win", plot.col[,2], "WinRain_GDFLCM3_2070")

  #- check
boxplot(t(frq.change.SumRain.GDFLCM3.2070$Current), horizontal = T, range = 0)
for( i in 1:nrow(frq.change.SumRain.GDFLCM3.2070$Current)){
points(frq.change.SumRain.GDFLCM3.2070$Current[i,"NCR"]+frq.change.SumRain.GDFLCM3.2070$FutureFrqChange[i,"NCR"], i, pch = 20, col = "red")
}

#-----------------------------------#
#    RANGE OF ALLELE FRQ CHANGES    #
#-----------------------------------#

#-- Boxplots of changes in frequency associated with future predicted climate

range.frq.change.by.var <- function(CCmodel){
  # create boxplots of allele frequency changes associated with individual variables for two time periods
  # dataIN = list created by "delta.MAT/Rain.increase" function containing two lists
  # 1 = "Current" = current allele freq
  # 2 = "FutureFrqChange" = change in wet-associated allele freq with climate change (includes + and - changes depending on direction of rainfall change)
  
  par(mar = c(3,4,0.5,0.5), fig = c(0,1,0,0.5))
  boxplot( as.vector( get(paste("frq.change.MAT.",CCmodel,".2050",sep=""))[["FutureFrqChange"]] ), at = 0.5, ylim = c(0,0.4), xlim = c(0,4), col = "blue") #, main = paste("Model:", CCmodel))
  boxplot( as.vector( get(paste("frq.change.MAT.",CCmodel,".2070",sep=''))[["FutureFrqChange"]] ), at = 1, add = T, axes = F, col = "red")
  boxplot( as.vector( abs(get(paste("frq.change.SumRain.",CCmodel,".2050",sep=''))[["FutureFrqChange"]]) ), at = 1.5, add = T, axes = F, col = "blue")
  boxplot( as.vector( abs(get(paste("frq.change.SumRain.",CCmodel,".2070",sep=''))[["FutureFrqChange"]]) ), at = 2, add=T, axes = F, col = "red")
  boxplot( as.vector( abs(get(paste("frq.change.WinRain.",CCmodel,".2050",sep=''))[["FutureFrqChange"]]) ), at = 2.5, add = T, axes = F, col = "blue")
  boxplot( as.vector( abs(get(paste("frq.change.WinRain.",CCmodel,".2070",sep=''))[["FutureFrqChange"]]) ), at = 3, add = T, axes = F, col = "red")
  mtext("Abs allele frq change", side = 2, line = 2)
  mtext(c("Temperature", "Summer Rain", "Winter Rain"), side = 1, at = c(0.75,1.75,2.75))
  legend(3,0.4, c("2050","2070"), col = c("blue","red"), pch = 15)
  
  par(mar = c(1,4,3,0.5), fig = c(0,1,0.5,1), new = T)
  SumRain.change <- ((SumRain.pred / enviro.data["rain_sum",])*100 ) - 100
  boxplot( as.vector( MAT.change[,paste("MAT_", CCmodel, "_2050", sep = "")] ), at = 0.5, ylim = c(0,6), xlim = c(0,4), col = "blue", main = paste("Model:", CCmodel))
  boxplot( as.vector( MAT.change[,paste("MAT_", CCmodel, "_2070", sep = "")]), at = 1, add = T, axes = F, col = "red")
  boxplot( as.vector( abs((SumRain.pred[,paste("SumRain_", CCmodel, "_2050", sep = "")] - enviro.data["rain_sum",])/10) ), at = 1.5, add = T, axes = F, col = "blue")
  boxplot( as.vector( abs((SumRain.pred[,paste("SumRain_", CCmodel, "_2070", sep = "")] - enviro.data["rain_sum",])/10) ), at = 2, add=T, axes = F, col = "red")
  boxplot( as.vector( abs((WinRain.pred[,paste("WinRain_", CCmodel, "_2050", sep = "")] - enviro.data["rain_win",])/10) ), at = 2.5, add = T, axes = F, col = "blue")
  boxplot( as.vector( abs((WinRain.pred[,paste("WinRain_", CCmodel, "_2070", sep = "")] - enviro.data["rain_win",])/10) ), at = 3, add = T, axes = F, col = "red")
  mtext(expression("Abs variable change ("*degree*"C, x10 mm)"), side = 2, line = 2)
  }

range.frq.change.by.var("Access1.0")
range.frq.change.by.var("HadGEM2.ES")
range.frq.change.by.var("GDFLCM3")

#-- Average change per degree temp or 10 mm rainfall based on linear models

  # create new output dataframes (empty)
lm.slope.temp.sumrain.winrain <- matrix(nrow = 62, ncol = 3)
colnames(lm.slope.temp.sumrain.winrain) <- c("bio01", "rain_sum", "rain_win")    

  # add slope data based on linear model
for( varX in c("bio01", "rain_sum", "rain_win")){   
  idx = 1
  # get loci with strong or very strong association for input environmental variable
    tmp <- outlier.frq[grep("strong", outlier.frq[,varX]),]
  # get slope of linear model for one unit of variable.
  for(snpX in rownames(tmp)){
    tmp.frq <- tmp[snpX,12:37]
    tmp.lm <- lm(unlist(tmp.frq) ~ enviro.data[varX,])
    lm.slope.temp.sumrain.winrain[idx,varX] <- summary(tmp.lm)$coefficients[2,1]
    idx = idx + 1
  }}

tmp.df <- abs(lm.slope.temp.sumrain.winrain)
tmp.df[,2] <- tmp.df[,2]*10
tmp.df[,3] <- tmp.df[,3]*10
boxplot(tmp.df)    
apply(tmp.df, 2, function(x) median(x,na.rm=T))
apply(tmp.df, 2, function(x) mean(x,na.rm=T))
apply(tmp.df, 2, function(x) sd(x,na.rm=T))

