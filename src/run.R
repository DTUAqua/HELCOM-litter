############################################################################################
## Script for calculating standardized indices of litter abundance in the Baltic Sea
## using DATRAS Litter Exchange data.
## Author: Casper W. Berg, DTU Aqua.
## June 2021
#############################################################################################

## remotes::install_github("DTUAqua/DATRAS/DATRAS")
## NOTE: "posprob" branch of surveyIndex package needed
remotes::install_github("casperwberg/surveyIndex/surveyIndex",ref="posprob")

library(DATRAS)
library(maps); library(mapdata)
library(surveyIndex)
library(marmap)
library(plot.matrix)
library(xtable)

source("readLitter.R")

litterTypes = c("Glass","Metal","Natural","Other","Plastic","Rubber")

litterTypesExt = c(litterTypes,"SUP","Fishing.related")

datafile = "../data/Litter Exchange Data_2022-08-26 12_41_43.zip"

d = readlitter(datafile,type="Weight")

d = subset(d,HaulDur>0 & HaulVal!="I")
d = subset(d,!is.na(d$lon) & !is.na(d$lat))


xtabs(~Year,d)
xtabs(~Year+Quarter,d)
d$Year = factor(d$Year)

## Make tables with data overview
makeTable<-function(x,fil,cap) cat(print(xtable( x, caption=cap, digits=0)),file=fil)

makeTable( xtabs(~Year+Quarter,d), fil="../report/yqtab.tex",cap="Weight: Number of hauls by year and quarter")
makeTable( xtabs(~Year+Country,d), fil="../report/yctab.tex",cap="Weight: Number of hauls by year and country")
makeTable( xtabs(~Gear+Country,d), fil="../report/gctab.tex",cap="Weight: Number of hauls by gear and country")
makeTable( xtabs(~Country+Quarter,d), fil="../report/cqtab.tex",cap="Weight: Number of hauls by country and quarter")



## data frame to DATRASraw object
df2dr<-function(x){
    x$haul.id = 1:nrow(x)
    dd = list()
    dd[[1]] = data.frame()
    dd[[2]] = x
    dd[[3]] = data.frame()
    class(dd)<-"DATRASraw"
    dd
}


drd = df2dr(d)

## Get prediction grid
bgrid <- getBathyGrid(drd,minDepth=1,maxDepth=1000,maxDist=Inf,resolution=3,shapefile="../shapefiles/ICES/ICES_areas.shp",select="ICES_SUB")
bgrid = subset(bgrid, ICES_SUB %in% as.character(21:32))

tmp = addSpatialData(df2dr(bgrid),shape="../shapefiles/EEZshape/EMODnet_HA_OtherManagementAreas_EEZ_v11_20210506.shp")
bgrid = tmp[[2]]
bgrid$Territory = factor(bgrid$Territory) ## drop empty factor levels

## Plot bathy grid
my.palette<-colorRampPalette(c("darkblue","mediumblue","lightblue1"))
my.palette.vec=my.palette(100);
png("../output/bathygrid.png",width=1200,height=800)
plot(bgrid$lon,bgrid$lat,col=rev(my.palette.vec)[cut(bgrid$Depth,100)],pch=15,cex=1.3)
maps::map("worldHires", fill = TRUE, plot = TRUE, 
                add = TRUE, col = grey(0.5))
points(d$lon,d$lat,col=2,pch=".",cex=3)
dev.off()

## Plot EEZ map
png("../output/EEZmap.png",width=1200,height=800)
plot(bgrid$lon,bgrid$lat,col=bgrid$Territory,pch=15,cex=1.3)
maps::map("worldHires", fill = TRUE, plot = TRUE, 
                add = TRUE, col = grey(0.9))
legend("bottomright",legend=levels(bgrid$Territory),col=1:nlevels(bgrid$Territory),pch=15,bg="white",cex=1.3)
dev.off()


### Standardized effort: 1 km^2 (1e6 m^2).
## Index is sum of grid points, so divide by number of grid points such that the index is
## the mean litter abundance. 
StdEffort = 1e6 / nrow(bgrid) 


## Define model formulas

NYEARS = length(unique(d$Year))


fm = paste0("s(ctime,k=",NYEARS,",bs='ds',m=c(1,0)) + s(lon,lat, bs='ds',m=c(1,0.5),k=128) + offset(log(EFFORT))")

## Year effect instead of spline
fm2 = paste0("Year + s(lon,lat, bs='ds',m=c(1,0.5),k=128) + offset(log(EFFORT))")

## Linear effect of time
fm3 = paste0("ctime + s(lon,lat, bs='ds',m=c(1,0.5),k=128) + offset(log(EFFORT))")

formulas = list(fm,fm2,fm3)


ages = 1

#################
## Fit models
##################

models = list()
models2 = list()
models3 = list()

for(lt in litterTypesExt){

    cat("Doing ",lt,"...\n")
    
    drd$Nage = matrix(d[,lt], nrow=nrow(d),ncol=1)
    colnames(drd$Nage)<-1
    
    system.time( models[[ lt ]]  <- getSurveyIdx(drd,ages,predD=bgrid,cutOff=0,fam="Tweedie",mc.cores=1,modelP=fm,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE,maxit=20) ) )

    system.time( models2[[ lt ]]  <- getSurveyIdx(drd,ages,predD=bgrid,cutOff=0,fam="Tweedie",mc.cores=1,modelP=fm2,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE,maxit=20) ) )

    system.time( models3[[ lt ]]  <- getSurveyIdx(drd,ages,predD=bgrid,cutOff=0,fam="Tweedie",mc.cores=1,modelP=fm3,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE,maxit=20) ) )

    
}


png("../output/allmodels.png",width=1200,height=800,pointsize=24)   
par(mfrow=c(2,4))
for(i in 1:length(models)){
    surveyIndex:::plot.SIlist(list(models[[i]],models2[[i]],models3[[i]]),main=names(models)[i])
}
dev.off()

### Plotting
maxBubble = 8

for(lt in litterTypesExt){
    png(paste0("../output/",lt,"%03d.png"),width=1200,height=800)
     
    ## Bubble plots - one per year
    myscale=maxBubble/max(sqrt(d[,lt]))
    
    xlims = range(d$lon)
    ylims = range(d$lat)
    op <- par(mfrow=n2mfrow(nlevels(d$Year)+1),mar=c(1,1,2,1),oma=c(0,0,2,0))
    for(yy in sort(unique(d$Year))){
        tmp = subset(d,Year==yy)
        mybubblePlot(tmp,response=lt,scale=myscale,rim=TRUE,xlim=xlims,ylim=ylims,axes=FALSE)
        box()
        title(yy)
    }
    title(lt,outer=TRUE,line=0,cex.main=2)

    legkg = c( 0, round((max(sqrt(d[,lt]))/maxBubble)^2,2),
              round((max(sqrt(d[,lt])))^2,2))   
    legend("bottomright",pch=c(3,16,16),col=c(2,1,1),pt.cex=c(1,1,maxBubble),legend=paste(legkg,"kg"),bg="white",x.intersp=3,y.intersp=3)

    ## All years together 
    mybubblePlot(d,response=lt,scale=myscale,rim=TRUE,xlim=xlims,ylim=ylims,axes=FALSE)
    title("All years")

    par(op)

    ## Fitted map
    surveyIdxPlots(models[[lt]],drd,myids=NULL,predD=bgrid,select="map",colors=rev(heat.colors(7)),legend=TRUE,legend.signif=2,map.cex=1.3,par=list(mfrow=c(1,1)),main=lt)


    ## Fitted total abundance
    surveyIndex:::plot.SIlist( list(models[[lt]]) ,main=paste(lt,"(kg / km^2)"))

    dev.off()
}


## Plot them all
png("../output/allidx.png",width=1200,height=800,pointsize=24)

allidxs = lapply(models[-3],function(x)x$idx)
maxY = max(sapply( allidxs,max))
par(mfrow=c(1,1),mar=c(4,4,3,3))
for(i in 1:length(allidxs)){
    ys = rownames(allidxs[[i]])

    if(i==1) plot(ys,allidxs[[i]],ylim=c(0,maxY),type="b",lwd=2,ylab="Index ( kg / km^2)") else lines(ys,allidxs[[i]],type="b",col=i,lwd=2)
    
}

legend("topright",col=1:length(allidxs),legend=names(allidxs),pch=1,lty=1,lwd=2)

dev.off()


## Export model summaries

sink("../output/summaries.txt")
lapply(models,function(x) { summary(x$pModels[[1]])  } )
cat("=====================\n")
lapply(models2,function(x) { summary(x$pModels[[1]])  } )
cat("=====================\n")
lapply(models3,function(x) { summary(x$pModels[[1]])  } )
sink()


################################
## Calculate indices by EEZ
################################
EEZmodels=list()

for(lt in litterTypesExt){
    cat("Doing ",lt,"...\n")
        
    EEZmodels[[lt]] = lapply(levels(bgrid$Territory),function(x){
        cat(x,"\n")
        pd = subset(bgrid,Territory==x)
        redoSurveyIndex(drd,models[[lt]],predD=pd,myids=NULL,predfix=list(EFFORT=1e6/nrow(pd)))
        })
}

## Arrange indices and uncertainties in matrix form
EEZmat = matrix(NA,nlevels(bgrid$Territory),length(litterTypesExt))
rownames(EEZmat)<-levels(bgrid$Territory)
colnames(EEZmat)<-litterTypesExt
EEZmatCV <- EEZmat

for(lt in litterTypesExt){
    for(ter in 1:nrow(EEZmat)){
        i = which(litterTypesExt == lt)
        EEZmat[ter,i] = tail(EEZmodels[[lt]][[ter]]$idx[,1],1)
        logSD = (tail(log(EEZmodels[[lt]][[ter]]$up[,1]),1) - tail(log(EEZmodels[[lt]][[ter]]$lo[,1]),1))/4 
        EEZmatCV[ter,i] = logSD
    }
}

col <- colorRampPalette(rev(c("red", "white", "blue")))
png("../output/EEZplot%03d.png",width=1200,height=800)

par(mfrow=c(1,1),mar=c(5,5,4,4))
plot(EEZmat,col=col,fmt.cell="%.2f",fmt.key="%.2f",las=1,xlab="",ylab="",main=paste("Litter density (kg / km^2) by EEZ",tail(levels(d$Year),1)))

col <- colorRampPalette(c("white","yellow","red"))
plot(EEZmatCV,col=col,fmt.cell="%.2f",fmt.key="%.2f",las=1,xlab="",ylab="",main=paste("Litter density uncertainty (CV) by EEZ",tail(levels(d$Year),1)))
dev.off()


#####################################
### Numbers insted of mass
#####################################

d = readlitter(datafile,type="Numbers")
## Note, different number of hauls when considering numbers!
d = subset(d,HaulDur>0 & HaulVal!="I")
d = subset(d,!is.na(d$lon) & !is.na(d$lat))
d$Year = factor(d$Year)
drd = df2dr(d)


makeTable( xtabs(~Year+Quarter,d), fil="../report/yqtab-n.tex",cap="Numbers: Number of hauls by year and quarter")
makeTable( xtabs(~Year+Country,d), fil="../report/yctab-n.tex",cap="Numbers: Number of hauls by year and country")
makeTable( xtabs(~Gear+Country,d), fil="../report/gctab-n.tex",cap="Numbers: Number of hauls by gear and country")
makeTable( xtabs(~Country+Quarter,d), fil="../report/cqtab-n.tex",cap="Numbers: Number of hauls by country and quarter")



StdEffort = 1e6 / nrow(bgrid)  ## Numbers pr km^2

nmodels = list()
nmodels2 = list()
nmodels3 = list()

for(lt in litterTypesExt){

    cat("Doing ",lt,"...\n")
    
    drd$Nage = matrix(d[,lt], nrow=nrow(d),ncol=1)
    colnames(drd$Nage)<-1
    
    system.time( nmodels[[ lt ]]  <- getSurveyIdx(drd,ages,predD=bgrid,cutOff=0,fam="negbin",mc.cores=1,modelP=fm,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE) ) )

    system.time( nmodels2[[ lt ]]  <- getSurveyIdx(drd,ages,predD=bgrid,cutOff=0,fam="negbin",mc.cores=1,modelP=fm2,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE,maxit=20) ) )

    system.time( nmodels3[[ lt ]]  <- getSurveyIdx(drd,ages,predD=bgrid,cutOff=0,fam="negbin",mc.cores=1,modelP=fm3,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE,maxit=20) ) )

    
}


png("../output/allmodels-numbers.png",width=1200,height=800,pointsize=24)   
par(mfrow=c(2,4))
for(i in 1:length(models)){
    surveyIndex:::plot.SIlist(list(nmodels[[i]],nmodels2[[i]],nmodels3[[i]]),main=names(models)[i])
}
dev.off()

## Maps
for(lt in litterTypesExt){
    png(paste0("../output/",lt,"-numbers.png"),width=1200,height=800)

        surveyIdxPlots(nmodels[[lt]],drd,myids=NULL,predD=bgrid,select="map",colors=rev(heat.colors(7)),legend=TRUE,legend.signif=2,map.cex=1.3,par=list(mfrow=c(1,1)),main=lt)

    dev.off()
}



## Output summaries
sink("../output/summaries-numbers.txt")
cat("============ Models on numbers ===============\n")
lapply(nmodels,function(x) { summary(x$pModels[[1]])  } )
cat("=====================\n")
lapply(nmodels2,function(x) { summary(x$pModels[[1]])  } )
cat("=====================\n")
lapply(nmodels3,function(x) { summary(x$pModels[[1]])  } )
sink()

##########################################
## Calculate indices by EEZ (numbers)
##########################################
EEZmodelsn=list()

for(lt in litterTypesExt){
    cat("Doing ",lt,"...\n")
        
    EEZmodelsn[[lt]] = lapply(levels(bgrid$Territory),function(x){
        cat(x,"\n")
        pd = subset(bgrid,Territory==x)
        redoSurveyIndex(drd,nmodels[[lt]],predD=pd,myids=NULL,predfix=list(EFFORT=1e6/nrow(pd)))
        })
}

## Arrange indices and uncertainties in matrix form
EEZmat = matrix(NA,nlevels(bgrid$Territory),length(litterTypesExt))
rownames(EEZmat)<-levels(bgrid$Territory)
colnames(EEZmat)<-litterTypesExt
EEZmatCV <- EEZmat

for(lt in litterTypesExt){
    for(ter in 1:nrow(EEZmat)){
        i = which(litterTypesExt == lt)
        EEZmat[ter,i] = tail(EEZmodelsn[[lt]][[ter]]$idx[,1],1)
        logSD = (tail(log(EEZmodelsn[[lt]][[ter]]$up[,1]),1) - tail(log(EEZmodelsn[[lt]][[ter]]$lo[,1]),1))/4 
        EEZmatCV[ter,i] = logSD
    }
}

col <- colorRampPalette(rev(c("red", "white", "blue")))
png("../output/EEZplotn%03d.png",width=1200,height=800)

par(mfrow=c(1,1),mar=c(5,5,4,4))
plot(EEZmat,col=col,fmt.cell="%.2f",fmt.key="%.2f",las=1,xlab="",ylab="",main=paste("Litter density (numbers / km^2) by EEZ",tail(levels(d$Year),1)))

col <- colorRampPalette(c("white","yellow","red"))
plot(EEZmatCV,col=col,fmt.cell="%.2f",fmt.key="%.2f",las=1,xlab="",ylab="",main=paste("Litter density uncertainty (CV) by EEZ (numbers)",tail(levels(d$Year),1)))
dev.off()



#############################
## Probability of encounter
#############################

pmodels = list()
pmodels2 = list()
pmodels3 = list()

## We want probability pr haul
StdEffort = TVS2TVL()[2]

for(lt in litterTypesExt){

    cat("Doing ",lt,"...\n")
    
    drd$Nage = matrix(d[,lt], nrow=nrow(d),ncol=1)
    colnames(drd$Nage)<-1
    
    system.time( pmodels[[ lt ]]  <- getSurveyIdx(drd,ages,predD=bgrid,cutOff=0,fam="negbin",mc.cores=1,modelP=fm,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE), doProbs=TRUE ) )

    system.time( pmodels2[[ lt ]]  <- getSurveyIdx(drd,ages,predD=bgrid,cutOff=0,fam="negbin",mc.cores=1,modelP=fm2,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE,maxit=20), doProbs=TRUE ) )

    system.time( pmodels3[[ lt ]]  <- getSurveyIdx(drd,ages,predD=bgrid,cutOff=0,fam="negbin",mc.cores=1,modelP=fm3,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE,maxit=20), doProbs=TRUE ) )

    
}


png("../output/allmodels-posprob.png",width=1200,height=800,pointsize=24)   
par(mfrow=c(2,4))
for(i in 1:length(pmodels)){
    surveyIndex:::plot.SIlist(list(pmodels[[i]],pmodels2[[i]],pmodels3[[i]]),main=names(pmodels)[i],posProb=TRUE)
}
dev.off()

## Fitted maps

lastyear = tail(levels(d$Year),1)

probcols = rev(hcl.colors(7,palette="Reds 3")) 

for(lt in litterTypesExt){
    png(paste0("../output/",lt,"-posprob.png"),width=1200,height=800)

    surveyIdxPlots(pmodels[[lt]],drd,myids=NULL,predD=bgrid,select="absolutemap",colors=probcols,legend=TRUE,legend.signif=2,map.cex=1.3,par=list(mfrow=c(1,1),oma=c(0,0,1,0)),year=lastyear,posProb=TRUE,scaleMap=FALSE)
    title(lt,outer=TRUE)
    dev.off()
    
}

##########################
## CSV output
##########################

allout = list()
## drop first year (2011) weight estimate in CSV output, because it is only available for weight
for(lt in litterTypesExt){
    allout[[lt]] = data.frame(Type=lt, Year = rownames(nmodels[[lt]]$idx), DensityMass=models[[lt]]$idx[-1,1], DensityMassLow = models[[lt]]$lo[-1,1], DensityMassHigh=models[[lt]]$up[-1,1],
                              DensityNumbers=nmodels[[lt]]$idx[,1], DensityNumbersLow = nmodels[[lt]]$lo[,1], DensityNumbersHigh=nmodels[[lt]]$up[,1],
                              DensityProb=pmodels[[lt]]$idx0[,1],    DensityProbLow = pmodels[[lt]]$lo0[,1], DensityProbHigh=pmodels[[lt]]$up0[,1])
                              
}



allout.df = do.call(rbind,allout)

write.csv(allout.df,file="../output/litterEstimates.csv",row.names=FALSE)


############################################################
## Trend models using only data from 2015 and onwards
############################################################

drd2 = subset(drd,!Year %in% as.character(2011:2014))
d2 = subset(d, !Year %in% as.character(2011:2014))

StdEffort = 1e6 / nrow(bgrid)  ## Numbers/Mass pr km^2

trend15models = list()
trend15modelsn = list()

for(lt in litterTypesExt){

    cat("Doing ",lt,"...\n")
    
    drd2$Nage = matrix(d2[,lt], nrow=nrow(d2),ncol=1)
    colnames(drd2$Nage)<-1

    system.time( trend15modelsn[[ lt ]]  <- getSurveyIdx(drd2,ages,predD=bgrid,cutOff=0,fam="negbin",mc.cores=1,modelP=fm3,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE,maxit=20) ) )

    system.time( trend15models[[ lt ]]  <- getSurveyIdx(drd2,ages,predD=bgrid,cutOff=0,fam="Tweedie",mc.cores=1,modelP=fm3,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE,maxit=20) ) )
   
}

## Output summaries
sink("../output/trend15summaries.txt")
cat("============ Trend models (2015 onwards)  ===============\n")
cat("============ Mass  ===============\n")
lapply(trend15models,function(x) { summary(x$pModels[[1]])  } )
cat("============ Numbers  ===============\n")
lapply(trend15modelsn,function(x) { summary(x$pModels[[1]])  } )
sink()
