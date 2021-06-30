############################################################################################
## Script for calculating standardized indices of litter abundance in the Baltic Sea
## using DATRAS Litter Exchange data.
## Author: Casper W. Berg, DTU Aqua.
## June 2021
#############################################################################################

## remotes::install_github("DTUAqua/DATRAS/DATRAS")
## remotes::install_github("casperwberg/surveyIndex/surveyIndex")

library(DATRAS)
library(maps); library(mapdata)
library(surveyIndex)
library(marmap)
library(plot.matrix)

source("readLitter.R")

litterTypes = c("Glass","Metal","Natural","Other","Plastic","Rubber")

litterTypesExt = c(litterTypes,"SUP","Fishing.related")

datafile = "../data/Litter Exchange Data_2021-06-08 12_39_40.csv"

d = readlitter(datafile,type="Weight")

d = subset(d,HaulDur>0 & HaulVal!="I")


xtabs(~Year,d)
xtabs(~Year+Quarter,d)
## Only few data points in 2011 -- drop

d = subset(d, Year != "2011")
d$Year = factor(d$Year)

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


png("../output/allmodels.png",width=1200,height=800)   
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
png("../output/allidx.png",width=1200,height=800)

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
d = subset(d,HaulDur>0 & HaulVal!="I")
d = subset(d, Year != "2011")
d$Year = factor(d$Year)
drd = df2dr(d)


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


png("../output/allmodels-numbers.png",width=1200,height=800)   
par(mfrow=c(2,4))
for(i in 1:length(models)){
    surveyIndex:::plot.SIlist(list(nmodels[[i]],nmodels2[[i]],nmodels3[[i]]),main=names(models)[i])
}
dev.off()
