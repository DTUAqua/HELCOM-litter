###################################################################
## Helper functions for handling DATRAS litter exchange data
## Author: Casper W. Berg, DTU Aqua.
###################################################################

TVS2TVL <- function(){
    ## Swept area (m^2) for 30 min hauls
    TVS_TVL = c("TVLobs"=82798,"TVLobspred"=87163,"TVSobs"=68764,"TVSobspred"=68184)
    TVS_TVL
}

readlitter <- function (file = "IBTS.csv", na.strings = c("-9", "-9.0", "-9.00", "-9.0000"), type="Weight"){
    ## Locating lines with headers
    lines <- readLines(file)
    i <- grep("RecordType", lines)
    headers <- lines[i]
    skip <- i - 1
    nrow <- c(diff(i), 0) - 1
    cat("Reading data files\n")
    d <- lapply(1:length(i), function(i) {
        emptyDF <- function(nm) {
            structure(rep(list(logical(0)), length(nm)), .Names = nm, 
                class = "data.frame")
        }
        if (nrow[i] == 0) {
            return(emptyDF(strsplit(headers[i], ",")[[1]]))
        }
        
        ans <- read.csv(file, nrow = nrow[i], 
                                          skip = skip[i], na.strings = na.strings, stringsAsFactors = TRUE)
        
        ans$StNo <- as.character(ans$StNo)
        ans <- renameDATRAS(ans)
        ans
    })
    names(d) <- sapply(d, function(x) {
        as.character(x$RecordType[1])
    })

    if (sum(sapply(d, nrow)) + length(i) != length(lines)) 
        stop("csv file appears to be corrupt.")
    if (is.null(d[[1]]$StatRec)) 
        d[[1]]$StatRec <- d[[1]]$AreaCode
    haul.id <- quote(factor(paste(Year, Quarter, Country, Ship, 
                                  Gear, StNo, HaulNo, sep = ":")))

    d[[1]]$haul.id = eval(haul.id,d[[1]])
    d[[2]]$haul.id = eval(haul.id,d[[2]])
    d[[2]]$haul.id = factor(d[[2]]$haul.id, levels=levels(d[[1]]$haul.id))
    
    LTnotHH = setdiff(levels(d[[2]]$haul.id),levels(d[[1]]$haul.id))
    HHnotLT = setdiff(levels(d[[1]]$haul.id),levels(d[[2]]$haul.id))
    if( length(LTnotHH) > 0 ){
        print(paste("Dropping",length(LTnotHH),"LT records with no HH info:"))
        print(LTnotHH)
        d[[2]] = subset(d[[2]],!haul.id %in% LTnotHH)
    }

    if( length(HHnotLT) > 0 ){
        print(paste("Dropping",length(HHnotLT),"HH records with no LT info:"))
        print(LTnotHH)
        d[[1]] = subset(d[[1]],!haul.id %in% HHnotLT)
    }

    
     ## Don't use RECO-LT (only total weight)
    badhauls = unique(d$LT$haul.id[ d$LT$LTREF=="RECO-LT" ] )
    cat("Dropping ", length(badhauls), " hauls of type RECO-LT\n")
    d[[2]] = subset(d[[2]], !haul.id %in% badhauls)
    d[[1]] = subset(d[[1]], !haul.id %in% badhauls)

    if(type=="Numbers"){
        badhauls = unique(d$LT$haul.id[ is.na(d$LT$LT_Items) ] )
        cat("Dropping ", length(badhauls), " hauls with only weight info\n")
        d[[2]] = subset(d[[2]], !haul.id %in% badhauls)
        d[[1]] = subset(d[[1]], !haul.id %in% badhauls)
    }

    if(type=="Weight"){
        badhauls = unique(d$LT$haul.id[ is.na(d$LT$UnitWgt) ] )
        cat("Dropping ", length(badhauls), " hauls with no weight unit info\n")
        d[[2]] = subset(d[[2]], !haul.id %in% badhauls)
        d[[1]] = subset(d[[1]], !haul.id %in% badhauls)

        d$LT$LT_Weight[ d$LT$UnitWgt=="g/haul" ] = d$LT$LT_Weight[ d$LT$UnitWgt=="g/haul" ] / 1000

    }
    
    if(nlevels(d[[1]]$haul.id) != nlevels(d[[2]]$haul.id)){
        print("Missing haul info:")
        cat("HH: " ,nlevels(d$HH$haul.id)," ,LT: ",nlevels(d$LT$haul.id),"\n","Setdiff:\n")
        print( setdiff(levels(d[[1]]$haul.id),levels(d[[2]]$haul.id)))
        stop()
    }
    
    conv = read.csv("../data/litter_conversion.csv")
    lookup <- function(x,table){
        x <- factor(x)
        levels(x) <- table[levels(x)]
        if(is.numeric(table)) as.numeric(as.character(x)) else x
    }

    Code2Type <- structure(as.character(conv$Type), names=conv$C.TS)
    Code2TypeRev <- structure(as.character(conv$Type), names=conv$C.TS.REV)

    d[[2]] = transform(d[[2]],Type= factor(
                      ifelse(as.character(LTREF)=="C-TS",
                             as.character( lookup(PARAM,Code2Type) ),
                             as.character( lookup(PARAM,Code2TypeRev) )
                      )))
                      

    if(type=="Weight"){
        agg = xtabs(LT_Weight ~ haul.id + Type,d[[2]])[ d[[1]]$haul.id, ]
        
    } else if(type=="Numbers"){
        agg = xtabs(LT_Items ~ haul.id + Type,d[[2]])[ d[[1]]$haul.id, ]
    }
    dn = dimnames(agg)[2][[1]]
    for(i in 1:length(dn))
        d[[1]][,dn[i]] = agg[,i]

    ## Extra categories:
    ## 1. Singe use plastic (SUP)
    Code2Type <- structure(as.character(conv$SUP), names=conv$C.TS)
    Code2TypeRev <- structure(as.character(conv$SUP), names=conv$C.TS.REV)
    d[[2]] = transform(d[[2]],SUP= factor(
                                  ifelse(as.character(LTREF)=="C-TS",
                                         as.character( lookup(PARAM,Code2Type) ),
                                         as.character( lookup(PARAM,Code2TypeRev) )
                                         )))
    if(type=="Weight"){
        agg = xtabs(LT_Weight ~ haul.id + SUP,d[[2]])[ d[[1]]$haul.id, ]
        
    } else if(type=="Numbers"){
        agg = xtabs(LT_Items ~ haul.id + SUP,d[[2]])[ d[[1]]$haul.id, ]
    }
    d[[1]]$SUP = agg[ ,"Yes"]  
         
    ## 2. Fishing related
    Code2Type <- structure(as.character(conv$Fishing.related), names=conv$C.TS)
    Code2TypeRev <- structure(as.character(conv$Fishing.related), names=conv$C.TS.REV)
    d[[2]] = transform(d[[2]],Fishing.related = factor(
                                  ifelse(as.character(LTREF)=="C-TS",
                                         as.character( lookup(PARAM,Code2Type) ),
                                         as.character( lookup(PARAM,Code2TypeRev) )
                                         )))
    if(type=="Weight"){
        agg = xtabs(LT_Weight ~ haul.id + Fishing.related,d[[2]])[ d[[1]]$haul.id, ]
        
    } else if(type=="Numbers"){
        agg = xtabs(LT_Items ~ haul.id + Fishing.related,d[[2]])[ d[[1]]$haul.id, ]
    }
    d[[1]]$Fishing.related = agg[,"Yes"]

    ## Extra variables (for surveyIndex package)
    d[[1]]$lon = (d[[1]]$ShootLong + d[[1]]$HaulLong)/2
    d[[1]]$lat = (d[[1]]$ShootLat + d[[1]]$HaulLat)/2
    
    d[[1]]$ctime = d[[1]]$Year + (d[[1]]$Month-1)/12 + (d[[1]]$Day-1)/365

    d[[1]]$timeOfYear = (d[[1]]$Month-1)/12 + (d[[1]]$Day-1)/365
    
    d[[1]]$Year = factor(d[[1]]$Year)

    d[[1]]$TimeShotHour = 12

    TVS_TVL = TVS2TVL()
    
    d[[1]]$EFFORT = TVS_TVL[2]/30 * d[[1]]$HaulDur
    d[[1]]$EFFORT[ d[[1]]$Gear == "TVS" ] = TVS_TVL[4]/30 * d[[1]]$HaulDur[ d[[1]]$Gear == "TVS" ]
    
    d[[1]]
}



mybubblePlot<-function (d, response = "HaulWgt", scale = NULL, col.zero = "red", 
    pch.zero = "+", rim = FALSE, ...) 
{
    d$resp.var <- d[,response]
    if (is.null(scale)) 
        scale = 4/max(sqrt(d$resp.var), na.rm = TRUE)
    plot(d$lon, d$lat, type = "n", xlab = "Longitude", ylab = "Latitude",...)
    map("worldHires", fill = TRUE, plot = TRUE, add = TRUE, col = grey(0.5))
    points(d$lon, d$lat, pch = 16, cex = scale * sqrt(d$resp.var))
    if (rim) 
        points(d$lon, d$lat, cex = scale * sqrt(d$resp.var), 
            pch = 1, lwd = 0.5, col = "blue")
    zero = subset(d, resp.var == 0)
    points(zero$lon, zero$lat, pch = pch.zero, col = col.zero)
}
