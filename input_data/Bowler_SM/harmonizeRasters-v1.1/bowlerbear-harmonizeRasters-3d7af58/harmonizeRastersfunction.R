# Function to read and harmonize global raster layers#
######################################################

#Main arguments is:

#x is a dataframe with 3 columns or list with 3 elements
#myfolder - folder containing original raster layer
#myproj - projection of original raster layer
#name - simple name to call raster layer in output (not necessary....)

#Additional optional arguments:
#newres = desired resolution in units of degrees or metres depending on type of projection
#newextent = desired extent in degrees
#timeperiod = desired time period if raster is a stack of many times points
#newproj = desired CRS geographic or projection
  #these are examples:  
  #newproj<-"+proj=eck4 +datum=WGS84"
  #newproj<-"+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
#plot = whether to plot the raster layers as part of the output
#summary = whether to summarise multiple input rasters as the "Mean" or the "Trend"
#aggregate = whether to aggregate the data by taking the mean of values ("Mean",default) or summing the values ("Sum")

#default is to harmonize and/all raster layers/bands to a global extent, 1 degree resolution, 
#in geographic projection

#Example use
#whan myrasters dataframe have just one row (i.e. one raster)
outout<-harmonizeRasters(myrasters)
output<-harmonizeRasters(myrasters,newres=500000,newproj="+proj=eck4 +datum=WGS84",plot=TRUE)

#apply to mutliple rasters and combine them as a stack
output<-lapply(1:nrow(myrasters),function(i)harmonizeRasters(myrasters[i,]))
output<-stack(output)

##############################################################################################

harmonizeRasters<-function(x, newres=1,newextent=extent(-180, 180, -90, 90),timeperiod=NULL, 
                           newproj="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs",
                           plot=FALSE,summary=NULL,aggregate="Mean"){
  
  #libraries we need
  require(raster)
  require(rgdal)
  require(gdalUtils)
  require(ncdf4)
  
  
  #Extract bits of x for later use in function
  name<-as.character(x["name"])
  myproj<-as.character(x["myproj"])  
  myfolder<-as.character(x["myfolder"]) 
  
  #Create reference raster layer that raster will be reprojected onto
  
  #set up a normal 1 degree degree grid
  ref<-newextent
  ref<-raster(ref)
  res(ref)<-1
  values(ref)<-1#dummy value
  projection(ref)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
  
  #project into whatever is the desired end projection resolution
  refProj<-projectRaster(ref, crs=newproj, res=newres)
  
  #Read in original raster layer
  
  #get file types in directory
  files<-list.files(myfolder)
  filetypes<-as.character(sapply(files,function(x)substr(x,nchar(x)-2,nchar(x))))
  
  #if file type is asci or tif
  if(any(filetypes%in%c("tif","asc","wcs",".nc"))){
    myfile<-files[which(filetypes%in%c("tif","asc","wcs",".nc"))]  
    r<-stack(paste(myfolder,myfile,sep="/"))
  }else {
    temp <- new("GDALReadOnlyDataset",myfolder)
    temp2<-asSGDF_GROD(temp)
    r <- stack(temp2)
  }
  
  #check whether we need to specific the crs projection of the raster (i.e., if R doesn't
  #automatically read it from the files)  
  if(is.na(crs(r))){
    projection(r)=myproj
  }
  
  #If the original raster has multiple time points, subset to time period of interest
  if(!is.null(timeperiod)){
    require(lubridate)
    dates<-as.character(sapply(names(r),function(x)substr(x,2,nchar(x))))
    dates<-gsub("\\.","/",dates)
    dates<-as.Date(dates)
    myyears=year(dates)
    r <- subset(r, names(r)[myyears%in%timeperiod]) 
  }
  
  #if in lat,lon reproject into the new proj that is in m
  if(grepl("longlat",myproj)){
    r <- projectRaster(r, crs=newproj)
  }
  
  #If the raster is on a finer resolution than that we are sampling it, we first need to aggregate it:
  #are we to aggregating by sum or mean?
  origRes<-res(r)
  agFactor<-floor(newres/origRes[1])
  if(agFactor>1){
    if(as.logical(aggregate=="Sum")){
      r<-aggregate(r,fact=agFactor,fun=sum,na.rm=T)
    }
    else if(as.logical(aggregate=="Mean")){
      r<-aggregate(r,fact=agFactor,fun=mean,na.rm=T)
    }else{
      r<-r
    }
  }
  
  #If the raster has multiple layers, do we just want a summary, i.e., the mean or trend over time
  #do we want to calculate a summary variable e.g. a mean or a trend
  if(!is.null(summary)){
    if(summary=="trend"){#trend, i.e. regression coefficient of a year effect
      myyears=myyears[myyears%in%timeperiod]
      
      require(spatial.tools)
      rStack<-brickstack_to_raster_list(r)
      
      #get average raster value per year
      rMean<-stack(llply(unique(myyears),function(x){
        rTemp <- stack(rStack[myyears%in%x])
        rTemp <- calc(rTemp,fun=function(x)mean(x))
        return(rTemp)
      }))
      
      uniqueYears<-unique(myyears)
      
      #get trend
      fun <- function(x) {
        m <- NA
        try( m <- lm(x ~uniqueYears)$coefficients[2] ,silent=T)
        m
      }
      
      r<- calc(rMean, fun)
      
    }else if(as.logical(summary=="T-stat")){ #t-statistic of the test
      myyears=myyears[myyears%in%timeperiod]

      require(spatial.tools)
      rStack<-brickstack_to_raster_list(r)
      
      #get average per year
      rMean<-stack(llply(unique(myyears),function(x){
        rTemp <- stack(rStack[myyears%in%x])
        rTemp <- calc(rTemp,fun=function(x)mean(x))
        return(rTemp)
      }))
      
      uniqueYears<-unique(myyears)
      
      getT<-function(x){
        if(all(!is.na(x))){
          summary(lm(x~uniqueYears))$coefficients[2,3]
        }else{
          NA
        }}
      
      r<- calc(rMean, getT)
      
    }else if(summary=="Mean"){
      r <- calc(r, mean)  
    }}
  
  
  #Clip raster to extent of reference grid
  extentGrid <- projectExtent(ref, crs(r))
  r<-crop(r,extentGrid)
  
  #Project raster onto the reference grid of raster
  rProj <- projectRaster(r, refProj) 
  rProj<-mask(rProj,refProj)

  #Convert raster into data frames
  rDF<-as.data.frame(rProj,xy=T)
  names(rDF)[3]<-"Value"
  rDF$Type<-name
  return(rDF)
}
