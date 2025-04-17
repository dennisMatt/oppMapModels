library(terra)
library(sf)
library(gdistance)
library(dplyr)


#build functions for exponential kernels (from https://rdrr.io/github/NINAnor/oneimpact/)

#####################Exp filter functions


filter_create <- function(r = NULL,
                          radius = NULL,
                          type = c("exp_decay", "bartlett", "circle", "threshold_decay",
                                   "gaussian_decay", "Gauss", "rectangle")[1],
                          zoi_limit = 0.05,
                          half_life = NULL,
                          zoi_hl_ratio = NULL,
                          sigma = NULL,
                          min_intensity = 0.01,
                          max_dist = 5000,
                          normalize = FALSE,
                          divisor = 1,
                          round_vals = NULL,
                          save_txt = FALSE,
                          save_format = c("GRASS_rmfilter", "raw")[1],
                          save_folder = NULL,
                          save_file = NULL,
                          parallel = TRUE,
                          ...) {
  
  # check the input data class of r
  if(class(r) %in% c("RasterLayer", "RasterBrick", "RasterStack", "SpatRaster")) {
    res <- terra::res(r)[1]
  } else {
    if(is.numeric(r) & r > 0) {
      res <- r
    } else
      stop("'r' must be either an input raster map or a numeric value corresponding to the resolution of a raster.")
    
  }
  
  # apply function
  if(type == "exp_decay") {
    parms <- set_filt_exp_decay(radius = radius,
                                zoi_limit = zoi_limit,
                                res = res,
                                half_life = half_life,
                                zoi_hl_ratio = zoi_hl_ratio,
                                min_intensity = min_intensity,
                                max_dist = max_dist)
  }
  
  if(type %in% c("step", "threshold", "circle", "threshold_decay", "step_decay")) {
    parms <- set_filt_step(radius = radius, res = res)
  }
  
  if(type %in% c("bartlett", "batlett_decay", "tent_decay", "linear_decay")) {
    parms <- set_filt_bartlett(radius = radius, res = res)
  }
  
  if(type %in% c("rectangle", "box")) {
    parms <- set_filt_rectangle(radius = radius, res = res)
  }
  
  if(type %in% c("Gauss", "gauss", "gaussian", "gaussian_decay")) {
    parms <- set_filt_gassian_decay(radius = radius,
                                    zoi_limit = zoi_limit,
                                    res = res,
                                    sigma = sigma,
                                    min_intensity = min_intensity,
                                    max_dist = max_dist)
  }
  
  # get parameters
  radius <- parms$radius
  radius_pix <- parms$radius_pix
  size_pix <- parms$size_pix
  
  # create distance matrix
  # distance in pixels to the central cell of the matrix
  dist_mat <- sqrt((matrix(c(1:size_pix), nrow = size_pix, ncol = size_pix, byrow = F) - (radius_pix + 1))^2+
                     (matrix(c(1:size_pix), nrow = size_pix, ncol = size_pix, byrow = T) - (radius_pix + 1))^2)
  # image(dist_mat)
  # plot(terra::rast(dist_mat))
  
  # apply function
  if(type == "exp_decay") {
    dist_mat <- exp(-parms$lambda * dist_mat)
  }
  
  if(type %in% c("step", "threshold", "circle", "threshold_decay", "step_decay")) {
    dist_mat <- 1 * (dist_mat*res <= radius)
  }
  
  if(type %in% c("bartlett", "batlett_decay", "tent_decay", "linear_decay")) {
    dist_mat <- pmax((1 + parms$lambda * dist_mat), 0)
  }
  
  if(type %in% c("rectangle", "box")) {
    dist_mat[] <- 1
  }
  
  if(type %in% c("Gauss", "gauss", "gaussian", "gaussian_decay")) {
    dist_mat <- exp(-parms$lambda * dist_mat**2)
  }
  # image(dist_mat)
  # plot(terra::rast(dist_mat))
  
  # normalize
  if(normalize)
    # dist_mat <- dist_mat/sum(dist_mat[1+radius_pix,])
    dist_mat <- dist_mat/sum(dist_mat)
  
  # image(dist_mat)
  # plot(terra::rast(dist_mat))
  
  # round decimals
  if(!is.null(round_vals))
    if(round_vals >= 0) dist_mat <- round(dist_mat, round_vals)
  # image(dist_mat)
  # plot(terra::rast(dist_mat))
  
  if(save_txt) {
    # save matrix outside R for use within GRASS GIS
    oneimpact::filter_save(filt = dist_mat, radius = radius, type = type,
                           save_format = save_format, save_folder = save_folder,
                           save_file = save_file, parallel = parallel,
                           divisor = divisor, separator = " ")
    
  }
  
  dist_mat
}

set_filt_exp_decay <- function(radius = NULL,
                               zoi_limit = 0.05,
                               half_life = NULL,
                               res = 100,
                               zoi_hl_ratio = NULL,
                               min_intensity = 0.01,
                               max_dist = 200){
  
  # define lambda depending on the input parameter
  if(!is.null(radius)) {
    
    # define radius in terms on number of pixels
    radius <- radius/res
    
    if(is.null(zoi_hl_ratio)) {
      lambda <- log(1/zoi_limit) / radius
    } else {
      half_life <- radius/zoi_hl_ratio
      lambda <- log(2)/half_life
    }
    
  } else {
    
    if(!is.null(half_life)) {
      # define radius or half life, depending on which is given as input
      half_life <- half_life/res
      lambda <- log(2)/half_life
    } else {
      stop("Either both 'radius' and 'zoi_limit' must be specified, or both 'half_life' and 'zoi_hl_ratio'.")
    }
  }
  
  # tmp <- exp(-lambda * c(0:round(half_life*6))/half_life)
  # define radius and size (diameter)
  tmp <- exp(-lambda * c(0:round(2*radius)))
  radius_pix <- min(which(tmp < min_intensity)[1], round(max_dist/res))
  size_pix <- 2*radius_pix + 1
  
  return(list(radius = radius, radius_pix = radius_pix, size_pix = size_pix, lambda = lambda))
}





#load in layers
riv<-vect("mergeSel.shp")
Gras<-vect("grasmere.shp")
riv<-crop(riv,Gras)
land<-rast("lcmGrasmere.tif")
ppgis<-rast("trees_25m_buffer_150m.tif")

rip<-vect("constraintsRip10m.shp")

demGras<-rast("demGrasmere.tif")


#classify LCM data to get broadleaf woodland class
lcmLevel<-data.frame(levels(as.factor(land$lcmGrasmere_1)))

classes<-c(0,0,0,0,1,1,1,1,0,0,0,0,0)

mat<-cbind(as.numeric(lcmLevel$lcmGrasmere_1),classes)

classesWood<-c(0,1,0,0,0,0,0,0,0,0,0,0,0)

matWood<-cbind(as.numeric(lcmLevel$lcmGrasmere_1),classesWood)

lcmWood<-classify(land$lcmGrasmere_1,matWood)


#rasterize river network
rivRast<-rasterize(riv,land,field=0,background=NA)


#get distance-to-woodland layer
lcmWoodDist<-lcmWood

lcmWoodDist[lcmWood==0]<-NA

woodDist<-distance(lcmWoodDist)


#reclassify for cost layer (values from Eycott et al.)
l<-data.frame(levels(as.factor(land$lcmGrasmere_1)))


costs<-c(NA,1,3,10,10,4.35,2.2,2.2,2.44,
         5,10,5,5)

costMat<-cbind(as.numeric(l$lcmGrasmere_1),
               as.numeric(costs))


costLCM<-classify(land$lcmGrasmere_1,costMat)


Q<-quantile(values(woodDist), probs = seq(0, 1, by = 0.05),na.rm=T)




#get 50 metre buffer of river network
riv50<-aggregate(buffer(riv,width=50))

#intersect with constrained riparian zone
riv50Rip<-st_intersection(st_as_sf(riv50),st_as_sf(rip))

riv50Rip<-vect(riv50Rip$geometry)


ppgis<-mask(ppgis,rip)


baselineEHI1<-5.27598601474247
baselineEHI5<-5.98392061874273
baselineEHI10<-7.29562094791212
baselinePropRip<-17.13881
baselineArea<-2757900




mcWood<-function(x){
  
  wood.i<-woodDist
  
  wood.i[wood.i>Q[x]]<-NA
  

  
  wood.i[!is.na(wood.i)]<-1
  wood.i[is.na(wood.i)]<-0
  
  sum(values(wood.i))/sum(values(lcmWood$lcmGrasmere_1))*100
  
  woodNew<-lcmWood+wood.i
  
  
 
  costNewWood<-woodNew-1
  
  costNewWood[costNewWood==-1]<-1 
  
  

  
  costNew<-costLCM*costNewWood
  
  costNew[costNew==0]<-1
  
  
  ####################Focal Edge Section
  
  levels(as.factor(land$lcmGrasmere_1))
  
  edges<-c(NA,0,20.01,49.97,29.97,15.8,15.33,15.35,15.34,
           5,5.01,75.55,75.54)
  
  
  edgeMat<-cbind(as.numeric(l$lcmGrasmere_1),
                 as.numeric(edges))
  
 
  
  edgeLCM<-classify(land$lcmGrasmere_1,edgeMat)
  

  
  edgeNew<-edgeLCM*costNewWood
  

  edges[edges==0]<-NA
  
  
  funData<-data.frame(cbind(edges,costs))
  
  
  #funData
  
  edgeFin<-land$lcmGrasmere_1
  
  costFin<-land$lcmGrasmere_1
  
  edgeFin<-0
  
  
  for(i in unique(edges[!is.na(edges)])){
    
    fun.i<-funData[funData$edges==i,]  
    
    
    edgeRast<-edgeLCM
    #plot(edgeNew)
    edgeRast[edgeRast!=i]<-0
    edgeRast[edgeRast!=0]<-1
    #plot(edgeRast)
    
    edgeRastZero<-edgeRast-1
    #plot(edgeRastZero)
    edgeRastZero[edgeRastZero==-1]<-1
    
    edgeRast<-extend(edgeRast,c(100,100),fill=0)
    edgeRast[is.na(edgeRast)]<-0
    # 
    filterEdge<-filter_create(r=edgeRast,type="exp_decay",radius=i, normalize = TRUE)
    
    
    focalD5Sum<-focal(edgeRast,w=filterEdge,fun="sum")
   
    
    focalD5Sum<-crop(focalD5Sum,ext(land))
    
    focalD5Sum<-focalD5Sum*edgeRastZero
    
    
    
    focalD5SumCost<-focalD5Sum*fun.i[!is.na(fun.i$costs),]$costs

    
    
    
    edgeFin<-edgeFin+focalD5Sum
    
    costFin<-costFin+focalD5SumCost
    
  }
  
  
  edgeFin[edgeFin>1]<-1
  

  
  ####################################
  woodNew<-woodNew-edgeFin
  

  
  woodNew[woodNew<0]<-0
  
  costNew<-costNew+costFin
  

  patches.i<-patches(woodNew,direction=8,zeroAsNA=TRUE)

  patches.i<-as.polygons(patches.i,dissolve=TRUE)
 
  patchesNew<-patches(wood.i,direction=8,zeroAsNA=TRUE)
  
  patchesNew<-as.polygons(patchesNew,dissolve=T)
  
  

  areaNew<-extract(wood.i,patchesNew,method="simple",fun="sum")
  
 
  areaNew<-sum(areaNew$lcmGrasmere_1,na.rm=T)*res(wood.i)[2]^2
  
  areaEdge<-extract(woodNew,patches.i,method="simple",fun="sum")
  
  areaPP<-extract(ppgis,patchesNew,method="simple",fun="max")
  

  
  areaPP<-areaPP$trees_25m
  
  ppProp<-mean(areaPP,na.rm=TRUE)
  

  areaRip<-st_intersection(st_as_sf(riv50),st_as_sf(patches.i))
  
 
  propRip<-sum(as.numeric(st_area(areaRip)))/sum(expanse(riv50))*100
  
  patches.i$Area<-areaEdge[2]*(res(woodNew)[2]^2)

  
  
  
  land_cost <- transition(1/raster(costNew), transitionFunction=mean, 8)
  
  #get least cost path between points 
  
  sites<-SpatialPoints(crds(centroids(patches.i)))
  
  costMat <- matrix(0, nrow=nrow(patches.i), ncol=nrow(patches.i))
  
  n.patch<-nrow(costMat)  
  
  for (i in 1:n.patch) {
    c <- gdistance::costDistance(land_cost, sites[i], sites) 
    costMat[i,] <- c*10

  }
  
  # this is a matrix of least cost distances between patches 
  
  #create basis for probability matrix
  distMat<-as.matrix(costMat, nrow=nrow(costMat), ncol=nrow(costMat)) 
  
  distMat <- apply(distMat, MARGIN=1, FUN=as.numeric) # make sure all elements are numeric here
  

  
  alpha05k= -log(0.05) / 1000
  alpha1k= -log(0.05) / 5000
  alpha5k= -log(0.05) / 10000
  
  
  

  # negative exponential of colonization kernel x distance get probability
  A.prob05k <- exp(-alpha05k * distMat) 
  A.prob1k <- exp(-alpha1k * distMat) 
  A.prob5k <- exp(-alpha5k * distMat) 
  
  # set diag to zero
  diag(A.prob1k) <- 1 
  diag(A.prob05k) <- 1 
  diag(A.prob5k) <- 1 
  
  # final matrix for connectivity graph
  A.prob05k <- as.matrix(A.prob05k)
  A.prob1k <- as.matrix(A.prob1k)
  A.prob5k <- as.matrix(A.prob5k)
  
  # final matrix for connectivity graph
  graph.Aprob05k <- graph_from_adjacency_matrix(A.prob05k, mode="undirected", weighted=T)
  graph.Aprob1k <- graph_from_adjacency_matrix(A.prob1k, mode="undirected", weighted=T)
  graph.Aprob5k <- graph_from_adjacency_matrix(A.prob5k, mode="undirected", weighted=T)
  
  ## calculate RHI
  
  # calculate all shortest paths between nodes
  pstar.mat05k <- distances(graph.Aprob05k, weights= -log(E(graph.Aprob05k)$weight))
  pstar.mat1k <- distances(graph.Aprob1k, weights= -log(E(graph.Aprob1k)$weight))
  pstar.mat5k <- distances(graph.Aprob5k, weights= -log(E(graph.Aprob5k)$weight))
  
  # back-transform to probabilities of connectedness
  pstar.mat05k <- exp(-pstar.mat05k)                                                   
  pstar.mat1k <- exp(-pstar.mat1k)  
  pstar.mat5k <- exp(-pstar.mat5k)  
  # get study area in m2
  AL <- ncell(land)*res(land)[1]^2
  
  # get area vector 
  areaN <- patches.i$Area
  # sum areas from the vector
  areaSum <- sum(areaN)
  areaMean<-mean(areaN,na.rm=T)
  # get product of all patch areas ij and multiply by probabilities above
  PCmat05k <- outer(areaN,areaN) * pstar.mat05k 
  PCmat1k <- outer(areaN,areaN) * pstar.mat1k 
  PCmat5k <- outer(areaN,areaN) * pstar.mat5k 
  
  pcMatSum05k <- sum(PCmat05k)
  pcMatSum1k <- sum(PCmat1k)
  pcMatSum5k <- sum(PCmat5k)
  
  # divide by total area of the study squared to get the PC metric  
  EHI05k <- sqrt(pcMatSum05k) / as.numeric(AL) 
  EHI1k <- sqrt(pcMatSum1k) / as.numeric(AL) 
  EHI5k <- sqrt(pcMatSum5k) / as.numeric(AL) 
  
  EHI05k<-EHI05k*100
  EHI1k<-EHI1k*100
  EHI5k<-EHI5k*100

  EHI1k
  
  df.i<-data.frame(EHI1=(EHI05k-baselineEHI1)/baselineEHI1*100,EHI5=(EHI1k-baselineEHI5)/baselineEHI5*100,EHI10=(EHI5k-baselineEHI10)/baselineEHI10*100,
                   meanDist=meanDist,meanArea=areaMean,sumArea=areaSum,npatch=n.patch,ripProp=propRip,ripPropPC=(propRip-baselinePropRip)/baselinePropRip*100,
                   sumAreaPC=areaNew/baselineArea*100,propStake=ppProp)
  
  
  
  print(x)

  print(df.i)
  return(df.i)
}  


nWood<-1:21



testWood<-lapply(nWood,mcWood)


dffin<-do.call("rbind",testWood)


dffin


write.csv(dffin,"woodData25_20_EdgeFinal.csv")



