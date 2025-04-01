library(terra)
library(sf)
library(gdistance)
library(dplyr)
library(doParallel)
setwd("E:/Cell")
riv<-vect("mergeSel.shp")
Gras<-vect("grasmere.shp")
riv<-crop(riv,Gras)
land<-rast("lcmGrasmere.tif")
ppgis<-rast("trees_25m_buffer_150m.tif")

rip<-vect("constraintsRip10m.shp")

ppgis

demGras<-rast("demGrasmere.tif")

lcmLevel<-data.frame(levels(as.factor(land$lcmGrasmere_1)))

classes<-c(0,0,0,0,1,1,1,1,0,0,0,0,0)

mat<-cbind(as.numeric(lcmLevel$lcmGrasmere_1),classes)

classesWood<-c(0,1,0,0,0,0,0,0,0,0,0,0,0)

matWood<-cbind(as.numeric(lcmLevel$lcmGrasmere_1),classesWood)

lcmWood<-classify(land$lcmGrasmere_1,matWood)

plot(lcmWood)

rivRast<-rasterize(riv,land,field=0,background=NA)

plot(rivRast)

lcmWoodDist<-lcmWood

lcmWoodDist[lcmWood==0]<-NA

woodDist<-distance(lcmWoodDist)

plot(woodDist)

alphaWood= -log(0.05)/50

# negative exponential of colonization kernel x distance get probability

probWood <- exp(-alphaWood * rivDist) 

probWood<-terra::mask(probWood,rip)

plot(probWood)

l<-data.frame(levels(as.factor(land$lcmGrasmere_1)))


costs<-c(NA,1,3,10,10,4.35,2.2,2.2,2.44,
         5,10,5,5)

costMat<-cbind(as.numeric(l$lcmGrasmere_1),
               as.numeric(costs))

costMat

costLCM<-classify(land$lcmGrasmere_1,costMat)


#writeRaster(costLCM,"costLCM_Test.tif", overwrite=TRUE)

gridWood = st_make_grid(rip, c(25, 25), what = "polygons", square = F)

gridWood<-st_intersection(gridWood,st_as_sf(rip))

gridWood<-vect(gridWood)

#plot(gridWood[1,])

gridWoodProb<-terra::extract(woodDist,gridWood,fun="mean",na.rm=T)

gridWood$woodProb<-gridWoodProb[,2]

nrow(gridWood)

hist(gridWood$woodProb)

#gridWood<-gridWood[!is.na(gridWood$woodProb),]

plot(gridWood)

ApplyQuintiles <- function(x) {
  cut(x, breaks=c(quantile(gridWood$woodProb, probs = seq(0, 1, by = 0.05))), 
      labels=as.character(1:20), include.lowest=TRUE)
}

gridWood$Q <- sapply(gridWood$woodProb, ApplyQuintiles)

woodDist<-mask(woodDist,rip)

plot(woodDist)

Q<-quantile(values(woodDist), probs = seq(0, 1, by = 0.05),na.rm=T)

Q#[21]


head(gridWood)

gridWood$Q<-as.numeric(gridWood$Q)

head(gridWood)

nrow(gridWood[gridWood$Q==5,])

#x=16


riv50<-aggregate(buffer(riv,width=50))

plot(riv50)
plot(rip,add=T)

riv50Rip<-st_intersection(st_as_sf(riv50),st_as_sf(rip))

riv50Rip<-vect(riv50Rip$geometry)
#writeVector(riv50Rip,"riv50Rip.shp")
plot(riv50Rip)

plot(rip,add=T)
x=20

head(gridWood)
tail(gridWood)

plot(woodDist)

#ppgis<-rast("trees_25m.tif")

ppgis<-mask(ppgis,rip)

ppgis
x

plot(lcmWood)
plot(woodRiv100)
sum(values(woodRiv100))
sum(values(lcmWood))

plot(woodRiv100)

Q

mcWood<-function(x){
  
  wood.i<-woodDist
  
  wood.i[wood.i>Q[x]]<-NA
  
  #plot(wood.i)
  
  wood.i[!is.na(wood.i)]<-1
  wood.i[is.na(wood.i)]<-0
  
  sum(values(wood.i))/sum(values(lcmWood$lcmGrasmere_1))*100
  
  woodNew<-lcmWood+wood.i
  
  #writeRaster(wood.i,"wood100_6.tif")
  #plot(lcmWood)
  plot(woodNew)
  
  #plot(costNewWood)
  
  costNewWood<-woodNew-1
  
  costNewWood[costNewWood==-1]<-1 
  
  
  #plot(costLCM)
  
  costNew<-costLCM*costNewWood
  
  costNew[costNew==0]<-1
  #plot(costNew)
  
  ####################Focal Edge Section
  
  levels(as.factor(land$lcmGrasmere_1))
  
  edges<-c(NA,0,20.01,49.97,29.97,15.8,15.33,15.35,15.34,
           5,5.01,75.55,75.54)
  
  
  edgeMat<-cbind(as.numeric(l$lcmGrasmere_1),
                 as.numeric(edges))
  
  #edgeMat
  
  edgeLCM<-classify(land$lcmGrasmere_1,edgeMat)
  
  #plot(edgeLCM)
  
  edgeNew<-edgeLCM*costNewWood
  
  #plot(edgeNew)
  
  #plot(land$lcmGrasmere_1)
  #plot(edgeNew)
  
  #edgeLCM<-edgeLCM*costNewWood
  #edges
  
  edges[edges==0]<-NA
  
  
  funData<-data.frame(cbind(edges,costs))
  
  
  #funData
  
  edgeFin<-land$lcmGrasmere_1
  
  costFin<-land$lcmGrasmere_1
  
  edgeFin<-0
  
  costFin<-0
  
  #i=15.80
  
  #edges
  
  for(i in unique(edges[!is.na(edges)])){
    
    fun.i<-funData[funData$edges==i,]  
    #fun.i
    
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
    
    #plot(focalD5Sum)
    
    #fun.i[!is.na(fun.i$costs),]$costs
    
    focalD5SumCost<-focalD5Sum*fun.i[!is.na(fun.i$costs),]$costs
    
    #focalD5SumCost[focalD5SumCost>fun.i[!is.na(fun.i$costs),]$costs]=0 
    
    #plot(focalD5SumCost)
    
    #print(i)
    
    
    
    edgeFin<-edgeFin+focalD5Sum
    
    costFin<-costFin+focalD5SumCost
    
  }
  
  
  edgeFin[edgeFin>1]<-1
  
  plot(edgeFin)
  
  plot(costFin)
  
  ####################################
  woodNew<-woodNew-edgeFin
  
  #plot(woodNew)
  
  woodNew[woodNew<0]<-0
  
  costNew<-costNew+costFin
  
  #plot(costNew)
  #writeRaster(costNew,"costNewTest.tif")
  patches.i<-patches(woodNew,direction=8,zeroAsNA=TRUE)
  #patches.i<-patches(lcmWood,direction=8,zeroAsNA=TRUE)
  patches.i<-as.polygons(patches.i,dissolve=TRUE)
  #plot(patches.i)
  
  patchesNew<-patches(wood.i,direction=8,zeroAsNA=TRUE)
  
  patchesNew<-as.polygons(patchesNew,dissolve=T)
  
  
  #plot(wood.i)
  #plot(patchesNew,add=T)
  areaNew<-extract(wood.i,patchesNew,method="simple",fun="sum")
  
  #plot(woodNew)
  areaNew<-sum(areaNew$lcmGrasmere_1,na.rm=T)*res(wood.i)[2]^2
  
  areaEdge<-extract(woodNew,patches.i,method="simple",fun="sum")
  
  areaPP<-extract(ppgis,patchesNew,method="simple",fun="max")
  
  #ppgis
  
  areaPP<-areaPP$trees_25m
  
  ppProp<-mean(areaPP,na.rm=TRUE)
  
  #plot(riv50Rip)
  #plot(patches.i,add=T,col="red")
  areaRip<-st_intersection(st_as_sf(riv50),st_as_sf(patches.i))
  
  #plot(areaRip)
  
  #plot(riv50,add=T)
  
  propRip<-sum(as.numeric(st_area(areaRip)))/sum(expanse(riv50))*100
  
  patches.i$Area<-areaEdge[2]*(res(woodNew)[2]^2)
  #sum(patches.i$Area)
  
  #print(propRip)
  #print(sum(patches.i$Area))/10000
  
  
  
  
  
  land_cost <- transition(1/raster(costNew), transitionFunction=mean, 8)
  
  #get least cost path between points 
  
  sites<-SpatialPoints(crds(centroids(patches.i)))
  
  costMat <- matrix(0, nrow=nrow(patches.i), ncol=nrow(patches.i))
  
  n.patch<-nrow(costMat)  
  
  for (i in 1:n.patch) {
    c <- gdistance::costDistance(land_cost, sites[i], sites) 
    costMat[i,] <- c*10
    #print(i)
  }
  
  # this is a matrix of least cost distances between patches 
  
  #create basis for probability matrix
  distMat<-as.matrix(costMat, nrow=nrow(costMat), ncol=nrow(costMat)) 
  
  distMat <- apply(distMat, MARGIN=1, FUN=as.numeric) # make sure all elements are numeric here
  
  
  
 
  
  meanDist<-extract(rivDist,patchesNew,method="simple",fun="max")
  
  meanDist<-max(meanDist[,2]) 
  
  alpha05k= -log(0.05) / 1000
  alpha1k= -log(0.05) / 5000
  alpha5k= -log(0.05) / 10000
  
  
  
  
  # init empty matrix for adjacency matrix 
  #A.prob <- matrix(0, nrow=nrow(distMat), ncol=ncol(distMat))
  
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
  
  
  #df.i
  print(x)
  #print(EHI1k)
  print(df.i)
  return(df.i)
}  

sum(values(lcmWood))*res(lcmWood)[2]^2

nRiv<-1:21



testWood<-mclapply(nRiv,mcWood,mc.cores=1)


baselineEHI1<-5.27598601474247
baselineEHI5<-5.98392061874273
baselineEHI10<-7.29562094791212
baselinePropRip<-17.13881
baselineArea<-2757900


dffin<-do.call("rbind",testWood)


dffin


write.csv(dffin,"woodData25_20_EdgeFinal.csv")



riv100<-rast("E:/Cell/riv100_4.tif")
wood100<-rast("E:/Cell/wood100_6.tif")

plot(riv100)
sum(values(riv100))
sum(values(wood100))
sum(values(rivWood))
plot(rip,add=T)
sum(values(lcmWood$lcmGrasmere_1))

rivWood<-riv100+wood100
rivWood[rivWood>1]<-1
plot(rivWood)
