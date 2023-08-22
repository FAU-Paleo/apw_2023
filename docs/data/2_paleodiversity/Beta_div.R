# Beta Diversity

pt <- "https://fau-paleo.github.io/apw_2023/data/2_paleodiversity/"
rads <- read.csv(paste0(pt,"Radiolarians.csv"), header=TRUE)

# Some manipulations for later
row.names(rads) <- rads$Sample
rads <- as.matrix(rads[,2:ncol(rads)])

# for manual computations
temp <- as.data.frame(t(rads))

####### Manual computation of similarities
 m.alpha <- length(which(temp >0))/ncol(temp) # mean alpha
  gamma <- nrow(temp) 
  beta <- gamma/m.alpha

####
  # Similarity and dissimilarity of two samples
  # Binary metrics
  test <- cbind(temp$Antarct.Tith.K201, temp$Antarct.Tith.K25) 

  a <- length(which(test[,1]>0 & test[,2]>0))
  b <- length(which(test[,1]>0 & test[,2]==0))
  c <- length(which(test[,1]==0 & test[,2]>0))
  
  
  jac <- 1-(a/(a+b+c))
  sor <- (b + c)/(2*a + b + c)
  sim <- min(b,c)/(a + min(b,c))
  nes <- sor-sim

 # Percent similarity
  # Transform to proportions
  p.test <- prop.table(test, 2)
  psim <- sum(pmin(p.test[,1], p.test[,2]))
  dis <- 1-psim



 # Use the dist function first (binary)
  dis <- dist(rads) # Euclidean
  dis <- dist(rads, method="binary") # This is Jaccard dissimilarity (check)
  
  # A simple cluster analysis
  clus <- hclust(dis) # "complete" is the default
  windows()
    plot(clus)

# Load the vegan package for Bray-Curtis distance 
library(vegan)

    
    ##### Cluster Analysis ####
    # Compute dissimilarity matrix
    rads <- prop.table(rads, 1) # Only if proportional data are being used, the bray method is equivalent to proportional dissimilary
    
    dis <- vegdist(rads) # default is Bray-Curtis
    r <- vegdist(t(rads))
    
    # Do the cluster analysis  
    clus <- hclust(dis, "single")
    plot(clus)
    
    
    X11()
    clus <- hclust(dis, "ward.D") # recommended method
    plot(clus, hang=-1) 
    
    
    # R mode and Q mode combined
    library(paleotree)
    twoWayEcologyCluster(r, dis, propAbund=t(rads), clustMethod="ward.D2")
    
####  Ordination ############################
  
  # Principal coordinate analysis
  eurodist
  loc <- cmdscale(eurodist)
  plot(loc, type="n")
  text(loc, rownames(loc), cex = 0.6)
  
  x <- loc[, 1] 
  y <- -loc[, 2]
  plot(x, y, type = "n", xlab = "", ylab = "", asp = 1, axes = FALSE, main = "cmdscale(eurodist)") 
  text(x, y, rownames(loc), cex = 0.6)
  
  
  
# Non-metric multidimensional scaling
  dis <- vegdist(rads)
  mMDS <- metaMDS(rads, k=2) # vegan 
  
  # How stressed is the analysis?
  mMDS$stress
  stressplot(mMDS)
  
  X11()
  plot(mMDS, display="sites", type="t")
   mtext(text=paste("Stress =", round(mMDS$stress, 3)), side=3)
  
   MDS <- cmdscale(dis)
   plot(MDS, type="n")
   text(MDS, rownames(MDS), cex = 0.6)

   MDS <- cmdscale(dis, eig=T)
   round(MDS$eig*100/sum(MDS$eig),1)

  
# Sites and species together   
   plot(mMDS)
   
 # Are groups identified by visual inspection statistically robust?
 # Define groups
   gr <- rep(1, nrow(rads)) # standard (Alps + Philippines)
   # Antarctic samples
   gr[grep("Antarc", row.names(rads))] <- 2
   # Add a third group (Oman)
   gr[grep("Oman", row.names(rads))] <- 3
   
   
   ordihull(mMDS, gr, display="sites", draw="polygon", 
            col=c("red", "blue", "green"))
   
   # Convex hull around groups in NMDS
   plot(mMDS, type="n")
   ordihull(mMDS, gr, display="sites", draw="polygon", 
            col=c("red", "blue", "green", "pink"))  
   

    # adonis testing
     adonis2(dis ~ gr)
     # summary(ado)
      
      groups <-  anosim(rads, gr)
      summary(groups)
      plot(groups)
      
      # parametric  
      mod <- betadisper(dis, gr) 
      plot(mod, main = "Metric MDS")
      
   
#####################
      #### Cassian Stuff ####
      cass <- read.csv(paste0(pt,"Cassian.csv"), header=TRUE)
      table(cass$environment)
      
      # Richness by site      
      table(cass$collection_no)
      # Shannon-wiener diversity 
      H <- tapply(cass$abund_value, cass$collection_no, diversity)
      
      
      gensp <- paste(cass$occurrence.genus_name, cass$occurrence.species_name)
      
      
      #create sites by species matrix
      ta.cas <- tapply(cass$abund_value, list(cass$collection_no, factor(gensp)), sum)
      ta.cas[is.na(ta.cas)] <- 0  
      
      # Richness
      table(cass$collection_no)
      # Shannon Diversity
      H <- diversity(ta.cas)
      
      #beta diversity
      prop.sp<-prop.table(ta.cas, margin=1) #proportions
      sp.beta<-vegdist(prop.sp) #relative bray-curtis

