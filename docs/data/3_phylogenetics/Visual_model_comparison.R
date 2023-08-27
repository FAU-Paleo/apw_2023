
library(dplyr)
library(viridisLite)
library(viridis)
library(phangorn)
library(permute)
library(lattice)
library(vegan)



## read in trees 
mk_one <- ape::read.tree(ADD PATH TO *posterior.trees)
mk_two <- ape::read.tree(ADD PATH TO *posterior.trees) 

### sample from posterior - here we are sample 100 trees but this number is up to you. 
### Increasing the number can massivly increase processing time. 
SAMPLE <- function(data){
  data[sample(length(data),100)]
}

mk_one.sample <- SAMPLE(mk_one)
mk_two.sample <- SAMPLE(mk_two)



#Prepare for graphic representation
col <- c(viridis::viridis(2))
pch <- c(0,1)
group<-rep(c(1,2),each=100)


#Robinson-Foulds ####  
#RF distances ####
Results <- function(df1,df2,ROWNAMES,COLNAMES){
  RF <- data.frame(matrix(ncol=length(df2),nrow=length(df1)))
  for (i in 1:length(df1)){
    for (j in 1:length(df2)){
      RF[i,j] <- phangorn::RF.dist(ape::unroot(df1[[i]]),ape::unroot(df2[[j]]),
                                   normalize=TRUE, check.labels=TRUE,rooted=FALSE)
      rownames(RF)[i] <- paste(ROWNAMES,i,sep="_")
      colnames(RF)[j] <- paste(COLNAMES,j,sep="_")
    }
  }
  return(RF)
}


# This part can take a while
RF1 <- Results(mk_one.sample,mk_one.sample,"mk_one", "mk_one")
RF2 <- Results(mk_two.sample,mk_two.sample,"mk_two","mk_two")
RF3 <- Results(mk_one.sample,mk_two.sample,"mk_one","mk_two")


join.tot <- cbind(rbind(RF1,t(RF3)),
                        rbind(t(RF2),t(RF3)))



dist.join <- stats::as.dist(join.tot, diag=TRUE)


join.group <- as.factor(rep(c("mk_one","mk_two"),
                            each=100))

# Your may get a warning message here but that's ok!
beta.test <- vegan::betadisper(dist.join, join.group)


plot(beta.test, axes = c(1,2), cex = 0.7, pch=pch, col=col,
     lty = "solid", lwd = 1, hull = FALSE, ellipse = FALSE,
     segments = TRUE, label = FALSE)
legend("bottomleft", legend=c("mk_one", "mk_two"), col=col, pch=pch, cex =0.5)









