# Alpha Diversity
# Script for course

library(vegan)


# PBDB collection 31618
# Bangtoupo F30, Qingyan, China: Pelsonian - Illyrian, China
pt <- "https://fau-paleo.github.io/apw_2022/data/2_paleodiversity/"
dat <- read.csv(paste0(pt,"Triassic_div.csv"), header=TRUE)

attach(dat)

# Compute diversity indices manually
# Immediate values
 S <- nrow(dat) # = number of species
 N <- sum(Individuals) # N
 cl <- length(levels(factor(Class))) # Number of classes
 gen <- length(levels(factor(Genus))) # Number of genera
 menh <- S/sqrt(N) # Menhinick Index
 marg <- (S-1)/log(N) # Margalef Index

# Prepare data for other indices
 psp <- Individuals/N  # Proportions of each individual in assemblage
 sha <- psp*log(psp)  # For Shannon-Index
 p2 <- psp^2          # For Hurlbert's PIE and Simpson's index

 # Calculate indices
      H <- -sum(sha)             ## Shannon H
      J <- H/log(S)              ## Shannon J
      E <- exp(H)/S              ## Equitability
      pD <- max(psp)             ## Berger-Parker Dominanz
      D <- sum(p2)               ## Simpson's D
      PIE <- S/(S-1)*(1-sum(p2)) ## Hurlbert's PIE
 
  # Fisher's alpha  S = a ln(1+N/a) 
     alpha <- 1
     F <- 1
     while(F <= S) {
           F <- alpha * log(1+N/alpha)
           alpha <- alpha+0.01
            }

  # Hill number
     # D=(SUM p_i^q)^1/(1-q)
     q = 0 # Modify here, Do not use q= 1 but rather q = 0.999
     sum(psp^q)^(1/(1-q))
     
    # Create a loop to compute Hill numbers for q values between 0 and 10
     
     
     
 ##### Subsampling ####
 # Get one rarefied diversity value for 100 individuals
    gensp <- paste(Genus, Species)
     abu <- rep(gensp, Individuals)
     
     trial <- 1000 # subsampling trials
     quota <- 100 # subsampling quota
     div <- numeric()
     
     for (i in 1:trial) {
       z <- sample(abu, quota) # subsampling without replacement
       div[i] <- length(levels(factor(z)))    
     }  
      
     
     # Expected number of species according to Good (1953)
     E.Sm <- S - sum((1-psp)^quota) 
     
   # Empirical rarefaction
     # Expected number of species if x individuals are collected

     trial <- 200 # subsampling trials
        rardiv <- numeric()
         erdiv <- numeric()
       count <- 1 
        sq <- seq(1, 651, by=10)
     for (j in sq)  {    # loop for the quota

      div <- numeric(trial)

     for (i in 1:trial) { # loop for the trial
       z <- sample(abu, j) # subsampling without replacement
       div[i] <- length(levels(factor(z)))    
           }  

      # Rarefied diversity
       rardiv[count] <- mean(div)
       erdiv[count] <- sd(div)
      count <- count+1
  }

        plot(sq, rardiv, type="l")
        segments(sq, rardiv-erdiv, sq, rardiv+erdiv)

op <- par(mfrow=c(2,1), lwd=2, mar=c(4,4,1,1))
 plot(sq, rardiv, type="l", xlab="", ylab="S")
 plot(sq, rardiv, log="xy", type="l", xlab="N", ylab="S")
  x <- seq(600); y <- x^0.8
  points(x,y, type="l", lty=2)
par(op)

#### Shareholder quorum subsampling ####
source(paste0(pt,"sqs_Holland.R"))
 SQ <- sqs(Individuals, 0.7, 100) # Shareholder quorum

# what is the estimated coverage of the full sample?
n1 <- length(which(Individuals==1))
u <- 1- n1/N # Good's u


# Summarize results
resname <- c("Species", "Individuals", "H", "J", "Simpson's D", "1-D", "PIE", "Berger-Parker", "Fisher's alpha",
   "Margalef", "Menhinick", "Rarefied.200", "SD_Rarefied.200", "SQS0.7")                       
resu <- c(S, N, H, J, D, 1-D, PIE, pD, alpha, marg, menh, rardiv[21], erdiv[21], SQ)
resround <- round(resu, 3)
x <- data.frame(cbind(resname, resround))

write.table(x, file="Result_div1.csv", sep=",", row.names=FALSE)
getwd()


##### Vegan package ####
library(vegan)
specnumber(Individuals)
diversity(Individuals)
diversity(Individuals, "simpson")
fisher.alpha(Individuals)
rarefy(Individuals, seq(5,600,by=5)) # rarefaction


####
#### Rank-abundance curves ####
 z <- sort(Individuals, decreasing=T)
 p.z <- sort(psp, decreasing = T)
 op <- par(cex=1.5, lwd=1.5, mar=c(5,4,2,1), mfrow=c(2,2))

 # Linear plot
 plot(z, type="b", xlab="Rank", ylab="Abundance", col="red", lwd=2)

 # Log plot
 plot(z, type="b", log="y", xlab="Rank", ylab="Abundance", col="red", lwd=2)

 # Cummulative plot
 plot(cumsum(z), type="l", log="y", xlab="Rank", ylab="Abundance", col="red", lwd=2)
 
 # Log proportions
 plot(p.z, type="b", log="y", xlab="Rank", ylab="Abundance", col="red", lwd=2)
 
par(op)


# Create model distributions given data
 # Normal distribution
   Nm <- mean(Individuals)
   Nsd <- sd(Individuals)
    nsa <- rnorm(S, Nm, Nsd)
    hist(nsa)
     z2 <- sort(nsa, decreasing=T)
     plot(z2, type="l", log="y", xlab="Rank", ylab="Abundance", col="red", lwd=2)
 # Lognormal distribution
    lN <- log(Individuals)
    lNm <- mean(lN)
    lNsd <- sd(lN)
    lsa <- rlnorm(S, lNm, lNsd) 
    hist(lsa)
     z3 <- sort(lsa, decreasing=T)
     plot(z3, type="l", log="y", xlab="Rank", ylab="Abundance", col="red", lwd=2)
     lsa <- sort(round(lsa), decreasing=T)+1
 # Geometric series
   # Determine rate
       zn <- z[1] # use dominant
       rt <- S/zn/min(Individuals) # overall slope
       z4 <- numeric()
       for (i in 1:S) z4 <- c(zn, z4/rt)
       # use regression model
       reg <- lm(log(z) ~ seq(1,length(z)))
       reg2 <- exp(predict(reg))
     plot(reg2, type="l", log="y", xlab="Rank", ylab="Abundance", col="red", lwd=2)

   

### Multiplot
plot(z, type="b", log="y", xlab="Rank", ylab="Abundance", col="red", lwd=2)
 points(z2, type="l", col="grey") # Broken Stick
 points(z3, type="l", col="blue") # 
 points(reg2, type="l", col="green")

library(vegan)
 fit <- radfit(Individuals)
 plot(fit)
 
fit
 
detach(dat)
 

 
   