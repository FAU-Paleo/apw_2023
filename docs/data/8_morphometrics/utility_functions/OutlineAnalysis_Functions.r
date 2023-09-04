## HEAD #####################################################################################################
#
# FUNCTION SET DESCRIPTION
#	Applying Outline Analyses in R
#
# DATASET DESCRIPTION
#	Morphometric datasets
#
# FURTHER READING
#	Claude, J. (2008) "Morphometrics with R". Gentleman, R., Hornik, K., and...
#		Parmigiani, G. (eds) "Use R!", vol. ii, 316 pp. (Springer).
#
#	Most of the functions in this sheet are at least based on functions developed by...
#		Claude (2008), and this publication should thus always be cited when using the code...
#		provided here. References to pages in that book are given in parentheses wherever necessary.
#	
# METADATA
#	Author: Manuel F. G. Weinkauf
#	E-mail: weinkauf.scientific@gmail.com
#	R-version: 4.2.1
#	RStudio-version: 2022.07.1
#	Code-version: 2.2
#	Date of last update: 31 October 2016

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

## BODY #####################################################################################################

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData/Outline")

#########################################################################
# Performing Zahn--Roskies Fourier Analysis (pp. 218f).                 #
# Necessary input variables:                                            #
#    Input: Digitized outline points coordinates.                       #
#           *matrix*                                                    #
#    Harmonics: Desired number of harmonics to be computed.             #
#               *numeric (integer)*                                     #
# Output data: List containing harmonics (ao, an, bn), variation of...	#
#              tangent angle (phi), distances along perimeter in...     #
#              radians (t), perimeter of outline, and first tangent...  #
#              angle thetao.                                            #
# Input dataset: Matrix of equally spaced x-y coordinates along...      #
#                digitized outline                                      #
#########################################################################

ZRfourier<-function (Input, Harmonics) {
	p<-dim(Input)[1]
	an<-numeric(Harmonics)
	bn<-numeric(Harmonics)
	tangvect<-Input-rbind(Input[p,], Input[-p,])
	perim<-sum(sqrt(apply((tangvect)^2, 1, sum)))
	v0<-(Input[1,]-Input[p,])
	tet1<-Arg(complex(real=tangvect[,1], imaginary=tangvect[,2]))
	tet0<-tet1[1]
	t1<-(seq(0, 2*pi, length=(p+1)))[1:p]
	phi<-(tet1-tet0-t1)%%(2*pi)
	ao<-2*sum(phi)/p
	for (i in 1:Harmonics) {
		an[i]<-(2/p)*sum(phi*cos(i*t1))
		bn[i]<-(2/p)*sum(phi*sin(i*t1))
	}
	list(ao=ao, an=an, bn=bn, phi=phi, t=t1, perimeter=perim, thetao=tet0)
}

#########################################################################
# Performing elliptic Fourier analysis (pp. 222f).                      #
# Necessary input variables:                                            #
#    Input: Digitized outline points coordinates.                       #
#           *matrix*                                                    #
#    Harmonics: Desired number of harmonics to be computed.             #
#               *numeric (integer)*                                     #
#               default=Number of outline points/2                      #
# Output data: List containing centroid coordinates (ao, co) and...     #
#              values for four coeffecients per harmonic.               #
# Input dataset: Matrix of equally spaced x-y coordinates along...      #
#                digitized outline.                                     #
#########################################################################

efourier<-function (Input, Harmonics=nrow(Input)/2) {
	#Consistency test
	if (Harmonics>nrow(Input)) {stop("There cannot be more harmonics than original data points!")}

	#Perform basic calculations
	p<-dim(Input)[1]
	Dx<-Input[,1]-Input[c(p, (1:p-1)),1]
	Dy<-Input[,2]-Input[c(p, (1:p-1)),2]
	Dt<-sqrt(Dx^2+Dy^2)
	t1<-cumsum(Dt)
	t1m1<-c(0, t1[-p])
	T<-sum(Dt)
	
	#Calculate harmonics
	an<-bn<-cn<-dn<-numeric(Harmonics)
	for (i in 1:Harmonics) {
		an[i]<-(T/(2*pi^2*i^2))*sum((Dx/Dt)*(cos(2*i*pi*t1/T)-cos(2*pi*i*t1m1/T)))
		bn[i]<-(T/(2*pi^2*i^2))*sum((Dx/Dt)*(sin(2*i*pi*t1/T)-sin(2*pi*i*t1m1/T)))
		cn[i]<-(T/(2*pi^2*i^2))*sum((Dy/Dt)*(cos(2*i*pi*t1/T)-cos(2*pi*i*t1m1/T)))
		dn[i]<-(T/(2*pi^2*i^2))*sum((Dy/Dt)*(sin(2*i*pi*t1/T)-sin(2*pi*i*t1m1/T)))
	}
	ao<-2*sum(Input[,1]*Dt/T)
	co<-2*sum(Input[,2]*Dt/T)
	
	#Resturn harmonic coefficients
	list(ao=ao, co=co, an=an, bn=bn, cn=cn, dn=dn)
}

#########################################################################
# Performing normalized elliptic Fourier analysis (p. 226).             #
# Necessary functions: efourier                                         #
# Necessary input variables:                                            #
#    Input: Digitized outline points coordinates.                       #
#           *matrix*                                                    #
#    Harmonics: Desired number of harmonics to be computed.             #
#               *numeric (integer)*                                     #
#               default=Number of outline points/2                      #
#    Rotation: Should rotation be normalized according to the longest...#
#              axis of the first harmonic (default) or a baseline.      #
#              *character*, either of "Ellipse" or "Baseline"           #
#              default="Ellipse"                                        #
#    Baseline: Vector containing the coordinates for the baseline (in...#
#              order x1, y1, x2, y2) if Rotation=="Baseline".           #
#              *vector* of length 4                                     #
#              default=NULL                                             #
#    start: Is the position of the starting point to be preserved?      #
#           *logical*                                                   #
#           TRUE: Starting point must be preserved                      #
#           FALSE: Starting point can not be preserved                  #
#           default=TRUE                                                #
# Output data: List containing values for four standardised...          #
#              coefficients (A-D) per harmonic; magnitude of semi-...   #
#              major axis of first harmonic (size); angle between...    #
#              starting point and semi-major axis of first harmonic...  #
#              (theta); orientation of first harmonic ellipse (psi);... #
#              harmonic coefficients ao and co.                         #
# Input dataset: Matrix of equally spaced x-y coordinates along...      #
#                digitized outline.                                     #
#########################################################################

NEF<-function (Input, Harmonics=nrow(Input)/2, Rotation="Ellipse", Baseline=NULL, start=TRUE) {
	#Consistency test
	if (Harmonics>nrow(Input)) {stop("There cannot be more harmonics than original data points!")}
	if (Rotation!="Ellipse" & Rotation!="Baseline") {stop("Rotation must be either 'Ellipse' or 'Baseline'!")}
	if (Rotation=="Baseline" & is.null(Baseline)) {stop("For Rotation='Baseline' a Baseline vector of length 4 must be provided!")}
	
	#Calculate elliptic Fourier reconstruction
	ef<-efourier(Input, Harmonics)
	A1<-ef$an[1]
	B1<-ef$bn[1]
	C1<-ef$cn[1]
	D1<-ef$dn[1]
	
	#Calculate rotation angle
	{if (Rotation=="Ellipse") {theta<-(0.5*atan(2*(A1*B1+C1*D1)/(A1^2+C1^2-B1^2-D1^2)))%%pi}
	else {
		m.base<-(Baseline[4]-Baseline[2])/(Baseline[3]-Baseline[1])
		m.null<-0
		theta<-as.numeric((0.5*atan((m.base-m.null)/(1+m.base*m.null)))%%pi)
	}
	}
	
	#Calulate rotation matrix
	phaseshift<-matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
	M2<-matrix(c(A1, C1, B1, D1), 2, 2)%*%phaseshift
	v<-apply(M2^2, 2, sum)
	if (v[1]<v[2]) {theta<-theta+pi/2} 
	theta<-(theta+pi/2)%%pi-pi/2
	Aa<-A1*cos(theta)+B1*sin(theta)
	Cc<-C1*cos(theta)+D1*sin(theta)
	scale<-sqrt(Aa^2+Cc^2)
	psi<-atan(Cc/Aa)
	if (Aa<0){psi<-psi+pi}
	size<-(1/scale)
	rotation<-matrix(c(cos(psi), -sin(psi), sin(psi), cos(psi)), 2, 2)
	
	#Apply rotation to shape
	A<-B<-C<-D<-numeric(Harmonics)
	if (start==TRUE) {theta<-0}
	for (i in 1:Harmonics) {
    		mat<-size*rotation%*%matrix(c(ef$an[i], ef$cn[i], ef$bn[i], ef$dn[i]), 2, 2)%*%matrix(c(cos(i*theta), sin(i*theta), -sin(i*theta), cos(i*theta)), 2, 2)
		A[i]<-mat[1,1]
		B[i]<-mat[1,2]
		C[i]<-mat[2,1]
		D[i]<-mat[2,2]
	}
	
	#Return normalizes Fourier reconstruction
	list(A=A, B=B, C=C, D=D, size=scale, theta=theta, psi=psi, ao=ef$ao, co=ef$co, Aa=Aa, Cc=Cc)
}

#########################################################################
# Reconstructing outline from Zahn--Roskies Fourier coefficients...     #
#    (p. 219).                                                          #
# Necessary input variables:                                            #
#    ao: Harmonic coefficient a0                                        #
#        *numeric (real)*                                               #
#    an, bn: Set of coefficients per harmonic.                          #
#            *list*                                                     #
#    Harmonics: Desired number of harmonics to be used for outline...   #
#               reconstruction.                                         #
#               *numeric (integer)*                                     #
#    Points: Desired number of outline points to be reconstructed.      #
#            IMPORTANT: Due to the applied mathematics value must be... #
#                       exactly the same as the number of outline...    #
#                       points originally extracted and used to...      #
#                       perform the Z--R Fourier Analysis.              #
#            *numeric (integer)*                                        #
#    thetao: First tangent angle thetao.                                #
#            *real*                                                     #
#            default=0                                                  #
# Output dataset: List of outline points reconstructed for harmonic...  #
#                 coefficients.                                         #
# Input dataset: List containing parameters from a Zahn--Roskies...     #
#                Fourier analysis.                                      #
#########################################################################

iZRfourier<-function (ao, an, bn, Harmonics, Points, thetao=0) {
	theta<-seq(0, 2*pi, length=Points+1)[-(Points+1)]
	harm<-matrix(NA, Harmonics, Points)
 	for (i in 1:Harmonics) {
		harm[i,]<-an[i]*cos(i*theta)+bn[i]*sin(i*theta)
	}
	phi<-(ao/2)+apply(harm, 2, sum)
	vect<-matrix(NA, 2, Points)
	Z<-complex(modulus=(2*pi)/Points, argument=phi+theta+thetao)
	Z1<-cumsum(Z)
	list(angle=theta, phi=phi, X=Re(Z1), Y=Im(Z1))
}

#########################################################################
# Reconstructing outline from elliptic Fourier coefficients (p. 223)    #
# Necessary input variables:                                            #
#    an, bn, cn, dn: Set of coefficients per harmonic.                  #
#                    *list*                                             #
#    Harmonics: Desired number of harmonics to be used for outline...   #
#               reconstruction.                                         #
#               *numeric (integer)*                                     #
#    Points: Desired number of outline points to be reconstructed.      #
#            *numeric (integer)*                                        #
# Output dataset: List of outline points reconstructed for harmonic...  #
#                 coefficients.                                         #
# Input dataset: List containing centroid coordinates (a0, co) and...   #
#                values for four coeffecients per harmonic; number of...#
#                harmonics to be used; number of points to be...        #
#                reconstructed.                                         #
#########################################################################

iefourier<-function (an, bn, cn, dn, Harmonics, Points, ao=0, co=0) {
	theta<-seq(0,2*pi, length=Points+1)[-(Points+1)]
	harmx<-matrix(NA, Harmonics, Points)
	harmy<-matrix(NA, Harmonics, Points)
	for (i in 1:Harmonics) {
		harmx[i,]<-an[i]*cos(i*theta)+bn[i]*sin(i*theta)
		harmy[i,]<-cn[i]*cos(i*theta)+dn[i]*sin(i*theta)
	}
	x<-(ao/2)+apply(harmx, 2, sum)
	y<-(co/2)+apply(harmy, 2, sum)
	list(x=x, y=y)
}

#########################################################################
# Calculating cumulative power of harmonics (p. 228f)                   #
# Necessary input variables:                                            #
#    an, bn, cn, dn: Set of coefficients per harmonic.                  #
#                    *list*                                             #
#    first: Should the first harmonic be included in the calculations?  #
#           *logical*                                                   #
#           FALSE: Do not include first harmonic.                       #
#           default=FALSE                                               #
#    last: Last harmonic for which power should be calculated.          #
#          *numeric (integer)*                                          #
#          default=20                                                   #
# Output data: List containing cumulative power of harmonics.           #
# Input dataset: List of harmonics from elliptic Fourier analysis...    #
#                (can, but do not have to be normalized).               #
#########################################################################

Power<-function(an, bn, cn, dn, first=FALSE, last=20) {
	P<-(an^2+bn^2+cn^3+dn^2)/2
	{if (first==FALSE){PowerSum<-(cumsum(P[-1])/sum(P[-1]))[1:last]}
		else{PowerSum<-(cumsum(P)/sum(P))[1:last]}
	}
	return(PowerSum)
}

#########################################################################
# Averaging data from several measurements and calculating...           #
#    measurement error (pp. 63--67)                                     #
# Necessary input variables:                                            #
#    Objects: List of data set names containing harmonics.              #
#             *vector* of *characters*                                  #
#    Centroids: Do the last two columns contain centroid coordinates?   #
#               *logical*                                               #
#               TRUE: Centroid coordinates contained in last two...     #
#                     columns.                                          #
#               FALSE: Data contain coefficients only.                  #
#               default=TRUE                                            #
#    Output: Name of output file with averaged Fourier coefficients...  #
#            to be exported.                                            #
#            *character*                                                #
# Output data: File with average fourier coefficients, error...         #
#              calculation.                                             #
# Input dataset: Files with elliptic Fourier coefficients, harmonics... #
#                in columns, specimens in rows, last two columns may... #
#                contain centroid coordinates.                          #
#########################################################################

OutlineAverage<-function (Objects, Centroids=TRUE, Output) {
	#Read data into array
	Read.Obj<-function (X) {eval(parse(text=X))}
	TT<-Read.Obj(Objects[1])
	Dat<-array(NA, dim=c(dim(TT)[1], dim(TT)[2], length(Objects)))
	for (i in 1:(length(Objects))) {
		Dat[,,i]<-as.matrix(Read.Obj(Objects[i]))
	}
	
	#Calculate mean harmonics
	EFA.Mean<-apply(Dat, c(1, 2), mean)
	colnames(EFA.Mean)<-colnames(TT)
	rownames(EFA.Mean)<-rownames(TT)
	write.table(EFA.Mean, paste(Output, "MeanCoefficients.txt", sep="_"))
	
	#Eliminate centroid columns
	if (Centroids==TRUE) {DatFull<-Dat; Dat<-Dat[,-c(dim(Dat)[2]-1, dim(Dat)[2]),]}
	
	#Define factors
	Session.factor<-gl(dim(Dat)[3], (dim(Dat)[1]))
	Individual.factor<-as.factor(rep((1:dim(Dat)[1]), dim(Dat)[3]))
	
	#Calculate relative errors of replications
	RelErr<-list()
	for (i in 1:(dim(Dat)[2])) {
		Hm<-vector(length=0)
		for (j in 1:(dim(Dat)[3])) {
			Hm<-append(Hm, as.vector(t(Dat[,i,j])))
		}
		SE<-summary(aov(Hm~Session.factor))
		{if (SE[[1]][2,3]>SE[[1]][1,3]) {pSE<-1} else {pSE<-0}}
		ME<-summary(aov(Hm~Individual.factor))
		{if (ME[[1]][1,3]>=ME[[1]][2,3]) {pME<-1} else {pSE<-0}}
		s2within<-MSwithin<-ME[[1]][2,3]
		MSamong<-ME[[1]][1,3]
		s2among<-(MSamong-MSwithin)/2
		Err<-s2within/(s2within+s2among)*100
		{if (i%%4==1) {RelErr$F1$Sessionfactor<-append(RelErr$F1$Sessionfactor,pSE);RelErr$F1$Individualfactor<-append(RelErr$F1$Individualfactor,pME);RelErr$F1$RelativeError<-append(RelErr$F1$RelativeError,Err)}
		else if (i%%4==2) {RelErr$F2$Sessionfactor<-append(RelErr$F2$Sessionfactor,pSE);RelErr$F2$Individualfactor<-append(RelErr$F2$Individualfactor,pME);RelErr$F2$RelativeError<-append(RelErr$F2$RelativeError,Err)}
		else if (i%%4==3) {RelErr$F3$Sessionfactor<-append(RelErr$F3$Sessionfactor,pSE);RelErr$F3$Individualfactor<-append(RelErr$F3$Individualfactor,pME);RelErr$F3$RelativeError<-append(RelErr$F3$RelativeError,Err)}
		else if (i%%4==0) {RelErr$F4$Sessionfactor<-append(RelErr$F4$Sessionfactor,pSE);RelErr$F4$Individualfactor<-append(RelErr$F4$Individualfactor,pME);RelErr$F4$RelativeError<-append(RelErr$F4$RelativeError,Err)}}
	}
	
	#Export measurement error report
	sink(paste(Output, "ErrorReport.txt", sep="_"), type="output")
	writeLines("ERROR REPORT \nLists per coefficient type (F1--F4) per harmonic: \n-Whether the session factor (influence of replicated outline extraction) introduced a significant error: 1---no systematic error, 0---systematic error, consider exclusion/re-evaluation \n-Whether the individuals showed larger between-individual variance than within-individual variance (error): 1---between-individual variance dominates, 0---within-individual variance dominates, consider exclusion/re-evaluation \n-The relative measurement error (in %) after Yezerinac et al. (1992) Syst. Biol. 41: 471-482 \n")
	print(RelErr)
	sink()
}


#--------------------------------------------

## EXAMPLES *************************************************************************************************

#Converting images
#ImageConversion("Sp", 1, 2, Scaling=100, ".jpg", ".ppm")

#Calculate replication error
#setwd("C:/R_TestData/Outline/Rep1")
#EFA1<-matrix(NA, 2, 42)
#colnames(EFA1)<-c(paste("A", 1:10, sep="."), paste("B", 1:10, sep="."), paste("C", 1:10, sep="."), paste("D", 1:10, sep="."), "ao", "co")
#for (i in 1:2) {
#	File<-paste("Sp", i, ".txt", sep="")
#	Dat<-read.table(File, header=TRUE, sep="\t")
#	ef<-NEF(Dat, Harmonics=10, start=TRUE)
#	EFA1[i,]<-c(ef$A, ef$B, ef$C, ef$D, ef$ao, ef$co)
#}
#setwd("C:/R_TestData/Outline/Rep2")
#EFA2<-matrix(NA, 2, 42)
#colnames(EFA2)<-c(paste("A", 1:10, sep="."), paste("B", 1:10, sep="."), paste("C", 1:10, sep="."), paste("D", 1:10, sep="."), "ao", "co")
#for (i in 1:2) {
#	File<-paste("Sp", i, ".txt", sep="")
#	Dat<-read.table(File, header=TRUE, sep="\t")
#	ef<-NEF(Dat, Harmonics=10, start=TRUE)
#	EFA2[i,]<-c(ef$A, ef$B, ef$C, ef$D, ef$ao, ef$co)
#}
#setwd("C:/R_TestData/Outline")
#OutlineAverage(Objects=c("EFA1", "EFA2"), Centroids=TRUE, "Sp1")

#Perform Zahn--Roskies Fourier analysis
#setwd("C:/R_TestData/Outline/Rep1")
#Dat<-read.table("Sp1.txt, header=TRUE, sep="\t")
#ZRF<-ZRfourier(Dat, 40)
#ZR.Recon<-iZRfourier(ZRF$ao, ZRF$an, ZRF$bn, Harmonics=40, Points=70, thetao=ZRF$thetao)
#ZR.Recon<-iZRfourier(ZRF$ao, ZRF$an, ZRF$bn, Harmonics=40, Points=100, thetao=ZRF$thetao)#this does not work, because the number of points must be the same as was originally extracted (70)
#plot(ZR.Recon$X, ZR.Recon$Y, type="l", asp=1)

#Perform elliptic Fourier analysis
#setwd("C:/R_TestData/Outline/Rep1")
#Dat<-read.table("Sp1.txt, header=TRUE, sep="\t")
#ef1<-efourier(Dat)
#ef2<-NEF(Dat)
#layout(matrix((1:9),3,3))
#par=(mar=c(2,2,2,2))
#for (i in 1:9) {
#	ief1<-iefourier(ef1$an, ef1$bn, ef1$cn, ef1$dn, i, 64, ef1$ao, ef1$co)
#	plot(Dat, type="l", asp=1, frame=FALSE, main=paste("Harmonics 0 to", i, sep=" "), col="grey")
#	polygon(Dat, col="grey", border=NA)
#	lines(ief1$x, ief1$y, type="l")
#}
#layout(matrix((1:9),3,3))
#par=(mar=c(2,2,2,2))
#for (i in 1:9) {
#	ief1<-iefourier(ef2$A, ef2$B, ef2$C, ef2$D, i, 64, ef2$ao, ef2$co)
#	plot(ief1$x, ief1$y, type="l", asp=1, frame=FALSE, main=paste("Harmonics 0 to", i, sep=" "), col="black", lwd=2)
#}

##***********************************************************************************************************

## FOOT #####################################################################################################

## VERSION HISTORY ******************************************************************************************
# Version 1.0
#	Date: XXX
#	Description of changes:
#		-Finished Program
#
# Version 1.1
#	Date: XXX
#	Description of changes:
#		-Added option in OutlineExtraction to export the raw (i.e. unsmoothed) outlines as well
#
# Version 1.2
#	Date: XXX
#	Description of changes:
#		-Fixed a problem that would occur in the success table of the outline extraction, if StartNum was not 1
#
# Version 1.3
#	Date: XXX
#	Description of changes:
#		-Included function to automatically average over replicates and calculate measurement errors
#
# Version 2.0
#	Date: XXX
#	Description of changes:
#		-Zahn--Roskies Fourier analysis functions included, fixed a variable naming issue in OutlineAverage
#
# Version 2.0.1
#	Date: XXX
#	Description of changes:
#		-Power function modified so that results are immediately returned
#
# Version 2.1
#	Date: XXX
#	Description of changes:
#		-Outline extraction modified, so that a NTS file is exported; Read.NTS function added
#		-NEF function modified, so that uniform orientation can also be achieved according to manually chosen baseline
#
# Version 2.1.1
#	Date: XXX
#	Description of changes:
#		-Removed Read.NTS and transferred into separate file MorphoFiles_Function.r
#
# Version 2.2
#	Date: 31 October 2016
#	Description of changes:
#		-Removed ImageConversion, Conte, EquiDist, smoothout, and OutlineExtraction to include in separate function
#			file MorphometricExtraction_Functions.r
##***********************************************************************************************************

#############################################################################################################

