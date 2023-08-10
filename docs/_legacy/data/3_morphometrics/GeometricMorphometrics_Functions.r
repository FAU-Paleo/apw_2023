## HEAD #####################################################################################################
#
# FUNCTION SET DESCRIPTION
#	Geometric Morphometrics in R
#
# DATASET DESCRIPTION
#	Morphometric datasets
#
# FURTHER READING
#	Claude, J. (2008) "Morphometrics with R". Gentleman, R., Hornik, K., and...
#		Parmigiani, G. (eds) "Use R!", vol. ii, 316 pp. (Springer).
#	Ezard, Th. H. G., Pearson, P. N., and Purvis, A. (2010) "Algorithmic approaches to aid...
#		species' delimitation in multidimensional morphospace" BMC Evolutionary Biology...
#		10: Article 175
#	Zelditch, M. L., Swiderski, D. L., and Sheets, H. D. (2012) "Geometric...
#		Morphometrics for Biologists"
#	
# METADATA
#	Author: Manuel F. G. Weinkauf
#	E-mail: weinkauf.scientific@gmail.com
#	R-version: 4.2.1
#	RStudio-version: 2022.07.1
#	Code-version: 1.11.2
#	Date of last update: 7 October 2019

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

## BODY #####################################################################################################

#Creation of test dataset
#t1<-rgamma(50, 1.5, .1)
#t2<-rnorm(50, 10, 2)
#t3<-rgamma(50, 2.8, .2)
#t4<-rnorm(50, 28, 6)
#t5<-rgamma(50, 20, .1)
#t6<-rnorm(50, 5, .4)
#tt<-matrix(c(t1, t2, t3, t4, t5, t6), ncol=6, byrow=FALSE)
#t1<-rgamma(50, 2.1, .1)
#t2<-rnorm(50, 18, 2)
#t3<-rgamma(50, 1.1, .2)
#t4<-rnorm(50, 45, 6)
#t5<-rgamma(50, 19, .1)
#t6<-rnorm(50, 12, .4)
#tt<-rbind(tt, matrix(c(t1, t2, t3, t4, t5, t6), ncol=6, byrow=FALSE))
#rownames(tt)<-c(rep("Spec1", 50), rep("Spec2", 50))

#**************************************************************************************
#Setting working dierctory
#setwd("C:/R_TestData/GeometricMorphometrics")

#########################################################################
# Mirror landmark configurations                                        #
# Necessary input variables:                                            #
#    Input: 2D landmark configuration to be mirrored.                   #
#           *matrix*                                                    #
#    Axis: Around which axis ("x"/"y") should the reflection be...      #
#          performed?                                                   #
#          *character*                                                  #
#    RevCoord: Revised order of coordinates, that arises from...        #
#              reflection. Should be a vector that gives the new...     #
#              position of coordinates 1:n in Input data.               #
#              NOTE: This only rearranges the rows of the reflected...  #
#              matrix, so that the landmarks in the rotated specimen... #
#              correspond to those in the unrotated specimen. Since...  #
#              rotation can have a large influence here, this must...   #
#              be worked out and checked by the user.                   #
# Output data: Mirrored landmark configuration.                         #
# Input dataset: Landmark configuration, landmarks in rows, x and y...	#
#                in columns.                                            #
#########################################################################

LM.Mirror<-function (Input, Axis=NULL, RevCoord=NULL) {
	#Check data consistency
	if (is.null(Axis)) {stop("Axis for reflection must be specified!")}
	if (Axis!="x" && Axis!="y") {stop("Axis must be either 'x' or 'y'!")}
	if (is.null(RevCoord)) {stop("RevCoord must be provided!")}
	if (length(RevCoord)!=nrow(Input)) {stop("RevCoord must be of same length as landmark configuration!")}
	
	#Read data into variable
	Data<-Input
	
	#Mirror data
	{if (Axis=="y") {Data[,1]<--Data[,1]}
	else {Data[,2]<--Data[,2]}}
	
	#Revise order of landmarks
	DataOut<-matrix(Data[RevCoord,], nrow(Data), ncol(Data))
	
	DataOut
}

#########################################################################
# Averaging data from several measurements and calculating...           #
#    measurement error                                                  #
# based on Claude (2008), pp. 63--67                                    #
# Necessary function files: MorphoFiles_Function.r                      #
# Necessary input variables:                                            #
#    Objects: List of morphometric data sets in shapes format. Must...  #
#             be provided as character vector containing the names of...#
#             the objects.                                              #
#             *character*                                               #
#    Meta: List of metadata (filenames and scales) as produced by...    #
#          Read.TPS. must be in the same as the Objects.                #
#          NOTE: If the second element (corresponding to Scale) is...   #
#          not null, it is assumed that the data were scaled upon...    #
#          import and will be scaled back before export!                #
#          *list* with two elements of type *vector*                    #
#    ErrorFile: Name of the file in which to write the results of the...#
#               measurement error calculations.                         #
#               *character*                                             #
#    Output: Name of TPS file with averaged landmarks to be exported.   #
#            *character*                                                #
# Output data: Morphometric data file in TPS format, and calculated...  #
#              replication errors.                                      #
# Input dataset: Morphometric data files in format as produced by...	#
#                Read.NTS, Read.TPS, Read.IMP or Read.PAST.             #
#########################################################################

LMAverage<-function (Objects, Meta, Error.File, Output) {
	#Evaluate parameters
	M<-length(Objects)
	N<-dim(eval(parse(text=Objects[1])))[3]
	LM<-dim(eval(parse(text=Objects[1])))[1]
	DM<-dim(eval(parse(text=Objects[1])))[2]
	
	#Convert complete landmark dataset into list
	Read.Obj<-function (X) {eval(parse(text=X))}
	Data<-lapply(Objects, Read.Obj)
	
	#Set up results array for mean landmark data
	LM.Mean<-array(NA, dim=c(LM, DM, N), dimnames=dimnames(eval(parse(text=Objects[1]))))
	Image.files<-vector(mode="character", length=N)
	Image.scales<-vector(mode="numeric", length=N)
	
	#Set up vectors for error calculations
	Values<-vector(mode="numeric", length=0)
	Session.factor<-vector(mode="numeric", length=0)
	Individual.factor<-vector(mode="numeric", length=0)
	
	#Calculate mean position for all landmarks
	for (i in 1:N) {
		Pos<-seq(from=i, to=M*N, by=N)
		if (!is.null(Meta[[1]])) {
			F<-Meta[[1]][Pos]
			{if (length(unique(F))!=1) {stop("Not all data from the same specimen, check Meta!")}
			else {Image.files[i]<-unique(F)}}
		}
		Image.scales[i]<-mean(Meta[[2]][Pos])
		Temp<-simplify2array(lapply(Data, "[", ,,i))
		LM.Mean[,,i]<-apply(Temp, c(1, 2), mean)
		Values<-append(Values, as.vector(Temp))
		Session.factor<-append(Session.factor, gl(dim(Temp)[3], dim(Temp)[1]*dim(Temp)[2]))
		Individual.factor<-append(Individual.factor, rep.int(seq(from=((DM*LM)*(i-1)+1), to=((DM*LM)*(i-1)+1)+((DM*LM)-1), by=1), M))
	}
	
	#Export TPS file
	{if (!is.null(Meta[[1]])) {FF<-Image.files} else {FF<-NULL}}
	{if (!is.null(Meta[[2]])) {SC<-Image.scales} else {SC<-NULL}}
	{if (is.null(SC)) {SCT<-FALSE} else {SCT<-TRUE}}
	Write.TPS(LM.Mean, Centroids=NULL, Filenames=FF, Scaling=SCT, Scale=SC, Output=Output)
	
	#Calculate errors
	Session.factor<-as.factor(Session.factor)
	Individual.factor<-as.factor(Individual.factor)
	Session.Error<-summary(aov(Values~Session.factor))
	Individual.Error<-summary(aov(Values~Individual.factor))
	
	s2within<-MSwithin<-Individual.Error[[1]][2,3]
	MSamong<-Individual.Error[[1]][1,3]
	s2among<-(MSamong-MSwithin)/M
	Rel.Error<-(s2within/(s2within+s2among))*100
	
	#Export measurement error report
	sink(Error.File, type="output")
	writeLines("---------- \nSESSION ERROR REPORT \nTests whether or not there is a strong discrepancy between replicated measurements \nShould in the ideal case be insignificant, with the Mean Sq of the Residuals (error variance) being larger than the Mean Sq of the Session.factor (systematic error of replications) \n")
	print(Session.Error)
	writeLines("\n\n---------- \nINDIVIDUAL ERROR REPORT \nTests whether or not there is a strong variance within individual measurements \nShould in the ideal case be significant, with Mean Sq of the Individual.factor (variance between individuals) being larger than Mean Sq of the Residuals (variance within individuals, i.e. measurement error) \n")
	print(Individual.Error)
	writeLines("\n\n---------- \nRELATIVE MEASUREMENT ERROR \nGives the percent measurement error as defined by Yezerinac et al. (1992) Syst. Biol. 41: 471-482 \n")
	print(paste(round(Rel.Error[1], digits=4), "%", sep=" "))
	sink()
}

#########################################################################
# Rotating configurations along first principal axis                    #
# based on Claude (2008), p. 165                                        #
# Necessary input variables:                                            #
#    Data: Array (shapes format) containing landmark data.              #
#          *array*                                                      #
# Output data: Morphometric data file in shapes format.                 #
# Input dataset: Morphometric data files in shapes format.              #
#########################################################################

aligne<-function(Data) {
	B<-Data
	n<-dim(Data)[3]
	k<-dim(Data)[2]
	for (i in 1:n) {
		Ms<-scale(Data[,,i], scale=FALSE)
		sv<-eigen(var(Ms))
		M<-Ms%*%sv$vectors
		B[,,i]<-M
	}
	B
}

#########################################################################
# Partial generalized Procrustes superimposition of landmark data       #
# based on Claude (2008), pp. 163f                                      #
# Necessary input variables:                                            #
#    Data: Array (shapes format) containing landmark data.              #
#          *array*                                                      #
# Output data: Partially Procrustes superimposed morphometric data in...#
#              shapes format.                                           #
# Input dataset: Morphometric data in shapes format.                    #
#########################################################################

pgPs<-function (Data) {
	#Read dimensions of dataset
	p<-dim(Data)[1]
	k<-dim(Data)[2]
	n<-dim(Data)[3]
	
	#Create array to store scaled and rotated configurations
	temp2<-temp1<-array(NA, dim=c(p, k, n))
	Siz<-numeric(n)

	#Translate and scale configurations
	##Function to translate set of landmarks
	trans1<-function(M) {scale(M, scale=FALSE)}
	
	##Function to calculate centroid size
	centsiz<-function (M) {
		p<-dim(M)[1]
		size<-sqrt(sum(apply(M, 2, var))*(p-1))
		list("centroid_size" = size,"scaled" = M/size)
	}
	
	#Perform alignment and scaling for all configurations
	for (i in 1:n) {
		Acs<-centsiz(Data[,,i])
		Siz[i]<-Acs[[1]]
		temp1[,,i]<-trans1(Acs[[2]])
	}

	#Define Procrustes distances quantity
	Qm1<-dist(t(matrix(temp1, k*p,n)))
	Q<-sum(Qm1)
	iter<-0
	
	#Loop minimization of Procrustes distances until convergence is reached
	##Function to calculate mean shape
	mshape<-function (A) {apply(A, c(1, 2), mean)}
	
	##Function to calculate interlandmark distances
	ild2<-function (M1, M2) {sqrt(apply((M1-M2)^2, 1, sum))}
	
	##Function to perform ordinary Procrustes superimposition
	pPsup<-function (M1, M2) {
		k<-ncol(M1)
		Z1<-trans1(centsiz(M1)[[2]])
		Z2<-trans1(centsiz(M2)[[2]])
		sv<-svd(t(Z2)%*%Z1)
		U<-sv$v
		V<-sv$u
		Delt<-sv$d
		sig<-sign(det(t(Z1)%*%Z2))
		Delt[k]<-sig*abs(Delt[k])
		V[,k]<-sig*V[,k]
		Gam<-U%*%t(V)
		beta<-sum(Delt)
		list(Mp1=Z1%*%Gam, Mp2=Z2, rotation=Gam, DP=sqrt(sum(ild2(Z1%*%Gam, Z2)^2)), rho=acos(beta))
	}
  
	##Induce loop
	while (abs(Q)>0.00001) {
		for (i in 1:n) {
			M<-mshape(temp1[,,-i])
			temp2[,,i]<-pPsup(temp1[,,i], M)[[1]]
		}
		Qm2<-dist(t(matrix(temp2, k*p, n)))
		Q<-sum(Qm1)-sum(Qm2)
		Qm1<-Qm2
		iter=iter+1
		temp1<-temp2
	}

	#Return list of results
	list("rotated"=temp2, "it.number"=iter, "Q"=Q, "intereucl.dist"=Qm2, "mshape"=centsiz(mshape(temp2))[[2]], "cent.size"=Siz)
}

#########################################################################
# Generalized resistant-fit superimposition of landmark data            #
# based on Claude (2008), pp. 178f                                      #
# Necessary functions: pgPs                                             #
# Necessary input variables:                                            #
#    Data: Array (shapes format) containing landmark data.              #
#          *array*                                                      #
#    AGoal: Anticipated goal of precision in iterative alignment...     #
#           (higher value means lower precisison).                      #
#           *numeric (real)*                                            #
#           default=0.0005                                              #
#    IterMax: Maximum number of iterations after which to terminate...  #
#             without result.                                           #
#             *numeric (integer)*                                       #
#             default=1000000                                           #
# Output data: Generalized resistant-fit superimposed morphometric...	#
#              data in shapes format.                                   #
# Input dataset: Morphometric data in shapes format.                    #
#########################################################################

GRF<-function (Data, AGoal=0.0005, IterMax=1000000) {
	#Read dimensions of dataset
	p<-dim(Data)[1]
	k<-dim(Data)[2]
	n<-dim(Data)[3]

	#Fit all objects as partial Procrustes superimposition
	Data<-pgPs(Data)$rotated
	D<-B<-array(NA, dim=c(p, k, n))

	#Scale all configurations
	##Function to compute median size of configurations
	medsize<-function (M) {
		mat<-as.matrix(dist(M))
		median(apply(mat, 2, median, na.rm=TRUE))
	}
	
	##Perform scaling
	for (i in 1:n) {
		B[,,i]<-Data[,,i]/medsize(Data[,,i])
	}

	#Compute and scale consensus
	Y<-apply(B, 1:2, median)
	Y<-Y/medsize(Y)
	A0<-10
	iter<-1

	#Superimpose configurations on consensus
	##Function to compute angles between homologous vectors
	argallvec<-function (X1, X2) {
		p<-dim(X1)[1]
		m<-m<-matrix(NA, p, p)
		for (i in 1:p){
			for (j in 1:p) {
				m[i,j]<-Arg(complex(1,X2[i,1],X2[i,2])-complex(1,X2[j,1],X2[j,2]))-Arg(complex(1,X1[i,1],X1[i,2])-complex(1,X1[j,1],X1[j,2]))  
			}
		}
		((m+pi)%%(2*pi))-pi
	}
	
	##Perform superimposition
	while(A0>AGoal){
		for (i in 1:n){
			M<-as.matrix(dist(Y))/as.matrix(dist(B[,,i]))
			beta<-median(apply(M, 2, median, na.rm=TRUE))
			ARG<-argallvec(B[,,i], Y)
			phi<-median(apply(ARG, 2, median, na.rm=TRUE))
			Gam<-matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), 2, 2)
			alpha<-Y-beta*B[,,i]%*%Gam
			D[,,i]<-beta*B[,,i]%*%Gam+matrix(apply(alpha, 2, median), p, k, byrow=TRUE)
		}
		Yb<-apply(D, 1:2, median)
		A0<-median(sqrt(apply((Yb-Y)^2, 1, sum)))
		Y<-Yb<-Yb/medsize(Yb)
		B<-D
		iter<-iter+1
		if (iter>IterMax) {stop("Max. number of iterations reached without solution!")}
	}

	#Return results
	list("rotated"=D, "limit"=A0, "iteration"=iter, "medshape"=Yb)
}

#########################################################################
# Perform orthogonal projection of superimposed coordinates             #
# based on Claude (2008), p. 169                                        #
# Necessary input variables:                                            #
#    Dat: Dataset of superimposed landmarks in shapes format.           #
#         *array*                                                       #
# Output data: Orthogonal projection of landmarks data.                 #
# Input dataset: Superimposed morphometric data in shapes format.       #
#########################################################################

ORP<-function (Dat) {
	#Get dimensions of dataset
	p<-dim(Dat)[1]
	k<-dim(Dat)[2]
	n<-dim(Dat)[3]
	
	#Calculate orthogonal projections
	##Function to calculate centroid size
	centsiz<-function (M) {
		p<-dim(M)[1]
		size<-sqrt(sum(apply(M, 2, var))*(p-1))
		list("centroid_size" = size,"scaled" = M/size)
	}
	##Function to calculate mean shape
	mshape<-function (A) {apply(A, c(1,2), mean)}
	##Perform calculations
	Y1<-as.vector(centsiz(mshape(Dat))[[2]])
	oo<-as.matrix(rep(1, n))%*%Y1
	I<-diag(1, k*p)
	mat<-matrix(NA, n, k*p)
	for (i in 1:n) {mat[i,]<-as.vector(Dat[,,i])}
	Xp<-mat%*%(I-(Y1%*%t(Y1)))
	Xp1<-Xp+oo
	array(t(Xp1), dim=c(p, k, n))
}

#########################################################################
# Analyze allometry in morphometric data                                #
# based on Zelditch et al. (2016), chapter 5 and Handbook pp. 143ff     #
# Required packages: vegan                                              #
# Necessary input variables:                                            #
#    Dat: Procrustes fitted landmark data as returned by shapes package.#
#         *list*                                                        #
#    DistType: Type of morphometric distance calculated for uni-...     #
#              variate allometry analysis. Either Riemannian or...      #
#              full Procrustes distances.                               #
#              *character*, either of "Riemannian" or "Procrustes"      #
#              default="Riemannian"                                     #
#    SL: How many of the landmarks are semi-landmarks?                  #
#        *integer*                                                      #
#        default=0                                                      #
#    Permutation: Should a permutation test (PERMANOVA) for the multi...#
#                 variate allometry analysis be performed?              #
#                 *logical*                                             #
#                 default=TRUE                                          #
# Output data: List containing raw data and statistical results of...   #
#              for allomatry analyses.                                  #
# Input dataset: Normalized landmarks with additional data as created...#
#                by shapes package (e.g. procGPA).                      #
#########################################################################

#Load packages
require(vegan)

Allometry<-function (Dat, DistType="Riemannian", SL=0, Permutation=TRUE) {
	#Test data for consistency
	if (DistType!="Riemannian" & DistType!="Procrustes") {stop("Distance type must be either of 'Riemannian' or 'Procrustes'!")}
	SL<-round(SL, digits=0)
	if (dim(Dat$rotated)[1]<SL) {stop("Number of semi-landmarks cannot be larger than number of landmarks!")}
	
	#Prepare results object
	Res<-list()

	#Univariate allometry analysis: Correlation between shape and size
	#Find smallest specimen as basis
	B<-which(Dat$size==min(Dat$size))
	Base.Shape<-Dat$rotated[,,B]
	
	#Calculate distances between smalles and remaining individuals
	Dist<-vector(length=0, mode="numeric")
	{if (DistType=="Riemannian") {
		for (i in 1:(dim(Dat$rotated)[3])) {
			if (i!=B) {Dist<-append(Dist, riemdist(Dat$rotated[,,i], Base.Shape))}
		}}
	else if (DistType=="Procrustes") {
		for (i in 1:(dim(Dat$rotated)[3])) {
			if (i!=B) {Dist<-append(Dist, procdist(Dat$rotated[,,i], Base.Shape, type="full"))}
		}}
	}

	#Calculate linear regression fit between specimen size and difference to smallest individual
	fit<-lm(Dist~Dat$size[-B])
	Res$Univariate$Statistics<-summary(fit)
	
	#Plot results
	win.graph(20, 10, 10)
	layout(matrix(c(1, 2), 1, 2))
	plot(Dat$size[-B], Dist, main="Univariate allometry", xlab="Size", ylab="Shape distance")
	curve(coef(fit)[2]*x+coef(fit)[1], add=TRUE, lwd=2, col="blue")
	Res$Univariate$Values<-cbind(Dat$size[-B], Dist)
	colnames(Res$Univariate$Values)<-c("Size", "Shape.Distance")
	
	#Univariate allometry analysis
	#Prepare data
	Size<-log(Dat$size)
	Shape<-Dat$rotated
	n<-dim(Shape)[3]
	m<-dim(Shape)[2]
	k<-dim(Shape)[1]
	Shape<-t(apply(Shape, 3, t))
	
	#Calculate model
	Allo.Mod<-lm(Shape~Size)
	Mean.Shape<-apply(Shape, 2, mean)
	Mean.Shape<-rep(1, n)%*%t(Mean.Shape)
	num<-sum(apply((Allo.Mod$fitted.values-Mean.Shape)^2, 1, sum))/1
	den<-sum(apply((Shape-Allo.Mod$fitted.values)^2, 1, sum))/(n-1-1)
	FL<-k-SL
	Fs<-num/den
	P<-1-pf(Fs, 1*((FL*m)+SL)-4, (n-1-1)*((FL*m)+SL)-4)
	vexp<-sum(diag(var(Allo.Mod$fitted.values)))/sum(diag(var(Shape)))
	Res$Multivariate$Statistics<-c(Fs, P, vexp)
	names(Res$Multivariate$Statistics)<-c("F-ratio", "p.value", "Var.explained")
	
	#Permuation F-ratio test
	if (Permutation==TRUE) {
		Perm<-adonis(Shape~Size, method="euclidean")
		Res$Multivariate$Permuation<-Perm
	}
	
	#Visualize shape change
	fit.val<-Allo.Mod$fitted.values
	matr<-Dat$rotated[,,which(Size==(max(Size)))]
	matt<-Dat$rotated[,,which(Size==(min(Size)))]
	XLIM<-c(min(c(matr[,1], matt[,1])), max(c(matr[,1], matt[,1])))
	YLIM<-c(min(c(matr[,2], matt[,2])), max(c(matr[,2], matt[,2])))
	plot(matr[,1], matr[,2], asp=1, type="n", xlim=XLIM, xlab="", ylim=YLIM, ylab="", axes=FALSE) 
	arrows(matt[,1], matt[,2], matr[,1], matr[,2], length=0.1, lwd=2)
	points(matt[,1], matt[,2], pch=21, col="black", bg="red", cex=1.3)
	points(matr[,1], matr[,2], pch=23, col="black", bg="blue", cex=0.8)
	title(main="Ontogenetic shape change", sub="Multivariate allometry")
	box()
	legend("topright", pch=c(21, 23), col="black", pt.bg=c("red", "blue"), legend=c("Small", "Large"), cex=1.5)
	Res$Multivariate$Values<-array(data=c(matr, matt), dim=c(nrow(matt), ncol(matt), 2), dimnames=list(NULL, c("x", "y"), c("Max.Size", "Min.Size")))
	
	#Return results
	return(Res)
}

#########################################################################
# Plot mean shapes with deformation arrows between two groups           #
# based on Claude (2008), pp. 181f                                      #
# Necessary input variables:                                            #
#    Dat: Dataset of superimposed landmarks in shapes format.           #
#         *array*                                                       #
#    Set1: Indices of matrices of group set 1.                          #
#          *numeric (integer)*                                          #
#    Set2: Indices of matrices of group set 2.                          #
#          *numeric (integer)*                                          #
#    Ampl: Should changes be amplified? Either FALSE or a numeric...    #
#          (real) value giving degree of amplification.                 #
#          default=FALSE                                                #
#    Points: Point symbol for plotting of outlines of the two sets.     #
#            *vector* of length 2                                       #
#            default=c(3, 4)                                            #
#    Arrows: Arrow parameters length, angle, and line width (compare... #
#            help(arrows)).                                             #
#            *vector* of length 3                                       #
#            default=c(0.1, 30, 2)                                      #
#    Col: Colour of the points for the outlines of the two sets.        #
#         *vector* of length 2                                          #
#         default=c("black", "black")                                   #
#    Arrow.Col: Colour of the deplacement arrows.                       #
#               *character*                                             #
#               default="blue"                                          #
#    Legend.Text: Text to put in legend to mark describe both sets.     #
#                 *vector* of length 2                                  #
#                 default=c("Set1", "Set2")                             #
# Output data: Plot showing the shape change from Set 1 to Set 2.       #
# Input dataset: Morphometric data in shapes format.                    #
#########################################################################

ShapeDeform<-function (Dat, Set1, Set2, Ampl=FALSE, Points=c(3, 4), Arrows=c(0.1, 30, 2), Col=c("black", "black"), Arrow.Col="blue", Legend.Text=c("Set1", "Set2")) {
	#Test sets for consistency
	if (any(Set1 %in% Set2)==TRUE) {stop("Some elements belong to more than one set!")}

	#Subset data
	D1<-Dat[,,Set1]
	D2<-Dat[,,Set2]
	
	#Calculate mean shape of sets
	mshape<-function (A) {apply(A, c(1, 2), mean)}
	MS1<-mshape(D1)
	MS2<-mshape(D2)
	
	#Plot shape deformation
	plot(MS1[,1], MS1[,2], asp=1, xlab="", ylab="", axes=FALSE, pch=Points[1], col=Col[1])
	points(MS2[,1], MS2[,2], pch=Points[2], col=Col[2])
	text(MS1, pos=2, labels=1:nrow(MS1), col="grey50", font=1, cex=0.7)
	{if (Ampl==FALSE) {arrows(MS1[,1], MS1[,2], MS2[,1], MS2[,2], length=Arrows[1], angle=Arrows[2], lwd=Arrows[3], col=Arrow.Col)}
	else {
		AR<-MS2+(MS2-MS1)*Ampl
		arrows(MS1[,1], MS1[,2], AR[,1], AR[,2], length=Arrows[1], angle=Arrows[2], lwd=Arrows[3], col=Arrow.Col)
	}}
	legend("bottomright", legend=Legend.Text, bty="n", pch=c(Points[1], Points[2]), col=c(Col[1], Col[2]))
}

#########################################################################
# Plot thin-plate splines for two groups                                #
# based on Claude (2008), pp. 184-186                                   #
# Necessary input variables:                                            #
#    Dat: Dataset of superimposed landmarks in shapes format.           #
#         *array*                                                       #
#    Set1: Indices of matrices of group set 1.                          #
#          *numeric (integer)*                                          #
#    Set2: Indices of matrices of group set 2.                          #
#          *numeric (integer)*                                          #
#    n: Number of grid cell columns.                                    #
#       *numeric (integer)*                                             #
#       default=20                                                      #
#    Ampl: Should changes be amplified? Either FALSE or a numeric...    #
#          (real) value giving degree of amplification.                 #
#    Lines: Lines to be plotted to indicate form. A vector containing...#
#           the points to be joint in consecutive order.                #
#           *vector*                                                    #
#           default=NULL                                                #
#    Points: Point symbol for plotting of outlines of the two sets.     #
#            *vector* of length 2                                       #
#            default=c(3, 4)                                            #
#    Col: Colour of the points for the outlines of the two sets.        #
#         *vector* of length 2                                          #
#         default=c("red", "blue")                                      #
#    Scale: Scaling parameter for point symbols (cex).                  #
#           *real*                                                      #
#           default=1                                                   #
#    Line.width: With of connecting lines for landmarks.                #
#                *integer*                                              #
#                default=2                                              #
#    Legend.Text: Text to put in legend to mark describe both sets.     #
#                 *vector* of length 2                                  #
#                 default=c("Set1", "Set2")                             #
# Output data: Plot showing the shape change from Set 1 to Set 2.       #
# Input dataset: Morphometric data in shapes format.                    #
#########################################################################

tps<-function (Dat, Set1, Set2, n=20, Ampl=FALSE, Lines=NULL, Points=c(3, 4), Col=c("red", "blue"), Scale=1, Line.width=2, Legend.Text=c("Set1", "Set2")) {
	#Test sets for consistency
	if (any(Set1 %in% Set2)==TRUE) {stop("Some elements belong to more than one set!")}

	#Subset data
	D1<-Dat[,,Set1]
	D2<-Dat[,,Set2]
	
	#Calculate mean shape of sets
	mshape<-function (A) {apply(A, c(1, 2), mean)}
	MS1<-mshape(D1)
	MS2<-mshape(D2)
	if (Ampl!=FALSE) {MS2Old<-MS2; MS2<-MS2+(MS2-MS1)*Ampl}
	
	#Estimate grid size
	xm<-min(MS1[,1])
	ym<-min(MS1[,2])
	xM<-max(MS1[,1])
	yM<-max(MS1[,2])
	rX<-xM-xm
	rY<-yM-ym
	
	#Calculate coordinates of line intersections
	n<-round(n, digits=0)
	a<-seq(from=xm-1/5*rX, to=xM+1/5*rX, length=n)
	b<-seq(from=ym-1/5*rX, to=yM+1/5*rX, by=(xM-xm)*7/(5*(n-1)))
	m<-round(0.5+(n-1)*(2/5*rX+yM-ym)/(2/5*rX+xM-xm))
	M<-as.matrix(expand.grid(a,b))
	
	#Calculate interpolated coordinates
	tps2d<-function (M, matr, matt) {
		p<-dim(matr)[1]
		q<-dim(M)[1]
		n1<-p+3
		P<-matrix(NA, p, p)
		for (i in 1:p) {
			for (j in 1:p) {
				r2<-sum((matr[i,]-matr[j,])^2)
				P[i,j]<-r2*log(r2)
			}
		}
		P[which(is.na(P))]<-0
		Q<-cbind(1, matr)
		L<-rbind(cbind(P, Q), cbind(t(Q), matrix(0, 3, 3)))
		m2<-rbind(matt, matrix(0, 3, 2))
		coefx<-solve(L)%*%m2[,1]
		coefy<-solve(L)%*%m2[,2]
		fx<-function (matr, M, coef) {
			Xn<-numeric(q)
			for (i in 1:q) {
				Z<-apply((matr-matrix(M[i,],p,2,byrow=TRUE))^2, 1, sum)
				Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))
			}
			Xn
		}
		matg<-matrix(NA, q, 2)
		matg[,1]<-fx(matr, M, coefx)
		matg[,2]<-fx(matr, M, coefy)
		matg
	}
	
	ngrid<-tps2d(M, MS1, MS2)
	
	#Plot results
	plot(ngrid[,1], ngrid[,2], cex=0.2, asp=1, axes=FALSE, xlab="", ylab="")
	for (i in 1:m) {lines(ngrid[(1:n)+(i-1)*n,])}
	for (i in 1:n) {lines(ngrid[(1:m)*n-i+1,])}
	if (Ampl!=FALSE) {MS1<-MS1+(MS2Old-MS1)*Ampl}
	points(MS1[,1], MS1[,2], pch=Points[1], col=Col[1], cex=Scale)
	points(MS2[,1], MS2[,2], pch=Points[2], col=Col[2], cex=Scale)
	if (!is.null(Lines)) {
		lines(MS1[Lines,], col=Col[1], lwd=Line.width)
		lines(MS2[Lines,], col=Col[2], lwd=Line.width)
	}
	legend("bottomright", legend=Legend.Text, bty="n", pch=c(Points[1], Points[2]), col=c(Col[1], Col[2]))
}

#########################################################################
# Plot vectorized thin-plate splines for two groups                     #
# Required packages: sp                                                 #
# based on Claude (2008), pp. 187-189                                   #
# Necessary input variables:                                            #
#    Dat: Dataset of superimposed landmarks in shapes format.           #
#         *array*                                                       #
#    Set1: Indices of matrices of group set 1.                          #
#          *numeric (integer)*                                          #
#    Set2: Indices of matrices of group set 2.                          #
#          *numeric (integer)*                                          #
#    n: Number of grid cell columns.                                    #
#       *numeric (integer)*                                             #
#       default=20                                                      #
#    Ampl: Should changes be amplified? Either FALSE or a numeric...    #
#          (real)  value giving degree of amplification.                #
#    Lines: Lines to be plotted to indicate form. A vector containing...#
#           the points to be joint in consecutive order.                #
#           *vector*                                                    #
#           default=NULL                                                #
#    Dens: Number of points to be interpolated.                         #
#          *numeric (integer)*                                          #
#          default=1000                                                 #
#    Plot: What type of plot should be produced? Options are:           #
#          "Vector": Pure vector plot                                   #
#          "Contour": Contour plot                                      #
#          "Surface": Surface plus vector plot (vectors 1/10 of Dens)   #
#          *character*                                                  #
# Output data: Plot showing the shape change from Set 1 to Set 2.       #
# Input dataset: Morphometric data in shapes format.                    #
#########################################################################

#Load packages
require(sp)

Vector.tps<-function (Dat, Set1, Set2, Lines=NULL, Dens=1000, Plot=NULL) {
	#Test sets for consistency
	if (any(Set1 %in% Set2)==TRUE) {stop("Some elements belong to more than one set!")}
	if (is.null(Lines)) {stop("Outline polygone (Lines) required!")}
	if (Lines[1]!=Lines[length(Lines)]) {stop("Closed outline polygone (Lines) required!")}
	if (Plot!="Vector" && Plot!="Contour" && Plot!="Surface") {stop("Plot must be one of either 'Vector', 'Contour', or 'Surface'")}
	Dens<-round(Dens, digits=0)

	#Subset data
	D1<-Dat[,,Set1]
	D2<-Dat[,,Set2]
	
	#Calculate mean shape of sets
	mshape<-function (A) {apply(A, c(1, 2), mean)}
	MS1<-mshape(D1)
	MS2<-mshape(D2)
	
	#Calculate interpolated coordinates
	tps2d<-function (M, matr, matt) {
		p<-dim(matr)[1]
		q<-dim(M)[1]
		n1<-p+3
		P<-matrix(NA, p, p)
		for (i in 1:p) {
			for (j in 1:p) {
				r2<-sum((matr[i,]-matr[j,])^2)
				P[i,j]<-r2*log(r2)
			}
		}
		P[which(is.na(P))]<-0
		Q<-cbind(1, matr)
		L<-rbind(cbind(P, Q), cbind(t(Q), matrix(0, 3, 3)))
		m2<-rbind(matt, matrix(0, 3, 2))
		coefx<-solve(L)%*%m2[,1]
		coefy<-solve(L)%*%m2[,2]
		fx<-function (matr, M, coef) {
			Xn<-numeric(q)
			for (i in 1:q) {
				Z<-apply((matr-matrix(M[i,],p,2,byrow=TRUE))^2, 1, sum)
				Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))
			}
			Xn
		}
		matg<-matrix(NA, q, 2)
		matg[,1]<-fx(matr, M, coefx)
		matg[,2]<-fx(matr, M, coefy)
		matg
	}
	
	#Calculate deformation vectors
	VDef<-spsample(Polygon(MS1[Lines,]), Dens, type="regular")
	VRef<-VDef@coords
	VTar<-tps2d(VRef, MS1, MS2)
	def<-sqrt(apply((VTar-VRef)^2, 1, sum))
	xl<-length(unique(VRef[,1]))
	yl<-length(unique(VRef[,2]))
	im<-matrix(NA, xl, yl)
	xind<-(1:xl)[as.factor(rank(VRef[,1]))]
	yind<-(1:yl)[as.factor(rank(VRef[,2]))]
	for (i in 1:(length(xind))) {
		im[xind[i], yind[i]]<-def[i]
	}
	
	#Plot results
	{if (Plot=="Vector") {
		plot(MS1, asp=1, xlab="", ylab="", axes=FALSE, pch=3)
		points(MS2, pch=4, col="grey50")
		lines(MS1[Lines,], col="black", lwd=2)
		lines(MS2[Lines,], col="grey50", lwd=2, lty=2)
		for (i in 1:(nrow(VRef))) {
			arrows(VRef[i,1], VRef[i,2], VTar[i,1], VTar[i,2], length=0.05)
		}
	}
	else if (Plot=="Contour") {
		plot(MS1, asp=1, xlab="", ylab="", axes=FALSE, pch=3)
		lines(MS1[Lines,], col="black", lwd=2)
		contour(sort(unique(VRef[,1])), sort(unique(VRef[,2])), im, axes=FALSE, frame=FALSE, add=TRUE)
	}
	else {
		image(sort(unique(VRef[,1])), sort(unique(VRef[,2])), im, col=heat.colors(40), asp=1, xlab="", ylab="", axes=FALSE, frame=FALSE)
		VDef<-spsample(Polygon(MS1[Lines,]), Dens/10, type="regular")
		VRef<-VDef@coords
		VTar<-tps2d(VRef, MS1, MS2)
		for (i in 1:(nrow(VRef))) {
			arrows(VRef[i,1], VRef[i,2], VTar[i,1], VTar[i,2], length=0.05)
		}
	}
	}
}

#########################################################################
# Analyze Euclidean distance matrix (mEDMA) and Euclidean variance-...  #
#    covariance matrix (vEDMA)                                          #
# based on Claude (2008), pp. 190-197                                   #
#    Dat: Dataset of superimposed landmarks in shapes format.           #
#         *array*                                                       #
# Output data: Matrix with mean form (M) and interlandmark distances... #
#              (FM) in Euclidean space (mEDMA) or variance-covariance...#
#              matrix of landmarks (vEDMA).                             #
# Input dataset: Morphometric data in shapes format.                    #
#########################################################################

mEDMA<-function (Dat) {
	#Calculate Euclidean average of configurations set
	##Function to calculate form matrix in vectorized form
	fm<-function (M) {
		mat<-as.matrix(dist(M))
		mat[col(mat)<row(mat)]
	}
	
	#Function to calculate average of configurations
	EDMA<-function (A) {
		#Extract dataset dimensions
		n<-dim(A)[3]
		p<-dim(A)[1]
		k<-dim(A)[2]
		
		#Calculate Euclidean distances matrix
		E<-matrix(NA, n, p*(p-1)/2)
		for (i in 1:n) {
			E[i,]<-(fm(A[,,i]))^2
		}
		
		#Calculate single mean squared Euclidean distances
		Em<-apply(E, 2, mean)
		
		#Compute sample variances of squared interlandmark distances
		S<-(apply(t((t(E)-Em)^2), 2, sum))/n
		
		#Compute Omega
		{if (k==2) {omega<-(Em^2-S)^0.25}
		else if (k==3) {omega<-(Em^2-1.5*S)^0.25}
		else {stop("Data must be 2D or 3D only!")}}
		Om<-diag(0, p)
		Om[row(Om)>col(Om)]<-omega
		Om<-t(Om)
		Om[row(Om)>col(Om)]<-omega
		if (any(is.nan(Om))) {
			Om[which(is.nan(Om))]<-0
			warning("Landmark variance>landmark distance in some cases, results may not be fully reliable")
		}
		Om
	}
	
	#Function to perform multidimensional scaling
	MDS<-function (mat, k) {
		p<-dim(mat)[1]
		C1<-diag(p)-1/p*matrix(1,p,p)
		B<--0.5*C1%*%mat^2%*%C1
		eC<-eigen(B)
		eve<-eC$vectors
		eva<-eC$values
		MD<-matrix(NA, p, k)
		for (i in 1:k) {MD[,i]<-sqrt(eva[i])*eve[,i]}
		MD
	}
	
	#Calculate estimate of mean form matrix
	k<-dim(Dat)[2]
	Eu<-EDMA(Dat)
	M<-MDS(Eu, k)
	list("M"=M, "FM"=as.matrix(dist(M)))
}

vEDMA<-function (Dat) {
	#Get dimensions of dataset
	p<-dim(Dat)[1]
	k<-dim(Dat)[2]
	n<-dim(Dat)[3]
	Bs<-array(NA, dim=c(p,p,n))
	for (i in 1:n){
		Cc<-apply(Dat[,,i], 2, mean)
		Ac<-t(t(Dat[,,i])-Cc)
		Bs[,,i]<-Ac%*%t(Ac)
	}
	B<-apply(Bs, 1:2, mean)
	M<-mEDMA(Dat)$M
	Ek<-(B-M%*%t(M))/k
	Ek
}

#########################################################################
# Identify most influential landmarks in configuration                  #
# Required functions: mEDMA                                             #
# based on Claude (2008), pp. 191f                                      #
# Necessary input variables:                                            #
#    Dat: Dataset of superimposed landmarks in shapes format.           #
#         *array*                                                       #
#    Set1: Indices of matrices of group set 1.                          #
#          *numeric (integer)*                                          #
#    Set2: indices of matrices of group set 2.                          #
#          *numeric (integer)*                                          #
# Output data: Plot showing the distribution of shape deformation...	#
#              over all landmarks and values of sum of divergence for...#
#              all landmarks.                                           #
# Input dataset: Morphometric data in shapes format.                    #
#########################################################################

FDM<-function (Dat, Set1, Set2) {
	#Test sets for consistency
	if (any(Set1 %in% Set2)==TRUE) {stop("Some elements belong to more than one set!")}

	#Subset data
	D1<-Dat[,,Set1]
	D2<-Dat[,,Set2]
	
	#Calculate mean form
	EM1<-mEDMA(D1)$M
	EM2<-mEDMA(D2)$M
	
	#Calculate sum of divergence
	FD<-as.matrix(dist(EM1))/as.matrix(dist(EM2))
	FD1<-abs(FD-median(FD, na.rm=TRUE))
	rownames(FD)<-rownames(FD1)<-1:nrow(D1)
	
	#Plot results
	I<-cbind(rep(1:nrow(D1), nrow(D1)), as.numeric(FD))
	plot(I[-which(is.na(FD)),], ylim=c(0, max(FD, na.rm=TRUE)), cex=0, xlab="Landmark", ylab="Form difference", font.lab=2)
	text(I[-which(is.na(FD)),], label=gl(nrow(D1), (nrow(D1)-1)), cex=0.7)
	abline(h=median(FD, na.rm=TRUE), lty=3)
	
	#Give output
	print("Sum of divergences")
	print(sort(round(apply(FD1, 2, sum, na.rm=TRUE), digits=2), decreasing=TRUE))
}

#########################################################################
# PCA of shape changes                                                  #
# Required functions: pgPs, aligne, ORP                                 #
# Required packages: shapes                                             #
# based on Claude (2008), pp. 234ff                                     #
# Necessary input variables:                                            #
#    Dat: Dataset of raw landmarks in shapes format.                    #
#         *array*                                                       #
#    Sym: Optional vector for symbols coding groups in PCA scatterplot. #
#         *vector*                                                      #
#    Col: Optional vector for colours coding groups in PCA scatterplot. #
#         *vector*                                                      #
#    Lines: Lines to be plotted to indicate form. A vector containing...#
#           the points to be joint in consecutive order.                #
#           *vector*                                                    #
#           default=NULL                                                #
#    Type: Type of plot created. Possible options are:                  #
#          "Explorative": Plot PCA, explained variance, and extreme...  #
#                         forms.                                        #
#          "Uniform": Plot UCA with uniform deformation components.     #
#          "Nonaffine": Plot relative warps with non-affine deformation.#
#          *character*                                                  #
#          default=NULL                                                 #
#    n: Number of grid cell columns.                                    #
#       *numeric (integer)*                                             #
#       default=20                                                      #
#    Ref: Index of two landmarks used for alignment along x-axis for... #
#         uniform and nonaffine deformation. If null, the principal...  #
#         axis of the configuration will be used.                       #
#         *numeric (integer)*, length=2                                 #
#         default=NULL                                                  #
# Output data: Plot showing the distribution of shape deformation...	#
#              over all landmarks and values of sum of divergence for...#
#              all landmarks.                                           #
# Input dataset: Morphometric data in shapes format.                    #
#########################################################################

#Load packages
require(shapes)

LM.PCA<-function (Dat, Sym=NULL, Col=NULL, Lines=NULL, Type=NULL, n=20, Ref=NULL) {
	#Test data for consistency
	{if (is.null(Sym)) {Sym<-rep(3, dim(Dat)[3])}
	else if (length(Sym)!=dim(Dat)[3]) {stop("Sym vector must be of same length as number of specimens!")}}
	{if (is.null(Col)) {Col<-rep("cornflowerblue", dim(Dat)[3])}
	else if (length(Col)!=dim(Dat)[3]) {stop("Col vector must be of same length as number of specimens!")}}
	if (is.null(Lines) && Type=="Explorative") {stop("Lines must be provided if Type=='Explorative'")}
	if (Type!="Explorative" && Type!="Uniform" && Type!="Nonaffine") {stop("Type must be one of either 'Explorative', 'Uniform', or 'Nonaffine'")}
	n<-round(n, digits=0)
	if (!is.null(Ref)) {
		if (length(Ref)!=2 | Ref[1]==Ref[2]) {stop("'Ref' mus have exactly two different values!")}
		if (Ref[1]<1 | Ref[2]<2 | Ref[1]>dim(Dat)[1] | Ref[2]>dim(Dat)[1]) {stop("'Ref' must contain valid landmark indices!")}
	}

	#Perform orthogonal rectangular projection
	#gos<-ORP(pgPs(aligne(Dat))$rotated)
	gos<-ORP(procGPA(Dat)$rotated)
	m<-t(matrix(gos, (dim(gos)[1]*2), dim(gos)[3]))
	pcs<-prcomp(m)
	
	#Functions
	##Function to calculate mean shape
	mshape<-function (A) {apply(A, c(1, 2), mean)}
	##Function to translate set of landmarks
	trans1<-function(M) {scale(M, scale=FALSE)}
	##Function to calculate centroid size
	centsiz<-function (M) {
		p<-dim(M)[1]
		size<-sqrt(sum(apply(M, 2, var))*(p-1))
		list("centroid_size" = size,"scaled" = M/size)
	}
	##Function to calculate interlandmark distances
	ild2<-function (M1, M2) {sqrt(apply((M1-M2)^2, 1, sum))}
	##Function to perform ordinary Procrustes superimposition
	pPsup<-function (M1, M2) {
		k<-ncol(M1)
		Z1<-trans1(centsiz(M1)[[2]])
		Z2<-trans1(centsiz(M2)[[2]])
		sv<-svd(t(Z2)%*%Z1)
		U<-sv$v
		V<-sv$u
		Delt<-sv$d
		sig<-sign(det(t(Z1)%*%Z2))
		Delt[k]<-sig*abs(Delt[k])
		V[,k]<-sig*V[,k]
		Gam<-U%*%t(V)
		beta<-sum(Delt)
		list(Mp1=Z1%*%Gam, Mp2=Z2, rotation=Gam, DP=sqrt(sum(ild2(Z1%*%Gam, Z2)^2)), rho=acos(beta))
	}
	##Function to produce alignment
	procalign<-function (A, X) {
		pA<-pgPs(A)
		n<-dim(A)[3]
		k<-dim(A)[2]
		A<-pA$rotated
		msh<-pA$mshape
		A1<-A
		{if (is.null(X)) {
			sv<-eigen(var(msh))
			V<-sv$vectors
			rotmsh<-msh%*%V
		}
		else {
			Coord<-msh[X, ]
			Slope<-(Coord[2,2]-Coord[1,2])/(Coord[2,1]-Coord[1,1])
			Deg<-atan(Slope)
			V<-matrix(c(cos(Deg), -sin(Deg), sin(Deg), cos(Deg)), 2, 2, byrow=TRUE)
			rotmsh<-msh%*%V
		}}
		for (i in 1:n) {A1[,,i]<-pPsup(A[,,i], rotmsh)$Mp1}
		list("rotated"=ORP(A1), "meansh"=rotmsh)
	}
	##Function to calculate uniform deformation
	uniform2D<-function (A, X) {
		#Get dimensions of data
		n<-dim(A)[3]
		kp<-dim(A)[1]*dim(A)[2]
		temp<-procalign(A, X)
		msh<-temp$meansh
		proc<-temp$rotated
		X<-t(matrix(proc, kp, n))
		V<-X-rep(1, n)%*%t(as.vector(msh))
		alph<-sum(msh[,1]^2)
		gam<-sum(msh[,2]^2)
		U1<-c(sqrt(alph/gam)*msh[,2], sqrt(gam/alph)*msh[,1])
		U2<-c(-sqrt(gam/alph)*msh[,1], sqrt(alph/gam)*msh[,2])
		score<-V%*%cbind(U1, U2)
		list("scores"=score, "uniform"=cbind(U1, U2), "meanshape"=msh, "rotated"=proc)
	}
	##Function to calculate thin-plate splines
	tps2d<-function (M, matr, matt) {
		p<-dim(matr)[1]
		q<-dim(M)[1]
		n1<-p+3
		P<-matrix(NA, p, p)
		for (i in 1:p) {
			for (j in 1:p) {
				r2<-sum((matr[i,]-matr[j,])^2)
				P[i,j]<-r2*log(r2)
			}
		}
		P[which(is.na(P))]<-0
		Q<-cbind(1, matr)
		L<-rbind(cbind(P, Q), cbind(t(Q), matrix(0, 3, 3)))
		m2<-rbind(matt, matrix(0, 3, 2))
		coefx<-solve(L)%*%m2[,1]
		coefy<-solve(L)%*%m2[,2]
		fx<-function(matr, M, coef) {
			Xn<-numeric(q)
			for (i in 1:q) {
				Z<-apply((matr-matrix(M[i,], p, 2, byrow=TRUE))^2, 1, sum)
				Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))
			}
			Xn
		}
		matg<-matrix(NA, q, 2)
		matg[,1]<-fx(matr, M, coefx)
		matg[,2]<-fx(matr, M, coefy)
		matg
	}
	tps<-function (matr, matt, n) {
		xm<-min(matt[,1])
		ym<-min(matt[,2])
		xM<-max(matt[,1])
		yM<-max(matt[,2])
		rX<-xM-xm; rY<-yM-ym
		a<-seq(xm-1/5*rX, xM+1/5*rX, length=n)
		b<-seq(ym-1/5*rX, yM+1/5*rX,by=(xM-xm)*7/(5*(n-1)))
		m<-round(0.5+(n-1)*(2/5*rX+ yM-ym)/(2/5*rX+ xM-xm))
		M<-as.matrix(expand.grid(a,b))
		ngrid<-tps2d(M, matr, matt)
		plot(ngrid[,1], ngrid[,2], cex=0.2, asp=1, axes=FALSE, xlab="", ylab="")
		for (i in 1:m) {lines(ngrid[(1:n)+(i-1)*n,])}
		for (i in 1:n) {lines(ngrid[(1:m)*n-i+1,])}
	}
	
	#Create plots
	{if (Type=="Explorative") {
		mesh<-apply(gos, 2, mean)
		mesh<-as.vector(mshape(Dat))
		max1<-matrix(mesh+max(pcs$x[,1])*pcs$rotation[,1], dim(gos)[1], dim(gos)[2])
		min1<-matrix(mesh+min(pcs$x[,1])*pcs$rotation[,1], dim(gos)[1], dim(gos)[2])
		max2<-matrix(mesh+max(pcs$x[,2])*pcs$rotation[,2], dim(gos)[1], dim(gos)[2])
		min2<-matrix(mesh+min(pcs$x[,2])*pcs$rotation[,2], dim(gos)[1], dim(gos)[2])
		XLIM1<-c(min(min1[,1], max1[,1]), max(min1[,1], max1[,1]))
		YLIM1<-c(min(min1[,2], max1[,2]), max(min1[,2], max1[,2]))
		XLIM2<-c(min(min2[,1], max2[,1]), max(min2[,1], max2[,1]))
		YLIM2<-c(min(min2[,2], max2[,2]), max(min2[,2], max2[,2]))
		##Plot results
		layout(matrix(c(1, 2, 3, 4), 2, 2))
		plot(pcs$x[,1], pcs$x[,2], pch=Sym, col=Col, xlab="PC 1", ylab="PC 2")
		barplot(pcs$sdev^2/sum(pcs$sdev^2), ylab="% of variance")
		title(sub="PC rank", mgp=c(0, 0, 0))
		plot(min1[,1], min1[,2], axes=FALSE, frame=FALSE, asp=1, xlim=XLIM1, ylim=YLIM1, xlab="", ylab="", pch=3)
		points(max1[,1], max1[,2], pch=4)
		title(sub="PC1", mgp=c(-4, 0, 0))
		lines(max1[Lines,], lty=1)
		lines(min1[Lines,], lty=2)
		legend("bottomright", pch=c(4, 3), lty=c(1, 2), legend=c("Max", "Min"), bty="n")
		plot(min2[,1], min2[,2], axes=FALSE, frame=FALSE, asp=1, xlim=XLIM2, ylim=YLIM2, xlab="", ylab="", pch=3)
		points(max2[,1], max2[,2], pch=4)
		title(sub="PC2", mgp=c(-4, 0, 0))
		lines(max2[Lines,], lty=1)
		lines(min2[Lines,], lty=2)
		legend("bottomright", pch=c(4, 3), lty=c(1, 2), legend=c("Max", "Min"), bty="n")
		return(list("PCA"=pcs, "PC1.Min"=min1, "PC1.Max"=max1, "PC2.Min"=min2, "PC2.Max"=max2))
	}
	else if (Type=="Uniform") {
		gou<-uniform2D(Dat, Ref)
		msh<-gou$meanshape
		Un<-gou$uniform
		##Plot results
		win.graph(16, 10, 10)
		layout(matrix(c(1, 1, 2, 3, 1, 1, 4, 5), 2, 4, byrow=TRUE))
		par(mar=c(5, 4, 4, 2))
		plot(gou$scores[,1], gou$scores[,2], pch=Sym, col=Col, xlab="U 1", ylab="U 2", asp=1)
		par(mar=c(4, 1, 4, 1))
		tps(msh, matrix(as.vector(msh)+Un[,1]*max(gou$scores[,1]), dim(Dat)[1], dim(Dat)[2]), n)
		points(matrix(as.vector(msh)+Un[,1]*max(gou$scores[,1]), dim(Dat)[1], dim(Dat)[2]), pch=16, col="blue")
		if (!is.null(Lines)) {
			lines(matrix(as.vector(msh)+Un[,1]*max(gou$scores[,1]), dim(Dat)[1], dim(Dat)[2])[Lines,], col="blue", lwd=2)
		}
		title("U 1: right")
		tps(msh, matrix(as.vector(msh)+Un[,1]*min(gou$scores[,1]), dim(Dat)[1], dim(Dat)[2]), n)
		points(matrix(as.vector(msh)+Un[,1]*min(gou$scores[,1]), dim(Dat)[1], dim(Dat)[2]), pch=16, col="blue")
		if (!is.null(Lines)) {
			lines(matrix(as.vector(msh)+Un[,1]*min(gou$scores[,1]), dim(Dat)[1], dim(Dat)[2])[Lines,], col="blue", lwd=2)
		}
		title("U 1: left")
		tps(msh, matrix(as.vector(msh)+Un[,2]*max(gou$scores[,2]), dim(Dat)[1], dim(Dat)[2]), n)
		points(matrix(as.vector(msh)+Un[,2]*max(gou$scores[,2]), dim(Dat)[1], dim(Dat)[2]), pch=16, col="blue")
		if (!is.null(Lines)) {
			lines(matrix(as.vector(msh)+Un[,2]*max(gou$scores[,2]), dim(Dat)[1], dim(Dat)[2])[Lines,], col="blue", lwd=2)
		}
		title("U 2: top")
		tps(msh, matrix(as.vector(msh)+Un[,2]*min(gou$scores[,2]), dim(Dat)[1], dim(Dat)[2]), n)
		points(matrix(as.vector(msh)+Un[,2]*min(gou$scores[,2]), dim(Dat)[1], dim(Dat)[2]), pch=16, col="blue")
		if (!is.null(Lines)) {
			lines(matrix(as.vector(msh)+Un[,2]*min(gou$scores[,2]), dim(Dat)[1], dim(Dat)[2])[Lines,], col="blue", lwd=2)
		}
		title("U 2: bottom")
		return(gou)
	}
	else {
		gou<-uniform2D(Dat, Ref)
		msh<-gou$meanshape
		Un<-gou$uniform
		kp<-dim(Dat)[2]*dim(Dat)[1]
		N<-dim(Dat)[3]
		X<-t(matrix(gou$rotated, kp, N))
		V<-X-t(t(rep(1, N)))%*%as.vector(msh)
		Ben<-diag(1, kp)-Un%*%solve(t(Un)%*%Un)%*%t(Un)
		LSR<-svd(V%*%Ben)
		##Calculate scores of nonaffine components
		score<-LSR$u%*%diag(LSR$d)
		NonUnif<-LSR$v
		##Plot results
		win.graph(16, 10, 10)
		layout(matrix(c(1, 1, 2, 3, 1, 1, 4, 5), 2, 4, byrow=TRUE))
		par(mar=c(5, 4, 4, 2))
		plot(score[,1], score[,2], pch=Sym, col=Col, xlab="RW 1", ylab="RW 2", asp=1)
		par(mar=c(4, 1, 4, 1))
		tps(msh, matrix(as.vector(msh)+NonUnif[,1]*max(score[,1]), dim(Dat)[1], dim(Dat)[2]), n)
		points(matrix(as.vector(msh)+NonUnif[,1]*max(score[,1]), dim(Dat)[1], dim(Dat)[2]), pch=16, col="blue")
		if (!is.null(Lines)) {
			lines(matrix(as.vector(msh)+NonUnif[,1]*max(score[,1]), dim(Dat)[1], dim(Dat)[2])[Lines,], col="blue", lwd=2)
		}
		title("RW 1: right")
		tps(msh, matrix(as.vector(msh)+NonUnif[,1]*min(score[,1]), dim(Dat)[1], dim(Dat)[2]), n)
		points(matrix(as.vector(msh)+NonUnif[,1]*min(score[,1]), dim(Dat)[1], dim(Dat)[2]), pch=16, col="blue")
		if (!is.null(Lines)) {
			lines(matrix(as.vector(msh)+NonUnif[,1]*min(score[,1]), dim(Dat)[1], dim(Dat)[2])[Lines,], col="blue", lwd=2)
		}
		title("RW 1: left")
		tps(msh, matrix(as.vector(msh)+NonUnif[,2]*max(score[,2]), dim(Dat)[1], dim(Dat)[2]), n)
		points(matrix(as.vector(msh)+NonUnif[,2]*max(score[,2]), dim(Dat)[1], dim(Dat)[2]), pch=16, col="blue")
		if (!is.null(Lines)) {
			lines(matrix(as.vector(msh)+NonUnif[,2]*max(score[,2]), dim(Dat)[1], dim(Dat)[2])[Lines,], col="blue", lwd=2)
		}
		title("RW 2: top")
		tps(msh, matrix(as.vector(msh)+NonUnif[,2]*min(score[,2]), dim(Dat)[1], dim(Dat)[2]), n)
		points(matrix(as.vector(msh)+NonUnif[,2]*min(score[,2]), dim(Dat)[1], dim(Dat)[2]), pch=16, col="blue")
		if (!is.null(Lines)) {
			lines(matrix(as.vector(msh)+NonUnif[,2]*min(score[,2]), dim(Dat)[1], dim(Dat)[2])[Lines,], col="blue", lwd=2)
		}
		title("RW 2: bottom")
		return(list("scores"=score, "Mean.Shape"=msh, "NonUniform.Deform"=NonUnif))
	}}
}

#########################################################################
# Thin-plate splines of PCA results                                     #
# based on Claude (2008), pp. 234ff                                     #
# Necessary input variables:                                            #
#    Mean.Shape: Dataset of mean shape landmarks.                       #
#                *matrix*                                               #
#    Deform: Deformation vectors along PCA axes.                        #
#            *matrix*                                                   #
#    Scores: PCA scores.                                                #
#            *matrix*                                                   #
#    n: Number of grid cell columns.                                    #
#       *numeric (integer)*                                             #
#       default=20                                                      #
#    Lines: Lines to be plotted to indicate form. A vector containing...#
#           the points to be joint in consecutive order.                #
#           *vector*                                                    #
#           default=NULL                                                #
#    Axis: Deformation at which end of the PCA is to plot?              #
#          Either of "PC1min", "PC1max", "PC2min", or "PC2max"          #
#          *character*                                                  #
#          default=NULL                                                 #
#    col: Colour of plot points.                                        #
#         *character*                                                   #
#         default="blue"                                                #
#    pch: Plot symbols.                                                 #
#         *integer*                                                     #
#         default=16                                                    #
#    Line.col: Colour of connecting lines.                              #
#              *character*                                              #
#              default="blue"                                           #
#    Line.width: Width of connecting lines.                             #
#                *integer*                                              #
#                default=2                                              #
#    Title: Optional title of the plot.                                 #
#           *character*                                                 #
#           default=NULL                                                #
# Output data: Thin-plate splines of deformation along PC axes.         #
# Input dataset: Results from LM.PCA.                                   #
#########################################################################

LM.PCA.tps<-function (Mean.Shape, Deform, Scores, n=20, Lines=NULL, Axis=NULL, col="blue", pch=16, Line.col="blue", Line.width=2, Title=NULL) {
	#Test data consistency
	if (is.null(Axis)) {stop("Axis must be either of 'PC1min', 'PC1max', 'PC2min', or 'PC2max'!")}
	if (Axis!="PC1min" & Axis!="PC1max" & Axis!="PC2min" & Axis!="PC2max") {stop("Axis must be either of 'PC1min', 'PC1max', 'PC2min', or 'PC2max'!")}
	
	#Function to calculate thin-plate splines
	tps2d<-function (M, matr, matt) {
		p<-dim(matr)[1]
		q<-dim(M)[1]
		n1<-p+3
		P<-matrix(NA, p, p)
		for (i in 1:p) {
			for (j in 1:p) {
				r2<-sum((matr[i,]-matr[j,])^2)
				P[i,j]<-r2*log(r2)
			}
		}
		P[which(is.na(P))]<-0
		Q<-cbind(1, matr)
		L<-rbind(cbind(P, Q), cbind(t(Q), matrix(0, 3, 3)))
		m2<-rbind(matt, matrix(0, 3, 2))
		coefx<-solve(L)%*%m2[,1]
		coefy<-solve(L)%*%m2[,2]
		fx<-function(matr, M, coef) {
			Xn<-numeric(q)
			for (i in 1:q) {
				Z<-apply((matr-matrix(M[i,], p, 2, byrow=TRUE))^2, 1, sum)
				Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))
			}
			Xn
		}
		matg<-matrix(NA, q, 2)
		matg[,1]<-fx(matr, M, coefx)
		matg[,2]<-fx(matr, M, coefy)
		matg
	}
	tps<-function (matr, matt, n) {
		xm<-min(matt[,1])
		ym<-min(matt[,2])
		xM<-max(matt[,1])
		yM<-max(matt[,2])
		rX<-xM-xm; rY<-yM-ym
		a<-seq(xm-1/5*rX, xM+1/5*rX, length=n)
		b<-seq(ym-1/5*rX, yM+1/5*rX,by=(xM-xm)*7/(5*(n-1)))
		m<-round(0.5+(n-1)*(2/5*rX+ yM-ym)/(2/5*rX+ xM-xm))
		M<-as.matrix(expand.grid(a,b))
		ngrid<-tps2d(M, matr, matt)
		plot(ngrid[,1], ngrid[,2], cex=0.2, asp=1, axes=FALSE, xlab="", ylab="")
		for (i in 1:m) {lines(ngrid[(1:n)+(i-1)*n,])}
		for (i in 1:n) {lines(ngrid[(1:m)*n-i+1,])}
	}
	
	#Calculate and plot thin-plate splines
	{if (Axis=="PC1min") {
		matt<-matrix(as.vector(Mean.Shape)+Deform[,1]*min(Scores[,1]), nrow(Mean.Shape), ncol(Mean.Shape))
	}
	else if (Axis=="PC1max") {
		matt<-matrix(as.vector(Mean.Shape)+Deform[,1]*max(Scores[,1]), nrow(Mean.Shape), ncol(Mean.Shape))
	}
	else if (Axis=="PC2min") {
		matt<-matrix(as.vector(Mean.Shape)+Deform[,2]*min(Scores[,2]), nrow(Mean.Shape), ncol(Mean.Shape))
	}
	else {
		matt<-matrix(as.vector(Mean.Shape)+Deform[,2]*max(Scores[,2]), nrow(Mean.Shape), ncol(Mean.Shape))
	}
	}
	tps(Mean.Shape, matt, n)
	points(matt[,1], matt[,2], pch=pch, col=col)
	if (!is.null(Lines)) {
		lines(matt[Lines,], col=Line.col, lwd=Line.width)
	}
	if (!is.null(title)) {title(Title)}
}

#########################################################################
# Hotelling-Lawley trace statistic of significance of difference...     #
#    between predefined groups                                          #
# Required functions: ORP                                               #
# Required packages: MASS                                               #
# based on Claude (2008), pp. 252ff                                     #
# Necessary input variables:                                            #
#    Dat: Dataset of raw landmarks in shapes format.                    #
#         *array*                                                       #
#    Groups: Coding of groups in Dat.                                   #
#            *vector*                                                   #
#    exact: calculate F-approximation following second moment estimate? #
#           *logical*                                                   #
#           TRUE: Use equation by McKeon (1974).                        #
#           FALSE: Use simple estimate.                                 #
#           default=FALSE                                               #
#    Lines: Lines to be plotted to indicate form, A vector containing...#
#           the points to be joint in consecutive order.                #
#           *vector*                                                    #
#           default=NULL                                                #
# Output data: Statistics for significance of difference between groups.#
# Input dataset: Morphometric data in shapes format.                    #
#########################################################################

#Load packages
require(MASS)
require(grDevices)

Hotellingsp<-function (Dat, Groups, exact=FALSE, Lines=NULL) {
	#Check data consistency
	if (length(Groups)!=dim(Dat)[3]) {stop("Groups must be of same length as Dat!")}
	Groups<-as.factor(Groups)
	
	#Decompose data into Euclidean space
	Dat.Eucl<-ORP(pgPs(aligne(Dat))$rotated)
	n<-dim(Dat.Eucl)[3]
	Eucl.Coord<-t(matrix(Dat.Eucl, (dim(Dat.Eucl)[1])*2, n))
	
	#Calculate linear model of dataset
	fit<-lm(Eucl.Coord~Groups)
	dfef<-length(levels(Groups))-1
	dfer<-n-length(levels(Groups))
	SSef<-(n-1)*var(fit$fitted.values)
	SSer<-(n-1)*var(fit$residuals)
	
	#Calculate shape space dimensions
	p<-qr(SSef+SSer)$rank
	
	#Calculate Hotelling's T2
	s<-min(dfef, p)
	m<-(dfer-p-1)/2
	t1<-(abs(p-dfef)-1)/2
	Ht<-sum(diag(SSef%*%ginv(SSer)))
	
	#Calculate F-statistics
	Fapprox<-Ht*(2*(s*m+1))/(s^2*(2*t1+s+1))
	ddfnum<-s*(2*t1+s+1)
	ddfden<-2*(s*m+1)
	pval<-1-pf(Fapprox, ddfnum, ddfden)
	
	#Calculate second-moment estimate for F-approximation
	if (exact==TRUE) {
		b<-(p+2*m)*(dfef+2*m)/((2*m+1)*(2*m-2))
		c1<-(2+(p*dfef+2)/(b-1))/(2*m)
		Fapprox<-((4+(p*dfef+2)/(b-1))/(p*dfef))*(Ht/c1)
		ddfnum<-p*dfef
		ddfden<-4+(p*dfef+2)/(b-1)
	}
	
	#Plot CVA results
	##Plot scores
	win.graph(16, 10, 10)
	layout(matrix(c(1, 1, 1, 1, 2, 3), 2, 3))
	CVA<-lda(Eucl.Coord, Groups)
	CVA2<-lda(Eucl.Coord, Groups, CV=TRUE)
	LDP<-predict(CVA)
	Class.Corr<-round(((sum(Groups==LDP$class)/(length(Groups)))*100), digits=3)
	Class.Corr2<-round(((sum(Groups==CVA2$class)/(length(Groups)))*100), digits=3)
	Score<-predict(CVA)$x
	{if (ncol(Score)==1) {
		YLIM<-c(0, max(c(max(hist(Score[Groups==levels(Groups)[1],1], plot=FALSE)$density), max(hist(Score[Groups==levels(Groups)[2],1], plot=FALSE)$density))))
		hist(Score[Groups==levels(Groups)[1],1], xlim=c(min(Score[,1]), max(Score[,1])), ylim=YLIM, xlab="LD 1", main="Discriminant scores", freq=FALSE, col=rgb(1, 0, 0, 0.5))
		hist(Score[Groups==levels(Groups)[2],1], xlim=c(min(Score[,1]), max(Score[,1])), freq=FALSE, col=rgb(0, 1, 0, 0.5), add=TRUE)
		legend("topright", legend=c(paste("Group", levels(Groups)[1], sep=" "), paste("Group", levels(Groups)[2], sep=" ")), fill=c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5)), bg="white", cex=0.7)
	}
	else {
		palette(rainbow(length(levels(Groups))))
		SymSet<-seq(1:length(levels(Groups)))
		Sym<-vector(mode="numeric", length=length(Groups))
		for (i in 1:length(Groups)) {
			Sym[i]<-SymSet[match(Groups[i], levels(Groups[i]))]
		}
		plot(Score[,1], Score[,2], pch=Sym, col=Sym, asp=1, xlab="LD 1", ylab="LD 2", main="Discriminant scores")
		for (i in 1:length(SymSet)) {
			Hull.Group<-which(Sym==SymSet[i])
			Score.Sub<-Score[Hull.Group,1:2]
			Conv.Hull<-chull(Score.Sub[,])
			polygon(Score.Sub[Conv.Hull,], border=palette()[i])
		}
		legend("topright", legend=paste("Group", levels(Groups)[1:length(levels(Groups))], sep=" "), pch=SymSet, col=palette(), bg="white")
	}
	}
	##Calculate average shapes of groups
	LD<-CVA$scaling
	msh<-apply(Eucl.Coord, 2, mean)
	n<-dim(Eucl.Coord)[1]
	dfw<-n-length(levels(Groups))
	SSw<-var(fit$residuals)*(n-1)
	VCVw<-SSw/dfw
	LDs<-VCVw%*%LD
	LD1M<-(matrix(msh+max(Score[,1])*LDs[,1], dim(Dat)[1], dim(Dat)[2]))
	LD1m<-(matrix(msh+min(Score[,1])*LDs[,1], dim(Dat)[1], dim(Dat)[2]))
	if (ncol(Score)>1) {
		LD2M<-(matrix(msh+max(Score[,2])*LDs[,2], dim(Dat)[1], dim(Dat)[2]))
		LD2m<-(matrix(msh+min(Score[,2])*LDs[,2], dim(Dat)[1], dim(Dat)[2]))
	}
	##Plot mean shapes
	plot(LD1M[,1], LD1M[,2], axes=FALSE, frame=FALSE, asp=1, xlab="", ylab="", pch=4, main="LD 1")
	points(LD1m[,1], LD1m[,2], pch=3)
	if (!is.null(Lines)) {lines(LD1M[Lines,], lty=1); lines(LD1m[Lines,], lty=2)}
	legend("bottomright", pch=c(4, 3), lty=c(1, 2), legend=c("Max", "Min"), bty="n")
	if (ncol(Score)>1) {
		plot(LD2M[,1], LD2M[,2], axes=FALSE, frame=FALSE, asp=1, xlab="", ylab="", pch=4, main="LD 2")
		points(LD2m[,1], LD2m[,2], pch=3)
		if (!is.null(Lines)) {lines(LD2M[Lines,], lty=1); lines(LD2m[Lines,], lty=2)}
		legend("bottomright", pch=c(4, 3), lty=c(1, 2), legend=c("Max", "Min"), bty="n")
	}
	legend("bottomright", pch=c(4, 3), lty=c(1, 2), legend=c("Max", "Min"), bty="n")
	
	#Return results
	return(list(unlist(list("df.effect"=dfef, "df.error"=dfer, "T2"=Ht, "Approx.F"=Fapprox, "df.1"=ddfnum, "df.2"=ddfden, "p-value"=pval, "Correct classification (%)"=Class.Corr)), "CVA.Model"=CVA))
}

#########################################################################
# Clustering of morphometric data                                       #
# Required packages: MASS, shapes, ape, cluster                         #
# based on Claude (2008), pp. 252ff                                     #
# Necessary input variables:                                            #
#    Dat: Dataset of superimposed landmarks in shapes format.           #
#         *array*                                                       #
#    Groups: Coding of groups in Dat.                                   #
#            *vector*                                                   #
#    K: Number of clusters for kmeans clustering.                       #
#       *numeric (integer)*                                             #
#       default=3                                                       #
#    Max.Groups: Maximum number of groups for elbow plot.               #
#                *integer*                                              #
#                default=round(dim(Dat)[3]/3, digits=0)                 #
# Output data: Plots showing UPGMA and complete clustering, elbow-...   #
#              plot, and partitional clustering plot alongside...       #
#              assiciated results values.                               #
# Input dataset: Morphometric data in shapes format.                    #
#########################################################################

#Load packages
require(MASS)
require(shapes)
require(ape)
require(cluster)

MorphoCluster<-function (Dat, Groups, K=3, Max.Groups=round(dim(Dat)[3]/3, digits=0)) {
	#Test data consistency
	K<-round(K, digits=0)
	if (length(Groups)!=dim(Dat)[3]){stop("Vector Groups must be of same length as Dat")}
	
	#Transform groups
	Groups<-as.factor(Groups)
	
	#Decompose data into Euclidean space
	Dat.Eucl<-ORP(Dat)
	Eucl.Coord<-t(matrix(Dat.Eucl, (dim(Dat.Eucl)[1])*2, dim(Dat.Eucl)[3]))
	
	#Set up results for export
	Res<-list()
	
	#Plot results
	win.graph(16, 14, 10)
	par(mar=c(4, 4, 3, 0), oma=c(0, 0, 0, 0))
	layout(matrix(1:4, 2, 2))
	##Define colour-coding
	ColSet<-rainbow(length(levels(Groups)))
	Col<-ColSet[match(Groups, levels(Groups))]
	##Prepare and plot clusters
	Dist.Mat<-dist(Eucl.Coord, method="euclidean")
	Clust<-hclust(Dist.Mat, method="average")
	Clust<-as.phylo(Clust)
	Clust$tip.label<-paste(as.character(Groups), Clust$tip.label, sep="-")
	plot(Clust, type="fan", tip.color=Col, main="UPGMA")
	Res$UPGMA$Dist<-Dist.Mat
	Res$UPGMA$Cluster<-Clust
	Clust<-hclust(Dist.Mat, method="complete")
	Clust<-as.phylo(Clust)
	Clust$tip.label<-paste(as.character(Groups), Clust$tip.label, sep="-")
	plot(Clust, type="fan", tip.color=Col, main="Complete")
	Res$Complete$Dist<-Dist.Mat
	Res$Complete$Cluster<-Clust
	
	#Perform partitional clustering
	d.f<-dim(Eucl.Coord)[1]-1
	SStot<-sum(diag(var(Eucl.Coord)))*d.f
	expl<-0
	for (i in 2:Max.Groups) {
		mod<-pam(dist(Eucl.Coord), k=i)
		mod1<-lm(Eucl.Coord~as.factor(mod$clustering))
		expl[i]<-sum(diag(var(mod1$fitted.values)))*d.f/SStot
	}
	plot(1:Max.Groups, expl, xlab="Number of clusters", ylab="Explained variance (fraction)", pch=16, type="b")
	Res$Elbow<-matrix(c(1:Max.Groups, expl), length(expl), 2)
	colnames(Res$Elbow)<-c("Groups", "Expl.variance")
	KMeans<-pam(dist(Eucl.Coord), k=K, keep.diss=TRUE)
	Class<-length(which(KMeans$clustering==as.numeric(Groups)))/length(Groups)
	plot(KMeans, which.plot=1, shade=TRUE, col.p=Col, col.txt=Col, col.clus="grey90", main=paste("Correct classification: ", round(Class*100, digits=2), "%", sep=""))
	Res$Partitional$Cluster<-KMeans
	Res$Partitional$Correct<-Class
	
	#Export results
	return(Res)
}

#########################################################################
# Two-block partial least squares for correlation between shape...      #
#    and other variables                                                #
# For significance test compare Zelditch et al. (2012), pp. 169ff       #
# Required packages: pls                                                #
# Necessary input variables:                                            #
#    Shape.Data: Dataset of superimposed landmarks in shapes format.    #
#                *array*                                                #
#    Parameter.Matrix: Matrix containing other set of parameters to...  #
#                      calculate regression. Variables in columns,...   #
#                      one row with set of variables per specimen....   #
#                      The first column is expected to contain the...   #
#                      sample-coding.                                   #
#                      *matrix*                                         #
# Output data: Test statistics and model for partial least squares...   #
#              regression between the two datasets.                     #
# Input dataset: Morphometric data in shapes format and other set of... #
#                parameters to correlate with. Parameter.Matrix can...  #
#                be another morphometric dataset, but also any other... #
#                data (e.g. environmental data), but must be numeric.   #
#########################################################################

#Load packages
require(pls)

PLS.Shape<-function (Shape.Data, Parameter.Matrix) {
	#Transpose shape data to matrix of form x1, y1, x2, y2,...
	Shape.Data<-t(apply(Shape.Data, 3, t))
	
	#Prepare parameter matrix and sample coding
	NA.rows<-complete.cases(Parameter.Matrix[,-1])
	Samples<-Parameter.Matrix[NA.rows,1]
	{if (is.null(colnames(Parameter.Matrix))) {Params<-paste("Param.", 1:(ncol(Parameter.Matrix)-1))}
	else {Params<-colnames(Parameter.Matrix)[-1]}}
	Parameter.Matrix<-as.matrix(Parameter.Matrix[,-1])
	colnames(Parameter.Matrix)<-Params
	Param.N<-ncol(Parameter.Matrix)
	
	#Test data for consistency
	if (nrow(Shape.Data)!=nrow(Parameter.Matrix)) {stop("Shape and parameter matrix do not have the same length!")}
	
	#Define function to calculate trace of square matrix
	Trace<-function (Mat) {
		Vals<-vector(mode="numeric", length=ncol(Mat))
		for (i in 1:(ncol(Mat))) {
			Vals[i]<-Mat[i,i]
		}
		Tr<-sum(Vals)
		return(Tr)
	}
	
	#Combine data into data frame
	PLS.Data<-data.frame(matrix(integer(0), nrow=nrow(Shape.Data)))
	PLS.Data$Shape<-Shape.Data
	PLS.Data<-cbind(PLS.Data, Parameter.Matrix)
	
	#Calculate significance of correlation between blocks
	R1<-var(Shape.Data, na.rm=TRUE)
	R2<-var(Parameter.Matrix, na.rm=TRUE)
	R12<-cov(Shape.Data, Parameter.Matrix, use="complete.obs")
	##Calculate Escoufiers coefficient
	RV<-Trace(R12%*%t(R12))/(sqrt(Trace(R1%*%t(R1))*Trace(R2%*%t(R2))))
	##Perform randomization test
	RV.Rand<-vector(mode="numeric", length=1000)
	for (i in 1:1000) {
		Rand.Vec<-sample(nrow(Shape.Data), replace=FALSE)
		Rand.Shape<-Shape.Data[Rand.Vec,]
		R1<-var(Rand.Shape, na.rm=TRUE)
		R12<-cov(Rand.Shape, Parameter.Matrix, use="complete.obs")
		RV.Rand[i]<-Trace(R12%*%t(R12))/(sqrt(Trace(R1%*%t(R1))*Trace(R2%*%t(R2))))
	}
	Large<-length(which(abs(RV.Rand)>=abs(RV)))
	Small<-length(which(abs(RV.Rand)<abs(RV)))
	P<-Large/(Large+Small)
	
	#Calculate partial least squares
	Indep<-paste(colnames(Parameter.Matrix), collapse="+")
	PLS<-plsr(as.formula(paste("Shape ~ ", Indep)), data=PLS.Data, na.action=na.exclude, validation="LOO", jackknife=TRUE)
	
	#Return results
	return(list("Significance"=P, "Number.Parameters"=Param.N, "Parameter.Names"=Params, "Sample.Names"=Samples, "Model"=PLS))
}

#########################################################################
# Two-block partial least squares visualization                         #
# Required packages: pls                                                #
# Necessary input variables:                                            #
#    Model: Output as produced by PLS.Shape.                            #
#    PC1: First principal component to use.                             #
#         *numeric (integer)*                                           #
#         default=1                                                     #
#    PC2: Second principal component to use.                            #
#         *numeric (integer)*                                           #
#         default=2                                                     #
#    Type: Which type of diagnostic plots should be produced? Options...#
#          are:                                                         #
#          "Scoreplot": Scores-and loadings-plot                        #
#          "Biplot": Residual variation and biplot                      #
#          "T.U.Plot": Correlation plot between T and U scores          #
#          *character*                                                  #
#          default=NULL                                                 #
#    Multiplic: Optional multiplicator for loadings vectors. Only...    #
#               meaningful for Type=="Biplot".                          #
#               *numeric (real)*                                        #
#    default=NULL                                                       #
# Output data: A series of diagnostic plots for the PLS.                #
# Input dataset: PLS Model list as produced by PLS.Shape.               #
#########################################################################

#Load packages
require(pls)

PLS.Visual<-function (Model, PC1=1, PC2=2, Type=NULL, Multiplic=NULL) {
	#Test data consistency
	PC1<-round(PC1, digits=0)
	PC2<-round(PC2, digits=0)
	if (is.null(Type)) {stop("'Type' must be provided!")}
	if (Type!="Scoreplot" & Type!="Biplot" & Type!="T.U.Plot") {stop("Type must be either of 'Scoreplot', 'Biplot', or 'T.U.Plot'!")}

	#Prepare data
	Mod<-Model$Model
	
	#Gather data
	PLS.scoresX<-scores(Mod)
	PLS.scoresY<-Yscores(Mod)
	PLS.loadingsX<-loadings(Mod)
	PLS.loadingsY<-Yloadings(Mod)
	PLS.residuals<-explvar(Mod)
	
	#PLS1 visualization
	if (Model$Number.Parameters==1) {
		if (PC1!=1) {PC1<-1; warning("PLS1 model, 'PC1' set to 1.")}
		#Scoreplot/loadingsplot
		if (Type=="Scoreplot") {
			win.graph(19, 10, 10)
			layout(matrix(c(1, 2), 1, 2))
			plot(max(PLS.scoresX[,PC1])*1000, 1000, xlim=range(PLS.scoresX[,PC1]), xlab=paste("PC ", PC1, sep=""), ylim=c(-1, 1), ylab="Meaningless", main="Scores")
			text(PLS.scoresX[,PC1], 0, labels=Model$Sample.Names)
			plot(max(c(PLS.loadingsX[,PC1], PLS.loadingsY[,PC1]))*1000, 1000, xlim=range(c(PLS.loadingsX[,PC1], PLS.loadingsY[,PC1])), xlab=paste("PC ", PC1, sep=""), ylim=c(1, 2), ylab="Meaningless", main="Loadings")
			text(PLS.loadingsX[,PC1], 1, labels=rownames(PLS.loadingsX), col="red")
			text(PLS.loadingsY[,PC1], 2, labels=rownames(PLS.loadingsY), col="blue")
		}
	
		#Residual variation/biplot
		if (Type=="Biplot") {
			hist(PLS.scoresX[,PC1], xlab=paste("PC ", PC1, " (", round(PLS.residuals[PC1], digits=1), "%)", sep=""), col="grey50", main="Biplot", sub=paste("Positive with ", rownames(PLS.loadingsX), sep=""))
		}
	
		#T-U-scoreplots
		if (Type=="T.U.Plot") {
			TU1<-lm(PLS.scoresY[,PC1]~PLS.scoresX[,PC1])
			R<-summary(TU1)$adj.r.squared
			P<-summary(TU1)$coefficients[2,"Pr(>|t|)"]
			plot(PLS.scoresX[,PC1], PLS.scoresY[,PC1], xlab="T scores", ylab="U scores", sub=paste("R2 = ", R, ", p = ", P, sep=""), main=paste("PC ", PC1, sep=""))
			curve(TU1$coefficients[2]*x+TU1$coefficients[1], add=TRUE, lwd=2)
		}
	}
	
	#PLS2 visualization
	if (Model$Number.Parameters>1) {
		#Scoreplot/loadingsplot
		if (Type=="Scoreplot") {
			win.graph(19, 10, 10)
			layout(matrix(c(1, 2), 1, 2))
			plot(max(PLS.scoresX[,PC1])*1000, max(PLS.scoresX[,PC2])*1000, xlim=range(PLS.scoresX[,PC1]), xlab=paste("PC ", PC1, sep=""), ylim=range(PLS.scoresX[,PC2]), ylab=paste("PC ", PC2, sep=""), main="Scores")
			text(PLS.scoresX[,PC1], PLS.scoresX[,PC2], labels=Model$Sample.Names)
			plot(max(c(PLS.loadingsX[,PC1], PLS.loadingsY[,PC1]))*1000, max(c(PLS.loadingsX[,PC2], PLS.loadingsY[,PC2]))*1000, xlim=range(c(PLS.loadingsX[,PC1], PLS.loadingsY[,PC1])), xlab=paste("PC ", PC1, sep=""), ylim=range(c(PLS.loadingsX[,PC2], PLS.loadingsY[,PC2])), ylab=paste("PC ", PC2, sep=""), main="Loadings")
			text(PLS.loadingsX[,PC1], PLS.loadingsX[,PC2], labels=rownames(PLS.loadingsX), col="red")
			text(PLS.loadingsY[,PC1], PLS.loadingsY[,PC2], labels=rownames(PLS.loadingsY), col="blue")
		}
	
		#Residual variation/biplot
		if (Type=="Biplot") {
			win.graph(19, 10, 10)
			layout(matrix(c(1, 2), 1, 2))
			barplot(PLS.residuals, col="blue", ylim=c(0, 100), ylab="Variation explained (%)", main="Variance explained")
			plot(PLS.scoresX[,PC1], PLS.scoresX[,PC2], pch=16, xlab=paste("PC ", PC1, " (", round(PLS.residuals[PC1], digits=1), "%)", sep=""), ylab=paste("PC ", PC2, " (", round(PLS.residuals[PC2], digits=1), "%)", sep=""), main="Biplot")
			{if (is.null(Multiplic)) {text(PLS.loadingsX[,PC1], PLS.loadingsX[,PC2], labels=rownames(PLS.loadingsX), col="red")}
			else {text(PLS.loadingsX[,PC1]*Multiplic, PLS.loadingsX[,PC2]*Multiplic, labels=rownames(PLS.loadingsX), col="red")}
			}
		}
	
		#T-U-scoreplots
		if (Type=="T.U.Plot") {
			win.graph(19, 10, 10)
			layout(matrix(c(1, 2), 1, 2))
			TU1<-lm(PLS.scoresY[,PC1]~PLS.scoresX[,PC1])
			R<-summary(TU1)$adj.r.squared
			P<-summary(TU1)$coefficients[2,"Pr(>|t|)"]
			plot(PLS.scoresX[,PC1], PLS.scoresY[,PC1], xlab="T scores", ylab="U scores", sub=paste("R2 = ", R, ", p = ", P, sep=""), main=paste("PC ", PC1, sep=""))
			curve(TU1$coefficients[2]*x+TU1$coefficients[1], add=TRUE, lwd=2)
			TU2<-lm(PLS.scoresY[,PC2]~PLS.scoresX[,PC2])
			R<-summary(TU2)$adj.r.squared
			P<-summary(TU2)$coefficients[2,"Pr(>|t|)"]
			plot(PLS.scoresX[,PC2], PLS.scoresY[,PC2], xlab="T scores", ylab="U scores", sub=paste("R2 = ", R, ", p = ", P, sep=""), main=paste("PC ", PC2, sep=""))
			curve(TU2$coefficients[2]*x+TU2$coefficients[1], add=TRUE, lwd=2)
		}
	}
}

#########################################################################
# Function to calculate and visualize optimal species delimitation.     #
# Required packages: pcaPP, mclust, mvoutlier, vegan                    #
# Necessary input variables:                                            #
#    Morpho.Data: Matrix containing multidimensional morphological data.#
#                 *matrix*                                              #
#    k: Desired number of total components to calculate in PCA...       #
#       (before dimension reduction).                                   #
#       *numeric (integer)*                                             #
#       default=number of parameters                                    #
#    Clust.N: Maximum number of clusters to check for in Bayesian...    #
#             clustering.                                               #
#             *numeric (integer)*                                       #
#             default=number of observations/10                         #
#    Ts, As: Scaling parameters for arrow and text positions, in case...#
#            automatic estimation yields unsatisfactory results.        #
#            *numeric (real)*                                           #
#            default=0.9 (Ts)/0.8 (As)                                  #
#    Group.Assignment: Desired group assignment of specimens. This is...#
#                      calculated internally but this value can be...   #
#                      overwritten.                                     #
#                      *vector*                                         #
#                      default=NULL                                     #
#    A.priori: A priori group assignments to be output alongside...     #
#              automated Group.Assignment for comparison.               #
#              *vector*                                                 #
# Output data: PCA plot and list containing number of dimensions,...    #
#              optimal number of clusters, species designation and...   #
#              outlier analysis, and robust PCA.                        #
# Input dataset: Matrix (variables in columns, observations in rows)... #
#                containing multidimensional morphological data.        #
#########################################################################

#Loading packages
require(pcaPP)
require(mclust)
require(mvoutlier)
require(vegan)

Species.Delimit<-function (Morpho.Data, k=nrow(Morpho.Data), Clust.N=round((nrow(Morpho.Data)/10), digits=0), Ts=0.9, As=0.8, Group.Assignment=NULL, A.priori=NULL)  {
	#Test data consistency
	k=round(k, digits=0)
	if (k>ncol(Morpho.Data)) {stop("k must be smaller or equal to the number of parameters!")}
	Clust.N<-round(Clust.N, digits=0)
	if (Clust.N>nrow(Morpho.Data)) {stop("Potential number of clusters cannot be larger than number of observations!")}
	if (is.null(colnames(Morpho.Data))) {colnames(Morpho.Data)<-paste("Param", 1:ncol(Morpho.Data), sep=".")}
	if (!is.null(Group.Assignment) & length(Group.Assignment)!=nrow(Morpho.Data)) {stop("Group.Assignment vector of wrong length! Each observation must be assigned to one group.")}
	if (!is.null(Group.Assignment) & !is.numeric(Group.Assignment)) {stop("Group.Assignment must be numeric vector!")}
	if (!is.null(A.priori) & length(A.priori)!=nrow(Morpho.Data)) {stop("A.priori vector of wrong length! Each observation must be assigned to one group.")}
	
	#Calculate robust PCA on Qn-method main direction and scale on median average deviation
	PCA<-PCAgrid(Morpho.Data, k=k, center=median, method="qn", scale=mad)
	
	#Dimension reduction of data
	k.Red<-max(which(PCA$sdev^2>bstick(PCA)))
	
	#Analyse PCA
	##Get PCA prediction values
	PCA.predict<-predict(PCA)[,1:k.Red]
	##Find optimal number of clusters
	PCA.clust<-Mclust(PCA.predict, G=1:Clust.N)
	
	#Find outliers in data
	{if (is.null(Group.Assignment)) {
		grps<-as.numeric(PCA.clust$classification)
	}
	else {grps<-as.numeric(Group.Assignment)}
	}
	out01<-sign2(as.data.frame(PCA.predict), qcrit=.975)
	{if (is.null(A.priori)) {out<-data.frame(Predicted.designation=grps, Outlier.0=out01$wfinal01)}
	else {out<-data.frame(A.priori.designation=as.character(A.priori), Predicted.designation=grps, Outlier.0=out01$wfinal01)}}
	if (!is.null(rownames(Morpho.Data))) {rownames(out)<-rownames(Morpho.Data)}

	#Visualize results
	##Plot observations
	XLIM=c(min(PCA$scores[,1]), max(PCA$scores[,1]))
	YLIM=c(min(PCA$scores[,2]), max(PCA$scores[,2]))
	plot(PCA$scores[,1], PCA$scores[,2], col=grps, pch=grps, xlab="rPC 1", xlim=XLIM, ylab="rPC 2", ylim=YLIM)
	##Plot parameters
	X.range<-c(min(PCA$loadings[1:k.Red,1]), max(PCA$loadings[1:k.Red,1]))
	Y.range<-c(min(PCA$loadings[1:k.Red,2]), max(PCA$loadings[1:k.Red,2]))
	Scale.range<-c(XLIM[1]/X.range[1], XLIM[2]/X.range[2], YLIM[1]/Y.range[1], YLIM[2]/Y.range[2])
	Scale.Text<-min(Scale.range)*Ts
	Scale.Arrow<-min(Scale.range)*As
	for (i in 1:k.Red) {
		arrows(0, 0, PCA$loadings[i,1]*Scale.Arrow, PCA$loadings[i,2]*Scale.Arrow, length=0.15, col="blue", lwd=2)
		text(PCA$loadings[i,1]*Scale.Text, PCA$loadings[i,2]*Scale.Text, labels=rownames(PCA$loadings)[i], col="blue")
	}
	
	#Return results
	return(list("Retained.Dimensions"=k.Red, "Cluster.Model"=PCA.clust, "Species.Analysis"=out, "rPCA"=PCA))
}

#--------------------------------------------

## EXAMPLES *************************************************************************************************

#Converting images
#ImageConversion("Foram", 1, 1, ImageType=".jpg", OutputType=".ppm")

#Extract landmarks
#LMExtract(Image="Foram", StartNum=1, StopNum=2, Output="Landmarks1", N=3, Export="NTS")
#LMExtract(Image="Foram", StartNum=1, StopNum=2, ScaleParam=c(5.84818294686911, 5.84818294686911), Output="Landmarks2", N=3, Export="NTS")
#LMExtract(Image="Foram", InputType="tif", StartNum=1, StopNum=2, Output="Landmarks3", N=3, Export="TPS")
#LMExtract(Image="Foram", StartNum=1, StopNum=2, Output="Landmarks4", Scale=FALSE, ScaleParam=c(0.2, 0.2), N=3, Export="TPS")

#Mirror landmarks
#LMExample<-matrix(c(2, 5, 6, 5, 5, 1, 4, 3, 2, 1), 5, 2, byrow=TRUE)
#LMExampleRef<-LM.Mirror(LMExample, Axis="y", RevCoord=c(2, 1, 5, 4, 3))
#LMExampleRef2<-LM.Mirror(LMExample, Axis="y", RevCoord=1:5)
#win.graph(20, 10, 10)
#layout(matrix(c(1, 2), 1, 2))
#plot(LMExample, asp=1)
#text(LMExample, pos=2, labels=(1:5))
#plot(LMExampleRef, asp=1)
#text(LMExampleRef, pos=2, labels=(1:5))
#points(LMExampleRef2, col="red")
#text(LMExampleRef2, pos=1, labels=(1:5))

#Read files
#Data<-Read.NTS("Landmarks.nts")
#Data<-Read.NTS("Landmarks.nts", na.remove=FALSE)
#Data<-Read.TPS("Landmarks3.tps")
#Data<-Read.TPS("Landmarks1.tps", na.remove=FALSE)
#Data2<-Read.TPS("Landmarks3.tps", Scale=FALSE, na.remove=FALSE)

#Write files
#Write.TPS(Data$LMData, Scaling=FALSE, Filenames=Data$Filenames, Scale=Data$Scale, Output="Test.tps")

#Average data
#Dat1<-Read.TPS("Landmarks3.tps")
#Dat2<-Read.TPS("Landmarks4.tps")
#LMAverage(Objects=c("Dat1$LMData", "Dat2$LMData"), Meta=list(c(Dat1$Filenames, Dat2$Filenames), c(Dat1$Scale, Dat2$Scale)), Error.File="ErrorTest.txt", Output="AverageTest.tps")

#Superimpose data
#require(shapes)
#Dat1<-gorm.dat
##Bookstein baseline registration
#Super1<-bookstein2d(Dat1, l1=1, l2=2)
##Partial Procrustes superimposition
#Super2<-pgPs(Dat1)
##Full Procrustes superimposition
#Super3<-procGPA(Dat1)
##Generalized resistant-fit superimposition
#Super4<-GRF(Dat1)

#win.graph(20, 20, 10)
#layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE))
#plotshapes(Super1$bshpv, joinline=c(1, 2))
#mtext("Bookstein registration", side=3, line=1)
#plotshapes(Super2$rotated)
#mtext("Partial procrustes", side=3, line=1)
#plotshapes(Super3$rotated)
#mtext("Full procrustes", side=3, line=1)
#plotshapes(Super4$rotated)
#mtext("Resistant-fit", side=3, line=1)

#Test allometry
#require(shapes)
#go<-procGPA(gorf.dat)
#Allo<-Allometry(go, Permutation=FALSE)
#Allo2<-Allometry(go, DistType="Procrustes")

#Plot shape deformation between groups
#require(shapes)
#gor<-array(c(gorf.dat, gorm.dat), dim=c(8, 2, 59))
#go<-procGPA(gor)
#LMDat<-go$rotated
##Vector deformation
#win.graph(20,10,10)
#layout(matrix(c(1, 2), 1, 2))
#ShapeDeform(LMDat, 1:30, 31:59)
#ShapeDeform(LMDat, 1:30, 31:59, Ampl=4)
##Thin plate splines
#win.graph(20,10,10)
#layout(matrix(c(1, 2), 1, 2))
#tps(LMDat, 1:30, 31:59, Lines=c(1, 6:8, 2:5, 1))
#tps(LMDat, 1:30, 31:59, n=16, Lines=c(1, 6:8, 2:5, 1), Ampl=4)
##Vectorized thin plate splines
#win.graph(26,10,10)
#layout(matrix(c(1, 2, 3), 1, 3))
#Vector.tps(LMDat, 1:30, 31:59, Lines=c(1, 6:8, 2:5, 1), Dens=600, Plot="Vector")
#Vector.tps(LMDat, 1:30, 31:59, Lines=c(1, 6:8, 2:5, 1), Plot="Contour")
#Vector.tps(LMDat, 1:30, 31:59, Lines=c(1, 6:8, 2:5, 1), Dens=10000, Plot="Surface")

#Explorative data mining
#require(shapes)
#gor<-array(c(gorf.dat, gorm.dat), dim=c(8, 2, 59))
#PCA1<-LM.PCA(gor, Sym=c(rep(16, 30), rep(17, 29)), Col=c(rep("blue", 30), rep("darkgreen", 29)), Lines=c(1, 6:8, 2:5, 1), Type="Explorative", n=20)
#PCA2<-LM.PCA(gor, Sym=c(rep(16, 30), rep(17, 29)), Col=c(rep("blue", 30), rep("darkgreen", 29)), Lines=c(1, 6:8, 2:5, 1), Type="Uniform", n=16)
#PCA3<-LM.PCA(gor, Sym=c(rep(16, 30), rep(17, 29)), Col=c(rep("blue", 30), rep("darkgreen", 29)), Lines=c(1, 6:8, 2:5, 1), Type="Nonaffine", n=16)
#LM.PCA.tps(PCA2$meanshape, PCA2$uniform, PCA2$scores, n=15, Lines=c(1, 6:8, 2:5, 1), Axis="PC1max", col="red", Line.col="darkgreen")
#LM.PCA.tps(PCA3$Mean.Shape, PCA3$NonUniform.Deform, PCA3$scores, n=20, Lines=c(1, 6:8, 2:5, 1), Axis="PC1max", col="red", pch=12, Line.col="darkgreen", Line.width=4, Title="Non-Uniform")

#Hotelling-Lawley trace
#Hotellingsp(gor, Groups=c(rep(1, 30), rep(2, 29)), exact=TRUE, Lines=c(1, 6:8, 2:5, 1))

#Automated species delimitation
#Spec<-Species.Delimit(tt, A.priori=rownames(tt))
#colnames(tt)<-paste("Value", 1:ncol(tt), sep="_")
#Spec2<-Species.Delimit(tt, Ts=0.5, As=0.4)
#Species<-c(rep(1, 50), rep(2, 50))
#Spec3<-Species.Delimit(tt, Group.Assignment=Species)

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
#		-Added a line to LMExtract that exports copy of image with landmarks
#
# Version 1.2
#	Date: XXX
#	Description of changes:
#		-Added the functionality to supply LMExtract with predefined scaling parameter
#
# Version 1.3
#	Date: XXX
#	Description of changes:
#		-Various little improvements in results presentation of LM.PCA and Hotellingsp
#
# Version 1.4
#	Date: XXX
#	Description of changes:
#		-Fixed a problem that would break FDM if landmark-distance variance was larger than mean landmark-distance
#
# Version 1.5
#	Date: XXX
#	Description of changes:
#		-Added function PLS.Shape
#
# Version 1.5.1
#	Date: XXX
#	Description of changes:
#		-Added function PLS.Visual
#
# Version 1.5.2
#	Date: XXX
#	Description of changes:
#		-Fixed some issues in PLS.Shape and PLS.Visual
#
# Version 1.5.3
#	Date: XXX
#	Description of changes:
#		-Modified Hotellingsp, so that the proportion of correct classification and the lda model is exported as well
#
# Version 1.5.4
#	Date: XXX
#	Description of changes:
#		-LM.PCA now allows to override alignment along principal axis by manual alignment for uniform and nonaffine deformation
#
# Version 1.5.5
#	Date: XXX
#	Description of changes:
#		-Improved visual output of function Hotellingsp
#
# Version 1.6
#	Date: XXX
#	Description of changes:
#		-Improved Read.NTS for missing value encoding and adapted NTS export of landmarks to use standard NA-encoding (-999)
#
# Version 1.6.1
#	Date: XXX
#	Description of changes:
#		-Removed Read.NTS, Read.TPS, and Write.TPS and transferred into separate file MorphoFiles_Function.r
#
# Version 1.7
#	Date: XXX
#	Description of changes:
#		-Added functionality to work with tif images in LM.extract
#
# Version 1.8
#	Date: XXX
#	Description of changes:
#		-Added function Allometry
#
# Version 1.8.1
#	Date: XXX
#	Description of changes:
#		-Corrected LM.PCA function for new rtiff handling of matrix plotting, LM.PCA now returns results
#
# Version 1.9
#	Date: XXX
#	Description of changes:
#		-Added function LM.PCA.tps
#		-Corrected other functions for new rtiff handling of matrix plotting, LM.PCA now returns results
#		-Added new plot customisation options to ShapeDeform and tps functions
#		-Removed ImageConversion and LMExtract functions to include in separate functions file MorphometricExtraction_Functions.r
#
# Version 1.10
#	Date: XXX
#	Description of changes:
#		-Added function Species.Delimit from former file MorphoSpeciesDelimitation_Function.r, which is now deprecated
#		-Added the functionality to Species.Delimit to export an a priori species delimitation
#
# Version 1.11
#	Date: XXX
#	Description of changes:
#		-Added calculation of leave-one-out cross-validated correct classification to Hotellingsp
#
# Version 1.11.1
#	Date: XXX
#	Description of changes:
#		-Fixed an error in Species.Delimit where k used ncol instead of nrow
#
# Version 1.11.2
#	Date: 7 October 2019
#	Description of changes:
#		-Added export of result values to MorphoCluster
##***********************************************************************************************************

#############################################################################################################







































