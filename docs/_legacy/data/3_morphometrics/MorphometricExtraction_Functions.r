## HEAD #####################################################################################################
#
# FUNCTION SET DESCRIPTION
#	Extracting morphometric data from images in R
#
# DATASET DESCRIPTION
#	Images, possibly binarized
#	
# METADATA
#	Author: Manuel F. G. Weinkauf
#	E-mail: weinkauf.scientific@gmail.com
#	R-version: 4.2.1
#	RStudio-version: 2022.07.1
#	Code-version: 2.0
#	Date of last update: 19 August 2022

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

## BODY #####################################################################################################

#**************************************************************************************
#Setting working dierctory
#setwd("C:/R_TestData/GeometricMorphometrics")

#########################################################################
# Image picking                                                         #
# Necessary input variables:                                            #
#    Dir: Directory path.                                               #
#         *character*                                                   #
#    pattern: Which file types to look for (i.e. the types of your)...  #
#             image files. By default scans for the common file types...#
#             "jpg", "png", "tif", "bmp", "gif", and the common type... #
#             used in moprhometrics "ppm".                              #
#             *character*                                               #
# Output data: Vector of all image files to use.                        #
# Input dataset: Computer directory.                                    #
#########################################################################

ImagePicking<-function (Dir, pattern=c("bmp", "gif", "jpg", "jpeg", "png", "tif", "tiff", "ppm")) {
	#Set up pattern parameter
	pattern<-paste("\\.", pattern, "$", sep="", collapse="|")

	#List all files
	File.list<-list.files(Dir, pattern=pattern, ignore.case=TRUE)
	names(File.list)<-1:length(File.list)
	print(File.list)
	
	#Provide the choice of files by number
	writeLines("Please input the numbers of the files you wish to use. \nSeveral numbers can be input separated by commas. \nDO NOT ENTER SPACES!")
	Image.No<-readline(prompt="Enter the image numbers: ")
	Image.No<-scan(text=Image.No, quiet=TRUE, sep=",")
	Images<-File.list[Image.No]
	
	#Output image names
	Images<-as.data.frame(Images)
	Images[]<-lapply(Images, as.character)
	return(Images)
}

#########################################################################
# Image conversion                                                      #
# Necessary programs: Image Magic                                       #
# Necessary input variables:                                            #
#    Images: A list of all the images to be converted.                  #
#            *data frame* with image names in first column.             #
#    Scaling: Relative scaling of the image during conversion.          #
#             *numeric (integer)*                                       #
#             default=100                                               #
#    OutputType: Image type for output as "xyz".                        #
#                *character*                                            #
# Output data: Images in OutputType format.                             #
# Input dataset: Images to convert.                                     #
#########################################################################

ImageConversion<-function(Images, Scaling=100, OutputType) {
	#Test data consistency
	if (!is.data.frame(Images)) {stop("'Images' must be a data frame!")}
	if (ncol(Images)!=1) {warning("'Images' has more than one column. Only first column will be used!")}
	
	#Read image list into vector
	Images<-as.vector(Images[,1])
	
	#Convert images
	for (i in 1:length(Images)){
		Name<-strsplit(Images[i], split=".", fixed=TRUE)[[1]][1]
		com<-paste("convert ", Images[i], " ", "-resize ", Scaling, "%", " ", Name, ".", OutputType, sep="")
		shell(com)
	}
}

#########################################################################
# Outline rastering                                                     #
# based on Claude (2008), p. 47                                         #
# Necessary input variables:                                            #
#    x: Manually chosen starting point                                  #
#       *integer*                                                       #
#    imagematrix: Name of image to use                                  #
#                 *string*                                              #
# Output data: List of x-y coordinates of all pixel along outline.      #
# Input dataset: Image in ppm format.                                   #
#########################################################################

Conte<-function (x, imagematrix) {
	#Reading image as image matrix and calculate dimensions
	I<-imagematrix
	x<-rev(x)
	x[1]<-dim(I)[1]-x[1]
	
	#Going from the starting point to the left, until reaching a pixel that significantly differs in grey value from the pixel next to it
	while (abs(I[x[1], x[2]]-I[x[1], (x[2]-1)])<0.1) {x[2]<-x[2]-1}
	#Defining that pixel as starting point of the outline
	a<-1
	
	#Writing indices of eight pixels surrounding starting position (0,0) into matrix
	M<-matrix(c(0, -1, -1, -1, 0, 1, 1, 1, 1, 1, 0, -1, -1, -1, 0, 1), 2, 8, byrow=TRUE)
	#Attaching first and last column of matrix as ninth and tenth column on matrix
	M<-cbind(M[,8], M, M[,1])
	
	#Initialising parameters
	X<-0
	Y<-0
	x1<-x[1]
	x2<-x[2]
	SS<-NA
	S<-6
	
	#Allocating all pixels that belong to the outline (running clockwise)
	while ((any(c(X[a], Y[a])!=c(x1, x2)) | length(X)<3)) {
		if(abs(I[x[1]+M[1,S+1], x[2]+M[2,S+1]]-I[x[1], x[2]])<0.1) {a<-a+1; X[a]<-x[1]; Y[a]<-x[2]; x<-x+M[,S+1]; SS[a]<-S+1; S<-(S+7)%%8}
		else if (abs(I[x[1]+M[1,S+2], x[2]+M[2,S+2]]-I[x[1], x[2]])<0.1) {a<-a+1; X[a]<-x[1]; Y[a]<-x[2]; x<-x+M[,S+2]; SS[a]<-S+2; S<-(S+7)%%8}
		else if (abs(I[x[1]+M[1,(S+3)], x[2]+M[2,(S+3)]]-I[x[1], x[2]])<0.1) {a<-a+1; X[a]<-x[1]; Y[a]<-x[2]; x<-x+M[,(S+3)]; SS[a]<-S+3; S<-(S+7)%%8}
		else S<-(S+1)%%8
	}
	list(X = (Y[-1]), Y = ((dim(I)[1]-X))[-1])
}

#########################################################################
# Curvelinear fitting                                                   #
# based on Claude (2008), p. 52                                         #
# Necessary input variables:                                            #
#    Input: Raw pixel coordinates of outline as extracted with...       #
#           function Conte.                                             #
#           *vector*                                                    #
#    n: Desired number of curvelinearily equally spaced points along... #
#       the  outline.                                                   #
#       *integer*                                                       #
# Output data: Matrix with equally spaced coordinates of outline.       #
# Input dataset: Vector of raw x and y coordinates.                     #
#########################################################################

EquiDist<-function (Input, n) {
	RcEqX<-(Input$X[seq(1, length(Input$X), length=(n+1))])[-1]
	RcEqY<-(Input$Y[seq(1, length(Input$Y), length=(n+1))])[-1]
	X_Temp<-as.matrix(RcEqX)
	Y_Temp<-as.matrix(RcEqY)
	M<-cbind(X_Temp, Y_Temp)
	colnames(M)<-c("X", "Y")
	rownames(M)<-paste("Coord", 1:n, sep="")
	M
}

#########################################################################
# Smoothing outline                                                     #
# based on Claude (2008), p. 55                                         #
# Necessary input variables:                                            #
#    M: Matrix of equally spaced outline points.                        #
#       *matrix*                                                        #
#    n: Number of iterations.                                           #
#       *integer*                                                       #
# Output data: List of smoothed x-y coordinates of outline.             #
# Input dataset: Curvelinear equally spaced x-y coordinates of outline. #
#########################################################################

smoothout<-function(M, n){
	p<-dim(M)[1]
	a<-0
	while(a<=n){
		a<-a+1
		Ms<-rbind(M[p,], M[-p,])
		Mi<-rbind(M[-1,], M[1,])
		M<-M/2+Ms/4+Mi/4
	}
	M
}

#########################################################################
# Outline extraction                                                    #
# Necessary packages: pixmap, tiff, splancs                             #
# Necessary functions: Conte, EquiDist, smoothout                       #
# Necessary input variables:                                            #
#    Images: A list of all the images to be converted.                  #
#            *vector* of image names.                                   #
#    Output: Name of the output file (excluding extension).             #
#            *string*                                                   #
#    Specimen.Labels: Names for the specimens. If NULL, numbers will... #
#                     be used. If "Names", the name part of the...      #
#                     images will be used. If vector, is must be a...   #
#                     a character vector with names of the same...      #
#                     as Images.                                        #
#                     either NULL or "Names" or *character*             #
#                     default=NULL                                      #
#    OutlinePoints: Desired number of equidistantly spaced points...    #
#                   along outline (see function EquiDist).              #
#                   *integer*                                           #
#                   default=100                                         #
#    Smoothing: Desired number of iterations for outline smoothing...   #
#               (see function smoothout).                               #
#               *integer*                                               #
#               default=1                                               #
#    Baseline: If True, you have to provide two points which can...     #
#              serve as baseline for uniform rotation during the NEF... #
#              (with option Rotation="Baseline"). Note that the two...  #
#              points have to be digitised in the same order in each... #
#              image.                                                   #
#              *logical*                                                #
#              TRUE=Extract baseline coordinates                        #
#              FALSE=Do not extract baseline coordinates (uniform...    #
#                    rotation will be achieved on the basis of the...   #
#                    longest axis of the first harmonic during NEF).    #
#                    default=FALSE                                      #
#    Scale: Do the images contain a scale bar that should should be...  #
#           used to extract the size of the object as cross-sectional...#
#           area?                                                       #
#           *logical*                                                   #
#           TRUE=Include step to extract size information               #
#           FALSE=Extract outline coordinates only                      #
#           default=FALSE                                               #
#    ScaleParam: Vector containing the known 1 px=y units conversion... #
#                for each image.                                        #
#                *numeric (real)*, length=number of images              #
#                default=NULL                                           #
#    RawVersion: Defines whether or not the raw (i.e. unsmoothed...     #
#                version of the outline coordinates should be ex-...    #
#                ported as well.                                        #
#                *logical*                                              #
#                TRUE=Export raw version                                #
#                FALSE=Do not export raw version                        #
#                default=FALSE                                          #
# Output data: .nts file x and y coordinates of specified number of...  #
#              equidistant points along outline.                        #
# Input dataset: Images to extract data from. Supported formats:...     #
#                .ppm, .tif. For .ppm files a rudimentary...            #
#                thresholding is performed, .tif images are read as...  #
#                they are. It is highly recommended for all images...   #
#                to be converted into high-contrast threshold images... #
#                before extraction: black or grey object of interest... #
#                on white background.                                   #
#########################################################################

#Loading packages
require(stringr)
require(sfsmisc)
require(pixmap)
require(tiff)
require(splancs)

OutlineExtraction<-function (Images, Output, Specimen.Labels=NULL, OutlinePoints=100, Smoothing=1, Baseline=FALSE, Scale=FALSE, ScaleParam=NULL, RawVersion=FALSE) {
	#Test data consistency
	if (!is.vector(Images)) {stop("'Images' must be a vector!")}
	{if (!is.null(Specimen.Labels)) {
		if (Specimen.Labels!="Names" && length(Specimen.Labels)!=nrow(Images)) {stop("Specimen.Labels must have same length as number of images!")}}
	}
	if (Scale==TRUE & !is.null(ScaleParam)) {Scale<-FALSE; warning("Cannot have scalebar and scale parameter at the same time, Scale has been set to FALSE")}
	if (!is.null(ScaleParam) & length(ScaleParam)!=length(Images)) {stop("Must supply one value of ScaleParam per image")}
	{if (!is.null(Specimen.Labels)) {
		if (Specimen.Labels!="Names" && length(Specimen.Labels)!=length(Images)) {stop("Specimen.Labels must have same length as number of images!")}}
	}
		
	#Read image list into vector
	InputType<-unlist(lapply(strsplit(Images, split=".", fixed=TRUE), tail, n=1L))
	InputType<-tolower(InputType)
	if (!all(InputType%in%c("tif", "tiff", "ppm"))) {stop("Only '.ppm' and '.tif' format images accepted!")}
	
	#Setup success report matrix
	ExtFail<-as.data.frame(matrix(NA, length(Images), 3))
	colnames(ExtFail)<-c("ID", "Image", "Success")
	
	#Setting up sizes matrix
	if (Scale==TRUE) {
		Sizes<-as.data.frame(matrix(NA, length(Images), 3))
		colnames(Sizes)<-c("ID", "Image", "Area")
	}
	
	#Setting up baseline matrix
	if (Baseline==TRUE) {
		Base<-as.data.frame(matrix(NA, length(Images), 6))
		colnames(Base)<-c("ID", "Image", "x1", "y1", "x2", "y2")
	}
	
	#Setting up results matrix
	if (RawVersion==TRUE) {CoordRes.Raw<-matrix(NA, length(Images), OutlinePoints*2)}
	CoordRes<-matrix(NA, length(Images), OutlinePoints*2)
	xc<-seq.int(from=1, to=OutlinePoints*2, by=2)
	yc<-seq.int(from=2, to=OutlinePoints*2, by=2)
	
	#Start data acquisition
	for (i in 1:length(Images)) {
		print(i)
		FileName<-Images[i]
		
		#Read and plot image
		{if (InputType[i]=="ppm") {
			y<-read.pnm(FileName, cellres=1)
		}
		else if (InputType[i]=="tif" | InputType[i]=="tiff") {
			y<-readTIFF(FileName, native=TRUE, convert=TRUE)
			{if (length(dim(y))>2) {
				if (dim(y)[3]==4) {Transparency<-TRUE} else {Transparency<-FALSE}
			}
				else {Transparency<-FALSE}
			}
			if (Transparency==TRUE) {transparent<-y[,,4]==0}
			{if (length(dim(y))==2) {y<-as.raster(y)}
			else {y<-as.raster(y[,,1:3])}}
			if (Transparency==TRUE) {y[transparent]<-NA}
		}
		}
		
		#Manually chosing starting point and decide whether or not outline is well fitted
		cont<-"n"
		while(cont=="n") {
			#Setting margins and plotting image
			{if (InputType[i]=="ppm") {
				#Converting image to grey scale
				y<-as(y, "pixmapGrey")
				y@grey[which(y@grey>=0.9)]<-1
				y@grey[which(y@grey<0.9)]<-0.7
				par(mar=(c(1, 1, 1, 1)))
				plot(y)
				}
			else if (InputType[i]=="tif" | InputType[i]=="tiff") {
				par(mar=(c(0.1, 0.1, 0.1, 0.1)))
				plot(x=c(0, dim(y)[2]), y=c(0, dim(y)[1]), type="n", asp=1, axes=FALSE)
				rasterImage(y, xleft=0, ybottom=0, xright=dim(y)[2], ytop=dim(y)[1], interpolate=FALSE)
				}
			}
		
			#Set scale
			if (Scale==TRUE) {
				writeLines("Please provide scale of the image. \nClick on two points on the scalebar.")
				flush.console()
				a<-locator(2, type="o", pch=8, lwd=2, col="grey60", lty=11)
				scale.px<-sqrt(sum(diff(a$x)**2+diff(a$y)**2))
				scale.length<-as.numeric(readline("How long is the implied line? "))
			}
		
			#Extract baseline
			if (Baseline==TRUE) {
				writeLines("Please provide baseline. \nClick on two landmarks of the object (order-sensitive).")
				flush.console()
				a<-locator(2, type="o", pch=8, lwd=2, col="grey60", lty=11)
				{if (is.null(Specimen.Labels)) {Base[i,"ID"]<-as.character(i)}
				else if (Specimen.Labels=="Names") {Base[i,"ID"]<-strsplit(FileName, split=".", fixed=TRUE)[[1]][1]}
				else {Base[i,"ID"]<-as.character(Specimen.Labels[i])}}
				Base[i,"Image"]<-FileName
				Base[i,"x1"]<-a$x[1]
				Base[i,"y1"]<-a$y[1]
				Base[i,"x2"]<-a$x[2]
				Base[i,"y2"]<-a$y[2]
			}

			#Extract outline
			print("Click within the object to the right of the starting point!")
			flush.console()
			start<-locator(1)
			
			{if (InputType[i]=="ppm") {
				Rc<-Conte(c(round(start$x), round(start$y)), y@grey)
			}
			else if (InputType[i]=="tif" | InputType[i]=="tiff") {
				rgb.vals<-col2rgb(y)
				grey.vals<-(0.299*rgb.vals[1,] + 0.587*rgb.vals[2,] + 0.114*rgb.vals[3,])/255
				y.grey<-matrix(grey.vals, nrow(y), ncol(y), byrow=TRUE)
				Rc<-Conte(c(round(start$x), round(start$y)), y.grey)
			}
			}
			
			lines(Rc$X, Rc$Y, lwd=3, col="deeppink1")
			arrows(0, Rc$Y[1], Rc$X[1], Rc$Y[1], length=0.3, col="red", lwd=4)
			
			cont<-readline("Is the outline correct? (y=proceed,n=try again,c=cancel and proceed)")
		}
		
		#Save copy of image with outlines for later comparison
		if (cont=="y") {
			dev.copy(png, filename=paste(strsplit(FileName, split=".", fixed=TRUE)[[1]][1], "_Outline.png", sep=""), width=5, height=5, units="in", res=200);
			dev.off();
		}
		
		#Write success report
		{if (is.null(Specimen.Labels)) {ExtFail[i,"ID"]<-as.character(i)}
		else if (Specimen.Labels=="Names") {ExtFail[i,"ID"]<-strsplit(FileName, split=".", fixed=TRUE)[[1]][1]}
		else {ExtFail[i,"ID"]<-as.character(Specimen.Labels[i])}}
		ExtFail[i,"Image"]<-FileName
		{if (cont=="y") {ExtFail[i,"Success"]<-1}
		else {ExtFail[i,"Success"]<-0}}

		#Normalising curve for curvelinear equally spaced number of points
		RcEqual<-EquiDist(Rc, OutlinePoints)
		
		#Smoothing outline
		RcSmooth<-smoothout(RcEqual, Smoothing)
		
		#Calculating sizes
		{if (Scale==TRUE) {
			{if (cont=="y") {
				RcScaled<-cbind(Rc$X,Rc$Y)
				for (j in 2:nrow(RcScaled)) {
					XD<-RcScaled[j,1]-RcScaled[1,1]
					YD<-RcScaled[j,2]-RcScaled[1,2]
					RcScaled[j,1]<-RcScaled[1,1]+(XD/scale.px)*scale.length
					RcScaled[j,2]<-RcScaled[1,2]+(YD/scale.px)*scale.length
				}
				RcScaled<-as.matrix(rbind(RcScaled, RcScaled[1,]))
				{if (is.null(Specimen.Labels)) {Sizes[i,"ID"]<-as.character(i)}
				else if (Specimen.Labels=="Names") {Sizes[i,"ID"]<-strsplit(FileName, split=".", fixed=TRUE)[[1]][1]}
				else {Sizes[i,"ID"]<-as.character(Specimen.Labels[i])}}
				Sizes[i,"Image"]<-FileName
				Sizes[i,"Area"]<-areapl(RcScaled)
			}
			else {
				{if (is.null(Specimen.Labels)) {Sizes[i,"ID"]<-as.character(i)}
					else if (Specimen.Labels=="Names") {Sizes[i,"ID"]<-strsplit(FileName, split=".", fixed=TRUE)[[1]][1]}
					else {Sizes[i,"ID"]<-as.character(Specimen.Labels[i])}}
					Sizes[i,"Image"]<-FileName
					Sizes[i,"Area"]<-NA
				}
			}
		}
		else if (!is.null(ScaleParam)) {
			{if (cont=="y") {
				RcScaled<-cbind(Rc$X,Rc$Y)
				for (j in 2:nrow(RcScaled)) {
					XD<-RcScaled[j,1]-RcScaled[1,1]
					YD<-RcScaled[j,2]-RcScaled[1,2]
					RcScaled[j,1]<-RcScaled[1,1]+XD*ScaleParam[i]
					RcScaled[j,2]<-RcScaled[1,2]+YD*ScaleParam[i]
				}
				RcScaled<-as.matrix(rbind(RcScaled, RcScaled[1,]))
				{if (is.null(Specimen.Labels)) {Sizes[i,"ID"]<-as.character(i)}
				else if (Specimen.Labels=="Names") {Sizes[i,"ID"]<-strsplit(FileName, split=".", fixed=TRUE)[[1]][1]}
				else {Sizes[i,"ID"]<-as.character(Specimen.Labels[i])}}
				Sizes[i,"Image"]<-FileName
				Sizes[i,"Area"]<-areapl(RcScaled)
			}
			else {
				{if (is.null(Specimen.Labels)) {Sizes[i,"ID"]<-as.character(i)}
					else if (Specimen.Labels=="Names") {Sizes[i,"ID"]<-strsplit(FileName, split=".", fixed=TRUE)[[1]][1]}
					else {Sizes[i,"ID"]<-as.character(Specimen.Labels[i])}}
					Sizes[i,"Image"]<-FileName
					Sizes[i,"Area"]<-NA
				}
			}
		}
		}
		
		#Writing x-y coordinates of outline into variable
		{if (cont=="y") {
			if (RawVersion==TRUE) {
				CoordRes.Raw[i,xc]<-RcEqual[,"X"]
				CoordRes.Raw[i,yc]<-RcEqual[,"Y"]
			}
			CoordRes[i,xc]<-RcSmooth[,"X"]
			CoordRes[i,yc]<-RcSmooth[,"Y"]
		}
		else if (cont=="c") {
			if (RawVersion==TRUE) {CoordRes.Raw[i,]<--999}
			CoordRes[i,]<--999
		}
		}
	}
	
	#Save outline coordinates as NTS file
	if (RawVersion==TRUE) {
		FileName<-paste(Output, "_Raw.nts", sep="")
		{if (any(ExtFail[,"Success"]==0)) {firstl<-paste(1, paste(length(Images), "L", sep=""), OutlinePoints*2, 1, -999, "dim=2", sep=" ")}
		else {firstl<-paste(1, paste(length(Images), "L", sep=""), OutlinePoints*2, 0, "dim=2", sep=" ")}}
		{if (is.null(Specimen.Labels)) {L<-1:length(Images)}
		else if (Specimen.Labels=="Names") {L<-unlist(lapply(strsplit(Images, split=".", fixed=TRUE), "[[", 1))}
		else {L<-Specimen.Labels}}
		secondl<-paste(L, sep="", collapse=" ")
		##Create file and write header
		cat(firstl, secondl, file=FileName, sep="\n", append=FALSE)
		##Create data body
		for (i in 1:nrow(CoordRes.Raw)) {
			B<-paste(CoordRes.Raw[i,], sep="", collapse=" ")
			cat(B, file=FileName, sep="\n", append=TRUE)
		}
	}
	FileName<-paste(Output, ".nts", sep="")
	{if (any(ExtFail[,"Success"]==0)) {firstl<-paste(1, paste(length(Images), "L", sep=""), OutlinePoints*2, 1, -999, "dim=2", sep=" ")}
	else {firstl<-paste(1, paste(length(Images), "L", sep=""), OutlinePoints*2, 0, "dim=2", sep=" ")}}
	{if (is.null(Specimen.Labels)) {L<-1:length(Images)}
	else if (Specimen.Labels=="Names") {L<-unlist(lapply(strsplit(Images, split=".", fixed=TRUE), "[[", 1))}
	else {L<-Specimen.Labels}}
	secondl<-paste(L, sep="", collapse=" ")
	##Create file and write header
	cat(firstl, secondl, file=FileName, sep="\n", append=FALSE)
	##Create data body
	for (i in 1:nrow(CoordRes)) {
		B<-paste(CoordRes[i,], sep="", collapse=" ")
		cat(B, file=FileName, sep="\n", append=TRUE)
	}
	
	#Save list of successes, i.e. in which specimens did outline extraction fail, size information, and baseline coordinates
	FileName<-paste(Output, "_Success.txt", sep="")
	write.table(ExtFail, FileName, sep="\t", row.names=FALSE)
	if (Scale==TRUE) {FileName<-paste(Output, "_Area.txt", sep=""); write.table(Sizes, FileName, sep="\t", row.names=FALSE)}
	if (Baseline==TRUE) {FileName<-paste(Output, "_Baseline.txt", sep=""); write.table(Base, FileName, sep="\t", row.names=FALSE)}
}

#########################################################################
# Landmark extraction                                                   #
# Necessary packages: pixmap, tiff, jpeg                                #
# Necessary input variables:                                            #
#    Images: A list of all the images to be converted.                  #
#            *vector* of image names.                                   #
#    Output: Name of the output file (excluding extension).             #
#            *character*                                                #
#    Specimen.Labels: Names for the specimens. If NULL, numbers will... #
#                     be used. If "Names", the name part of the...      #
#                     images will be used. If vector, is must be a...   #
#                     a character vector with names of the same...      #
#                     as Images.                                        #
#                     either NULL or "Names" or *character*             #
#                     default=NULL                                      #
#    Scale: Does the image contain a scale for which the the point...   #
#           coordinates should be normalized?                           #
#           NOTE: if Scale==TRUE and Export=="NTS" the coordinates...   #
#           will be exported in scaled format. However, if...           #
#           Export=="TPS" the coordinates will be exported as pixel...  #
#           coordinates and the SCALE= parameter will be included for...#
#           later conversion.                                           #
#           *logical*                                                   #
#           TRUE: Scale bar present                                     #
#           FALSE: Scale bar missing or not used                        #
#           default=TRUE                                                #
#    ScaleParam: Vector containing the known 1 px=y units conversion... #
#                for each image.                                        #
#                *numeric (real)*, length=number of images              #
#                default=NULL                                           #
#    N: Number of landmark points to extract.                           #
#       *numeric (integer)*                                             #
#    Export: In which format should results be exported? It can be...   #
#            one of the following: "NTS", "TPS"                         #
#            *character*                                                #
#            default="TPS"                                              #
# Output data: List of x-y coordinates of landmarks in .nts or .tps...  #
#              file.                                                    #
# Input dataset: Images to extract data from. Supported formats:...     #
#                .ppm, .tif, .jpg.                                      #
#########################################################################

#Load packages
require(stringr)
require(sfsmisc)
require(pixmap)
require(tiff)
require(jpeg)

LMExtract<-function(Images, Output, Specimen.Labels=NULL, Scale=TRUE, ScaleParam=NULL, N, Export="TPS") {
	#Test data for consistency
	if (!is.vector(Images)) {stop("'Images' must be a vector!")}
	if (Export!="NTS" & Export!="TPS") {stop("Export format must be either NTS or TPS!")}
	if (Scale==TRUE & !is.null(ScaleParam)) {Scale<-FALSE; warning("Cannot have scalebar and scale parameter at the same time, Scale has been set to FALSE")}
	if (!is.null(ScaleParam) & length(ScaleParam)!=length(Images)) {stop("Must supply one value of ScaleParam per image")}
	{if (!is.null(Specimen.Labels)) {
		if (Specimen.Labels!="Names" && length(Specimen.Labels)!=length(Images)) {stop("Specimen.Labels must have same length as number of images!")}}
	}
	
	#Read image list into vector
	InputType<-unlist(lapply(strsplit(Images, split=".", fixed=TRUE), tail, n=1L))
	InputType<-tolower(InputType)
	if (!all(InputType%in%c("tif", "tiff", "ppm", "jpg", "jpeg"))) {stop("Only '.ppm', '.tif', and '.jpg' format images accepted!")}
	
	#Set up results matrices
	{if (Export=="NTS") {
		CoordRes<-matrix(NA, length(Images), N*2)
		xc<-seq.int(from=1, to=N*2, by=2)
		yc<-seq.int(from=2, to=N*2, by=2)
	}
	else {
		CoordRes<-array(NA, dim=c(N, 2, length(Images)))
		Meta<-list()
		Meta$Files<-vector(mode="character", length=length(Images))
		Meta$ID<-vector(mode="character", length=length(Images))
		Meta$Scale<-vector(mode="numeric", length=length(Images))
	}}
	
	#Setup success report matrix
	ExtFail<-as.data.frame(matrix(NA, length(Images), 3))
	colnames(ExtFail)<-c("ID", "Image", "Success")
	
	#Extract landmarks
	for (i in 1:length(Images)) {
		print(i)
		FileName<-Images[i]

		#Read and plot image
		{if (InputType[i]=="ppm") {
			y<-read.pnm(FileName)
		}
		else if (InputType[i]=="tif" | InputType[i]=="tiff") {
			y<-readTIFF(FileName, native=TRUE, convert=TRUE)
			if (dim(y)[3]==4) {Transparency<-TRUE} else {Transparency<-FALSE}
			if (Transparency==TRUE) {transparent<-y[,,4]==0}
			y<-as.raster(y[,,1:3])
			if (Transparency==TRUE) {y[transparent]<-NA}
		}
		else if (InputType[i]=="jpg" | InputType[i]=="jpeg") {
			y<-readJPEG(FileName, native=TRUE)
		}}
		
		cont<-"n"
		while(cont=="n") {
			{if (InputType[i]=="ppm") {
				par(mar=(c(1, 1, 1, 1)))
				plot(y)
			}
			else if (InputType[i]=="tif" | InputType[i]=="tiff") {
				par(mar=(c(0.1, 0.1, 0.1, 0.1)))
				plot(x=c(0, dim(y)[2]), y=c(0, dim(y)[1]), type="n", asp=1, axes=FALSE)
				rasterImage(y, xleft=0, ybottom=0, xright=dim(y)[2], ytop=dim(y)[1], interpolate=FALSE)
			}
			else if (InputType[i]=="jpg" | InputType[i]=="jpeg") {
				par(mar=(c(0.1, 0.1, 0.1, 0.1)))
				plot(x=c(0, dim(y)[2]), y=c(0, dim(y)[1]), type="n", asp=1, axes=FALSE)
				rasterImage(y, xleft=0, ybottom=0, xright=dim(y)[2], ytop=dim(y)[1], interpolate=FALSE)
			}
			}
	
			#Set scale
			if (Scale==TRUE) {
				writeLines("Please provide scale of the image. \nClick on two points on the scalebar.")
				flush.console()
				a<-locator(2, type="o", pch=8, lwd=2, col="grey60", lty=11)
				scale.px<-sqrt(sum(diff(a$x)**2+diff(a$y)**2))
				scale.length<-as.numeric(readline("How long is the implied line? "))
			}
	
			#Extract and label landmarks
			print(paste("Digitize ", N, " landmarks by clicking into the image!"))
			flush.console()
			LM<-locator(n=N, type="p", pch=3, col="red", lwd=2)
			text(LM, pos=2, labels=1:N, col="red", font=2)
			
			cont<-readline("Are the landmarks correct? (y=proceed, n=try again, c=cancel and proceed)")
		}
		
		#Save copy of image with points for later comparison
		if (cont=="y") {
			dev.copy(png, filename=paste(strsplit(FileName, split=".", fixed=TRUE)[[1]][1], "_Landmarks.png", sep=""), width=5, height=5, units="in", res=200);
			dev.off();
		}
		
		#Write success report
		{if (is.null(Specimen.Labels)) {ExtFail[i,"ID"]<-as.character(i)}
		else if (Specimen.Labels=="Names") {ExtFail[i,"ID"]<-strsplit(FileName, split=".", fixed=TRUE)[[1]][1]}
		else {ExtFail[i,"ID"]<-as.character(Specimen.Labels[i])}}
		ExtFail[i,"Image"]<-FileName
		{if (cont=="y") {ExtFail[i,"Success"]<-1}
		else {ExtFail[i,"Success"]<-0}}
	
		#Write landmark coordinates into table
		{if (Scale==TRUE & Export=="NTS") {LM$x<-LM$x*(scale.length/scale.px); LM$y<-LM$y*(scale.length/scale.px)}
		else if (Scale==FALSE & !is.null(ScaleParam) & Export=="NTS") {LM$x<-LM$x*ScaleParam[i]; LM$y<-LM$y*ScaleParam[i]}
		else if (Scale==TRUE & Export=="TPS") {ScaleTPS=scale.length/scale.px}
		else if (Scale==FALSE & !is.null(ScaleParam) && Export=="TPS") {ScaleTPS=ScaleParam[i]}}
	
		{if (Export=="NTS") {
			{if (cont=="y") {
				CoordRes[i,xc]<-LM$x
				CoordRes[i,yc]<-LM$y
			}
			else if (cont=="c") {
				CoordRes[i,]<--999
			}
			}
		}
		else {
			if (cont=="y") {
				CoordRes[,,i]<-cbind(LM$x, LM$y)
				if (Scale==TRUE | !is.null(ScaleParam)) {Meta$Scale[i]<-ScaleTPS}
			}
			Meta$Files[i]<-FileName
			{if (is.null(Specimen.Labels)) {Meta$ID[i]<-as.character(i)}
			else if (Specimen.Labels=="Names") {Meta$ID[i]<-strsplit(FileName, split=".", fixed=TRUE)[[1]][1]}
			else {Meta$ID[i]<-as.character(Specimen.Labels[i])}}
		}}
	}
	
	#Export coordinates
	{if (Export=="NTS") {
		FileName<-paste(Output, ".nts", sep="")
		{if (any(ExtFail[,1]==0)) {firstl<-paste(1, paste(length(Images), "L", sep=""), N*2, 1, -999, "dim=2", sep=" ")}
		else {firstl<-paste(1, paste(length(Images), "L", sep=""), N*2, 0, "dim=2", sep=" ")}}
		{if (is.null(Specimen.Labels)) {L<-1:length(Images)}
		else if (Specimen.Labels=="Names") {L<-unlist(lapply(strsplit(Images, split=".", fixed=TRUE), "[[", 1))}
		else {L<-Specimen.Labels}}
		secondl<-paste(L, sep="", collapse=" ")
		##Create file and write header
		cat(firstl, secondl, file=FileName, sep="\n", append=FALSE)
		##Create data body
		for (i in 1:nrow(CoordRes)) {
			B<-paste(CoordRes[i,], sep="", collapse=" ")
			cat(B, file=FileName, sep="\n", append=TRUE)
		}
	}
	else {
		FileName<-paste(Output, ".tps", sep="")
		firstl<-paste("LM=", N, sep="")
		##Create file
		for (j in 1:dim(CoordRes)[3]) {
			cat(firstl, file=FileName, sep="\n", append=TRUE)
			for (i in 1:dim(CoordRes)[1]) {
				B<-paste(CoordRes[i,,j], sep="", collapse=" ")
				cat(B, file=FileName, sep="\n", append=TRUE)
			}
			cat(paste("IMAGE=", Meta$Files[j], sep=""), file=FileName, sep="\n", append=TRUE)
			cat(paste("ID=", Meta$ID[j], sep=""), file=FileName, sep="\n", append=TRUE)
			if (Scale==TRUE | !is.null(ScaleParam)) {cat(paste("SCALE=", Meta$Scale[j], sep=""), file=FileName, sep="\n", append=TRUE)}
		}
	}
	}
	
	#Export failure report
	write.table(ExtFail, paste(Output, "_SuccessReport.txt", sep=""), sep="\t", row.names=FALSE)
}

#########################################################################
# Plotting image, digitizing points along spiral                        #
# Necessary packages: pixmap, rtiff                                     #
# Necessary input variables:                                            #
#    Images: A list of all the images to be converted.                  #
#            *data frame* with image names in first column.             #
#    Output: Name of the output file (excluding extension).             #
#            *character*                                                #
#    Specimen.Labels: Names for the specimens. If NULL, numbers will... #
#                     be used.                                          #
#                     default=NULL                                      #
#    Guidelines: Shall guiders be drawn to enable to take take points...#
#                at equidistant spaces? Not recommended if measuring... #
#                structures that have natural guides (like chambers).   #
#                *logical*                                              #
#                TRUE=Guidelines will be drawn                          #
#                FALSE=No guidelines will be drawn                      #
#                default=TRUE                                           #
#    Density: Sampling density (degrees). Only meaningful if...         #
#             Guidelines=TRUE.                                          #
#             *numeric (integer)*                                       #
#             default=45                                                #
#    Equidistant: Shall points along spiral be equidistant?             #
#                 *logical*                                             #
#                 TRUE=points must be equidistant                       #
#                 FALSE=points do not need to be equidistant            #
#                 default=TRUE                                          #
#                 IMPORTANT: If FALSE is chosen points along the...     #
#                            spiral outline are allowed to be not...    #
#                            equidistant (whether they actually are...  #
#                            depends on the value of Density). This...  #
#                            could likely bias the results of the...    #
#                            further analysis!                          #
#    Normalize: Should the spirals be normalized for unit size before...#
#               coordinate export?                                      #
#               *logical*                                               #
#               TRUE: Normalize size to radius 1                        #
#               FALSE: Leave size as is                                 #
#               default=TRUE                                            #
#    Double: Should a double spiral be extracted, for instance on...    #
#            the external and internal side of the same shell?          #
#            *logical*                                                  #
#            default=FALSE                                              #
# Output data: File of type .spiral (effectively a .csv file without... #
#              row names) containing x- and y-coordinates and polar...  #
#              coordinates (distance from center t, angle theta in...   #
#              radians) for digitized points along spiral.              #
# Input dataset: Images in InputType format.                            #
#########################################################################

#Loading packages
require(stringr)
require(sfsmisc)
require(pixmap)
require(tiff)
require(jpeg)

SpiralExtraction<-function(Images, Output, Specimen.Labels=NULL, Guidelines=TRUE, Density=45, Equidistant=TRUE, Normalize=TRUE, Double=FALSE){
	#Test data for consistency
	if(Equidistant==TRUE){
		if(180%%Density!=0){stop(paste("180 cannot be divided by ", Density, " without remainder. Choose another value to get equidistant points along the spiral!", sep=""))}
	}
	{if (!is.null(Specimen.Labels)) {
		if (Specimen.Labels!="Names" && length(Specimen.Labels)!=nrow(Images)) {stop("Specimen.Labels must have same length as number of images!")}}
	}
	
	#Read image list into vector
	Images<-as.vector(Images[,1])
	InputType<-unlist(lapply(strsplit(Images, split=".", fixed=TRUE), tail, n=1L))
	InputType<-str_to_lower(InputType)
	if (!all(InputType%in%c("tif", "tiff", "ppm"))) {stop("Only '.ppm' and '.tif' format images accepted!")}
	
	#Set up temporary results objects
	Res<-list()
	if (Double==TRUE) {Res2<-list()}
	
	#Setup success report matrix
	ExtFail<-matrix(NA, length(Images), 1)
	colnames(ExtFail)<-c("Success")
	{if (is.null(Specimen.Labels)) {rownames(ExtFail)<-1:length(Images)}
	else if (Specimen.Labels=="Names") {rownames(ExtFail)<-Images}
	else {rownames(ExtFail)<-Specimen.Labels}}
	
	#Extract spiral
	for (k in 1:length(Images)) {
		print(k)
		Image<-Images[k]
		
		#Read image
		{if (InputType[k]=="ppm") {
			y<-read.pnm(Image)
		}
		else if (InputType[k]=="tif" | InputType[k]=="tiff") {
			y<-readTiff(Image)
		}}
		
		cont<-NA
		while(any(is.na(cont), cont=="n", (cont!="y" && cont!="c"))) {
			#Plot image
			par(mar=c(1, 1, 1, 1))
			plot(y)
			
			#Set midpoint for spiral
			writeLines("Please choose the center of the spiral. \nClick within the image.")
			flush.console()
			start<-locator(1)
			
			#Plot reference lines
			{if (Guidelines==TRUE) {
				abline(h=start$y, col="pink", lwd=2)
				if(90%%Density==0){abline(v=start$x, col="pink", lwd=2)}
				Ang<-0
				while(Ang<(180-Density)){
					Ang<-Ang+Density
					{if(Ang!=90){abline(a=start$y-start$x*tan(Ang*(pi/180)), b=tan(Ang*(pi/180)), col="pink", lwd=2)}
					else{}}
				}
			}
			else {
				points(start$x, start$y, pch=16, col="black")
			}}
			
			#Digitize spiral outline
			writeLines("Please digitize the points along the spiral. \nWhen you are finished right-click and choose stop.")
			flush.console()
			Coord<-locator(n=1000, type="o", lwd=2, pch=3, col="red")
			if (Double==TRUE) {
				writeLines("Please digitize the points along the other side of the spiral. \nWhen you are finished right-click and choose stop.")
				flush.console()
				Coord2<-locator(n=1000, type="o", lwd=2, pch=3, col="yellow")
			}
			
			bringToTop(-1)
			cont<-readline("Are the landmarks correct? (y=proceed, n=try again, c=cancel and proceed)")
		}
		
		#Save image with spiral on for later comparison
		dev.copy(png, filename=paste(strsplit(Image, split=".", fixed=TRUE)[[1]][1], "_Spiral.png", sep=""), width=5, height=5, units="in", res=200);
		dev.off();
		
		#Write success report
		{if (cont=="y") {ExtFail[k,1]<-1}
		else {ExtFail[k,1]<-0}}
		
		#Calculate lengths and angles of data
		Res[[k]]<-matrix(NA, length(Coord$x), 4)
		colnames(Res[[k]])<-c("x", "y", "t", "theta")
		if (cont=="y") {
			Res[[k]][,"x"]<-Coord$x
			Res[[k]][,"y"]<-Coord$y
			for (i in 1:length(Coord$x)) {
				Res[[k]][i,"t"]<-sqrt((Coord$x[i]-start$x)^2+(Coord$y[i]-start$y)^2)
				{if((Coord$x[i]-start$x)>=0 && (Coord$y[i]-start$y)>=0){Res[[k]][i,"theta"]<-atan(abs((Coord$y[i]-start$y))/abs((Coord$x[i]-start$x)))}
				else if((Coord$x[i]-start$x)<0 && (Coord$y[i]-start$y)>=0){Res[[k]][i,"theta"]<-pi-atan(abs((Coord$y[i]-start$y))/abs((Coord$x[i]-start$x)))}
				else if((Coord$x[i]-start$x)<0 && (Coord$y[i]-start$y)<0){Res[[k]][i,"theta"]<-pi+atan(abs((Coord$y[i]-start$y))/abs((Coord$x[i]-start$x)))}
				else {Res[[k]][i,"theta"]<-(2*pi)-atan(abs((Coord$y[i]-start$y))/abs((Coord$x[i]-start$x)))}
				}
			}
		}
		if (Double==TRUE) {
			Res2[[k]]<-matrix(NA, length(Coord2$x), 4)
			colnames(Res2[[k]])<-c("x", "y", "t", "theta")
			if (cont=="y") {
				Res2[[k]][,"x"]<-Coord2$x
				Res2[[k]][,"y"]<-Coord2$y
				for (i in 1:length(Coord2$x)) {
					Res2[[k]][i,"t"]<-sqrt((Coord2$x[i]-start$x)^2+(Coord2$y[i]-start$y)^2)
					{if((Coord2$x[i]-start$x)>=0 && (Coord2$y[i]-start$y)>=0){Res2[[k]][i,"theta"]<-atan(abs((Coord2$y[i]-start$y))/abs((Coord2$x[i]-start$x)))}
					else if((Coord2$x[i]-start$x)<0 && (Coord2$y[i]-start$y)>=0){Res2[[k]][i,"theta"]<-pi-atan(abs((Coord2$y[i]-start$y))/abs((Coord2$x[i]-start$x)))}
					else if((Coord2$x[i]-start$x)<0 && (Coord2$y[i]-start$y)<0){Res2[[k]][i,"theta"]<-pi+atan(abs((Coord2$y[i]-start$y))/abs((Coord2$x[i]-start$x)))}
					else {Res2[[k]][i,"theta"]<-(2*pi)-atan(abs((Coord2$y[i]-start$y))/abs((Coord2$x[i]-start$x)))}
					}
				}
			}
		}
		
		#Normalize spiral outline
		##Normalize size for radius=1
		if (Normalize==TRUE) {
			{if (Double==TRUE) {
				Size<-max(c(Res[[k]][,"t"], Res2[[k]][,"t"]))
				}
			else {Size<-max(Res[[k]][,"t"])}
			}
			Res[[k]][,"t"]<-Res[[k]][,"t"]/Size
			if (Double==TRUE) {
				Res2[[k]][,"t"]<-Res2[[k]][,"t"]/Size
			}
		}
		##Normalize rotation for start-point at radian=0
		Rotation<-Res[[k]][1,"theta"]
		for (i in 1:nrow(Res[[k]])) {
			{if (!is.na(Res[[k]][i,"theta"]) & Res[[k]][i,"theta"]>=Rotation) {Res[[k]][i,"theta"]<-Res[[k]][i,"theta"]-Rotation}
			else {Res[[k]][i,"theta"]<-(2*pi)+(Res[[k]][i,"theta"]-Rotation)}}
		}
		if (Double==TRUE) {
			for (i in 1:nrow(Res2[[k]])) {
				{if (!is.na(Res2[[k]][i,"theta"]) & Res2[[k]][i,"theta"]>=Rotation) {Res2[[k]][i,"theta"]<-Res2[[k]][i,"theta"]-Rotation}
				else {Res2[[k]][i,"theta"]<-(2*pi)+(Res2[[k]][i,"theta"]-Rotation)}}
			}
		}
	}
	
	#Prepare output file
	Res.Final<-matrix(NA, length(Res)*4, max(unlist(lapply(lapply(Res, dim), '[[', 1))))
	{if (is.null(Specimen.Labels)) {L<-1:length(Images)}
	else if (Specimen.Labels=="Names") {L<-unlist(lapply(strsplit(Images, split=".", fixed=TRUE), "[[", 1))}
	else {L<-Specimen.Labels}}
	rownames(Res.Final)<-paste(c("x", "y", "t", "theta"), rep(L, each=4), sep=".")
	for (i in 1:length(Res)) {
		Start.Line<-i+((i-1)*3)
		L<-nrow(Res[[i]])
		Res.Final[Start.Line,1:L]<-Res[[i]][,"x"]
		Res.Final[Start.Line+1,1:L]<-Res[[i]][,"y"]
		Res.Final[Start.Line+2,1:L]<-Res[[i]][,"t"]
		Res.Final[Start.Line+3,1:L]<-Res[[i]][,"theta"]
	}
	if (Double==TRUE) {
		Res.Final2<-matrix(NA, length(Res2)*4, max(unlist(lapply(lapply(Res2, dim), '[[', 1))))
		{if (is.null(Specimen.Labels)) {L<-1:length(Images)}
		else if (Specimen.Labels=="Names") {L<-unlist(lapply(strsplit(Images, split=".", fixed=TRUE), "[[", 1))}
		else {L<-Specimen.Labels}}
		rownames(Res.Final2)<-paste(c("x", "y", "t", "theta"), rep(L, each=4), sep=".")
		for (i in 1:length(Res2)) {
			Start.Line<-i+((i-1)*3)
			L<-nrow(Res2[[i]])
			Res.Final2[Start.Line,1:L]<-Res2[[i]][,"x"]
			Res.Final2[Start.Line+1,1:L]<-Res2[[i]][,"y"]
			Res.Final2[Start.Line+2,1:L]<-Res2[[i]][,"t"]
			Res.Final2[Start.Line+3,1:L]<-Res2[[i]][,"theta"]
		}
	}
	
	#Export results
	{if (Double==FALSE) {write.table(Res.Final, paste(Output, ".spiral", sep=""), sep=",", col.names=FALSE)}
	else {
		write.table(Res.Final, paste(Output, "_SideRed.spiral", sep=""), sep=",", col.names=FALSE)
		write.table(Res.Final2, paste(Output, "_SideYellow.spiral", sep=""), sep=",", col.names=FALSE)
	}
	}
	write.table(ExtFail, paste(Output, "_SuccessReport.txt", sep=""), sep="\t")
}

## EXAMPLES *************************************************************************************************

#Picking images
#MorphoImages<-ImagePicking(Dir="C:/R_TestData")

#Converting images
#ImageConversion(Images=MorphoImages, OutputType="ppm")
#ImageConversion(Images=MorphoImages, Scaling=50, OutputType="tif")

#Extract outlines
#setwd("C:/R_TestData/Outlines")
#MorphoImages<-ImagePicking(Dir=getwd())
#OutlineExtraction(Images=MorphoImages, Output="Spirula_Rep1", Specimen.Labels="Names", OutlinePoints=60, Smoothing=1, Scale=FALSE, RawVersion=TRUE)
#OutlineExtraction(Images=MorphoImages, Output="Spirula_Rep2", Specimen.Labels=paste("Spec", 1:3, sep="."), OutlinePoints=70, Smoothing=1, Scale=FALSE, RawVersion=FALSE)

#Extract landmarks
#setwd("C:/R_TestData/Landmarks")
#MorphoImages<-ImagePicking(Dir=getwd())
#LMExtract(Images=MorphoImages, Output="Landmarks1", N=5, Export="NTS")
#LMExtract(Images=MorphoImages, ScaleParam=rep(5.84818294686911, 3), Output="Landmarks2", Specimen.Labels=paste("Spec", 1:3, sep="_"), N=5, Export="NTS")
#LMExtract(Images=MorphoImages, Output="Landmarks3", Specimen.Labels="Names", Scale=FALSE, N=5, Export="TPS")

#Extract spiral morphology
#setwd("C:/R_TestData/Spirals")
#MorphoImages<-ImagePicking(Dir=getwd())
#SpiralExtraction(Images=MorphoImages, Output="SpiralForm1", Guidelines=TRUE, Density=90, Equidistant=TRUE, Normalize=TRUE)
#SpiralExtraction(Images=MorphoImages, Output="SpiralForm2", Specimen.Labels="Names", Guidelines=TRUE, Density=90, Equidistant=TRUE, Normalize=TRUE, Double=TRUE)

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
#		-Added function SpiralExtraction
#
# Version 1.1.1
#	Date: XXX
#	Description of changes:
#		-Added possibility to provide specimen labels manually in OutlineExtraction and LMExtract
#
# Version 1.1.2
#	Date: XXX
#	Description of changes:
#		-Numbering of specimens in Spiral.Extraction now based on start and stop number
#
# Version 1.2
#	Date: XXX
#	Description of changes:
#		-Added functionality to Spiral.Extraction to extract two parallel spirals and enhanced visuals
#
# Version 1.2.1
#	Date: XXX
#	Description of changes:
#		-Fixed an error where Spiral.Extraction would fail if StartNum was different from 1
#
# Version 1.2.2
#	Date: XXX
#	Description of changes:
#		-Fixed an error where OutlineExtraction would fail if StartNum was different from 1
#
# Version 1.3
#	Date: XXX
#	Description of changes:
#		-All extraction functions select images via an image-name list, eliminating the need to rename image files
#
# Version 2.0
#	Date: 19 August 2022
#	Description of changes:
#		-Complete rehaul of tiff and jpeg functionality due to discontinuation of old packages
##***********************************************************************************************************

#############################################################################################################






































