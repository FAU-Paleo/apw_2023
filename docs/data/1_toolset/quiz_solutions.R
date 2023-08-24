####################################################
# excerpt from Entry Quiz - Fossilworks Course, 2014
# by Michal Kowalewski, 2014.06.12

# solutions and related...
 
####################################################
# Part 1 - operation on vectors
# x is a vector with following values:
	#parentheses will always print
(  x<- c(4,3,5,6,6,2,-1)  )
 
#1.What is the value (or values) of y1?
(  y1<- x[3]  )
 
#2.What is the value (or values) of y2?
(  y2<- x[c(1,7)]  )
 
#3.What is the value (or values) of y3?
(  y3<- length(rep(x[2:3],x[1]))  )
 
#4.What is the value (or values) of y4?
(  y4<- x[which(x==6)]  )

# numeric subscripts
x[which(x==6)]

# logical subscript
x[x==6]


# with missing values they are different!
xNA <- x
xNA[1] <- NA

# numeric subscripts
xNA[which(xNA==6)]

# logical subscript
xNA[xNA==6]

 
#5.What is the value (or values) of y5?
(  y5<- sort(x)[1:2]  )
 
#6.What is the value (or values) of y6?
(  y6<- rank(-x, ties="first")[length(x)]  )

rank(-x)

noTies <- x[-5]
multiple <- c(x, 6)

rank(multiple)
rank(multiple, ties="first")


#7.What is the value (or values) of y7?
x1<- x[order(x,decreasing=TRUE)]
x2<- sort(x,decreasing=TRUE)
(  y7<- x1-x2  )
 
#8.What is the value (or values) of y8?
(  y8<- dim(cbind(x,x))  )
 
#9.What is the value (or values) of y9?
(  ifelse (sum(x>2)==5, y9<-2, y9<-1)  )
 
# How many true values do you have in a vector
sum(x>2)

# The same
any(x>2)
sum(x>2)>0

# are all the values TRUE?  
all(x>2)
sum(x>2)==length(x)

x>2

as.numeric(x>2)

#10. What is the value (or values) of y10?
(  y10<- length(sample(x,2,replace=TRUE))  )

sample(x,40,replace=TRUE)
sample(x,40)

# shuffling
sample(x,length(x))
a# same as
sample(x)

####################################################
# Part 2 - operation on matrices
# m1 is a matrix that looks like this:
# 1 3 4 0
# 3 1 4 0
# 1 3 4 2
# 3 1 4 2
 
typeof()
class()

typeof(1)
typeof(1.0)
typeof(1L)
typeof(TRUE)
typeof("a")
typeof(12+4i)


class(1)
class(1L)
typeof(m1)
class(m1)


(  m1<- cbind(c(1,3,1,3),c(3,1,3,1),c(4,4,4,4),c(0,0,2,2))  )
 
#11. What is the value (or values) of y11
(  y11<- m1[2:3,4]  )

is.vector(m1)
is.matrix(m1)
is.matrix(y11)
is.vector(y11)


(  y11<- m1[2:3,4, drop=FALSE]  )

a <- 2:3
b <- 4 
d <- m1[a, b]
d[,1 ]


a <- 2:3
b <- 4 
d <- m1[a, b, drop=FALSE]
d[,1 ]
 

#12. What is the value (or values) of y12
(  y12<- dim(m1[which(m1[,4]>0),])  )

#
a <- m1[,4]
b <- a>0
c <- which(b)
d <- m1[c, ]
dim(d)

 
#13. What is the value (or values) of y14
(  y13<- apply(m1,2,mean)  )

(  y13<- apply(X=m1,MARGIN=2,FUN=mean)  )

# set up a container with accurate structure
result <- rep(NA, ncol(m1))

for(i in 1:ncol(m1)){
	# focal column
	y <- m1[,i]
	# storing the current result
	result[i] <- mean(y)
}

apply(X=m1,MARGIN=2,FUN=min)  

#### Custom function
# Inline function definitions
apply(X=m1,MARGIN=2,FUN=function(x){
	length(unique(x))
})  

# define and store function beforehand
UniqueEntries <- function(x){
	length(unique(x))
}

# call to it
apply(X=m1,MARGIN=2,FUN=UniqueEntries)

1-> Row-wise
2-> column-wise


 
#14. What is the value (or values) of y15 (Hint: you do not need a calculator to solve this)
(  y14<- 1 + m1[1,4] * (sqrt(m1[2,3])+log(m1[3,3])+sum(m1)-exp(m1[2,1])/m1[3,4]**m1[1,1])  )
 

4*((((((((3+2))))))))

#15. What is the value (or values) of y16
(  y15<- sqrt(nrow(m1)*ncol(m1))  )


#16. What is the value (or values) of y19
	#[,2] is a categorical variable
(  y16<- tapply(X=m1[,1],INDEX=m1[,2],FUN=sum) )

table(m1[,2])
class(table(m1[,2]))
names(table(m1[,2]))
as.numeric(table(m1[,2]))

diff(table(m1[,2]))

# verbose description
X=m1[,1]
INDEX=m1[,2]
FUN=sum


X=c(5, 4, 3, 2, 1, 2, 3)
INDEX=c("a", "b", "a", "a", "c", "b", "b")
FUN=mean


# entries <- names(table(INDEX))
entries <- sort(unique(INDEX))
result <- rep(NA, length(entries))

for(i in 1:length(entries)){
	focal <- entries[i]
	ind <- which(INDEX == focal)
	result[i] <- FUN(X[ind])
}

names(result) <- entries


 
#17. What is the value (or values) of y20
m2<- m1
m2[m2==0] <- 1
(  y17<- sum(m2)-sum(m1)  )
 
 
####################################################
# Part 3 - programing statements
# matrix m1 and vector x will be used here again
 
#18. What is the value (or values) of y21
for (i in 1:100) {
  z<- mean(m1)
}
(  y18<- length(z)  )
 
# Known structure before iteration
set.seed(0)
trials <- 100
container <- rep(NA, trials )
for (i in 1:trials) {
	rand <- sample(1:10, 1)
	if(rand%%2==0){
	  container[i] <- rand 
	}
}

container[!is.na(container)]

# Unkown structure before iteration
set.seed(0)

container <- NULL
for (i in 1:trials) {
	rand <- sample(1:10, 1)
	if(rand%%2==0){
	  container <- c(container, rand)
	}
}



#19. What is the value (or values) of y19
v2<- NULL
for (j in 1:100) {
  v1<- cbind(j,nrow(m1))
  v2<- rbind(v2,v1)
}
(  y19<- v2[67,]  )


#20. What is the value (or values) of y20
(  y20<- x[-length(x)]-x[-1]  )


################################################################################
# define a single list
li <- list(a=c(TRUE, FALSE), b=1:4)

# list 'wrapper' still around
li[1]

# omitting the list 'wrapper'
li[[1]]

# acessing elements with names
li[["a"]]
# the same as
li$a

# a second list
li2 <- list(d=letters[1:5], e=4:8)

# combining lists on the same level
# c(<list>, <list>)
c(li, li2)

# a more complex, multi-level list
complex <- c(li, list(li2))

# digging deeper to access more 'buried' elements
complex[[1]]
complex[[3]][["d"]]


# example for digging values out of a complex
a<- rnorm(1:10)
b <- 1:10

mod <- lm(a ~ b )

str(mod)
mod$qr$pivot

var <- "qr"
mod[[var]]$pivot

