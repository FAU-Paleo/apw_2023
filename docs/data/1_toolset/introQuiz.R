####################################################
# excerpt from Entry Quiz - Fossilworks Course, 2014
# by Michal Kowalewski, 2014.06.12

 
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
 
#5.What is the value (or values) of y5?
(  y5<- sort(x)[1:2]  )
 
#6.What is the value (or values) of y6?
(  y6<- rank(-x, ties="first")[length(x)]  )
 
#7.What is the value (or values) of y7?
x1<- x[order(x,decreasing=TRUE)]
x2<- sort(x,decreasing=TRUE)
(  y7<- x1-x2  )
 
#8.What is the value (or values) of y8?
(  y8<- dim(cbind(x,x))  )
 
#9.What is the value (or values) of y9?
(  ifelse (sum(x>2)==5, y9<-2, y9<-1)  )
 
#10. What is the value (or values) of y10?
(  y10<- length(sample(x,2,replace=TRUE))  )
 
####################################################
# Part 2 - operation on matrices
# m1 is a matrix that looks like this:
# 1 3 4 0
# 3 1 4 0
# 1 3 4 2
# 3 1 4 2
 
(  m1<- cbind(c(1,3,1,3),c(3,1,3,1),c(4,4,4,4),c(0,0,2,2))  )
 
#11. What is the value (or values) of y11
(  y11<- m1[2:3,4]  )
 
#12. What is the value (or values) of y12
(  y12<- dim(m1[which(m1[,4]>0),])  )
 
#13. What is the value (or values) of y14
(  y13<- apply(m1,2,mean)  )
 
#14. What is the value (or values) of y15 (Hint: you do not need a calculator to solve this)
(  y14<- 1 + m1[1,4] * (sqrt(m1[2,3])+log(m1[3,3])+sum(m1)-exp(m1[2,1])/m1[3,4]**m1[1,1])  )
 
#15. What is the value (or values) of y16
(  y15<- sqrt(nrow(m1)*ncol(m1))  )

#16. What is the value (or values) of y19
	#[,2] is a categorical variable
(  y16<- tapply(m1[,1],INDEX=m1[,2],FUN=sum) )
 
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
 
#19. What is the value (or values) of y22
v2<- NULL
for (j in 1:100) {
  v1<- cbind(j,nrow(m1))
  v2<- rbind(v2,v1)
}
(  y19<- v2[67,]  )

#20. What is the value (or values) of y25
(  y20<- x[-length(x)]-x[-1]  )

