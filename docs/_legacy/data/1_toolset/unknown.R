TheThing <- function(one, two, three, da, n=1000){
	base<- runif(n, 0, 2*22/7)

	ps <- data.frame(one=cos(base), two=sin(base))*three
	ps$one <- ps$one+one
	ps$two <- ps$two+two

	# add scatter
	jitter <- rnorm(n, mean=0, sd=three*da)
	ps$one <- ps$one + jitter
	ps$two <- ps$two + sample(jitter)

	return(ps)

}
