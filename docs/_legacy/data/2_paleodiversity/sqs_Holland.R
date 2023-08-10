# If you use this in a published work, please cite it as an internet resource. The format will vary
# depending on your journal, but it will look something like this:
# Holland, S. M. 2015. Estimating diversity with SQS. Accessed: <insert date you downloaded it>. 
#   URL: http://strata.uga.edu/8370/rtips/shareholderQuorumSubsampling.html
# Also, be sure to cite John Alroy's paper that introduced the method:
# Alroy, J. 2010. Geographical, environmental and intrinsic biotic controls on Phanerozoic 
#   marine diversification. Palaeontology 53:1211-1235.

sqs <-function(abundance, quota=0.9, trials=100, ignore.singletons=FALSE, exclude.dominant=FALSE) {
	# abundance is a vector of integers representing the abundance of every species

	if ((quota <= 0 || quota >= 1)) {
		stop('The SQS quota must be greater than 0.0 and less than 1.0')
	}

	# compute basic statistics
	specimens <- sum(abundance)
	numTaxa <- length(abundance)
	singletons <- sum(abundance==1)
	doubletons <- sum(abundance==2)
	highest <- max(abundance)
	mostFrequent <- which(abundance==highest)[1]

	if (exclude.dominant == FALSE) {
		highest <- 0
		mostFrequent <- 0
	}

	# compute Good's u
	u <- 0
	if (exclude.dominant == TRUE) {
		u <- 1 - singletons / (specimens - highest)
	} else {
		u <- 1 - singletons / specimens
	}

	if (u == 0) {
		stop('Coverage is zero because all taxa are singletons')
	}

	# re-compute taxon frequencies for SQS
	frequencyInitial <- abundance - (singletons + doubletons / 2) / numTaxa
	frequency <- frequencyInitial / (specimens - highest)

	# return if the quorum target is higher than estimated coverage
	if ((quota > sum(frequency)) || (quota >= sum(abundance))) {
		stop('SQS quota is too large, relative to the estimated coverage')
	}

	# create a vector, length equal to total number of specimens,
	#  each value is the index of that species in the abundance array
	ids <- unlist(mapply(rep, 1:numTaxa, abundance))

	# subsampling trial loop
	richness <- rep(0, trials)                 # subsampled taxon richness
	for (trial in 1:trials) {
		pool <- ids                            # pool from which specimens will be sampled
		specimensRemaining <- length(pool)     # number of specimens remaining to be sampled
		seen <- rep(0, numTaxa)                # keeps track of whether taxa have been sampled
		subsampledFrequency <- rep(0, numTaxa) # subsampled frequencies of the taxa
		coverage <- 0

		while (coverage < quota) {
			# draw a specimen
			drawnSpecimen <- sample(1:specimensRemaining, size=1)
			drawnTaxon <- pool[drawnSpecimen]

			# increment frequency for this taxon
			subsampledFrequency[drawnTaxon] <- subsampledFrequency[drawnTaxon] + 1

			# if taxon has not yet been found, increment the coverage
			if (seen[drawnTaxon] == 0) {
				if (drawnTaxon != mostFrequent && (ignore.singletons == 0 || abundance[drawnTaxon] > 1)) {
					coverage <- coverage + frequency[drawnTaxon]
				}
				seen[drawnTaxon] <- 1

				# increment the richness if the quota hasn't been exceeded,
				# and randomly throw back some draws that put the coverage over quota
				if (coverage < quota || runif(1) <= frequency[drawnTaxon]) {
					richness[trial] <- richness[trial] + 1
				} else {
					subsampledFrequency[drawnTaxon] <- subsampledFrequency[drawnTaxon] - 1
				}
			}

			# decrease pool of specimens not yet drawn
			pool[drawnSpecimen] <- pool[specimensRemaining]
			specimensRemaining <- specimensRemaining - 1
		}
	}

	# compute subsampled richness
	s2  <- richness[richness>0]
	subsampledRichness <- exp(mean(log(s2))) * length(s2)/length(richness)
	return(round(subsampledRichness, 1))
}
