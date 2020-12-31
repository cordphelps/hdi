


visualizeBeta <- function(mean, mode, shape1, shape2) {

	# https://bookdown.org/content/3686/inferring-a-binomial-probability-via-exact-mathematical-analysis.html#specifying-a-beta-prior.
	mu <- mean
	omega <- mode 
	kappa <- shape1 + shape2


  	data.tbl <- tibble(shape1 = shape1, shape2 = shape2) %>% 

  	mutate(	a        = str_c("a = ", shape1),
         	b        = str_c("b = ", shape2),
         	kappa    = kappa,
         	mu_omega = omega) %>% 

  	expand(nesting(shape1, shape2, a, b, kappa, mu_omega), 
         x = seq(from = 0, to = 1, length.out = 100)) 



	gg <- ggplot(data=data.tbl, aes(x = x)) +

  		geom_vline(xintercept = .8, color = "white") +
  		geom_line(aes(y = dbeta(x, shape1 = shape1, shape2 = shape2)),
            color = "grey50", size = 1.25) +
  
  
  		scale_x_continuous(expression(theta), breaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1)) +
  		ylab(expression(p(theta*"|"*a*", "*b))) +
  		coord_cartesian(ylim = c(0, 5))

  	rtn <- list()
  	rtn[[1]] <- gg
  	rtn[[2]] <- data.tbl 

  	return(rtn)

}

bernoulli_likelihood <- function(theta, data) {

	# https://bookdown.org/content/3686/inferring-a-binomial-probability-via-exact-mathematical-analysis.html#the-likelihood-function-the-bernoulli-distribution
  
  	# `theta` = success probability parameter ranging from 0 to 1
  	# `data`  = the vector of data (i.e., a series of 0s and 1s)

  	n <- length(data)
  	z <- sum(data)
  
  	return(theta^z * (1 - theta)^(n - z))
  
}

#########################################
# Shape parameters from central tendency and scale:
#########################################

betaABfromMeanSD <- function(mean, sd) {
	# https://bookdown.org/content/3686/inferring-a-binomial-probability-via-exact-mathematical-analysis.html#the-likelihood-function-the-bernoulli-distribution
  if (mean <= 0 | mean >= 1) stop("must have 0 < mean < 1")
  if (sd <= 0) stop("sd must be > 0")
  kappa <- mean * (1 - mean)/sd^2 - 1
  if (kappa <= 0) stop("invalid combination of mean and sd")
  a <- mean * kappa
  b <- (1.0 - mean) * kappa
  return(list(a = a, b = b))
}


betaABfromMeanKappa <- function(mean, kappa) {
  if (mean <= 0 | mean >= 1) stop("must have 0 < mean < 1")
  if (kappa <= 0) stop("kappa must be > 0")
  a <- mean * kappa
  b <- (1.0 - mean) * kappa
  return(list(a = a, b = b))
}

betaABfromModeKappa <- function(mode, kappa) {
	# https://bookdown.org/content/3686/inferring-a-binomial-probability-via-exact-mathematical-analysis.html

	# "effective sample size" kappa = shape1  + shape2

  if (mode <= 0 | mode >= 1) stop("must have 0 < mode < 1")
  if (kappa <= 2) stop("kappa must be > 2 for mode parameterization")
  a <- mode * (kappa - 2) + 1
  b <- (1.0 - mode) * (kappa - 2) + 1
  return(list(a = a, b = b))
}

hdi_of_icdf <- function(name, width = .95, tol = 1e-8, ... ) {

  # https://bookdown.org/content/3686/inferring-a-binomial-probability-via-exact-mathematical-analysis.html

  # Arguments:
  #   `name` is R's name for the inverse cumulative density function
  #   of the distribution.
  #   `width` is the desired mass of the HDI region.
  #   `tol` is passed to R's optimize function.
  # Return value:
  #   Highest density iterval (HDI) limits in a vector.
  # Example of use: For determining HDI of a beta(30, 12) distribution, type
  #   `hdi_of_icdf(qbeta, shape1 = 30, shape2 = 12)`
  #   Notice that the parameters of the `name` must be explicitly stated;
  #   e.g., `hdi_of_icdf(qbeta, 30, 12)` does not work.
  # Adapted and corrected from Greg Snow's TeachingDemos package.
  
  incredible_mass <- 1.0 - width
  interval_width <- function(low_tail_prob, name, width, ...) {
    name(width + low_tail_prob, ...) - name(low_tail_prob, ...)
  }
  
  opt_info <- optimize(interval_width, c(0, incredible_mass), 
                       name = name, width = width, 
                       tol = tol, ...)
  
  hdi_lower_tail_prob <- opt_info$minimum
  
  return(c(name(hdi_lower_tail_prob, ...),
           name(width + hdi_lower_tail_prob, ...)))
  
}


overlap <- function(start1, end1, start2, end2) {

	# determine if there is overlap between two numeric intervals

	# https://stackoverflow.com/questions/36035074/how-can-i-find-an-overlap-between-two-given-ranges

	if (start1 > end2 || end1 < start2) {
    	# // no overlap
    	return(FALSE)
	}
	else {
    	# // overlap
    	return(TRUE)
	}

}

hdi.num <- function(theta, prior.a, prior.b, n, z, index) {

	# HDIs for the posterior trial 
	hdiPostA <- hdi_of_icdf(name = qbeta, 
                        	shape1 = z + prior.a,
                        	shape2 = n - z + prior.b, 
                        	width = theta)

	if (index == 1) {
		return(hdiPostA[[1]])
	} else {
		return(hdiPostA[[2]])
	}

}



compareObservations <- function(prior.data, trialA.data, trialB.data, beginTheta) {

# https://web.stanford.edu/class/cs109/reader/11%20Parameter%20Estimation.pdf

	preObsTraps <- length(prior.data)            # 30
	preObsSpiders <- sum(prior.data)              #  5
	preObsMean <- preObsSpiders / preObsTraps
	#ci <- .95

  
	# define the prior
	essPrior <- betaABfromMeanKappa(mean = preObsMean, kappa = preObsTraps)


	data <- tibble(theta = seq(from = beginTheta, to = 0.99, by=0.01)) %>%

								 add_column(priorParamA = 0) %>% add_column(priorParamB = 0) %>% 
								 add_column(nTrialA = 0) %>% add_column(zTrialA = 0) %>% 
								 add_column(hdiLowerTrialA = 0) %>% add_column(hdiUpperTrialA = 0) %>% 
								 add_column(nTrialB = 0) %>% add_column(zTrialB = 0) %>% 
								 add_column(hdiLowerTrialB = 0) %>% add_column(hdiUpperTrialB = 0) %>%
								 add_column(overlap=TRUE) %>% add_column(postTrialA=0, postTrialB=0) 

	# mutate() will not work with non-vectorized hdi.num() and overlap()

	for (i in 1:nrow(data)) {


		data$priorParamA[[i]] <- essPrior$a 

		data$priorParamB[[i]] <- essPrior$b

		data$nTrialA[[i]] <- length(trialA.data)

		data$zTrialA[[i]] <- sum(trialA.data)

		data$hdiLowerTrialA[[i]] <- hdi.num(theta = data$theta[[i]],
						 							prior.a = data$priorParamA[[i]], 
													prior.b = data$priorParamB[[i]], 
													n = data$nTrialA[[i]],
													z = data$zTrialA[[i]],
													index = 1)

		data$hdiUpperTrialA[[i]] <- hdi.num(theta = data$theta[[i]],
						 							prior.a = data$priorParamA[[i]], 
													prior.b = data$priorParamB[[i]], 
													n = data$nTrialA[[i]],
													z = data$zTrialA[[i]],
													index = 2)

		data$nTrialB[[i]] <- length(trialB.data)

		data$zTrialB[[i]] <- sum(trialB.data)

		data$hdiLowerTrialB[[i]] <- hdi.num(theta = data$theta[[i]],
						 							prior.a = data$priorParamA[[i]], 
													prior.b = data$priorParamB[[i]], 
													n = data$nTrialB[[i]],
													z = data$zTrialB[[i]],
													index = 1)

		data$hdiUpperTrialB[[i]] <- hdi.num(theta = data$theta[[i]],
						 							prior.a = data$priorParamA[[i]], 
													prior.b = data$priorParamB[[i]], 
													n = data$nTrialB[[i]],
													z = data$zTrialB[[i]],
													index = 2)

		data$overlap[[i]] <- overlap(data$hdiLowerTrialA[[i]], data$hdiUpperTrialA[[i]],
									 data$hdiLowerTrialB[[i]], data$hdiUpperTrialB[[i]])


		data$postTrialA[[i]] <- dbeta( data$theta[[i]], 
                                       shape1 = data$zTrialA[[i]] + data$priorParamA[[i]], 
                                       shape2 = data$nTrialA[[i]] - data$zTrialA[[i]] + data$priorParamB[[i]])

		data$postTrialB[[i]] <- dbeta( data$theta[[i]], 
                                       shape1 = data$zTrialB[[i]] + data$priorParamA[[i]], 
                                       shape2 = data$nTrialB[[i]] - data$zTrialB[[i]] + data$priorParamB[[i]])



	}




	return(data)

}


chopObservations <- function(sourceData, chopFactor) {

	# chop a 270 observation dataset into 10 gradually increasing subsets
	# chopFactor = 27 

	# https://www.statmethods.net/management/subset.html

	rtn <- list()

	 for (i in 1:(length(sourceData)/chopFactor)) {

		rtn[[i]]  <- sourceData[c(1:(chopFactor*i))]

	}

	return(rtn)

}

buildRandomPriorData <- function(sourceVector, returnSampledRows) {

	# sourceVector <- c(1,2,3,4,5,6,7,8,9,0)
	# returnSampledRows <- 31

	columns <- 1

	set.seed(999)

	tibble.tbl <- tibble( sampledData = replicate(columns, sample(sourceVector, replace=TRUE, size=returnSampledRows)) )

	return(tibble.tbl)
}

makeDistributions <- function(label, prior.data, trial.data) {

	#=============================================================
	# build prior, likelihood, and posterior

	# Determining the Effective Sample Size of a Parametric Prior
	# # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3081791/
	# effective sample size (ESS) = a + b = 'kappa'


	## --- READ READ READ https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3081791/ --- ##
	#      "If the prior is elicited from a domain expert, then an informative prior is desirable"
	# for the prior, examine the first 30 observations; 5 positives
	# 5/30 = 0.167

	kappa <- 6   # it just does from spider observations (binomial distribution shape factor)

	essPrior <- betaABfromMeanKappa((mean = sum(prior.data) / length(prior.data)), kappa = kappa)

	# https://bookdown.org/content/3686/inferring-a-binomial-probability-via-exact-mathematical-analysis.html
	# paragraph 6.2.1 
	modePrior <- (essPrior$a -1) / (essPrior$a + essPrior$b -2)

	distributions.tbl <- tibble(theta = seq(from = 0, to = 1, length.out = 100)) %>% 

                  			mutate(prior = dbeta(theta, 
                                       shape1 = essPrior$a, 
                                       shape2 = essPrior$b)) %>% 

                 			mutate(like = bernoulli_likelihood(
  										theta = theta, data = trial.data)) %>%

                  			mutate(post = dbeta(theta, 
                                       shape1 = sum(trial.data) + essPrior$a, 
                                       shape2 = length(trial.data) - sum(trial.data) + essPrior$b))


    modePost <- distributions.tbl$theta[which.max(distributions.tbl$post)]

    if (label=="SNH") {
    	greenColor <- '#66a182' # green-ish
    	postLabel <- 'green-ish'
    } else {
    	greenColor <- '#d1495d' # red-ish
    	postLabel <- 'red-ish'
    }

	gg <- ggplot(data=distributions.tbl, aes(x = theta)) +
  

  			geom_ribbon(aes(ymin = 0, ymax = as.numeric(prior)), fill = "blue", alpha=0.4) +

  			geom_ribbon(aes(ymin = 0, ymax = as.numeric(like)), fill = "red", alpha=0.4) +

			geom_ribbon(aes(ymin = 0, ymax = as.numeric(post)), fill = greenColor, alpha=0.6) +

			geom_vline(aes(xintercept = modePrior), linetype = 2, color = "blue") +
  			geom_vline(aes(xintercept = modePost), linetype = 2, color = "black") +

  			theme_bw() +

  			#scale_y_continuous(limits=c(0, 10), breaks=c(0, 4, 8)) +
  			scale_x_continuous(limits=c(0, 1), breaks=c(0, .25, .5, .75)) +

  			ylim(0, 15) +

  			labs(	x = expression(theta),
         		 	y = paste( expression(p(theta)), sep=""), 
 					title = paste("beta distribution, ", label, " trial", 
 								  "\nprior (blue): ", length(prior.data), " sampled observations, mode= ", round(modePrior, digits=2),
 								  "\nposterior (", postLabel, "): ", length(trial.data), 
 								  " actual observations, mode= ", round(modePost, digits=2), sep=""))


    return(list(distributions.tbl, gg))
}


ggPair <- function(prior, trial, labelText) {

	# plot prior and posterior distributions for control and SNH

	 rtnControl.lst <- makeDistributions(
            label = labelText[[1]],
            prior.data = prior[[1]],
            trial.data = trial[[1]])
  
  	rtnSNH.lst <- makeDistributions(
            label = labelText[[2]],
            prior.data = prior[[2]],
            trial.data = trial[[2]])

  	return(list(rtnControl.lst[[2]],rtnSNH.lst[[2]]))

}


findMode <- function(tibble, prior, trial, labelText, index) {

	# plot prior and posterior distributions for control and SNH

	# from makeDistributions()



	kappa <- 6   # it just does from spider observations (binomial distribution shape factor)

	essControlPrior <- betaABfromMeanKappa((mean = sum(prior[[1]]) / length(prior[[1]])), kappa = kappa)
	essSNHPrior <- betaABfromMeanKappa((mean = sum(prior[[2]]) / length(prior[[2]])), kappa = kappa)

	# https://bookdown.org/content/3686/inferring-a-binomial-probability-via-exact-mathematical-analysis.html
	# paragraph 6.2.1 
	modeControlPrior <- (essControlPrior$a -1) / (essControlPrior$a + essControlPrior$b -2)
	modeSNHPrior <- (essSNHPrior$a -1) / (essSNHPrior$a + essSNHPrior$b -2)

	temp.tbl <- tibble(theta = seq(from = 0, to = 1, length.out = 100)) %>% 

                  			mutate(postControl = dbeta(theta, 
                                       shape1 = sum(trial[[1]]) + essControlPrior$a, 
                                       shape2 = length(trial[[1]]) - sum(trial[[1]]) + essControlPrior$b)) %>%

							mutate(postSNH = dbeta(theta, 
                                       shape1 = sum(trial[[2]]) + essSNHPrior$a, 
                                       shape2 = length(trial[[2]]) - sum(trial[[2]]) + essSNHPrior$b))

    modeControlPost <- temp.tbl$theta[which.max(temp.tbl$postControl)]
    modeSNHPost <- temp.tbl$theta[which.max(temp.tbl$postSNH)]

    rm(temp.tbl)

    # https://stackoverflow.com/questions/61444678/replace-values-in-tibble-in-r-4-0
    tibble[index, 2:5] <- as.list(c(modeControlPrior, modeSNHPrior, modeControlPost, modeSNHPost))

  	return(tibble)

}





