##################### Question 2
# function definitions:
inv_logit = function(x){exp(x)/(1+exp(x))}

getLogPost = function(X,m,y,b,beta.0,Sigm.0.inv){
	# returns a value that is proportional to the log posterior probability
	# requires data design matrix X, # trials m, response y, parameter vector or matrix to test b, and multivariate normal prior on the parameters with hyperparameters beta.0 and Sigma.0.inv
	b = matrix(b)
	pr = sapply(1:nrow(X),function(i){
		i = inv_logit(X[i,]%*%b)
		if(is.na(i)){return(1)} # without this, it returned NaN for very large i
		i
	})
	sum(log(pr**y*(1-pr)**(m-y)))-.5*t(b-beta.0)%*%Sigma.0.inv%*%(b-beta.0)
	
}

"bayes.logreg" <- function(m,y,X,beta.0,Sigma.0.inv,niter=10000,burnin=1000,print.every=1000,retune=100,verbose=TRUE){
	
	# create matrix to store samples. this is what function returns 
	p = ncol(X)
	sample_matrix = matrix(nrow=niter,ncol=p)
	
	# pick starting values:
	start = 1/sapply(1:p, function(p){
		mean(abs(X[,p]))
	})
	sample_matrix[1,] = start
	post_current = getLogPost(X,m,y,start,beta.0,Sigma.0.inv) 
	
	# find indices where tuning will happen
	tune_indices = seq(from=retune,to=burnin,by=retune)
	print(tune_indices)
	n_tunes = length(tune_indices)
	
	# create matrices for storing acceptance rates and proposal std. devs that we tried	
	prop_sd_matrix = matrix(nrow=(n_tunes+1),ncol=p)  # add one row for starting std. devs
	prop_sd_matrix[1,] = abs(start)/25
	accept_rate_matrix = matrix(nrow=(n_tunes),ncol=p)
	
	# the tuning iteration we're on 
	tune_iter = 1
	
	# set acceptances to 0 so we can count them and retune as needed
	accept = rep(0,p)
	
	# a logical vector identifying whether or not we've found a propsoal sd that gives 30-60% acceptance...if we have we skip tuning 
	is_done = rep(FALSE,times=p)
	
	# find indices where we'll print updates
	print_indices = floor(niter/print.every)
	
	for(i in 2:niter){
		# print updates?
		if(i %in% print_indices){
			print(paste('the chain is at index',i))
		}
		
		# tuning:
		if(i %in% tune_indices){
			accept_probs = accept/retune
			print(paste('with proposal sds:',prop_sd_matrix[tune_iter,]))
			print(paste('acceptance rate was:',accept_probs))
			
			# store acceptance rates 
			accept_rate_matrix[tune_iter,] = accept_probs
			for(par in 1:p){
				if(is_done[par]){ # don't bother if acceptance rate already OK
					prop_sd_matrix[(tune_iter+1),par] =prop_sd_matrix[tune_iter,par]
				} else{
					# if it's the first tuning, we don't have any previous sds we've tried
					if(tune_iter==1){
					low_max = 0
					high_min = 20 
					}else{
						accept_too_low = which(accept_rate_matrix[1:(tune_iter-1),par]<.3)
						accept_too_high = which(accept_rate_matrix[1:(tune_iter-1),par]>.6)
						sd_too_low = prop_sd_matrix[1:(tune_iter-1),par][accept_too_high]
						sd_too_high = prop_sd_matrix[1:(tune_iter-1),par][accept_too_low]
						low_max = max(0,sd_too_low)
						high_min = min(20,sd_too_high)
					}
					
					if(accept_probs[par]<.3){
						# acceptance rate too low, so decrease sd
						prop_sd_matrix[(tune_iter+1),par] =(prop_sd_matrix[tune_iter,par]+low_max)/2					
					}else if(accept_probs[par]>.6){
						# acceptance rate too high, so increase sd
						prop_sd_matrix[(tune_iter+1),par] =(prop_sd_matrix[tune_iter,par]+high_min)/2					
					}else{
						# acceptance rate just right, so don't change sd ...
						prop_sd_matrix[(tune_iter+1),par] =prop_sd_matrix[tune_iter,par]
						# ...and note that we already found a good sd for this parameter so we don't bother fiddling again
						is_done[par] = TRUE
					}
				}
				
						
			}
			
			# reset acceptance counters 
			accept = rep(0,p) 
			tune_iter = tune_iter+1
		} # tuning done 
		# iterate through parameters and update one at a time:
		for(j in 1:p){
			# get new values for parameters we've already changed
			new_beta = sample_matrix[i,0:(j-1)] 
			
			# proposed new value of parameter j
			propose = rnorm(n=1,mean=sample_matrix[(i-1),j],sd=prop_sd_matrix[tune_iter,j])
			
			# for the parameters we haven't changed yet, get previous values
			if(j==p){
				old_beta = numeric(0)
			}else{
				old_beta = sample_matrix[(i-1),(j+1):p] 
			}
			
			# find log posterior of proposal
			proposed = c(new_beta,propose,old_beta)
			post_proposed = getLogPost(X,m,y,proposed,beta.0,Sigma.0.inv)
			
			# find alpha...the 0 is there since log(1) = 0 and we're on log scale
			log_accept_prob = min(0,post_proposed-post_current)
			
			# generate random number and decide whether or not to accept
			rand = log(runif(n=1,min=0,max=1))
			if(rand<log_accept_prob){
				sample_matrix[i,j] = propose
				post_current = post_proposed
				accept[j] = accept[j]+1
			}else{
				sample_matrix[i,j] = sample_matrix[(i-1),j]
			}
		}
	
		
	}
	print(paste('after burnin, acceptance rate was:',accept/(niter-burnin)))
	sample_matrix
}

# actually fitting problem 2 on my test dataset: 
d = read.csv('/Users/matthewmeisner/Documents/R/STA250-Fall2013/blr_data_1100.csv')
X = as.matrix(d[,c('X1','X2')])
m = d$n
y = d$y
p = ncol(X)
beta.0 <- matrix(c(0,0))
Sigma.0.inv <- diag(rep(1.0,p))
samples_problem2 = bayes.logreg(m,y,X,beta.0,Sigma.0.inv,niter=10000)



######################### Question 3
# reading data - a bit of a pain
breast_cancer = read.table('/Users/matthewmeisner/Documents/Stuff/HW1/BayesLogit/breast_cancer.txt')
bc_names = breast_cancer[1,]
breast_cancer = breast_cancer[-1,]

# convert to numerics from factors
for(i in 1:10){
	breast_cancer[,i] = as.numeric(as.character(breast_cancer[,i]))
	names(breast_cancer)[i] = as.character(bc_names[1,i])
}
names(breast_cancer)[11] = as.character(bc_names[1,11])
head(breast_cancer)

X = as.matrix(breast_cancer[,1:10])
p = ncol(X)
y = ifelse(breast_cancer$diagnosis=='M',1,0)
beta.0 <- matrix(rep(0,times=p))
Sigma.0.inv <- diag(rep(1/1000,p))
niter <- 10000
m = rep(1,times=nrow(X))

# standardize 
X_std = X
for(i in 1:ncol(X_std)){
	X_std[,i] = (X_std[,i]-mean(X_std[,i]))/sd(X_std[,i])
}

# add column of 1s for the intercept
X_std = cbind(rep(1,nrow(X_std)),X_std)

# sample
samples = bayes.logreg(m,y,X_std,beta.0,Sigma.0.inv,niter=100000)

# remove burnin 
samples_no_burnin = samples[-(1:1000),]

# compute lag 1 autocorrelations
ac = sapply(1:ncol(samples_no_burnin), function(par){
	t = samples_no_burnin[2:nrow(samples_no_burnin),par]
	t_minus1 = samples_no_burnin[1:(nrow(samples_no_burnin)-1),par]
	cor(t,t_minus1)
})
xtable(matrix(round(ac,2)))

# see which covariates related with response
ci = sapply(1:ncol(samples_no_burnin), function(par){
	quantile(samples_no_burnin[,par],c(.025,.975))
})
xtable(round(t(ci),2))

# posterior predictive check
postPredCheck = function(n,sample_matrix,X,fxns,true_data,debug=FALSE){
	# this function performs posterior predictive checks
	# n:number of paramater samples to take from posterior
	# sample_matrix: is matrix of posterior samples; columns corrspond to number of parameters
	# X: design matrix of original raw data
	# fxns: list of functions (as character strings) to be applied to simulated and true data for comparison. separate plot will be generated for each function applied
	# true_data: numeric of observed response varialbe
	samples = sapply(1:n, function(i){
		row = sample(1:nrow(sample_matrix),size=1)
		b = matrix(sample_matrix[row,])
		logit_p = X%*%b
		p = inv_logit(logit_p)
		mal = rbinom(n=length(p),size=1,prob=p)
		sapply(fxns, function(fxn){
			get(fxn)(mal)
		})
	
	})
	if(debug){print(samples)}
	
	trues = lapply(fxns, function(fxn){
		get(fxn)(true_data)
	})
	
	for(i in 1:length(fxns)){
		dev.new()
		title = paste(fxns[[i]],'of Simulated Datasets from Posterior')
		if(length(fxns)==1){
			hist(samples,main=title)
		}else{
			hist(samples[i,],main=title)
		}
		abline(v=trues[[i]],col='red',lwd=3)
		legend('topright',legend='True Value',lty=1,col='red',lwd=2)
	}	
	
}


postPredCheck(1000,samples_no_burnin,X_std,list('mean','sd'),y,debug=F)
