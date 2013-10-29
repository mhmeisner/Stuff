
##
#
# Logistic regression
# 
# Y_{i} | \beta \sim \textrm{Bin}\left(n_{i},e^{x_{i}^{T}\beta}/(1+e^{x_{i}^{T}\beta})\right)
# \beta \sim N\left(\beta_{0},\Sigma_{0}\right)
#
##

library(mvtnorm)
library(coda)

########################################################################################
########################################################################################
## Handle batch job arguments:

# 1-indexed version is used now.
args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)

####
# sim_start ==> Lowest simulation number to be analyzed by this particular batch job
###

#######################
sim_start <- 1000
length.datasets <- 200
#######################

if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}

# Simulation datasets numbered 1001-1200

########################################################################################
########################################################################################



"bayes.logreg" <- function(m,y,X,beta.0,Sigma.0.inv,niter=10000,burnin=1000,print.every=1000,retune=100,verbose=TRUE){
	p = ncol(X)
	sample_matrix = matrix(nrow=niter,ncol=p)
	
	# randomly pick starting values:
	start = 1/sapply(1:p, function(p){
		mean(abs(X[,p]))
	})
	#start = mvrnorm(n=1,mu=beta.0,Sigma = solve(Sigma.0.inv)) # /100 worked well when not standardized 
	sample_matrix[1,] = start
	post_current = getLogPost(X,m,y,start,beta.0) 
	
	tune_indices = seq(from=retune,to=burnin,by=retune)
	print(tune_indices)
	n_tunes = length(tune_indices)
		
	prop_sd_matrix = matrix(nrow=(n_tunes+1),ncol=p)  # add one row for starting std. devs
	prop_sd_matrix[1,] = abs(start)/25
	
	accept_rate_matrix = matrix(nrow=(n_tunes),ncol=p)
	
	tune_iter = 1
	
	#prop_sd_matrix[1,] = c(.0002,.0025)**.5
	
	accept = rep(0,p)
	is_done = rep(FALSE,times=p)
	
	for(i in 2:niter){
		# tuning:
		if(i %in% tune_indices){
			accept_probs = accept/retune
			print(paste('with proposal sds:',prop_sd_matrix[tune_iter,]))
			print(paste('acceptance rate was:',accept_probs))
			accept_rate_matrix[tune_iter,] = accept_probs
			for(par in 1:p){
				if(is_done[par]){
					prop_sd_matrix[(tune_iter+1),par] =prop_sd_matrix[tune_iter,par]
				} else{
					if(tune_iter==1){
					low_max = 0
					high_min = 10 # uhhh, not sure what to pick here. depends a lot on parameters. 
					}else{
						accept_too_low = which(accept_rate_matrix[1:(tune_iter-1),par]<.3)
						accept_too_high = which(accept_rate_matrix[1:(tune_iter-1),par]>.6)
						sd_too_low = prop_sd_matrix[1:(tune_iter-1),par][accept_too_high]
						sd_too_high = prop_sd_matrix[1:(tune_iter-1),par][accept_too_low]
						low_max = max(0,sd_too_low)
						high_min = min(10,sd_too_high)
					}
					
					if(accept_probs[par]<.3){
						prop_sd_matrix[(tune_iter+1),par] =(prop_sd_matrix[tune_iter,par]+low_max)/2					
					}else if(accept_probs[par]>.6){
						prop_sd_matrix[(tune_iter+1),par] =(prop_sd_matrix[tune_iter,par]+high_min)/2					
					}else{
						# if we're already in range, we don't fiddle anymore 
						prop_sd_matrix[(tune_iter+1),par] =prop_sd_matrix[tune_iter,par]
						is_done[par] = TRUE
					}
				}
				
						
			}
			
			accept = rep(0,p) # reset acceptance counters 
			tune_iter = tune_iter+1
		}
		for(j in 1:p){
			new_beta = sample_matrix[i,0:(j-1)] 
			propose = rnorm(n=1,mean=sample_matrix[(i-1),j],sd=prop_sd_matrix[tune_iter,j])
			if(j==p){
				old_beta = numeric(0)
			}else{
				old_beta = sample_matrix[(i-1),(j+1):p] 
			}
				
			proposed = c(new_beta,propose,old_beta)
			
			post_proposed = getLogPost(X,m,y,proposed,beta.0)
			log_accept_prob = min(0,post_proposed-post_current)
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
	print(paste('acceptance rate was:',accept/(niter-burnin)))
	sample_matrix
}


getLogPost = function(X,m,y,b,beta.0){
	b = matrix(b)
	pr = sapply(1:nrow(X),function(i){
		inv_logit(X[i,]%*%b) 
	})
	sum(log(pr**y*(1-pr)**(m-y)))-.5*t(b-beta.0)%*%Sigma.0.inv%*%(b-beta.0)
	
}


#################################################
# Set up the specifications:
beta.0 <- matrix(c(0,0))
Sigma.0.inv <- diag(rep(1.0,2))
niter <- 10000
# etc... (more needed here)
#################################################

# Read data corresponding to appropriate sim_num:
filepath = paste('~/Stuff/HW1/BayesLogit/data/blr_data_',sim_num,'.csv',sep='')
d = read.csv(filepath)

# Extract X and y:
X = as.matrix(d[,c('X1','X2')])
m = d$n
y = d$y

# Fit the Bayesian model:
samples = bayes.logreg(m,y,X,beta.0,Sigma.0.inv)

# Extract posterior quantiles...
quantiles = sapply(1:ncol(samples), function(par){
	quantile(samples[,par],(1:99)/100)
})


# Write results to a (99 x p) csv file...
out_filepath = paste('~/Stuff/HW1/BayesLogit/data/blr_res_',sim_num,'.csv',sep='')
write.csv(quantiles,file=out_filepath)

# Go celebrate.
 
cat("done. :)\n")






