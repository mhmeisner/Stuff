############ Load Libraries
library(RCUDA)
library(MASS)
library(xtable)

############ Load Kernel
m = loadModule("rtruncnorm2.ptx")
rtk = m$rtruncnorm_kernel
source('~/utility.R')

############ CPU rtruncnorm function.  
rtruncnormCPU = function(n,mu,sigma,lo,hi,maxtries=100,verbose=FALSE){
	# n: integer
	# mu, sigma, lo, and hi all need to be vectors of length n 
	sapply(1:n,function(i){
		ntries = 0
		while(ntries < maxtries){ # no need for accepted = T/F; just return the current sample when we get one we like 
			ntries <- ntries+1
			ran = rnorm(1,mu[i],sigma[i])
			if(ran>=lo[i]&ran<=hi[i]){
				#print(ntries)
				return(ran)
			}
		}
		#### Use the Robert sampling method if that didn't work:
		# note: Robert method as implemented here won't work for bounded intervals...truncation region needs to be from -Inf to a value or from a value to Inf
		if(verbose){cat('Using Robert Sampling Method \n')}
		# we will only sample from z~N(0,1,lower,Inf), so have to find c,k,and lower s.t. cz+k ~ N(mu,sigma**2,lo[i],Inf)...possibly with sign reversed if we're sampling from N(mu,sigma_sq,-Inf,hi[i]), since my code requires lower to be positive 
		# cz+k ~ N(k,c**2,k+c*lower,Inf)
		# c = 1, since sigma**2 always equals 1 for our purposes 
		# k = mu[i]
		# if we need to sample from (-Inf,b)
		reverse_sign = FALSE
		if(hi[i]==Inf){
			lower = lo[i]-mu[i]
		}else{
			lower = mu[i]-hi[i] # -(hi[i]-mu[i])
			reverse_sign = TRUE
		}
		#print(lower)
		alpha = (lower+sqrt(lower**2+4))/2
		accepted = FALSE
		while(!accepted){
			z = lower+rexp(1,alpha)
			if(lower<alpha){
				psi = exp(-(alpha-z)**2/2)
			}else{
				psi = exp((lower-alpha)**2/2)*exp(-(alpha-z)**2/2)
			}
			u = runif(1,0,1)
			if(u<psi){accepted=TRUE}
		}
		if(reverse_sign){
			sample = mu[i]-z
		}else{
			sample = mu[i]+z
		}
		return(sample)			
	})
}

############ GPU rtruncnorm function.  
rtruncnormGPU = function(n,mu,sigma,lo,hi,maxtries=2000L,curand_seq = 1L,debug=F){
	N = as.integer(n)
	bg <- compute_grid(n)
	grid_dims <- bg$grid_dims
	block_dims <- bg$block_dims
	if(debug){cat('grid dims are \n',grid_dims,'\n block dims are \n',block_dims,'\n')}
	x = numeric(n)
	if(debug){cat('Calling .cuda \n')}
	res = .cuda(rtk,'x'=x,N,mu,sigma,lo,hi,maxtries,1L,1L,curand_seq, blockDim = block_dims,gridDim= grid_dims,outputs='x')
	if(debug){cat('.cuda is done! \n')}
	res

}

############ see if empirical means from both functions match expected values (they do!):
mu = 2
sigma = 1
a = 0
b = 1.5
cpu_samples = rtruncnormCPU(n=10000,mu=rep(mu,10000),sigma=rep(sigma,10000),lo=rep(a,10000),hi=rep(b,10000))
gpu_samples = rtruncnormGPU(n=10000,mu=rep(mu,10000),sigma=rep(sigma,10000),lo=rep(a,10000),hi=rep(b,10000),debug=T)
mean(cpu_samples)
mean(gpu_samples)
expected = mu+sigma*(dnorm((a-mu)/sigma)-dnorm((b-mu)/sigma))/(pnorm((b-mu)/sigma)-pnorm((a-mu)/sigma))
expected

############ make sure both functions work for -Inf and Inf 
cpu_neg_inf = rtruncnormCPU(1000,rep(0,1000),rep(1,1000),rep(-Inf,1000),rep(-1000,1000))
gpu_neg_inf = rtruncnormGPU(1000,rep(0,1000),rep(1,1000),rep(-Inf,1000),rep(-1000,1000))
head(cpu_neg_inf)
head(gpu_neg_inf)
range(cpu_neg_inf)
range(gpu_neg_inf)
cpu_inf = rtruncnormCPU(1000,rep(0,1000),rep(1,1000),rep(1000,1000),rep(Inf,1000))
gpu_inf = rtruncnormGPU(1000,rep(0,1000),rep(1,1000),rep(1000,1000),rep(Inf,1000))
head(cpu_inf)
head(gpu_inf)

############ make sure both function work for tail regions: 
cpu_tail = rtruncnormCPU(1000,rep(0,1000),rep(1,1000),rep(-Inf,1000),rep(-10,1000))
gpu_tail = rtruncnormGPU(1000,rep(0,1000),rep(1,1000),rep(-Inf,1000),rep(-10,1000))
head(cpu_tail)
head(gpu_tail)
cpu_tail = rtruncnormCPU(1000,rep(0,1000),rep(1,1000),rep(-Inf,1000),rep(-100000,1000))
gpu_tail = rtruncnormGPU(1000,rep(0,1000),rep(1,1000),rep(-Inf,1000),rep(-100000,1000))
head(cpu_tail)
head(gpu_tail)


############ function to compare CPU and GPU runtimes for rtruncnrom generation 
compareCPUGPU = function(k,mu,sigma,lo,hi){
	# k should be numeric; will time generation 10*kk samples for each kk in k 
	# mu, sigma, lo and hi should all be numerics (floats); the same value for each will be used for generating each random numbers (not the best way to do it, in retrospect, but works for what we wanted to use this function for!)
	# returns runtimes of CPU and GPU for each k provided
	sapply(k,function(kk){
		cat('current k: ',kk,'\n')
		n = 10**kk
		CPU_time = system.time({res = rtruncnormCPU(n,rep(mu,n),rep(sigma,n),rep(lo,n),rep(hi,n))})
		GPU_time = system.time({res = rtruncnormGPU(n,rep(mu,n),rep(sigma,n),rep(lo,n),rep(hi,n))})
		c(CPU_time['elapsed'],GPU_time['elapsed'])
	})
}

############ compare runtimes and make plots
times = compareCPUGPU(1:7,2,1,0,1.5)
# on AWS, write outfile of runtimes:
write.table(times,'~/compare_times.csv') 
# to make graphs on my computer:
system('scp -i mhmeisner.pem ec2-user@ec2-54-203-145-38.us-west-2.compute.amazonaws.com:~/compare_times.csv ~/Documents')
t = read.table('~/Documents/compare_times.csv',row.names=NULL)[,-1]
par(mfrow=c(1,2))
plot(1:7,t[1,],type='o',xlab='Log_10(n)',ylab='Time (seconds)',main='CPU and GPU Runtimes for TN Sampling')
lines(1:7,t[2,],type='o',col='red',lty=2,pch=2)
legend(1.5,450,legend=c('CPU','GPU'),lty=c(1,2),pch=c(1,2),col=c('black','red'))
plot(1:4,t[1,1:4],type='o',xlab='Log_10(n)',ylab='Time (seconds)',main='CPU and GPU Runtimes for TN Sampling')
lines(1:4,t[2,1:4],type='o',col='red',lty=2,pch=2)
legend(1.25,.24,legend=c('CPU','GPU'),lty=c(1,2),pch=c(1,2),col=c('black','red'))

# output table:
xtable(t)

############ Probit MCMC functon (has an optin to use CPU or GPU)
probitMCMC = function(y,X,beta_0,sigma_0_inv,niter,burnin,useGPU=FALSE,debug=FALSE){
	# niter is total draws, including the burn in 
	n = as.numeric(length(y))
	p = ncol(X)
	sample_matrix = matrix(nrow=niter,ncol=p)
	
	# cov of conditional dist. of beta; find once since it's constant 
	cond_b_cov = solve(sigma_0_inv+t(X)%*%X) 
	
	# create upper/lower regions for TN:
	lower = ifelse(y==1,0,-Inf)
	upper = ifelse(y==1,Inf,0)
	
	#beta_start = mvrnorm(n=1,mu=beta_0,Sigma=solve(sigma_0_inv))
	current_b = beta_0
	sample_matrix[1,] = current_b
	for(i in 2:niter){
		print(i)
		if(debug){print(i)}
		# sample the z's:
		mu = as.numeric(X%*%matrix(current_b))
		if(debug){print(sum(mu>lower&upper==Inf))}
		if(useGPU){
			if(debug){
				cat('n is:',n,'\n')
				cat('mu is:',mu,'\n')
				cat('lower is:',lower,'\n')
				cat('upper is:',upper,'\n')
			}
			current_z = rtruncnormGPU(n,mu,rep(1,n),lower,upper,curand_seq =as.integer(i),debug=debug)
			
			# sometimes, my rtruncnormGPU function wasn't able to generate a sample in a reasonable number of tries (even using the Robert method). This probably means that I've done something wrong. But, if failed *very* rarely, and rather than have the whole code fail due to the NaN returned from rtruncnormGPU, I just reverted to the CPU function, which always seemed to work.  
			if(NaN %in% current_z){
				cat('GPU sampler returned NA; reverting to CPU method for this iteration \n')
				current_z = rtruncnormCPU(n,mu,rep(1,n),lower,upper)
			}
		}else{
			current_z = rtruncnormCPU(n,mu,rep(1,n),lower,upper)
		}
		
		current_z = matrix(current_z)
		
		# sample the b's
		current_b = mvrnorm(1,cond_b_cov%*%(sigma_0_inv%*%beta_0+t(X)%*%current_z),cond_b_cov)
		
		# store the b samples:  
		sample_matrix[i,] = current_b 
	}
	return(sample_matrix[(1+burnin):niter,])
}

############ function to run probit MCMC using both CPU and GPU methods, for various datasets, and to make a graph with posterior samples and true parameter values overlaid. Number of iterations can differ for each dataset. Runtimes of CPU and GPU are the outputs.
fitProbit = function(data_filepaths,par_filepaths,figure_outfilepaths,iter,burnin,debug=FALSE){
	runtimes = matrix(nrow=2,ncol=length(data_filepaths))
	for(i in 1:length(data_filepaths)){
		data_filepath = data_filepaths[i]
		par_filepath = par_filepaths[i]
		figure_outfilepath = figure_outfilepaths[i]
		d = read.table(data_filepath,header=T)	
		y = d[,1]
		X = as.matrix(d[,2:ncol(d)])
		beta_0 = matrix(rep(0,ncol(X)))
		sigma_0_inv = matrix(0,nrow=ncol(X),ncol=ncol(X))
		if(debug){cat('running CPU MCMC')}
		cpu_time = system.time({cpu_samples = probitMCMC(y,X,beta_0,sigma_0_inv,iter[i],burnin[i],debug=debug)})		
		if(debug){cat('running GPU MCMC')}
		gpu_time = system.time({gpu_samples = probitMCMC(y,X,beta_0,sigma_0_inv,iter[i],burnin[i],useGPU=T,debug=debug)})		
		runtimes[1,i] = cpu_time['elapsed']
		runtimes[2,i] = gpu_time['elapsed']
		correct_pars = read.table(par_filepath,header=T)[,1]
		print(correct_pars)
		pdf(figure_outfilepath)
		par(mfrow=c(4,2),mar=c(1,2,4,1))
		for(p in 1:8){
			plot(density(cpu_samples[,p],na.rm=T),main=paste('Parameter',p),xlab='')
			print(gpu_samples[,p])
			lines(density(gpu_samples[,p],na.rm=T),lty=2,col='red')
			abline(v=correct_pars[p])
			legend('topleft',legend=c('CPU Posterior Samples','GPU Posterior Samples'),lty=c(1,2),col=c('black','red'),cex=.5)
		}
		dev.off()
	}	
	return(runtimes)
}

############  make plots comparing runtimes of CPU and GPU for probitMCMC:
t = fitProbit(c('mini_data.txt','data_01.txt','data_02.txt','data_03.txt','data_04.txt','data_05.txt'),c('mini_pars.txt','pars_01.txt','pars_02.txt','pars_03.txt','pars_04.txt','pars_05.txt'),c('mini_plot.txt','plot_01.txt','plot_02.txt','plot_03.txt','plot_04.txt','plot_05.txt'),c(2000,2000,2000,500,50,25),c(500,500,500,50,5,3))
par(mfrow=c(1,2))
plot(2:7,t[1,],type='o',xlab='Log_10(n)',ylab='Time (seconds)',main='CPU and GPU Runtimes for Probit MCMC')
lines(2:7,t[2,],type='o',col='red',lty=2,pch=2)
legend('topleft',legend=c('CPU','GPU'),lty=c(1,2),pch=c(1,2),col=c('black','red'))
plot(2:5,t[1,1:4],type='o',xlab='Log_10(n)',ylab='Time (seconds)',main='CPU and GPU Runtimes for Probit MCMC')
lines(2:5,t[2,1:4],type='o',col='red',lty=2,pch=2)
legend('topleft',legend=c('CPU','GPU'),lty=c(1,2),pch=c(1,2),col=c('black','red'))
xtable(t)
