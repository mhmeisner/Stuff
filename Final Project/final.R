##############
library(lme4)
library(blme)
library(gplots)
library(rethinking)

############## Gibbs Sampler for Random Effects Only Model 
updateTheta = function(n_grps,grp_sizes,grp_means,mu,lambda_theta,lambda_e){
	sapply(1:n_grps,function(grp){
		n_i = grp_sizes[grp]
		y_bar_i = grp_means[grp]
		
		mean = (lambda_theta*mu+n_i*lambda_e*y_bar_i)/(lambda_theta+n_i*lambda_e)
		sd = 1/sqrt(lambda_theta+n_i*lambda_e)
				
		rnorm(1,mean,sd)
	})	
}

updateMu = function(n_grps,lambda_theta,s0,m0,thetas){
	mn = (s0*m0+n_grps*lambda_theta*mean(thetas))/(s0+n_grps*lambda_theta)
	sd = 1/sqrt(s0+n_grps*lambda_theta)
	rnorm(1,mn,sd)
}

updateLambdaTheta = function(n_grps,a1,b1,thetas,mu){
	shape = n_grps/2+a1
	rate = sum((thetas-mu)**2)/2+b1
	rgamma(1,shape,rate)
}

updateLambdaE = function(n,n_grps,a2,b2,y,grps,thetas){
	ss_theta_by_grp = sapply(1:n_grps,function(i){
		y_i = y[grps==i]
		theta_i = thetas[i]
		sum((y_i-theta_i)**2)
	})
 	ss_theta_devs = sum(ss_theta_by_grp)
	shape = n/2+a2
	rate = ss_theta_devs/2+b2
	rgamma(1,shape,rate)
}

justRanefGibbsSampler = function(y,grps,a1,b1,a2,b2,m0,s0,niter,burnin=0){
	# grps need to be 1-n_grps; grps should be vector of same length as y that contains group number 
	# a1, b1, a2, b2 are priors for data-level and group-level precision, respectively
	# m0, s0 are prior mean and precision for normal prior on mu (mean of random effects distribution)
	n = length(y)
	n_grps = max(grps)
	grp_sizes = table(grps)
	grp_means = aggregate(y,list(grps),mean)[,2]
	
	# initialize parameters in "reasonable" starting ranges: 
	mu = rnorm(1,mean(y),sd(y))
	lambda_e = runif(1,0,sd(y))
	lambda_theta = runif(1,0,sd(y))
	
	# allocate matrix to store samples in
	sample_matrix = matrix(nrow=niter,ncol=n_grps+3)
	
	# run the sampler! 
	for(i in 1:niter){
		thetas = updateTheta(n_grps,grp_sizes,grp_means,mu,lambda_theta,lambda_e)
		mu = updateMu(n_grps,lambda_theta,s0,m0,thetas)
		lambda_theta = updateLambdaTheta(n_grps,a1,b1,thetas,mu)
		lambda_e = updateLambdaE(n,n_grps,a2,b2,grp_sizes,y,grps,thetas)
		sample_matrix[i,] = c(thetas,mu,lambda_theta,lambda_e)
 	}
	sample_matrix[(1:burnin):niter,]
}

### compare to expected results from blmer
# generate random, clustered data:
grps = sample(1:20,1000,replace=T)
y = rnorm(1000,grps,1)

# fit with my code
s = justRanefGibbsSampler(y,grps,0.01,0.01,0.01,0.01,0,0,10000)
m = colMeans(s[5000:10000,])

# fit with blmer, which we know (hopefully) should work:
m1 = blmer(y~(1|grps))
coef(m1)

# compare 
plot(m[1:20],coef(m1)$grps[,1],xlab='Random Effects: Posterior Means from My GS',ylab='Correct blmer Estimates',main='Ran. Ef. Estimates from My Code vs. Correct Values')

############## Gibbs Sampler for Mixed Model With One Random Effect
makeRanefIncidenceMatrix = function(grps){
	#grps is a numeric of length n indicating which group each data point is in.  Groups MUST be numbered 1 thru n_grps
	# returns a n x n_grps incidence matrix indicating which group each data point is in
	n = length(grps)
	n_grps = max(grps)
	Z = matrix(0,nrow=n,ncol=n_grps)
	for(i in 1:n){
		Z[i,grps[i]] = 1
	}
	Z
}
updateTheta = function(X,Z,D,beta,lambda_e){
	M = solve(t(Z)%*%Z+solve(D)/lambda_e)
	mn = M%*%t(Z)%*%(y-X%*%beta)
	cov = M/lambda_e
	as.matrix(mvrnorm(1,mn,cov))
}

updateBeta = function(X,Z,theta,y,lambda_e){
	mn = solve(t(X)%*%X)%*%t(X)%*%(y-Z%*%theta)
	cov = solve(t(X)%*%X)/lambda_e 
	as.matrix(mvrnorm(1,mn,cov))
}

updateLambdaTheta = function(n_grps,a1,b1,theta){
	shape = n_grps/2+a1
	rate = b1+.5*t(theta)%*%theta
	rgamma(1,shape,rate)
}

updateLambdaE = function(n,a2,b2,y,X,beta,Z,theta){
	shape = n/2+a2
	rate = b2+.5*t(y-X%*%beta-Z%*%theta)%*%(y-X%*%beta-Z%*%theta)
	rgamma(1,shape,rate)
}

ranefGibbsSampler = function(y,grps,a1,b1,a2,b2,X,niter,burnin=0){
	# grps need to be 1-n_grps
	n = length(y)
	p = ncol(X)
	n_grps = max(grps)
	
	# make incidence matrix
	Z = makeRanefIncidenceMatrix(grps)
	
	# pick starting values
	lambda_e = runif(1,0,sd(y[,1]))
	lambda_theta = runif(1,0,sd(y[,1]))
	beta = matrix(mvrnorm(1,matrix(0,nrow=p),diag(1,p)))
	D = diag(1,n_grps)/lambda_theta # D is covariance matrix for the random effects 
	
	# run sampler! 
	sample_matrix = matrix(nrow=niter,ncol=n_grps+p+2) # ranef, fixef, and both variances
	for(i in 1:niter){
		theta = updateTheta(X,Z,D,beta,lambda_e)
		beta = updateBeta(X,Z,theta,y,lambda_e)
		lambda_theta = updateLambdaTheta(n_grps,a1,b1,theta)
		D = diag(1,n_grps)/lambda_theta
		lambda_e = updateLambdaE(n,a2,b2,y,X,beta,Z,theta)
		sample_matrix[i,] = c(t(beta),t(theta),lambda_theta,lambda_e)
 	}
	sample_matrix[(burnin+1):niter,]
}

# check:
n = 1000
p = 3
grps = sample(1:5,n,replace=T)
X = cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
b = matrix(c(3,4,6,2))
y = matrix(rnorm(n,X%*%b+matrix(grps),1))
s = ranefGibbsSampler(y,grps,0.01,0.01,0.01,0.01,X,10000)
m = colMeans(s[5000:1000,])

plot(s[,1])
acf(s[,5])
cor(s[1000:2000,2],s[1001:2001,2]) # var and ranef highly corrr, but not others 

blmer_test = blmer(y~X[,1]+X[,2]+X[,3]+X[,4]+(1|grps)-1)
summary(blmer_test)
ranef(blmer_test)


############## Gibbs Sampler for Mixed Model With Two Random Effects
updateTheta1 = function(y,X,Z1,Z2,D1,beta,lambda_e,theta_2){
	M = solve(t(Z1)%*%Z1+solve(D1)/lambda_e)
	mn = M%*%t(Z1)%*%(y-X%*%beta-Z2%*%theta_2)
	cov = M/lambda_e
	#print(cov)
	as.matrix(mvrnorm(1,mn,cov))
}

updateTheta2 = function(y,X,Z1,Z2,D2,beta,lambda_e,theta_1){
	M = solve(t(Z2)%*%Z2+solve(D2)/lambda_e)
	mn = M%*%t(Z2)%*%(y-X%*%beta-Z1%*%theta_1)
	cov = M/lambda_e
	as.matrix(mvrnorm(1,mn,cov))
}

updateBeta = function(X,Z1,Z2,theta_1,theta_2,y,lambda_e){
	mn = solve(t(X)%*%X)%*%t(X)%*%(y-Z1%*%theta_1-Z2%*%theta_2)
	cov = solve(t(X)%*%X)/lambda_e 
	as.matrix(mvrnorm(1,mn,cov))
}

updateLambdaTheta1 = function(n_grps_1,a1,b1,theta_1){
	shape = n_grps_1/2+a1
	rate = b1+.5*t(theta_1)%*%theta_1
	rgamma(1,shape,rate)
}

updateLambdaTheta2 = function(n_grps_2,a2,b2,theta_2){
	shape = n_grps_2/2+a2
	rate = b2+.5*t(theta_2)%*%theta_2
	rgamma(1,shape,rate)
}

#r = rgamma(10000,2,3)
#mean(r)

updateLambdaE = function(n,a3,b3,y,X,beta,Z1,Z2,theta_1,theta_2){
	shape = n/2+a3
	rate = b3+.5*t(y-X%*%beta-Z1%*%theta_1-Z2%*%theta_2)%*%(y-X%*%beta-Z1%*%theta_1-Z2%*%theta_2)
	rgamma(1,shape,rate)
}



twoRanefGibbsSampler = function(y,grps1,grps2,X,a1,b1,a2,b2,a3,b3,niter,burnin=0,debug=F){
	# grps1 and grps2 contain IDs for two sources of clustering, need to be 1-n_grps1 and 1-n_grps2; all other parameters as before 
	n = length(y)
	if(debug){print(n)}
	p = ncol(X)
	if(debug){print(dim(X))}
	n_grps_1 = max(grps1)
	n_grps_2 = max(grps2)
		
	Z1 = makeRanefIncidenceMatrix(grps1)
	Z2 = makeRanefIncidenceMatrix(grps2)
	
	lambda_e = runif(1,0,sd(y[,1]))
	lambda_theta_1 = runif(1,0,sd(y[,1]))
	lambda_theta_2 = runif(1,0,sd(y[,1]))
	beta = matrix(mvrnorm(1,matrix(0,nrow=p),diag(1,p)))
	theta_2 = matrix(rnorm(n_grps_2))
	
	D1 = diag(1,n_grps_1)/lambda_theta_1 # cov matrixes for random effects
	D2 = diag(1,n_grps_2)/lambda_theta_2
	
	sample_matrix = matrix(nrow=niter,ncol=n_grps_1+n_grps_2+p+3) # ranef, fixef, and all three variances
	for(i in 1:niter){
		theta_1 = updateTheta1(y,X,Z1,Z2,D1,beta,lambda_e,theta_2)
		theta_2 = updateTheta2(y,X,Z1,Z2,D2,beta,lambda_e,theta_1)
		beta = updateBeta(X,Z1,Z2,theta_1,theta_2,y,lambda_e)
		lambda_theta_1 = updateLambdaTheta1(n_grps_1,a1,b1,theta_1)
		lambda_theta_2 = updateLambdaTheta2(n_grps_2,a2,b2,theta_2)
		D1 = diag(1,n_grps_1)/lambda_theta_1
		D2 = diag(1,n_grps_2)/lambda_theta_2
		lambda_e = updateLambdaE(n,a3,b3,y,X,beta,Z1,Z2,theta_1,theta_2)
		sample_matrix[i,] = c(t(beta),t(theta_1),t(theta_2),lambda_theta_1,lambda_theta_2,lambda_e)
 	}
	sample_matrix[(burnin+1):niter,]
}

### fit to my data, predicting yield 
cc = complete.cases(d[,c('actual_yield',matrix_names)])

# convert weird labels to 1-n_grps
fields = as.numeric(droplevels(d$field_num_ucd[cc]))
years = as.numeric(as.factor(d$crop_year[cc]))

# and make sure that it worked! 
table(years,d$crop_year[cc])

# set up matrices and fit mdodel 
X = cbind(rep(1,sum(cc)),as.matrix(d[cc,matrix_names]))
y = as.matrix(d$actual_yield[cc])
s = twoRanefGibbsSampler(y,fields,years,X,0.001,0.001,0.001,0.001,0.001,0.001,2000,0)

# make trace plots:
par(mfrow = c(3,3))
for(i in 2:10){
	plot(s[,i],type='l')

# now for pests
cc = complete.cases(d[,c('may_june_total_insects',matrix_names)])
sum(cc)
fields = as.numeric(droplevels(d$field_num_ucd[cc]))
years = as.numeric(as.factor(d$crop_year[cc]))
X1 = X = cbind(rep(1,sum(cc)),as.matrix(d[cc,matrix_names]))
y1 = as.matrix(d$may_june_total_insects[cc])
s2 = twoRanefGibbsSampler(y1,fields,years,X=X1,0.001,0.001,0.001,0.001,0.001,0.001,2000,0,debug=T)
par(mfrow = c(3,3),mar=c(2,2,2,2))
for(i in 2:10){
	plot(s2[,i],type='l')
}

# remove burnin (not done in fxn so that we could see convergence time in trace plot)
s2s = s2[-(1:500),]
ss = s[-(1:500),]

# compare results to "true" values:
m = blmer(actual_yield~(1|field_num_ucd)+(1|crop_year)+alfalfa_matrix+barley_matrix+carrots_matrix+corn_matrix+cotton_matrix+gbeans_matrix+garlic_matrix+lettuce_matrix+melons_matrix+onions_matrix+potatoes_matrix+safflower_matrix+sugarbeets_matrix+tomatoes_matrix+wheat_matrix, data=d) 
precis(m)
plot(colMeans(ss)[2:16],fixef(m)[2:16],xlab='Posterior Means from My Code',ylab='Estimates from blmer',main='My MCMC Estimates vs blmer Estimates for Fixed Effects')

# extract means, medians, and 95% HPDIs for parameters of interest: 
yield_effects = sapply(c(2:16),function(i){
	diffs = ss[,i]
	c(mean(diffs),median(diffs),HPDI(diffs))
})

pest_effects = sapply(c(2:16),function(i){
	diffs = s2s[,i]
	c(mean(diffs),median(diffs),HPDI(diffs))
})

### make plots!
# convert to better units:
1 bales/acre * (217 kg/bale)*(1/.404686 acres/ha) = 
con = 217.7243/.404686 # units are (kg*)/(bale*)

ints = yield_effects
colnames(ints) =c('Alfalfa','Barley','Carrots','Corn','Cotton','Garb. Beans','Garlic','Lettuce','Melons',"Onions",'Potatoes','Safflower','Sugarbeets','Tomatoes','Wheat')
o = order(ints[1,], decreasing=T)
ints1 = ints[,o]
up = ints1[4,]
down = ints1[3,]

cps = c('Alfalfa','Barley','Carrots','Corn','Cotton','Garb. Beans','Garlic','Lettuce','Melons',"Onions",'Potatoes','Safflower','Sugarbeets','Tomatoes','Wheat')[o]

par(mar=c(4,6,4,1),mfrow=c(1,2))
plotCI(y=1:15,x=ints1[1,]*con,ui=up*con,li=down*con,pch=19,yaxt='n',ylab='',err='x',xlab='',xlim=c(-400,200),cex=.7,xaxt='n')
abline(v=0,lty=2)
axis(side=2,at = 1:15,labels = cps,las=2,xpd=T)
mn_yield = mean(d$actual_yield[w],na.rm=T)*con
tickk= seq(from=-400,to=200,by=200)
axis(side=3,at = tickk,labels = round(tickk/mn_yield,digits =2))
axis(side=1,at = tickk,labels = round(tickk))
text(-100,17.5,labels='Proportion Yield Difference from Cotton',xpd=T)
text(-100,-1.5, labels = 'Yield Difference from Cotton (kg/ha)',xpd=T)
text(-400,17.5, labels = 'A',xpd=T,cex=2)
sigs = c(1,2,3,7,13) # from bottom
xsig = up[sigs]*con
text(xsig+21,sigs,labels = '*',cex=2,)

colnames(pest_effects) =c('Alfalfa','Barley','Carrots','Corn','Cotton','Garb. Beans','Garlic','Lettuce','Melons',"Onions",'Potatoes','Safflower','Sugarbeets','Tomatoes','Wheat')
ints2 = pest_effects[,o]

up2 = ints2[4,]
down2 = ints2[3,]
par(mar=c(4,0,4,8))
plotCI(y=1:15,x=ints2[1,],ui=up2,li=down2,pch=19,yaxt='n',ylab='',err='x',xlab='',cex=.7,xlim=c(-.5,1),xaxt='n')
abline(v=0,lty=2)
mn_pest = mean(d$may_june_total_insects[w],na.rm=T)
tickk= seq(from=-.5,to=1,by=.5)
axis(side=1,at = tickk,labels = tickk)
axis(side=3,at = tickk,labels = round(tickk/mn_pest,digits =2))
text(0.25,17.5,labels='Proportion Pest Difference from Cotton',xpd=T)
text(0.25,-1.5, labels = 'Pest Difference from Cotton (insects/sweep)',xpd=T)
sigs = c(1,2,5,6,7,10,12,14)
xsig = up2[sigs]
text(xsig+.05,sigs,labels = '*',cex=2,)
text(-.5,17.5, labels = 'B',xpd=T,cex=2)


############## Bayesian Linear Regression 
updateBeta = function(sigma_0_inv,beta_0,X,y,lambda){
	cov = solve(sigma_0_inv+lambda*t(X)%*%X)
	mn = cov%*%(sigma_0_inv%*%beta_0+lambda*t(X)%*%y)	
	return(as.matrix(mvrnorm(1,mn,cov)))
}

updateLambda = function(a,b,n,X,beta,y){
	shape = a+n/2
	rate = b+.5*t(y-X%*%beta)%*%(y-X%*%beta)
	rgamma(1,shape,rate)
}

bayesianLinearRegression = function(y,X,beta_0,sigma_0_inv,a,b,niter,burnin=0,debug=F){
	# beta_0, sigma_0_inv are priors for beta
	# a, b are priors for gamma prior for precision
	n = nrow(X)
	p = ncol(X)
	lambda = runif(1,0,1)
	sample_matrix = matrix(nrow=niter,ncol=p+1) # add one for variance 
	for(i in 1:niter){
		if(debug){cat('\n current iteration is: ',i)}
		beta = updateBeta(sigma_0_inv,beta_0,X,y,lambda)
		lambda = updateLambda(a,b,n,X,beta,y)
		sample_matrix[i,] = c(beta,lambda)
	}
	sample_matrix[(burnin+1):niter,]	
}

# fit model:
y = as.matrix(yield_effects[1,])
X = cbind(rep(1,ncol(yield_effects)),as.matrix(pest_effects[1,]))
beta_0 = matrix(1,nrow=2)
sigma_0_inv = diag(.001,2)
blr = bayesianLinearRegression(y,X,beta_0,sigma_0_inv,0.001,0.001,10000,5000)
colMeans(blr)
tail(blr)

# get whole posterior, integrating over uncertainty in effects from previous model: 
beta_0 = matrix(0,nrow=2)
sigma_0_inv = diag(.001,2)

samples = lapply(1:1000,function(i){
	print(i)
	ye = matrix(ss[i,2:16])
	pe = matrix(s2s[i,2:16])
	X = cbind(rep(1,nrow(pe)),pe)
	blr = bayesianLinearRegression(ye,X,beta_0,sigma_0_inv,0.001,0.001,10000,5000)
	blr[,2]
})

all = unlist(samples)
mean(all)
median(all)
HPDI(all)
sum(all<=0)/length(all)
plot(density(all),xlab='Slope',ylab='Posterior Probability')
dens(all,show.HPDI=0.95,xlab='Slope',ylab='Posterior Probability',main='Posterior Distribution of Slope of Yield Effects vs. Pest Effects')

### make animated posterior samples: 
for(i in 1:100){
	filepath = paste("/Users/matthewmeisner/Documents/LaTeX/Homework/STA 250/figures/samp_effects",i,".pdf", sep = '')
	pdf(file=filepath)
	plot(x=s2s[i,2:16],y=ss[i,2:16], col ="black",pch=19,xlab='Effect on Pests',ylab='Effect on Yield',main='Effects of Surrounding Crops on Yield and Pest Densities',xlim=c(-.3,.7),ylim=c(-.2,.3))
	dev.off()
}

beta_0 = matrix(1,nrow=2)
sigma_0_inv = diag(.001,2)
samples2 = lapply(1:100,function(i){
	print(i)
	ye = matrix(ss[i,2:16])
	pe = matrix(s2s[i,2:16])
	X = cbind(rep(1,nrow(pe)),pe)
	blr = bayesianLinearRegression(ye,X,beta_0,sigma_0_inv,0.001,0.001,1000,500)
	list(blr[,1],blr[,2])
})

intercept_samples = unlist(lapply(samples2,function(i){i[[1]]}))
slope_samples = unlist(lapply(samples2,function(i){i[[2]]}))
HPDI(slope_samples)
HPDI(intercept_samples)
b0 = mean(intercept_samples)
b1 = mean(slope_samples)
pest_seq = seq(from = -.4, to = .8, length.out = 100)
int = sapply(pest_seq, function(i){
	HPDI(intercept_samples+i*slope_samples)
})


for(i in 1:100){
	filepath = paste("/Users/matthewmeisner/Documents/LaTeX/Homework/STA 250/figures/slopes",i,".pdf", sep = '')
	pdf(file=filepath)
	curve(b0+b1*x,lwd=5,xlab='Effect on Pests (insects/sweep)',ylab='Effect on Yield (bales/ac)',main='Effects of Surrounding Crops on Yield and Pest Densities',xlim=c(-.3,.7),ylim=c(-.2,.3),cex.lab = 1.5, col = "#08306B")
	lines(pest_seq,int[1,],lty=2,col="#9ECAE1",lwd=4)
	lines(pest_seq,int[2,],lty=2,col="#9ECAE1",lwd=4)
	points(x=pest_effects[1,],y=yield_effects[1,], col ="black",pch=19)
	b2 = intercept_samples[i]
	b3 = slope_samples[i]
	curve(b2+b3*x,lwd=4,lty=2,add = T,col="#08519C")
	dev.off()
}
