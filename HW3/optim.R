bisect = function(f,initial_interval,tolerance=.1,max_iter=1000,debug=FALSE){
	tolerance = abs(tolerance) # in case someone gave a negative number 
	lower = initial_interval[1]
	upper = initial_interval[2]
	if(f(lower)*f(upper)>0){stop('Not Necessarily a Root in Initial Interval')}
	if(f(lower)*f(upper)==0){
		root = ifelse(f(lower)==0,lower,upper)
		stop(root,' is already a root!')
	}
	
	iter = 1
	converged = FALSE
	
	while(!converged&(iter<=max_iter)){
		if(debug){cat('Starting Iteration',iter,'\n','Current Interval:',lower,upper,'\n')}
		mid = (lower+upper)/2
		value_at_mid = f(mid)
		if(debug){
			cat('Currently checking for root at:',mid,'\n','The value of f here is:',value_at_mid,'\n')
		}
		if(abs(value_at_mid)<=tolerance){
			root = mid
			converged =TRUE
		}else{
			if(f(lower)*f(mid)>0){
				lower = mid
			}else{
				upper = mid
			}
		}
		iter <- iter+1
	}
 	if(converged){
 		cat('Converged!\n Converged in',iter-1,'iterations \n Root is:',root,'\n f at the root is:',value_at_mid,'\n')
 		return(root)
 	}else{
 		cat(':( :( :( :( :( \n','Failed to converge in',iter-1,'iterations \n Last interval:',c(lower,upper),'\n')
 		return(NULL)
 	}
}

newtonRaphson = function(f,f_prime,start,tolerance=.1,max_iter=1000,debug=FALSE){
	tolerance = abs(tolerance) # in case someone gave a negative number 
	if(abs(f(start))<tolerance){stop(start,' is already a root (within specified tolerance)!')}
	
	iter = 1
	converged = FALSE
	x = start
	
	while(!converged&(iter<=max_iter)){
		if(debug){cat('Starting Iteration',iter,'\n')}
		x <- x-f(x)/f_prime(x)
		val = f(x)
		if(debug){
			cat('Currently checking for root at:',x,'\n','The value of f here is:',val,'\n')
		}
		if(abs(val)<=tolerance){
			root = x
			converged = TRUE
		}
		iter <- iter+1
	}
	
	if(converged){
 		cat('Converged!\n Converged in',iter-1,'iterations \n Root is:',x,'\n f at the root is:',val,'\n')
 		return(root)
 	}else{
 		cat(':( :( :( :( :( \n','Failed to converge in',iter-1,'iterations \n Last x tried was:',x,'\n And the value here was:',val,'\n')
 		return(NULL)
 	}
}

l_prime = function(x){
	125/(2+x)-38/(1-x)+34/x
}
l_double_prime = function(x){
	-125/(2+x)**2-38/(1-x)**2-34/x**2
}


# optimize with bisect
bisect(l_prime,c(0,1),tolerance=0.00001)

# optimizie with NR
# all return same value (try 100 times to help avoid potential local min): 
sapply(1:100,function(i){
	newtonRaphson(l_prime,l_double_prime,runif(1),tolerance=0.0001)
})
