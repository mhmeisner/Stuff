
mini <- FALSE

#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 1000
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  sim_num <- sim_start + as.numeric(args[1])
  sim_seed <- (762*(sim_num-1) + 121231)
}

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))

############ Find r and s indices:
sim_num_1_thru_rs = sim_num-sim_start # get number in 1-250 range:
# define r and s; calculate their product rs:
r = 50 
s = 5
rs = r*s
# create 1, 51, 101, 151, 201, 251 vector, so that I can then see which of these bins the current index falls in (this will determine s)
s_indices = seq(from=1,to=rs+1,by=r)
# determine which of the bins the current sim_num falls in:
in_bin = sapply(1:s, function(i){
	sim_num_1_thru_rs >= s_indices[i] & sim_num_1_thru_rs < s_indices[i+1]
})
# and then set the value of s accordingly
s_index = (1:s)[in_bin]
# find r by seeing how many iterations have passed since the start of the current index of s 
r_index = sim_num_1_thru_rs-r*(s_index-1)

#============================== Run the simulation study ==============================#

# Load packages:
library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)

# I/O specifications:
datapath <- "/home/pdbaines/data"
outpath <- "output/"

# mini or full?
if (mini){
	rootfilename <- "blb_lin_reg_mini"
} else {
	rootfilename <- "blb_lin_reg_data"
}

# Filenames:
# desc_file = '~/Documents/R/blb_lin_reg_mini.desc' # for testing
desc_file = paste0(datapath,'/',rootfilename,'.desc')

# Attach big.matrix :
bm = attach.big.matrix(desc_file)

# Remaining BLB specs:
gamma = 0.7
n = nrow(bm)
b = round(n**gamma) # need an integer, and b only needs to be approximately n**gamma

# Extract the subset:
set.seed(s_index) # make sure that jobs with the same s value use the same rows!
sub_indices = sample(1:n,b,replace=F) # get random indices
bm_sub = bm[sub_indices,] # create matrix with only the rows in this subset

# Reset simulation seed:
set.seed(sim_seed) # now that the appropriate rows for the current s index have been extracted, we want to use a different boostrap sample each time that we sample n points from that subset

# Bootstrap dataset:
reps =rmultinom(n=1,size=n,prob =rep(1,b)) # we want to sample n data points from this subset. the probs are all specified to be the same, since we will use sample n data points with replacment from our subset 

# Fit lm:
y_colname = paste('V',ncol(bm),sep='') # this is the column name of the last column, which is the response variable...we need this to properly specify the formula to lm

bm_df = as.data.frame(bm_sub) # convert to data frame; this is easier for lm 

f = formula(paste0(y_colname,'~.')) # create formula object using the character describing the column name for the response variable. the period means we use ALL other columns as predictors 

m = lm(f,data=bm_df,weights=reps) # fit the linear model 

# Output file:
outfile = paste0("output/","coef_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt")

# Save estimates to file:
write.table(coef(m),outfile,col.names=F,row.names=F)