# Function for Neighborhood Proportion Error (NPE)
# Reference: Konstorum et al. (2018) https://doi.org/10.1101/273862
# Authors: Anna Konstorum, Nathan Jekel (2019)
# Contact: anna.konstorum@yale.edu


# NPE_ERROR
#
# This function uses neighborhood homogeneity preservation to determine the accuracy
# of a dimension reduction (designed for gated data)
#
# INPUT ARGUMENTS:
# 1. sample_data: Must be a data frame with the rows being the events and the
# columns being the variables. Additionally, the last 2 columns
# must contain the integer label and the population name of each cell (in that order).
#
# 2. reduced_data: the data input must be a matrix or data.frame where the second to last
# column contains integers corresponding to cell population and the last
# column contains the names of cell populations
#
# 3. population_key: two column matrix with first row being integer labels and the second
# row being the names of the cell populations
#
# 4. k: NPE neighborhood size

npe_error<- function(sample_data, reduced_data, population_key, k){
  library(entropy)
  library(FNN)
  library(som)
  library(topicmodels)
  library(distrEx)
  library(distr)
  
  sd_cols <- ncol(sample_data)
  rd_cols <- ncol(reduced_data)
  # Normalize the data by column to have mean 0 and variance 1
  labels <- as.data.frame(sample_data[,(sd_cols-1):(sd_cols)])
  norm <- function(x) {(x-mean(x))/sd(x)}
  saved_sample_n <- lapply(sample_data[,1:(sd_cols-2)], norm)
  saved_sample_n <- cbind(saved_sample_n, labels)
  neighbors <- knn.index(saved_sample_n[,1:(sd_cols-2)], k)
  
  # For loop to replace all instances of 'NA' with 0 
  freq_high <- rep(0,nrow(neighbors))
  zero_saved_sample <- saved_sample_n
  zero_saved_sample[is.na(zero_saved_sample)] <- 0 
  
  # For loop to create a vector with the number of like neigbors each data point has
  for (i in 1:nrow(neighbors)){
    for (j in 1:ncol(neighbors)){
      if ((zero_saved_sample[i,ncol(zero_saved_sample)-1])==(zero_saved_sample[(neighbors[i,j]),(ncol(zero_saved_sample)-1)])){
        freq_high[i]=freq_high[i]+1
      }
    }
  }
  
  # Turn frequencies into labeled frequencies
  labeled_freq <- cbind(zero_saved_sample[,ncol(zero_saved_sample)-1],freq_high)
  
  # Split into groups
  freq_list_high <- split(labeled_freq[,2], labeled_freq[,1], drop = FALSE)
  freq_list_high$`0` <- NULL
  
  get_freq <- function(x){
    count_vec <- rep(0,k)
    for(i in 1:k){
      count <- length(which(x==i))
      count_vec[i] = count
    }
    if (sum(count_vec, na.rm = FALSE)==0){
      count_vec[1]<-1
    }
    return(count_vec)
  }
  
  get_dense <- function(x){
    sum=sum(as.vector(x), na.rm = FALSE)
    for(i in 1:k){
      x[i]<-x[i]/sum
    }
    return(x)
  }
  
  freqfreq_high <- sapply(freq_list_high, get_freq)
  freqdist_high <- apply(freqfreq_high,2, get_dense)
  
  
  
  ################################################################################
  # Repeat the process on the low dimensional data
  # Normalize the data by column to have mean 0 and variance 1
  labels <- as.data.frame(reduced_data[,(rd_cols-1):(rd_cols)])
  saved_sample_n <- lapply(reduced_data[,1:(rd_cols-2)], norm)
  saved_sample_n <- cbind(saved_sample_n, labels)
  neighbors <- knn.index(saved_sample_n[,1:(rd_cols-2)], k)
  
  
  # For loop to replace all instances of 'NA' with 0 
  freq_low <- rep(0,nrow(neighbors))
  zero_saved_sample <- saved_sample_n
  zero_saved_sample[is.na(zero_saved_sample)] <- 0 
  
  # For loop to create a vector with the number of like neighbors each data point has
  for (i in 1:nrow(neighbors)){
    for (j in 1:ncol(neighbors)){
      if ((zero_saved_sample[i,ncol(zero_saved_sample)-1])==(zero_saved_sample[(neighbors[i,j]),(ncol(zero_saved_sample)-1)])){
        freq_low[i]=freq_low[i]+1
      }
    }
  }
  
  # Turn frequencies into labeled frequencies
  labeled_freq <- cbind(zero_saved_sample[,ncol(zero_saved_sample)-1],freq_low)
  
  # Split into groups
  freq_list_low <- split(labeled_freq[,2], labeled_freq[,1], drop = FALSE)
  freq_list_low$`0` <- NULL
  
  
  freqfreq_low <- sapply(freq_list_low, get_freq)
  freqdist_low <- apply(freqfreq_low,2, get_dense)
  
  # calculate error
  num_pops <- nrow(population_key)
  
  # Convert lowD counts to discrete distributions
  bins <- seq(1,k,by=1)
  tv_dist <- rep(0,num_pops)
  for(i in 1:num_pops){
    ith_label<-as.character(i)
    if (ith_label %in% colnames(freqdist_low)){
      dist_low <- DiscreteDistribution(bins,freqdist_low[,ith_label])
      dist_high <- DiscreteDistribution(bins,freqdist_high[,ith_label])
      tv_dist[i] <- TotalVarDist(dist_high, dist_low) }
  }
  tv_error <- sum(tv_dist)
  return(tv_error)
}
