#all packages needed for this code
install.packages("flexclust")
install.packages("grid")
install.packages("lattice")
install.packages("modeltools")
install.packages("stats4")
library(flexclust)
library(grid)
library(lattice)
library(stats4)
install.packages("DescTools")
library(DescTools)
install.packages("RSKC")
install.packages("dplyr")
library(RSKC)
library(dplyr)
install.packages("sparcl")
library(sparcl)
library(data.table)
library('devtools')
install.packages("NbClust")
library(NbClust)
install.packages("cluster")
library(cluster)

#Initialize
set.seed(123)
CER_s_after_k_LVR <- numeric(20)
CER_s_after_k_SDR <- numeric(20)
CER_AGap <- numeric(20)
CER_CV_LVR <- numeric(20)
CER_CV_SDR <- numeric(20)
jointLVR_CER <- numeric(20)
jointSDR_CER <- numeric(20)
#create subset of indices to use from NbClust package
bro <- c("kl", "ch", "hartigan", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "mcclain", "gamma", "frey", "gplus", "tau", "dunn", "sdindex", "sdbw")
i = 1
#define the number of features and the mu 
p <- 
mu <- 

#run 20 simulations
while (i <= 20) {
  #create a matrix with the generated data of 60 rows and p columns
  x <- matrix(rnorm(60*p),ncol=p)
  #create the difference in the first 50 features for all classes
  x[1:20, 1:50] <- x[1:20, 1:50] + mu
  x[21:40, 1:50] <- x[21:40, 1:50] - mu
  #standardize as is done in the example for the sparcl package
  x <- scale(x, TRUE, TRUE)
  #create a vector of labels
  c <- rep(1:3, each = 20)
  #put labels together with data for later stage with CER
  dataused <- cbind(c, x)
  #generate row names
  rownames(dataused) <- seq(1,60,1)
  #seperate data for clustering
  dataused2 <- x

  k <- numeric(21)
  
  #for the subset of indices from the NbClust function calculate the optimal number of clusters in the data between 2 and 5
  for (j in 1:21) {
    method_used <- bro[j]
    output <- NbClust(data = dataused2, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 5, method = "kmeans", index = method_used)
    k[j] <- output$Best.nc[1]
  }
  
  #by majority rule select the optimal number of clusters from the selected in the for loop above
  Number_of_clusters <- Mode(as.vector(k))[1]
  
  best_s_LVR <- numeric(25)
  best_s_SDR <- numeric(25)
  best_s_AGap <- numeric(25)
  
  #generate sequence of possible values for s in between 1.1 and the square root of the number of features of dataused2
  s_values <- seq(1.1, sqrt(p), 0.5)
  j = 1
  
  #execute 25 times the method for selecting the optimal s with the largest value and standard deviation rule in addition to the Altered gap statistic
  while (j <= 25) {
    
    #shuffle data to change ordering
    datashuffled <- dataused2[sample(1:nrow(dataused2)),]
    
    #permute over range of possible s values to find the optimal one
    s <- KMeansSparseCluster.permute(datashuffled, Number_of_clusters, 5, wbounds = s_values, silent=TRUE)
    
    #best s by largest value rule
    best_s_LVR[j] <- s$bestw
    
    #determine best s by standard deviation rule
    id_largest <- which.max(s$gaps)
    id_st_dev <- which(s$gaps[1:id_largest] > s$gaps[id_largest] - s$sdgaps[id_largest])[1]
    best_s_SDR[j] <- s$wbounds[id_st_dev]
    
    #calculate AGap values and find optimal s as the s with the largest AGap
    AGap <- s$gaps - log(s$nnonzerows)
    id_optimal_AGap <- which.max(AGap)
    best_s_AGap[j] <- s$wbounds[id_optimal_AGap]
    
    j = j+1
  }

  #prepare for cross validation
  best_s_CV_LVR <- numeric(25)
  best_s_CV_SDR <- numeric(25)
  j = 1
  
  #loop for cross validation
  while (j <= 25) {
    #shuffle data to change ordering
    datashuffling <- dataused[sample(1:nrow(dataused)),]
    
    #give data a number to indicate the fold id
    datasplit <- cut(seq(1,nrow(datashuffling)),breaks=5,labels=FALSE)
    
    #prepare for steps of cross validation
    CER_list_LVR <- numeric(5)
    selected_s_LVR <- numeric(5)
    CER_list_SDR <- numeric(5)
    selected_s_SDR <- numeric(5)
    
    #perform cross validation with each fold as validation fold once
    for (t in 1:5) {
      #select data that belongs to a certain fold to separate into the training data and validation data
      validateIndex <- which(datasplit == t)
      datavalidate <- datashuffling[validateIndex,]
      datatrain <- datashuffling[-validateIndex,]
      
      #remove labels
      datavalidating <- datavalidate[, -c(1)]
      datatraining <- datatrain[, -c(1)]
      
      #find optimal s over range of s values in training set
      Train_s <- KMeansSparseCluster.permute(datatraining, Number_of_clusters, 3, wbounds = s_values, silent = TRUE)
      
      #best s in training data by largest value rule
      selected_s_LVR[t] <- Train_s$bestw
      
      #determine best s in training data by standard deviation rule
      id_max <- which.max(Train_s$gaps)
      id_st_dev <- which(Train_s$gaps[1:id_max] > Train_s$gaps[id_max] - Train_s$sdgaps[id_max])[1]
      selected_s_SDR[t] <- Train_s$wbounds[id_st_dev]
      
      #cluster the validation sample and calculate the CER
      clustering_LVR <- KMeansSparseCluster(datavalidating, Number_of_clusters, wbounds = Train_s$bestw, nstart = 20, silent = TRUE, maxiter = 10)
      labels_validate_LVR <- sapply(clustering_LVR,"[[",2)
      CER_list_LVR[t] <- CER(labels_validate_LVR, datavalidate[,1], nrow(datavalidate))
      
      clustering_SDR <- KMeansSparseCluster(datavalidating, Number_of_clusters, wbounds = Train_s$wbounds[id_st_dev], nstart = 20, silent = TRUE, maxiter = 10)
      labels_validate_SDR <- sapply(clustering_SDR, "[[",2)
      CER_list_SDR[t] <- CER(labels_validate_SDR, datavalidate[,1], nrow(datavalidate))
    }
    
    #find optimal (lowest value) CER with corresponding optimal s value with both rules
    id_CV_LVR <- which.min(CER_list_LVR)
    id_CV_SDR <- which.min(CER_list_SDR)
    best_s_CV_LVR[j] <- selected_s_LVR[id_CV_LVR]
    best_s_CV_SDR[j] <- selected_s_SDR[id_CV_SDR]
    j = j + 1
  }
  
  #determine optimal s by taking mean and execute clustering. save value of CER and the nonzero weights for all methods but the jointly determination of K and s
  optimal_s_LVR <- mean(best_s_LVR)
  LVR <- KMeansSparseCluster(dataused2, Number_of_clusters, optimal_s_LVR, 20, silent = TRUE)
  labels_LVR <- sapply(LVR,"[[",2)
  #if needed the code below can be used to save the number of nonzero weights
  weights_LVR <- sapply(LVR,"[[",1)
  nonzeroes_LVR <- length(which(weights_LVR !=0))
  #calculate the CER for the specific rule
  CER_s_after_k_LVR[i] <- (CER(labels_LVR, dataused[,1], nrow(dataused)))
  
  
  #below the steps are the same as above, but with different rules, in the order: standard deviation rule, altered gap statistic
  #cross-validation by largest value rule and cross-validation by standard deviation rule
  optimal_s_SDR <- mean(best_s_SDR)
  SDR <- KMeansSparseCluster(dataused2, Number_of_clusters, optimal_s_SDR, 20, silent = TRUE)
  labels_SDR <- sapply(SDR,"[[",2)
  weights_SDR <- sapply(SDR,"[[",1)
  nonzeroes_SDR <- length(which(weights_SDR!=0))
  CER_s_after_k_SDR[i] <- (CER(labels_SDR, dataused[,1], nrow(dataused)))
  
  optimal_s_AGap <- mean(best_s_AGap)
  AGap <- KMeansSparseCluster(dataused2, Number_of_clusters, optimal_s_AGap, 20, silent = TRUE)
  labels_AGap <- sapply(AGap,"[[",2)
  weights_AGap <- sapply(AGap,"[[",1)
  nonzeroes_AGap <- length(which(weights_AGap!=0))
  CER_AGap[i] <- (CER(labels_AGap, dataused[,1], nrow(dataused)))
  
  optimal_s_CV_LVR <- mean(best_s_CV_LVR)
  CV_LVR <- KMeansSparseCluster(dataused2, Number_of_clusters,  optimal_s_CV_LVR, 20, silent = TRUE)
  labels_CV_LVR <- sapply(CV_LVR,"[[",2)
  weights_CV_LVR <- sapply(CV_LVR,"[[",1)
  nonzeroes_CV_LVR <- length(which(weights_CV_LVR!=0))
  CER_CV_LVR[i] <- (CER(labels_CV_LVR, dataused[,1], nrow(dataused)))
  
  optimal_s_CV_SDR <- mean(best_s_CV_SDR)
  CV_SDR <- KMeansSparseCluster(dataused2, Number_of_clusters, optimal_s_CV_SDR, 20, silent = TRUE)
  labels_CV_SDR <- sapply(CV_SDR,"[[",2)
  weights_CV_SDR <- sapply(CV_SDR,"[[",1)
  nonzeroes_CV_SDR <- length(which(weights_CV_SDR!=0))
  CER_CV_SDR[i] <- CER(labels_CV_SDR, dataused[,1], nrow(dataused))
  
  
  #set range of clusters
  cluster_options <- seq(2, 5, 1)
  
  #construct matrices to save all values for the jointly algorithm
  amount <- length(cluster_options)*25
  JointKandS_LVR <- numeric(amount)
  dim(JointKandS_LVR) <- c(length(cluster_options), 25)
  JointKandS_SDR <- numeric(amount)
  dim(JointKandS_SDR) <- c(length(cluster_options), 25)
  t = 1
  
  #loop 25 times each time shuffling the data and then determining for all possible values of K the optimal value of s according to both rules
  while (t <= 25) {
    datashuffled <- dataused2[sample(1:nrow(dataused2)),]
    for (j in 1:length(cluster_options)) {
      s <- KMeansSparseCluster.permute(datashuffled, cluster_options[j], 5, wbounds = seq(1.1, sqrt(p), 0.5), silent=TRUE)
      JointKandS_LVR[j,t] <- s$bestw
      id_largest <- which.max(s$gaps)
      id_st_dev <- which(s$gaps[1:id_largest] > s$gaps[id_largest] - s$sdgaps[id_largest])[1]
      JointKandS_SDR[j,t] <- s$wbounds[id_st_dev]
    }
    t = t + 1
  }
  
  best_k_and_s_LVR <- numeric(length(cluster_options))
  best_k_and_s_SDR <- numeric(length(cluster_options))
  
  #calculate optimal value for s for each number of clusters by taking the mean over the rows corresponding to the number of clusters
  for (j in 1:length(cluster_options)) {
    best_k_and_s_LVR[j] <- mean(JointKandS_LVR[j,])
    best_k_and_s_SDR[j] <- mean(JointKandS_SDR[j,])
  }
  
  Gaps_LVR <- numeric(length(cluster_options))
  Gaps_SDR <- numeric(length(cluster_options))
  
  #find the gap statistic for each number of cluster and optimal s value combination
  for (j in 1:length(cluster_options)) {
    LVR_check <- KMeansSparseCluster.permute(dataused2, cluster_options[j], 5, best_k_and_s_LVR, silent=TRUE)
    Gaps_LVR[j] <- LVR_check$gaps[j]
    SDR_check <- KMeansSparseCluster.permute(dataused2, cluster_options[j], 5, best_k_and_s_SDR, silent=TRUE)
    Gaps_SDR[j] <- SDR_check$gaps[j]
  }
  
  #get the optimal value of s and K jointly using both rules
  id_joint_LVR_s <- which.max(Gaps_LVR)
  joint_LVR_s <- best_k_and_s_LVR[id_joint_LVR_s]
  joint_LVR_k <- cluster_options[id_joint_LVR_s]
  
  id_joint_SDR_s <- which.max(Gaps_SDR)
  joint_SDR_s <- best_k_and_s_SDR[id_joint_SDR_s]
  joint_SDR_k <- id_joint_SDR_s[id_joint_SDR_s]
  
  #execute clustering with optimal clusters and value for s for both rules
  jointLVR <- KMeansSparseCluster(dataused2, joint_LVR_k, joint_LVR_s, 20, silent = TRUE)
  jointLVR_labels <- sapply(jointLVR,"[[",2)
  weightsJointLVR <- sapply(jointLVR,"[[",1)
  nonzeroesJointLVR <- length(which(weightsJointLVR!=0))
  jointSDR <- KMeansSparseCluster(dataused2, joint_SDR_k, joint_SDR_s, 20, silent = TRUE)
  jointSDR_labels <- sapply(jointSDR,"[[",2)
  weightsJointSDR <- sapply(jointSDR,"[[",1)
  nonzeroesJointSDR <- length(which(weightsJointSDR!=0))
  
  #calculate the CER for jointly determining with both rules
  jointLVR_CER[i] <- CER(jointLVR_labels, dataused[,1], nrow(dataused))
  jointSDR_CER[i] <- CER(jointSDR_labels, dataused[,1], nrow(dataused))
  
  i = i + 1
}


