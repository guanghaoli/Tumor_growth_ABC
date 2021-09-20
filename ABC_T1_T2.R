
library(abc)
library(ggplot2)
library(plyr)

t1_mut_number <- c(25735,25677,25968,26130,25962,45049,68149,26321,25959,26387,25931,26116,21345,65884,50067,42174)
t1_sample <- c("t1f2","t1f4","t1f9","t1f11","t1f14","t1f23","t1f24","t1l1","t1l3","t1l6","t1l8","t1l10","t1l13","t1z1","t1z3","t1z5")
t1_number <- c(max(t1_mut_number),min(t1_mut_number),mean(t1_mut_number),median(t1_mut_number),sd(t1_mut_number))
t1_mut_poly <- t1_mut_number - 6869
t1_number <- c(max(t1_mut_poly),min(t1_mut_poly),mean(t1_mut_poly),median(t1_mut_poly),sd(t1_mut_poly))

#####************** Four intervals ***************#####
partition_t1 <- c(0.02,0.25,0.4,0.9,1.0)
sum_stat_interval <- function(ccf_dis,partition){
  stat_par <- c()
  par_len <- length(partition)
  for(ii in 1:(par_len-1)){
    par_value <- length(which(ccf_dis>partition[ii] & ccf_dis<partition[ii+1]))
    stat_par <- c(stat_par,par_value)
  }
  return(stat_par)
}

t1_ccf <- read.table("t1_MAFA_CCF.txt")

### threshold is 0.02
hist(t1_ccf[t1_ccf$V2>0.02,2],breaks = 100)
t1_ccf_obs <- sum_stat_interval(t1_ccf$V2,partition_t1)

#### fst and ks distance between samples
t1_fst <- read.table("t1_fst_result.txt",header = T,stringsAsFactors = F)
t1_ksd <- read.table("t1_ksd_result.txt",header = T,stringsAsFactors = F)
t1_ksd_fst <- c(median(t1_ksd$ksd.value),mean(t1_ksd$ksd.value),sd(t1_ksd$ksd.value),
                median(t1_fst$fst.value),mean(t1_fst$fst.value),sd(t1_fst$fst.value))
t1_ccf_obs <- c(t1_ccf_obs,t1_ksd_fst)

t1_ccf_obs <- c(t1_ccf_obs,t1_number)

############## The simulated Data ##################

t1_ccf_simu <- read.table("T1_sum_stat_partition.txt",header = T,stringsAsFactors = F,sep="\t")

param.value <- t1_ccf_simu[,c(1,2)]
sumstat.value <- t1_ccf_simu[,5:19]

#### cross validation for parameters
cv.rej.t1 <- cv4abc(param = param.value,sumstat = sumstat.value,nval = 10,tols = c(0.005,0.01,0.05),method = "rejection")
summary(cv.rej.t1)

#### chose 0.01
t1_rej <- abc(target = t1_ccf_obs,param = param.value,sumstat = sumstat.value,tol = 0.01,method = "rejection")

t1_rej_data <- data.frame(t1_rej$unadj.values)
write.table(t1_rej_data,"T1_posterior_U_comp_partition_0.01.txt",sep="\t",col.names = T,row.names = F,quote = F)

##### posterior of mutation rate and number of mutation #####
barplot(table(t1_rej_data$U))
barplot(table(t1_rej_data$N))



##### For T2 ########

t2_mut_number <- c(25356,25429,25579,25605,70093,29778,25229,24370,30193)
t2_sample <- c("T2F2","T2F5","T2F9","T2F13","T2Z1","T2Z6","T2Z9","T2Z11","T2Z13")
t2_number <- c(max(t2_mut_number),min(t2_mut_number),mean(t2_mut_number),median(t2_mut_number),sd(t2_mut_number))

t2_mut_poly <- t2_mut_number - 3721
t2_number <- c(max(t2_mut_poly),min(t2_mut_poly),mean(t2_mut_poly),median(t2_mut_poly),sd(t2_mut_poly))


#####************** Four intervals ***************#####
partition_t2 <- c(0.02,0.25,0.4,0.9,1.0)

t2_ccf <- read.table("t2_MAFA_CCF.txt")
### threshold is 0.02
hist(t2_ccf[t2_ccf$V2>0.02,2],breaks = 100)
t2_ccf_obs <- sum_stat_interval(t2_ccf$V2,partition_t2)
## fst and ks distance
t2_fst <- read.table("t2_fst_result.txt",header = T,stringsAsFactors = F)
t2_ksd <- read.table("t2_ksd_result.txt",header = T,stringsAsFactors = F)
t2_ksd_fst <- c(median(t2_ksd$ksd.value),mean(t2_ksd$ksd.value),sd(t2_ksd$ksd.value),
                median(t2_fst$fst.value),mean(t2_fst$fst.value),sd(t2_fst$fst.value))
t2_ccf_obs <- c(t2_ccf_obs,t2_ksd_fst)
t2_ccf_obs <- c(t2_ccf_obs,t2_number)

############## The simulated Data ##################

t2_ccf_simu <- read.table("T2_sum_stat_partition.txt",header = T,stringsAsFactors = F,sep="\t")
#### no weight
param.value <- t2_ccf_simu[,c(1,2)]
sumstat.value <- t2_ccf_simu[,5:19]

#### cross validation for parameters
cv.rej.t2 <- cv4abc(param = param.value,sumstat = sumstat.value,nval = 10,tols = c(0.005,0.01,0.05),method = "rejection")
summary(cv.rej.t2)

t2_rej <- abc(target = t2_ccf_obs,param = param.value,sumstat = sumstat.value,tol = 0.01,method = "rejection")

# t2_rej <- abc(target = t2_ccf_obs,param = param.value,sumstat = sumstat.value,tol = 0.05,method = "rejection")
t2_rej_data <- data.frame(t2_rej$unadj.values)
write.table(t2_rej_data,"T2_posterior_U_comp_partition_0.01.txt",sep="\t",col.names = T,row.names = F,quote = F)

##### posterior of mutation rate #####
barplot(table(t2_rej_data$U))
barplot(table(t2_rej_data$N))


##### Analysis each parameter's posterior distribution #####
t1_rej_data$label <- rep("T1",dim(t1_rej_data)[1])
t2_rej_data$label <- rep("T2",dim(t2_rej_data)[1])


All_data <- rbind(t1_rej_data,t2_rej_data)


data_U <- ddply(All_data,.(U,label),summarise,num=length(U))
data_U$freq <- data_U$num/dim(t1_rej_data)[1]
write.table(data_U,"Posterior_U.txt",sep="\t",row.names = F,col.names = T,quote = F)


data_N <- ddply(All_data,.(N,label),summarise,num=length(N))
data_N$freq <- data_N$num/dim(t1_rej_data)[1]

write.table(data_N,"Posterior_N.txt",sep="\t",row.names = F,col.names = T,quote = F)

