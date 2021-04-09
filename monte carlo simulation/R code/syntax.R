#----------------------------------------------------------------------------------------------#
# Author: Joseph Kush (jkush1@jhu.edu) (jmk7cj@virginia.edu)
# 
# Title: Statistical power for randomized controlled trials with clusters of varying size 
#        R code for Monte Carlo simulations
#
# Date: 4/20/2021
#
# Purpose: Master .R file to set up and run a Monte Carlo simulation study for 2-level CRT
#          designs with varying cluster size, comparing simulation-based power to calculations
#          of power based on either the arithmetic average number of L1 units per L2 cluster or
#          the harmonic mean number of L1 units per L2 cluster
#          Step 1: Generate facets for Monte Carlo study
#          Step 2: Generate data, estimate treatment effect, calculate power
#          Step 3: Store and view results
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Step 1: Generate facets for Monte Carlo study
#----------------------------------------------------------------------------------------------#
# Load necessary packages, remove all objects
library("lme4"); library("lmerTest"); require("psych"); library("plyr")
library("foreach"); library("doParallel"); library("doSNOW")
rm(list = ls())

# Facets to vary: 4 x 3 x 3 x 5 x 7 = 1,260 conditions
set.seed(123)
sims <- 1:1 # 1 replication as example
n_schools <- c(20, 30, 50, 60) # L2 sample size
rho <- c(0.05, 0.10, 0.20) # ICCs
effect_sizes <- c(0.2, 0.3, 0.4) # Standardized effect sizes 
min_students <- c(5, 10, 15, 20, 30) # Min i/j
max_students <- c(10, 15, 20, 25, 30, 40, 50) # Max i/j

# Create progress bar for simulation based on all conditions + replications
n_unique_sim_cells = sum(table(count(sims)))*sum(table(count(n_schools)))*
sum(table(count(rho)))*sum(table(count(effect_sizes)))*
sum(table(count(min_students)))*sum(table(count(max_students)))

processors <- makeCluster(detectCores()[1]-1) 
registerDoSNOW(processors)
pb <- txtProgressBar(max=n_unique_sim_cells, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Step 2: Generate data, estimate treatment effect, calculate power
#----------------------------------------------------------------------------------------------#
simulate_function = 
foreach(i=sims, .combine=rbind) %:%
foreach(j=n_schools, .combine=rbind) %:%  
foreach(icc=rho, .combine=rbind) %:%  
foreach(beta=effect_sizes, .combine=rbind) %:%  
foreach(min_ij=min_students, .combine=rbind) %:%  
foreach(max_ij=max_students, .combine=rbind, .options.snow=opts) %dopar% {

library(lme4); library(lmerTest); require(psych)

# Create clusters of different size
sch1 <- rep(min_ij, j/2)
sch2 <- rep(max_ij, j/2)
schools <- 1:j 
students <- c(sch1, sch2)
schid <- rep(schools, times=students) 
studid <- unlist(lapply(students, seq))
data <- data.frame(schid, studid)

# Random error at both L1 and L2
l2_sd <- sqrt(icc) 
l1_sd <- sqrt(1-icc)
l2var_a <- rnorm(n=(j/2), mean=0, sd=l2_sd)
l2var_b <- rnorm(n=(j/2), mean=0, sd=l2_sd)
l2var_a <- rep(l2var_a, each=min_ij) 
l2var_b <- rep(l2var_b, each=max_ij)
l2var <- c(l2var_a, l2var_b)
l1var <- rnorm(nrow(data), mean=0, sd=l1_sd)
data <- data.frame(schid, studid, l1var, l1var)

# Randomly assign treatment status to each school  
treat_a <- rbinom(n=j/2, size=1, prob=0.5) 
treat_b <- rbinom(n=j/2, size=1, prob=0.5)
treat_a <- rep(treat_a, each=min_ij)
treat_b <- rep(treat_b, each=max_ij)
treat <- c(treat_a, treat_b)

# Generate L1 outcome
y <- (beta*treat + l2var + l1var)
data <- data.frame(schid, studid, l1var, l2var, treat, y)

# Estimate treatment effect
model <- lmer(y ~ treat + (1|schid), data=data)

# Calculate simulation-based power using % of significant replications
p <- summary(model)[[10]][2,5] 
power <- ifelse(p<0.05, 1, 0)

# Calculate power using arithmetic average i/j
prob_treat <- round(mean(data[,5]), digit=2)
avg_iperj <- nrow(data) / j

aa_lambda <- beta * sqrt(prob_treat * (1-prob_treat) *j) / 
  sqrt(icc * (1-0) + (1-icc) * (1-0) / avg_iperj)

aa_power <- 1 - pt(qt(1-.05/2, df=j-0-2), df=j-0-2, aa_lambda)


# Calculate power using harmonic mean i/j
harm_avg_iperj = round( j / (((j/2)/min_ij) + ((j/2)/max_ij)) , digits=2)

hm_lambda <- beta * sqrt(prob_treat * (1-prob_treat) *j) / 
  sqrt(icc * (1-0) + (1-icc) * (1-0) / harm_avg_iperj)

hm_power <- 1 - pt(qt(1-.05/2, df=j-0-2), df=j-0-2, hm_lambda)

# Store values
v1 <- i; v2 <- beta; v3 <- prob_treat; v4 <- j; v5 <- min_ij; v6 <- max_ij
v7 <- avg_iperj; v8 <- harm_avg_iperj; v9 <- icc; v10 <- power
v11 <- aa_power; v12 <- hm_power

# Return list
list(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12)
}
#----------------------------------------------------------------------------------------------#  


#----------------------------------------------------------------------------------------------#
# Step 3: Store and view results
#----------------------------------------------------------------------------------------------#
results <- as.data.frame(simulate_function)
list_cols <- sapply(results, is.list)
results <- cbind(results[!list_cols], t(apply(results[list_cols], 1, unlist)))
colnames(results) <- c("iteration", "effect_size", "prob_treat", 
                          "schools", "min_ij", "max_ij", "avg_ij", 
                          "harmonic_avg_iperj", "icc", "sim_power", 
                          "aa_power", "hm_power")
results <- results[with(results, order(iteration, effect_size, schools, min_ij, max_ij, icc)),]
# Collapse across all replications for average values
results <- aggregate(. ~ effect_size + prob_treat + schools + min_ij + max_ij, data=results, FUN=mean)
results$iteration <- length(sims)
summary(results)
#----------------------------------------------------------------------------------------------#  