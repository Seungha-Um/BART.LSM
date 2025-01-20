library(statnet)
library(MASS)
library(reshape2)
library(dplyr)
library(HLSM)
library(zeallot) 
library(SoftBart)
library(gdata)
library(rlist)
library(ggplot2)
library(ggpubr)
library(invgamma)
library(fastDummies)


source("~/quantile_normalize.R")

expit <- function(x) (1/(1 + exp(-x)))
'%<-%' <- zeallot::'%<-%'
num_iter <- num_burn <- 5000
sim_num <- num_iter + num_burn
num_tree <- 200

in_vec <- seq(1, 10, 1)

p <- 2
N <- 100
num_grp <- 4
num_ind <- rep(N/num_grp, num_grp)
grp_vec <- rep(1:num_grp, num_ind)
rep_iter <- 10

set.seed(2023)
for(rr in 1:rep_iter){
  
  set.seed(rr)
  grp_vec <- sample(grp_vec)
  
  dummy_grp <- fastDummies::dummy_cols(grp_vec)[,-1]
  mu_min <- -1.5    
  mu_max <- 1.5   #3.5   
  mu_grp <- rbind(c(mu_max, mu_max), c(mu_min,mu_min), c(mu_max,mu_min), c(mu_min,mu_max))
  
  degree_sum <- 1
  while(degree_sum != 0) {  
    mu_mat <- mu_grp[grp_vec,]
    X <- t(apply(mu_mat, 1, function(x) mvrnorm(1, mu = x, Sigma = diag(0.3, p))))
    dist_X <- as.matrix(dist(X[,1:2], method = 'euclidean', diag = FALSE, upper = TRUE))
    temp <- 6 -  5 * dist_X
    A_p <- expit(temp)
    A <- sapply(1:N, function(i) rbinom(N, 1, A_p[i,]))
    # Make undirect/symmetric adjacency matrix
    lowerTriangle(A) <- upperTriangle(A, byrow=TRUE)
    diag(A) <- 0
    degree <- rowSums(A)
    degree_sum <- sum(degree == 0)
    if(sum(degree == 0)){print("isolated node exists")}
  }
  
  X <- scale(X) 
  # Between tie 
  nodes_in <- N/num_grp;  total_bet <- (nodes_in * (N - nodes_in)/2) * num_grp
  ties <- rep(NA, num_grp)
  
  for(k in 1:num_grp){
    id <- grp_vec == k
    ties[k] <- sum(rowSums(A[id,!id]) )
  }
  bet_ties <- (sum(ties/2))
  
  # Within
  total_with <- nodes_in * (nodes_in-1)/2 * num_grp
  ties2 <- rep(NA, num_grp)
  for(k in 1:num_grp){
    id <- grp_vec == k
    ties2[k] <- sum(rowSums(A[id,id]))
  }
  with_ties <- (sum(ties2/2))
  df_net <- data.frame(with_ties, bet_ties, total_with, total_bet)
  
  ############  LSM
  priors = list(MuBeta = c(0,0), VarBeta = c(100,100), MuZ = c(0,0),
                VarZ = c(5,5), PriorA =  100, PriorB = 150)
  
  mod_complete = LSM(Y = A,
                     tuneIn = TRUE, dd = 2, estimate.intercept = TRUE, niter = sim_num,
                     priors = priors)
  
  ZZlist = getLS(mod_complete, burnin = 0, thin = 1)[[1]]
  est_median <- apply(ZZlist, c(1,2), median)
  
  set.seed(rr)
  for(tt in 1:2){
    if(tt == 1){
      Treat_fun <-  function(xx){5 + 0 * xx[,1] + 0 * xx[,2] +
          rnorm(length(xx[,1]), mean = 0, sd = 2)}
    }else{
      Treat_fun <-  function(xx){5 + 0.5 * xx[,1] + 0.5 * xx[,2] +
          rnorm(length(xx[,1]), mean = 0, sd = 1)}
    }
    
    Treat <- Treat_fun(X)
    nei_treat <- sapply(1:N, function(i) sum(Treat[A[i,] == 1])/degree[i])
    c(nei_min, nei_max) %<-% quantile(unlist(nei_treat), c(0.025, 0.975))
    nei_vec_true <- seq(nei_min, nei_max, length.out = 10)
    nei_vec <- seq(4, 6, length = 10)
    
    for(q in 1:length(in_vec)){
      

        mu_fun <- function(Treat, G, X){
          Y_mu <- 10 + 3 * Treat + 10 * expit( 10 * (G - 5) ) + in_vec[q] * X[,1] + in_vec[q]* X[,2]
          return(Y_mu)
        }
        true_potential <- 10 + 3 * 5 + 10 * expit( 10 * (nei_vec - 5) )
    
        Y_mu <- mu_fun(Treat = Treat, G = nei_treat, X = X)
        Y <-  rnorm(N, Y_mu, 1)
        
        step_size <- 1000
        
        Y_scaled <- scale(Y)
        c(Y_dims, center_Y, scale_Y) %<-% attributes(Y_scaled)
        df_train <- cbind(nei_treat, Treat, est_median)
        colnames(df_train) <- c("nei_treat", "Treat", "est_Z1", "est_Z2")
        df_train <- quantile_normalize_bart(df_train)
        hypers <- Hypers(X = df_train, Y = Y_scaled, normalize_Y = FALSE)
        
        hypers$num_tree <- num_tree
        opts <- Opts(update_s = FALSE)
        forest <- MakeForest(hypers, opts, warn = FALSE)
        step_size <- 1000
        
        mu_hat <- matrix(NA, sim_num, N)
        mu_test_hat <- matrix(NA, sim_num/step_size, length(nei_vec)*N)
        
        ave_mu_test <- list()
        ACPI_df <- list()
        for(s in 1:sim_num){
          est_Z <- ZZlist[,,s]
          ################### LSM-BART
          df_train <- cbind(nei_treat, Treat, est_Z)
          colnames(df_train) <- c("nei_treat", "Treat", "est_Z1", "est_Z2")
          #df_train <- quantile_normalize_bart(df_train)
          
          temp_mat <- cbind(Treat, est_Z)
          rownames(temp_mat) <- NULL
          df_test <- data.frame(rep(nei_vec, each = N), temp_mat)
          colnames(df_test) <- c("nei_treat", "Treat", "est_Z1", "est_Z2")
          #df_test <- quantile_normalize_bart(df_test)
          XXtest <- quantile_normalize(df_train, df_test)
          df_train <- XXtest$X
          df_test <- XXtest$test_X
          
          df_test <- as.matrix(df_test)
          mu_hat[s,] <- forest$do_gibbs(df_train, Y_scaled, df_train, 1)
          mu_test_temp <- forest$do_predict(df_test) * scale_Y + center_Y
          
          if(s %% step_size == 0) mu_test_hat[(s/step_size),] <- mu_test_temp
          
          # To save estimates by group, creat a new df
          df_temp <- data.frame(nei_treat = rep(nei_vec, each = N), est_Y = mu_test_temp)
          
          ave_mu_test[[s]] <- data.frame(df_temp %>% group_by(nei_treat) %>% 
                                           summarize(mean = mean(est_Y)) %>%  mutate(iter = s))
          ACPI_df[[s]] <-  ave_mu_test[[s]]  %>%  
                            mutate(foo = lead(mean),
                                   IE = foo - mean,
                                   iter = s)
          if(s %% 100 == 0) cat("MCMC iter :", s,"coef", in_vec[q], "treat", tt, "iter",  rr,"\n") 
        }
        
        mu_hat_list <- mu_hat * scale_Y + center_Y
        ave_mu_test_list <- list.rbind(ave_mu_test)
        mu_test_hat_list <- mu_test_hat
        
        # save all outputs
        # out <- list()
        # out$mu_hat_list <- mu_hat_list
        # out$mu_test_hat_list <- mu_test_hat_list
        # out$ave_mu_test_list <- ave_mu_test_list 
        # out$true_potential <- true_potential
        # out$Y <- Y
        # out$Treat_list <-  Treat
        # out$nei_list <- nei_treat
        # out$nei_vec_true <- nei_vec_true
        # out$num_grp <- num_grp 
        # out$N <- N
        # out$nei_vec <- nei_vec
        # out$df_net <- df_net
        out1 <- ave_mu_test_list %>% filter(iter > num_burn) %>% group_by(nei_treat) %>%
          summarize(est_mean = mean(mean),
                    LI = quantile(mean, 0.025),
                    UI = quantile(mean, 0.975)) %>%
          mutate(true_po = true_potential,
                 leng = UI - LI,
                 cov = case_when(
                   true_po <= UI & LI <=true_po ~ 1,
                   TRUE ~ 0
                 ),
                 foo = lead(est_mean), 
                 IE = foo - est_mean,
                 foo2 = lead(true_po), 
                 IE_true = foo2 - true_po,
                 SE = (IE - IE_true)^2
          )
        
        true_PI <- data.frame(nei_treat = nei_vec, true_po = true_potential) %>% 
          mutate(foo2 = lead(true_po), 
                 IE_true = foo2 - true_po) %>% 
          filter(!is.na(foo2))
        
        
        out3 <- list.rbind(ACPI_df) %>% 
                filter(iter > num_burn, !is.na(IE)) %>% group_by(nei_treat) %>% 
                summarize(PI_mean = mean(IE),
                          PI_LI = quantile(IE, 0.025),
                          PI_UI = quantile(IE, 0.975)) %>%
                left_join(true_PI, by = "nei_treat") %>%
                mutate(leng = PI_UI - PI_LI,
                       cov = PI_LI <= IE_true & IE_true <= PI_UI,
                       SE = (IE_true - PI_mean)^2)
              
        out <- list()
        out$out1 <- out1
        out$out3 <- out3
        
        saveRDS(out, paste0("BART_LSM","coef", in_vec[q], "treat", tt,  "iter", rr, ".rds"))
      }  
    }  
  }
