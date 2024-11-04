###############load required packages
library(MASS)
library(dplyr)
library(rockchalk)
library(rBeta2009)
library(Matrix)
library(psych)
library(bbricks)
library(Rlab)
library(mvtnorm)
library(abind)
library(pROC)
library(caret)
library(lme4)
library(coxme)
library(Rcpp)
library(prodlim)
sourceCpp("./optimizeArray.cpp")

############### helper function
log_h2<-function(l,inverse_block_diag_matrices){
  result <- sapply(1:length(inverse_block_diag_matrices), function(j) {
    -0.5 * sum((t(l)[j, ] %*% inverse_block_diag_matrices[[j]] %*% l[, j]))
  })
  return(result)
}

rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}

################data generation
num_visit <- 3 #T
n=1000 #number of subjects
N=num_visit*n #number of total observations
t=rep(c(1:num_visit),n)
V <- 100
Q=5
conn_block <- ((Q+1)*Q)/2 

#generate covariates
s_1<-sample(1:3,n,replace=TRUE)
s_1<-scale(s_1)
s=c(rep(s_1,each=3))
X=cbind(s,t,s*t)

#generate delta
sparsity=0.9
ssum=ssum1=0
while(ssum!=ceiling((1-sparsity)*(conn_block)*3 )|| ssum1==0 ||ssum2==0){
  delta_sim=array(rbinom((conn_block)*3,1,(1-sparsity)),dim=c(3,conn_block))
  delta_sim[3,which(apply(delta_sim, 2, function(col) all(col == c(0, 0, 1)|col== c(1, 0, 1)|col == c(0, 1, 1)))==T)]=0
  delta_sim[3,which(apply(delta_sim, 2, function(col) all(col == c(1, 1, 0)))==T)]=sample(0:1,1)
  ssum<-sum(delta_sim==1)
  ssum1=sum(delta_sim[1,]==1)
  ssum2=sum(delta_sim[3,]==1)
}
delta_sim<-t(delta_sim)

#generate alpha and beta
alpha_sim=delta_sim*0.15
beta_sim=alpha_sim*delta_sim

#generate D and b
D<-array(c(1,0.3*1,0.3*1,1),dim=c(2,2))
b_sim=array(0,dim=c(conn_block,n,2))
for(i in 1:conn_block){
  b_sim[i,,]=rmvnorm(n, sigma=D)
}

#generate mit
Z_3=array(0,dim=c(n,3,2))
for(i in 1:n){
  Z_3[i,,]<-cbind(1,c(1,2,3))
}
mit=array(0, dim = c(N, conn_block))
re=array(0,dim=c(3,n))
rel=array(0,dim=c(N, conn_block))
for(m in 1:conn_block){
  for(k in 1:n){
    re[,k]=Z_3[i,,]%*%(b_sim[m,k,])
  }
  rel[,m]<-c(re)
  mit[,m]=X%*%beta_sim[m,]+c(re)+(rnorm(N,0,1))
}

#generate mit_qr
mit_qr <- array(0, dim=c(n, num_visit, conn_block))
for(q in 1:n){
  for(w in 1:num_visit){
    for (e in 1: conn_block ){
      kcal<-num_visit*(q-1)+w
      mit_qr[q,w,e]=mit[kcal,e]
    }
  }
}

#generate V allocation
c_latent_sim=array(0, dim = c(V, Q))
for(g in 1:V){
  k=sample(1:Q,1)
  c_latent_sim[g,k]=1
}

#generate A
noise=3
A <- optimizeArray(n, num_visit, V, Q, 
                   c_latent_sim, c(mit_qr), noise)
dim(A) <- c(n, num_visit, V, V)


################ main function
model<-function( A, X, mcmc_samples, Q){
  
  # X: design matrix for covariates(N*3:SNP,visit,SNP*visit) 
  # A: array of connectivity matrices (n*T*V*V)
  # mcmc_samples: number of MCMC iterations
  # Q: number of subnetworks
  
  N=dim(A)[1]*dim(A)[2] 
  pd <- dim(X)[2] 
  n=dim(A)[1] 
  num_visit<- dim(A)[2]
  V <- dim(A)[3]
  conn_block1 <- ((Q+1)*Q)/2
  A_f<-c(A)
  
  # initialization
  alpha <- array(NA, dim = c(mcmc_samples, conn_block1,pd))
  #alpha[1,,]<-array(rnorm(conn_block1*3),dim=c(conn_block1,3))
  alpha[1,,]<-alpha_sim
  
  delta<-array(NA, dim = c(mcmc_samples, conn_block1,pd))
  #delta[1,,]<-array(rep(1,conn_block1*3),dim=c(conn_block1,3))
  delta[1,,]<-delta_sim
  
  beta <- array(NA, dim = c(mcmc_samples, conn_block1,pd))
  beta[1,,]<-alpha[1,,]*delta[1,,]
  
  sigma2beta <- 1000
  sigma2alpha=1000
  
  b <- array(0, dim = c(mcmc_samples, conn_block1, n,2))

  D_ini<-array(c(1,0.3*1,0.3*1,1),dim=c(2,2))
  D<-array(0,dim=c(mcmc_samples,conn_block1,2,2))
  for(k in 1:conn_block1){
    D[1,k,,]<-D_ini
  }
  
  for(i in 1:conn_block1){
    b[1,i,,]=rmvnorm(n, sigma=D_ini)
  }
  
  sigma2 <- array(0.1, dim = c(mcmc_samples, conn_block1))
  
  sigma2qr_matrix <- array(1, dim = c(mcmc_samples, Q, Q)) 
  sigma2qr_matrix[1,,] <- matrix(1, nrow = Q, ncol = Q)
  a0 <- 1
  b0 <- 1
    
  paic <- matrix(0, nrow = mcmc_samples, ncol = Q)
  paic_sp_pos <- rep(0, times = Q)
  dir_vec <-  rep(3, times = Q)
  dir_vec_pos<- rep(0, times = Q)
  paic[1,]<- rdirichlet(1, dir_vec)
  
  c_latent <- array(0, dim = c(mcmc_samples, V, Q))
  # for(g in 1:V){
  #     k=sample(1:Q,1)
  #     c_latent[1,g,k]<-1
  #   }
  c_latent[1,,]<-c_latent_sim
  
  mik_qr_ini <- array(0, dim=c(n, num_visit, conn_block1))
  for (i in 1:n){
    for (j in 1:num_visit){
      for (s_subnetwork in 1:Q){
        for (t_subnetwork in s_subnetwork:Q){
          s_region<- which(c_latent[1,,s_subnetwork]==1)
          t_region<- which(c_latent[1,,t_subnetwork]==1)
          min_l_p<- min(s_subnetwork,t_subnetwork)
          max_l_p<- max(s_subnetwork,t_subnetwork)
          k_cal<- (min_l_p-1)*Q+ max_l_p- (min_l_p*(min_l_p-1)/2)
          mik_qr_ini[i,j,k_cal]<-  mean(A[i,j,s_region,t_region], na.rm = T)
        }
      }
    }
  }
  
  mik_qr<- array(0, dim = c(mcmc_samples, n, num_visit, conn_block1))
  mik_qr[1,,,] <- mik_qr_ini
  
  mik_ini <- array(0, dim = c(N, conn_block1))
  for (k in 1:conn_block1){
    v=c(t(mik_qr_ini[,,k]))
    mik_ini[,k]=v[v!=0]
  }
  
  mik <- array(0, dim = c(mcmc_samples, N, conn_block1))
  for(k in 1:mcmc_samples){
    mik[k, , ] <- mik_ini
  }
  
  Z=array(0,dim=c(n,num_visit,2))
  for(i in 1:n){
    Z[i,,]<-cbind(1,c(1:num_visit))
  }
  
  #MCMC
  for (i in 2:mcmc_samples){
    print(paste0("Iteration: ", i))
    matrix_list_z1 <- lapply(1:conn_block1, function(k) {
        lapply(1:n, function(e) {
          Z[e,,] %*% D[i-1,k,,] %*% t(Z[e,,]) 
        })
      })
      sigmare <- lapply(1:conn_block1, function(k) {
        lapply(1:n, function(e) {
          matrix_list_z1[[k]][[e]]+ diag(sigma2[(i-1),k], num_visit)
        })
      })
      
      inverse_block_diag_matrices <- lapply(sigmare, function(sublist) solve(bdiag(sublist)))
      
      # update alpha and delta
      t_veca0=mik[(i-1),,]-X[,2:3]%*%t(beta[i-1,,2:3])
      t_veca1=t_veca0-outer(X[,1],alpha[i-1,,1],"*")
      aaa=log_h2(t_veca1,inverse_block_diag_matrices)
      bbb=log_h2(t_veca0,inverse_block_diag_matrices)
      pa=exp(aaa-(log(exp(aaa)+exp(bbb))))
      delta[i,,1]<-ifelse(is.na(rbern(n=conn_block1,prob=pa)),delta[i-1,,1],rbern(n=conn_block1,prob=pa))
      
      t_vecb0=mik[(i-1),,]-X[,c(1,3)]%*%t(beta[i-1,,c(1,3)])
      t_vecb1=t_vecb0-outer(X[,2],alpha[i-1,,2],"*")
      aaa=log_h2(t_vecb1,inverse_block_diag_matrices)
      bbb=log_h2(t_vecb0,inverse_block_diag_matrices)
      pb=exp(aaa-(log(exp(aaa)+exp(bbb))))
      delta[i,,2]<-ifelse(is.na(rbern(n=conn_block1,prob=pb)),delta[i-1,,2],rbern(n=conn_block1,prob=pb))
      
      t_vecc0=mik[(i-1),,]-X[,c(1,2)]%*%t(beta[i-1,,c(1,2)])
      t_vecc1=mik[(i-1),,]-X %*% t(alpha[i-1,,])
      aaa=log_h2(t_vecc1,inverse_block_diag_matrices)
      bbb=log_h2(t_vecc0,inverse_block_diag_matrices)
      pc=exp(aaa-(log(exp(aaa)+exp(bbb))))
      delta[i,,3]<-ifelse(is.na(rbern(n=conn_block1,prob=pc)),delta[i-1,,3],rbern(n=conn_block1,prob=pc))
      delta[i,unique(c(which(delta[i,,2]==0),which(delta[i,,1]==0))),3]<-0
      
      p=sapply(1:conn_block1, function(k) which(delta[i,k,]==1))
      q=sapply(1:conn_block1, function(k) which(delta[i,k,]==0))
      for (k in 1:conn_block1){
        alpha[i,k,unlist(q[k])]=0
        if(length(unlist(p[k]))!=0){
          sigma_pos_inv<-t(X[,unlist(p[k])])%*%inverse_block_diag_matrices[[k]]%*%X[,unlist(p[k])]+solve(diag(sigma2alpha,length(unlist(p[k]))))
          sigma_pos<- solve(sigma_pos_inv)
          t_vec<- mik[(i-1),,k]
          mu_pos<- (sigma_pos%*%t(X[,unlist(p[k])])%*%inverse_block_diag_matrices[[k]]%*%t_vec)
          alpha[i,k,unlist(p[k])]<- as.vector(mvrnorm(n=1,mu=mu_pos,Sigma=sigma_pos))
        }
      }
      
      # update beta
      beta[i,,]<-delta[i,,]*alpha[i,,]
      if(i%%10==0){print(beta[i,,])}
      
      #update D
      for (k in 1:conn_block1){
        D[i,k,,]<-solve(rwish(1,2+n,solve(matrix(diag(2),2,2)+crossprod(b[(i-1),k,,]))))
      }
      
      #update b
      tvecc<-array(0,dim=c(conn_block1,N))
      for (k in 1:conn_block1){
        tvecc[k,]<- mik[(i-1),,k]- X%*%beta[i,k,]
      }
    
      for(k in 1:conn_block1){
        re=array(tvecc[k,],dim=c(num_visit,n))
        for(r in 1:n){
          sigma_pos_inv<-crossprod(Z[r,,])/sigma2[(i-1),k]+solve(D[i,k,,])
          sigma_pos<- solve(sigma_pos_inv)
          mu_pos<- (sigma_pos%*%t(Z[r,,])%*%re[,r]/sigma2[(i-1),k])
          b[i,k,r,]<- as.vector(mvrnorm(n=1,mu=mu_pos,Sigma=sigma_pos))
        }
      }
      
      #update sigma2_qr
      rel=array(0,dim=c(N, conn_block1))
      for(m in 1:conn_block1){
        ree=array(0,dim=c(num_visit, n))
        for(k in 1:n){
          ree[,k]=Z[k,,]%*%(b[i,m,k,])
        }
        rel[,m]<-c(ree)
      }
      
      for (k in 1:conn_block1){
        rate<- sum((tvecc[k,]-rel[,k])^2)
        sigma2[i,k]<- 1/rgamma(1,(0.001+(N/2)),(0.001+(rate/2)))
      }
      
      #update pi
      dir_vec_pos <- colSums(c_latent[(i-1),,] + dir_vec)
      paic[i,]<- rdirichlet(1, dir_vec_pos)
      
      #update V allocation
      c_latent[i,,]<- c_latent[(i-1),,]
      l_subnetwork <- max.col(c_latent[i, , ] == 1, "first")
      min_l_p<-outer(1:Q,l_subnetwork, FUN = pmin)
      max_l_p<-outer(1:Q,l_subnetwork, FUN = pmax)
      k_cal <- (min_l_p - 1) * Q + max_l_p- (min_l_p* (min_l_p - 1) / 2)
      
      mik_qr_f<-c(mik_qr[i-1,,,])
      sigma2qr_matrix_m<-sigma2qr_matrix[i-1,,]
      paic_v<-paic[i,]
      c_latent_old<-c_latent[(i-1),,]
      min_l_p_s<-min_l_p-1
      max_l_p_s<-max_l_p-1
      k_cal1 <- k_cal-1
      
      ccc<-c_latent_new( V, Q, max_l_p_s,min_l_p_s,
                         k_cal1, A_f,mik_qr_f,  sigma2qr_matrix_m,
                         paic_v, num_visit, n,c_latent_old) 
      for (s in 1:V){
        c_latent[i,s,]<- t(rmultinom(1, size = 1, prob = ccc[s,]))
      }
      
      #update mit_qr
      for (p in 1:Q){
        for (q in p:Q){
          if (p!=q){
            latent_p_indices <- which(c_latent[i,,p] == 1)
            latent_q_indices <- which(c_latent[i,,q] == 1)
            sum_mpq <- rowSums(A[, , latent_p_indices, latent_q_indices], dims = 2)
            
            num_mpq <- rep(length(latent_p_indices) * length(latent_q_indices),num_visit)
            k_cal<- (p-1)*Q+ q- (p*(p-1)/2)
            sigma_pos_inv<- (num_mpq/sigma2qr_matrix[(i-1),p,q])+(1/sigma2[i,k_cal])
            sigma_pos<- solve(diag(sigma_pos_inv,num_visit))
            re=array(0,dim=c(num_visit,n))
            re<- (sapply(1:n, function(m) Z[m,,] %*% b[i,k_cal, m,]))
            v=X%*%beta[i,k_cal,]+c(re)
            v=array(v,dim=c(num_visit,n))
            mu_pos<- ((sum_mpq/sigma2qr_matrix[(i-1),p,q])+(t(v)/sigma2[i,k_cal]))%*%sigma_pos
            mik_qr[i,,,k_cal]<-t(sapply(1:n, function(u) rmvnorm(1, mean = mu_pos[u,], sigma = sqrt(sigma_pos))))
          } else {
            sum_mpq<-array(0,dim=c(n,num_visit))
            num_mpq<-rep(0,num_visit)
            latent_p_indices <- which(c_latent[i,,p] == 1)
            latent_p_indices_f <- lapply(latent_p_indices, function(e) latent_p_indices[latent_p_indices >= e])
            sum_mpq <- t(sapply(1:n, function(w) rowSums(sapply(latent_p_indices, function(u) 
              sapply(1:num_visit, function(v) sum(A[w, v, u, latent_p_indices_f[[which(latent_p_indices == u)]]]))
            ))))
            num_mpq <- num_mpq + sum(lengths(latent_p_indices_f))
            k_cal<- (p-1)*Q+ q- (p*(p-1)/2)
            sigma_pos_inv<- (num_mpq/sigma2qr_matrix[(i-1),p,q])+(1/sigma2[i,k_cal])
            sigma_pos<- solve(diag(sigma_pos_inv,num_visit))
            re=array(0,dim=c(num_visit,n))
            re<- (sapply(1:n, function(m) Z[m,,] %*% b[i,k_cal, m,]))
            v=X%*%beta[i,k_cal,]+c(re)
            v=array(v,dim=c(num_visit,n))
            mu_pos<- ((sum_mpq/sigma2qr_matrix[(i-1),p,q])+(t(v)/sigma2[i,k_cal]))%*%sigma_pos
            mik_qr[i,,,k_cal]<-t(sapply(1:n, function(u) rmvnorm(1, mean = mu_pos[u,], sigma = sqrt(sigma_pos))))
          }
        }
      }
      
      mik[i,,1:conn_block1] <- apply(mik_qr[i,,,1:conn_block1], c(1, 3), function(x) x[x != 0])
      
      # update sigma2_epsilon_qr
      for (p in 1:Q){
        for (q in p:Q){
          e1 <- which(c_latent[i,,p] == 1)
          f1<- which(c_latent[i,,q] == 1)
          k_cal<- (p-1)*Q+ q- (p*(p-1)/2)
          sum_sigm2pq<- 0
          num_sigm2pq<- 0
          if (p!=q){
            sum_sigm2pq <- sum(
              sapply(e1, function(e) {
                rowSums(sapply(1:num_visit, function(s) {
                  rowSums(sapply(f1, function(f) {
                    (A[, s, e, f] - mik_qr[i, , s, k_cal])^2
                  }))
                }))
              })
            )
            num_sigm2pq <- length(e1) * length(f1) * n * num_visit
            sigma2qr_matrix[i,p,q]<- 1/rgamma(1,(a0+(num_sigm2pq/2)),(b0+(sum_sigm2pq/2)))
          } else {
            k_cal<- (p-1)*Q+ q- (p*(p-1)/2)
            latent_q_indices <- which(c_latent[i,,q] == 1)
            latent_q_indices_f <- lapply(latent_q_indices, function(e) latent_q_indices[latent_q_indices >= e])
            num_sigm2pq <- sum(lengths(latent_q_indices_f)) * n * num_visit
            sum_sigm2pq <- sum(
              sapply(latent_q_indices, function(e) {
                rowSums(sapply(1:num_visit, function(s) {
                  rowSums(sapply(latent_q_indices_f[[which(latent_q_indices == e)]], function(f) {
                    (A[, s, e, f] - mik_qr[i, , s, k_cal])^2
                  }))
                }))
              })
            )
            sigma2qr_matrix[i,p,q]<- 1/rgamma(1,(a0+(num_sigm2pq/2)),(b0+(sum_sigm2pq/2)))
          }
        }
      }
      
  }
  return(list(mik,beta,c_latent,delta,alpha,D,sigma2,sigma2qr_matrix,mik_qr,b))
}

######## run model
mcmc_samples=5000
r1=model( A, X, mcmc_samples, Q)


