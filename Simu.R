library(ald)
library(parallel)
library(doParallel)
library(optimParallel)
require(xtable)
require(magrittr)
library(nloptr)



#setwd("C:/Users/bdeketelaere/OneDrive - UCL/Documents/Sauvegardessimulations/Contrainte41")
#setwd("C:/Users/bdeketelaere/OneDrive - UCL/Documents/Code thèse/Thèse")
setwd("C:/Users/bdeketelaere/OneDrive - UCL/Documents/Sauvegardessimulations/Papier1apres11nov2023/5comp")



opts <- list( "algorithm"= "NLOPT_LN_COBYLA",
              "ftol_rel"= 0.005,
              "maxeval"= 5000,
              "print_level" = 0,
              maxtime=100)


Simulation5.3 <- function(M,n,lc,tau,K=4,C2=0.1,mmax=4,no_cores,seed=1){
set.seed(seed)
  loglik.IC.beta.general2.sigma.subfun <- function(param2){
    m=length(theta)-1;m_tilde=length(theta_tilde)-1

    lenbeta=ncol(as.matrix(X))+1
    beta=param2[1:lenbeta]
    sigma=param2[lenbeta+1]

    theta=c(1,param2[(lenbeta+2):(lenbeta+2+m-1)])
    theta_tilde=c(1,param2[(lenbeta+2+m):(lenbeta+2+m+m_tilde-1)])
    return(-1*loglik.IC.beta.general.sigma(beta,sigma,theta,theta_tilde,Y,tau,X))
  }


  double_constraint_fun_sigma.subfun <- function(x){
    if (m != 0 & m_tilde != 0) {
      part1 <- ((1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]))^2)/(1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]^2))
      part2 <- ((1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]))^2)/(1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]^2))
      common_term <- part1 - part2
      return(c(common_term,- common_term))
    } else if (m == 0 & m_tilde == 0) {
      return(c(0,0))
    } else if (m == 0 & m_tilde != 0) {
      part1 <- 1
      part2 <- ((1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]))^2)/(1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]^2))
      common_term <- part1 - part2
      return(c(common_term,- common_term))
    } else if (m != 0 & m_tilde == 0) {
      part1 <- ((1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]))^2)/(1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]^2))
      part2 <- 1
      common_term <- part1 - part2
      return(c(common_term,- common_term))
    }
  }
  
#Test pour 2 covar.

beta.CV4.01=matrix(NA,nrow=M,ncol=4)#derniere colonne de beta est sigma

eval_f=loglik.IC.beta.general2.sigma.subfun

# gamma=c(1);lengamma=length(gamma)
# beta=c(0,1,1);lenbeta=length(beta)


# opts <- list( "algorithm"= "NLOPT_LN_COBYLA",
#               "ftol_rel"= 1.0e-3,
#               "maxeval"= 5000,
#               #              "local_opts" = local_opts,
#               "print_level" = 3,
#               maxtime=100)


# tau=0.5
# n=100
# lc=0.5
# K=4
nK=n/K
# C2=0.1 #Constante CV

mtildemax=mmax

CV.crit=rep(NA,mmax+1)

count4.01=matrix(0,nrow=mmax+1)

BESTM4.01=rep(NA,M)


cl <- makeCluster(no_cores-1)  
registerDoParallel(cl)

#Déterminer quel est une valeur un peu optimale pour le maximum d'itération
clusterExport(cl=cl,list("Norm",
                         "check_function",
                         "ELD.cdf2.cont.sigma",
                         "loglik.IC.beta.general.sigma",
                         "loglik.IC.beta.general2.sigma.subfun",
                         "nloptr",
                         "double_constraint_fun_sigma.subfun",
                         "opts"))

for (kk in 1:M) {
  
  tic()
  set.seed(kk)
  
  eps=log(rweibull(n,shape=1.5,scale=2))-log(qweibull(p = tau,shape=1.5,scale=2))
  
  X1full=rbinom(n=n,1,prob=0.5)
  X2full=rbinom(n=n,1,prob=0.5)
  Xfull=cbind(X1full,X2full)
  Yfulluncen=X1full+X2full+eps
  Yfull=tranform.IC(Yfulluncen,l=0.75*lc,seed=kk)
  
  sum(Yfull[,2] > Yfull[,1]); mean(Yfull[,2]-Yfull[,1]);sum(Yfull[,2]-Yfull[,1]>0.1)
  
  #nK est le nombre d'observation dans un fold. Il se calcule comme n/K arrondi vers le haut (ceiling)
  #La k-ième ligne de fold_indices contient les indices des observations qui appartiennent au k-ième fold
  #Si K \not| n, alors dans le dernier fold il y a moins que nk observations. 
  #fold_indices laissera les valeurs pour ces observations à 0.  
  fold_indices <- matrix(0,ncol=nK,nrow=K)
  ind <- sample(1:n,n)
  for(i in 1:(K-1)) {
    fold_indices[i,1:nK] <- ind[((i-1)*nK+1):(i*nK)] #rempli les indices pour les K-1 premiers fold
  }
  fold_indices[K,1:(n-(K-1)*nK)] <- ind[((K-1)*nK+1):n] #rempli les indices pour le dernier fold
  
  
  simuler_points_contrainte_beta2V2.subfun <- function(n,m,m_tilde,sigmastart) {
    # Initialiser une matrice pour stocker les points
    points <- matrix(NA, nrow=n, ncol=length(beta)+1+m+m_tilde)
    
    # Compter le nombre de points satisfaisant la contrainte
    nb_points_satisfaisants <- 0
    
    #sera utilisé apres pour les theta qui sont généralement de plu en plus petits
    suite_geom=numeric(m)
    for (j in 1:m) {
      suite_geom[j] <- 1 * 0.75^(j-1)
    }
    suite_geom_2=numeric(m_tilde)
    for (j in 1:m_tilde) {
      suite_geom_2[j] <- 1 * 0.75^(j-1)
    }
    
    # Tant que nous n'avons pas suffisamment de points
    while (nb_points_satisfaisants < n) {
      # Générer aléatoirement des valeurs pour x1, x2 et x3
      x1 <- runif(1, min=-0.25, max=0.25)
      x2 <- runif(2, min=0.75, max=1.25)
      x3 <- sigmastart#Si sigma trop petit probleme de minimisation
      x4 <- runif(m, min=-0.75, max=1)* suite_geom
      x5 <- runif(m_tilde, min=-0.75, max=1) * suite_geom_2
      x=c(x1,x2,x3,x4,x5)
      
      # Vérifier si les valeurs satisfont la contrainte
      part1 <- ((1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]))^2)/(1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]^2))
      part2 <- 1
      common_term <- part1 - part2
      
      if ( abs(common_term) < 0.05) {
        # Si oui, stocker les valeurs dans la matrice de points
        nb_points_satisfaisants <- nb_points_satisfaisants + 1
        points[nb_points_satisfaisants,] <- c(x1, x2, x3,x4,x5)
      }
    }
    
    return(points)
  }
  

  for (m in 2:mmax) {
    theta=c(1,rep(1,m))
    m_tilde=m
    theta_tilde=c(1,rep(1,m_tilde))
    
    
    STARTS1=simuler_points_contrainte_beta2V2.subfun(1,m,m_tilde,sigmastart=0.2) 
    STARTS2=simuler_points_contrainte_beta2V2.subfun(1,m,m_tilde,sigmastart=0.3) 
    STARTS3=simuler_points_contrainte_beta2V2.subfun(1,m,m_tilde,sigmastart=0.4) 
    STARTS4=simuler_points_contrainte_beta2V2.subfun(1,m,m_tilde,sigmastart=0.5) 
    STARTS=rbind(STARTS1,STARTS2,STARTS3,STARTS4)
    
    nrowstarts=nrow(STARTS)
    
    # opts <- list( "algorithm"= "NLOPT_LN_COBYLA",
    #               "ftol_rel"= 0.005,
    #               "maxeval"= 5000,
    #               "print_level" = 3,
    #               maxtime=100)
    
    clusterExport(cl=cl,list("Yfull",
                             "Xfull",#ajouter Xfull ici si covariable
                             "beta",
                             "lenbeta",
                             "lengamma",
                             "m",
                             "m_tilde",
                             "STARTS",
                             "theta",
                             "theta_tilde",
                             "tau",
                             "opts"))
    beta0.best=rep(NA,K)
    beta1.best=rep(NA,K)
    beta2.best=rep(NA,K)
    gamma0.best=rep(NA,K)
    pred.error=rep(NA,K)
    
    for (k in 1:K) {
      #current_indices contient les indices des observations qui vont servir à estimer le modèle
      current_indices <- setdiff(1:n,fold_indices[k,]) 
      Y=Yfull[current_indices,] 
      X=Xfull[current_indices,] 
      
      Ytest=Yfull[-current_indices,]
      Xtest=Xfull[-current_indices,]
      lengthYtest=nrow(Ytest)
      
      clusterExport(cl=cl,list("Y","X","nloptr","eval_f","double_constraint_fun_sigma.subfun","opts"
      ))
      res.par=unlist(
        parLapply(cl,X=1:nrow(STARTS),fun=function(j){
          optim.temp=unlist(
            nloptr ( x0 = STARTS[j,],
                     eval_f = eval_f,
                     eval_g_ineq = double_constraint_fun_sigma.subfun,
                     opts = opts,
                     lb=c(rep(-Inf,length(beta)) ,0.05,rep(-2,m),rep(-2,m_tilde)),
                     ub=c(rep(+Inf,length(beta)) ,0.8,rep(4,m),rep(4,m_tilde))
            )[c("objective","solution","status","iterations")]
          )
          val.temp.par=optim.temp[1]
          coef.temp.par=optim.temp[2:(m+m_tilde+2+length(beta))]
          status.temp.par=optim.temp[(m+m_tilde+3+length(beta))]
          iterations.temp.par=optim.temp[(m+m_tilde+4+length(beta))]
          return(c(val.temp.par,coef.temp.par,status.temp.par,iterations.temp.par))
        }
        )
      )
      indexval=seq(from=1,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
      indexbeta0=seq(from=2,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
      indexbeta1=seq(from=3,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
      indexbeta2=seq(from=4,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
      indexgamma0=seq(from=5,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
      indexstatus=seq(from=2+lenbeta+lengamma+m+m_tilde,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
      indexiterations=seq(from=3+lenbeta+lengamma+m+m_tilde,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
      index.min.temp=which.min(res.par[indexval])
      beta0.best[k]=res.par[indexbeta0][index.min.temp]
      beta1.best[k]=res.par[indexbeta1][index.min.temp]
      beta2.best[k]=res.par[indexbeta2][index.min.temp]
      gamma0.best[k]=res.par[indexgamma0][index.min.temp]
      # res.par[indexstatus]
      # res.par[indexiterations]
      ## Compute sum of the evaluations of the check function (see equation 5)
      pred.error[k]=mean( 1/  apply(cbind(Ytest[,2]-Ytest[,1],rep(C2,lengthYtest)),MARGIN = 1,FUN=max)  *  (check_function((Ytest[,1]/2 + Ytest[,2]/2) - beta0.best[k] - beta1.best[k]%*%Xtest[,1] - beta2.best[k]%*%Xtest[,2] , tau )))
    }
    CV.crit[m+1]=mean(pred.error) 
    if(kk==1){print(c(m,m_tilde))}
  }
  
  BESTM4.01[kk]=which.min(CV.crit[3:(mmax+1)])+1
  
  
  m=BESTM4.01[kk]
  m_tilde=m
  
  theta=c(1,rep(1,m))
  theta_tilde=c(1,rep(1,m_tilde))
  
  Y=Yfull
  X=Xfull
  
  STARTS1=simuler_points_contrainte_beta2V2.subfun(2,m,m_tilde,sigmastart=0.2) 
  STARTS2=simuler_points_contrainte_beta2V2.subfun(2,m,m_tilde,sigmastart=0.3) 
  STARTS3=simuler_points_contrainte_beta2V2.subfun(2,m,m_tilde,sigmastart=0.4) 
  STARTS4=simuler_points_contrainte_beta2V2.subfun(2,m,m_tilde,sigmastart=0.5) 
  STARTS=rbind(STARTS1,STARTS2,STARTS3,STARTS4)
  nrowstarts=nrow(STARTS)
  
  # opts <- list( "algorithm"= "NLOPT_LN_COBYLA",
  #               "ftol_rel"= 0.005,
  #               "maxeval"= 5000,
  #               "print_level" = 3,
  #               maxtime=100)
  
  clusterExport(cl=cl,list("Y",#ajouter X ici si cov
                           "X",
                           "beta",
                           "m",
                           "m_tilde",
                           "STARTS",
                           "theta",
                           "theta_tilde",
                           "tau",
                           "opts"))
  
  res.par=unlist(
    parLapply(cl,X=1:nrow(STARTS),fun=function(j){
      optim.temp=unlist(
        nloptr ( x0 = STARTS[j,],
                 eval_f = eval_f,
                 eval_g_ineq = double_constraint_fun_sigma.subfun,
                 opts = opts,
                 lb=c(rep(-Inf,length(beta)) ,0.05,rep(-3,m),rep(-3,m_tilde)),
                 ub=c(rep(+Inf,length(beta)) ,1,rep(3,m),rep(3,m_tilde))
        )[c("objective","solution","status","iterations")]
      )
      val.temp.par=optim.temp[1]
      coef.temp.par=optim.temp[2:(m+m_tilde+2+length(beta))]
      status.temp.par=optim.temp[(m+m_tilde+3+length(beta))]
      iterations.temp.par=optim.temp[(m+m_tilde+4+length(beta))]
      return(c(val.temp.par,coef.temp.par,status.temp.par,iterations.temp.par))
    }
    )
  )
  indexval=seq(from=1,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
  indexbeta0=seq(from=2,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
  indexbeta1=seq(from=3,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
  indexbeta2=seq(from=4,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
  indexgamma0=seq(from=5,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
  indexstatus=seq(from=2+lenbeta+lengamma+m+m_tilde,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
  indexiterations=seq(from=3+lenbeta+lengamma+m+m_tilde,by=3+lenbeta+lengamma+m+m_tilde,length.out = nrowstarts)
  index.min.temp=which.min(res.par[indexval])
  beta.CV4.01[kk,1] = res.par[indexbeta0][index.min.temp]
  beta.CV4.01[kk,2] = res.par[indexbeta1][index.min.temp]
  beta.CV4.01[kk,3] = res.par[indexbeta2][index.min.temp]
  beta.CV4.01[kk,4] = res.par[indexgamma0][index.min.temp]
  
  #toc()
  
  count4.01[which.min(CV.crit)]=count4.01[which.min(CV.crit)]+1
  
  
  if(kk == 1){toc();print(beta.CV4.01[1,])}
  if(kk == 2){toc();print(beta.CV4.01[2,])}
  if(kk == 3){toc();print(beta.CV4.01[3,])}
  if(kk%%10==0){print(kk);toc()}
  
  nomfichiertemp=(c("Simu53Laguerre", "tau=", 100*tau,"n=",n,"lc=",10*lc,"M=",M,"seed=",seed))
  if(kk%%50==0){write.table(cbind(beta.CV4.01,BESTM4.01),paste(nomfichiertemp,collapse = ""))}
  
}
  return(list(count4.01,beta.CV4.01))
}



################################################################################
################################################################################
M=5
n=100
lc=1
tau=0.5
K=4
C2=0.1 #Constante CV
no_cores <- detectCores(logical = TRUE)
seed=1

A1=Simulation5.3(M,n,lc,tau,K=4,C2=0.1,mmax=4,no_cores,seed)
################################################################################
M=500
n=200
lc=1
tau=0.5
K=4
C2=0.1 #Constante CV
no_cores <- detectCores(logical = TRUE)
seed=1

A2=Simulation5.3(M,n,lc,tau,K=4,C2=0.1,mmax=4,no_cores,seed)
################################################################################
M=500
n=400
lc=1
tau=0.5
K=4
C2=0.1 #Constante CV
no_cores <- detectCores(logical = TRUE)
seed=1

A3=Simulation5.3(M,n,lc,tau,K=4,C2=0.1,mmax=4,no_cores,seed)
################################################################################
M=500
n=100
lc=2
tau=0.5
K=4
C2=0.1 #Constante CV
no_cores <- detectCores(logical = TRUE)
seed=1

A4=Simulation5.3(M,n,lc,tau,K=4,C2=0.1,mmax=4,no_cores,seed)
################################################################################
M=500
n=200
lc=2
tau=0.5
K=4
C2=0.1 #Constante CV
no_cores <- detectCores(logical = TRUE)
seed=1

A5=Simulation5.3(M,n,lc,tau,K=4,C2=0.1,mmax=4,no_cores,seed)
################################################################################
M=500
n=400
lc=2
tau=0.5
K=4
C2=0.1 #Constante CV
no_cores <- detectCores(logical = TRUE)
seed=1

A6=Simulation5.3(M,n,lc,tau,K=4,C2=0.1,mmax=4,no_cores,seed)
################################################################################
