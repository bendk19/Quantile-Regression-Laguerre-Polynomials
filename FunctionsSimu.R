
#install.packages("pracma")
library(pracma)#utilis? pour calcul? la norme d'un vecteur
#install.packages("qlcMatrix")
library(qlcMatrix)
#install.packages("xtable")
library(xtable)
#install.packages("GoFKernel")
library(GoFKernel)










## Implementation of the check function
check_function <- function(z,tau) {
  return(z*(tau-as.numeric(z<=0)))
}






##Transforme des observations d'une variable aleatoire (univariee) en des donnees
##censurees par intervalles. Les parametres sont adaptés a une Weibull

tranform.IC <- function(Time,l=1,seed){
  set.seed(seed)
  n=length(Time)
  L=R=rep(NA,n)
  for (i in 1:n) {
    set.seed(n*seed+i)
    incr=runif(5000,min=0,max=2*l)
    insptimetemp=rep(NA,5000)
    insptimetemp[1]=-15+incr[1]
    dd=1
    while(insptimetemp[dd] < 50){
      insptimetemp[dd+1]=insptimetemp[dd]+incr[dd]
      dd = dd+1
      if(dd == 4999){insptimetemp[5000] = 50}
    }
    insptime=insptimetemp[1:dd]
    R[i]=insptime[which.max(Time[i]<insptime)]
    L[i]=insptime[which.max2(Time[i]>insptime)]
  } 
  return(cbind(L,R))
}




transform.IC.extrval.V3 <- function(Time, lmin = 0.1, lmax = 1.9, seed, ximin = -10, ximax = 10) {
  set.seed(seed)
  n <- length(Time)
  L <- R <- numeric(n)
  
  interval_length <- ximax - ximin
  max_points <- floor(interval_length / lmin)#il faut l'arrondir vers le haut
  
  for (i in seq_len(n)) {
    set.seed(n * seed + i)
    
    # Increments=xi-i
    increments <- runif(max_points, min = lmin, max = lmax)
    
    # Yk,i
    insptime <- cumsum(c(ximin, increments))
    
    # il y aura qulques inspections time plus grandes que ximax, on les retire
    insptime <- insptime[insptime <= ximax]
    
    # L[i]=Première valeur plus grande que Time[i] et R[i]=Dernière valeur plus petite que Time[i]
    R[i] <- insptime[which(Time[i] < insptime)[1]]
    L[i] <- insptime[tail(which(Time[i] > insptime), 1)]
  }
  
  return(cbind(L, R))
}


which.min2 <- function(x, last.index = T, ...){
  if(last.index) max(which(x == min(x, ...))) else which.min(x)
}
# which.min2(c(1,1,1))
# which.min2(c(1,0,1))

which.max2 <- function(x, last.index = T, ...){
  if(last.index) max(which(x == max(x, ...))) else which.max(x)
}






ELD.cdf2.cont.sigma <- function(y, tau, mu, sigma, theta, theta_tilde) {
  # m and m_tilde are derived from the lengths of theta and theta_tilde
  m <- length(theta) - 1
  m_tilde <- length(theta_tilde) - 1
  n <- length(y)  # Total number of elements in y

  # Error handling for incorrect input dimensions
  if (length(mu) == 1) {
    print("Consider using ELD.cdf2.sigma instead")
  } else if (length(mu) != length(y)) {
    print("Dimension mismatch in mu or y"); return()
  }

  # Separate finite and infinite values in y
  yfinite <- y[is.finite(y)]
  yinfinite <- y[!is.finite(y)]

  # Use only finite values of y
  y <- yfinite
  n <- length(y)

  ################################### Start of Part 1 ###################################
  # Compute various coefficient matrices and arrays.
  # This section has a higher computational cost but is executed once.
  # Matrices Au and Ad store specific combinations of coefficients.

  # Generate matrix Au for upper bounds
  Au <- matrix(0, ncol = m + 1, nrow = m + 1)
  for (k in 1:(m + 1)) {
    Au[k, 1:k] <- choose(k - 1, 0:(k - 1)) * (-1)^(0:(k - 1)) / factorial(0:(k - 1))
  }

  # Generate matrix Ad for lower bounds
  Ad <- matrix(0, ncol = m_tilde + 1, nrow = m_tilde + 1)
  for (k in 1:(m_tilde + 1)) {
    Ad[k, 1:k] <- choose(k - 1, 0:(k - 1)) * (-1)^(0:(k - 1)) / factorial(0:(k - 1))
  }

  # Generate Bu and Bd matrices for factorial-based coefficients
  Bu <- matrix(0, ncol = 2 * m + 1, nrow = 2 * m + 1)
  for (j.j_prime in 0:(2 * m)) {
    Bu[j.j_prime + 1, (1:(j.j_prime + 1))] <- factorial(j.j_prime) / factorial(0:j.j_prime)
  }

  Bd <- matrix(0, ncol = 2 * m_tilde + 1, nrow = 2 * m_tilde + 1)
  for (j.j_prime in 0:(2 * m_tilde)) {
    Bd[j.j_prime + 1, (1:(j.j_prime + 1))] <- factorial(j.j_prime) / factorial(0:j.j_prime)
  }

  # Cu and Cd store combinations of Au and Ad coefficients
  Cu <- array(0, dim = c(m + 1, m + 1, m + 1, m + 1))
  for (k in 0:m) {
    for (k_prime in 0:m) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cu[j + 1, j_prime + 1, k + 1, k_prime + 1] <- Au[k + 1, j + 1] * Au[k_prime + 1, j_prime + 1]
        }
      }
    }
  }

  Cd <- array(0, dim = c(m_tilde + 1, m_tilde + 1, m_tilde + 1, m_tilde + 1))
  for (k in 0:m_tilde) {
    for (k_prime in 0:m_tilde) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cd[j + 1, j_prime + 1, k + 1, k_prime + 1] <- Ad[k + 1, j + 1] * Ad[k_prime + 1, j_prime + 1]
        }
      }
    }
  }

  # Du and Dd arrays store sums of coefficients with specific conditions
  Du <- array(0, dim = c(m + 1, m + 1, 2 * m + 1))
  for (k in 1:(m + 1)) {
    for (k_prime in 1:(m + 1)) {
      for (j_star in 1:(k + k_prime - 1)) {
        Du[k, k_prime, j_star] <-
          sum(Cu[,,k,k_prime][(row(as.matrix(Cu[,,k,k_prime])) + col(as.matrix(Cu[,,k,k_prime])) - 1) == j_star])
      }
    }
  }

  Dd <- array(0, dim = c(m_tilde + 1, m_tilde + 1, 2 * m_tilde + 1))
  for (k in 1:(m_tilde + 1)) {
    for (k_prime in 1:(m_tilde + 1)) {
      for (j_star in 1:(k + k_prime - 1)) {
        Dd[k, k_prime, j_star] <-
          sum(Cd[,,k,k_prime][(row(as.matrix(Cd[,,k,k_prime])) + col(as.matrix(Cd[,,k,k_prime])) - 1) == j_star])
      }
    }
  }

  ################################### Start of Part 2 ###################################

  # Prepare input data as a data frame to facilitate computation and maintain order
  data.temp <- data.frame(num = 1:n, y = y, mu = mu[1:length(y)])

  # Subset data for y >= mu
  data.up <- subset(data.temp, y >= mu[1:length(y)])
  num.up <- data.up$num
  yu <- data.up$y
  muu <- data.up$mu
  su <- (yu - muu) * tau / sigma
  lenu <- length(su)
  resu <- rep(NA, lenu)  # Result vector for y >= mu

  if (lenu != 0) {
    # Compute coefficients Eu, Gu, and Hu for y >= mu
    Eu <- matrix(NA, 2 * m + 1, lenu)
    for (l in 1:(2 * m + 1)) {
      Eu[l, ] <- (su)^(l - 1) * exp(-su)
    }

    Gu <- array(0, dim = c(2 * m + 1, 2 * m + 1, lenu))
    for (j.j_prime in 1:(2 * m + 1)) {
      for (l in 1:(j.j_prime)) {
        Gu[j.j_prime, l, ] <- Bu[j.j_prime, l] * Eu[l, ]
      }
    }

    Hu <- array(0, dim = c(2 * m + 1, lenu))
    for (i in 1:lenu) {
      Hu[, i] <- rowSums(Gu[, , i])
    }

    # Compute the result resu using precomputed arrays
    Ju <- array(0, dim = c(m + 1, m + 1, lenu))
    for (k in 1:(m + 1)) {
      for (k_prime in 1:(m + 1)) {
        Ju[k, k_prime, ] <- Du[k, k_prime, ] %*% Hu
      }
    }

    THETA <- (theta) %*% t(theta)
    Ku <- matrix(Ju[, , 1:lenu], ncol = (m + 1)^2, nrow = lenu, byrow = TRUE) *
      matrix(rep(THETA), byrow = TRUE, ncol = (m + 1)^2, nrow = lenu)

    if (m != 0) {
      resu <- 1 + (1 - tau) * -rowSums(Ku) / Norm(theta)^2
    } else if (m == 0) {
      resu <- 1 + ((1 - tau) * -Ku / Norm(theta)^2)
    }
  }

  # Similar process for y < mu
  data.down <- subset(data.temp, y < mu[1:length(y)])
  num.down <- data.down$num
  yd <- data.down$y
  mud <- data.down$mu
  sd <- -(yd - mud) * (1 - tau) / sigma
  lend <- length(sd)

  if (lend != 0) {
    # Compute coefficients Ed, Gd, Hd for y < mu
    resd <- rep(NA, lend)
    Ed <- matrix(NA, 2 * m_tilde + 1, lend)
    for (l in 1:(2 * m_tilde + 1)) {
      Ed[l, ] <- (sd)^(l - 1) * exp(-sd)
    }

    Gd <- array(0, dim = c(2 * m_tilde + 1, 2 * m_tilde + 1, lend))
    for (j.j_prime in 1:(2 * m_tilde + 1)) {
      for (l in 1:(j.j_prime)) {
        Gd[j.j_prime, l, ] <- Bd[j.j_prime, l] * Ed[l, ]
      }
    }

    Hd <- array(0, dim = c(2 * m_tilde + 1, lend))
    for (i in 1:lend) {
      Hd[, i] <- rowSums(Gd[, , i])
    }

    Jd <- array(0, dim = c(m_tilde + 1, m_tilde + 1, lend))
    for (k in 1:(m_tilde + 1)) {
      for (k_prime in 1:(m_tilde + 1)) {
        Jd[k, k_prime, ] <- Dd[k, k_prime, ] %*% Hd
      }
    }

    THETA_tilde <- (theta_tilde) %*% t(theta_tilde)
    Kd <- matrix(Jd[, , 1:lend], ncol = (m_tilde + 1)^2, nrow = lend, byrow = TRUE) *
      matrix(rep(THETA_tilde), byrow = TRUE, ncol = (m_tilde + 1)^2, nrow = lend)

    if (m_tilde != 0) {
      resd <- -tau * -rowSums(Kd) / Norm(theta_tilde)^2
    } else if (m_tilde == 0) {
      resd <- -tau * -Kd / Norm(theta_tilde)^2
    }
  }

  # Combine results and add 1s for infinite y values
  data.temp2 <- data.frame(num = c(num.down, num.up), res = c(resd, resu))
  c(data.temp2[order(data.temp2$num), ]$res, rep(1, length(yinfinite)))
}










loglik.IC.beta.general.sigma<-function(beta, sigma, theta,theta_tilde,Y,tau,X)
{
  n <- nrow(Y)
  X=cbind(rep(1,n),X)
  
  longY <- c(Y[,1],Y[,2])
  
  res = ELD.cdf2.cont.sigma(longY,tau,rep(X%*%beta,2),sigma,theta,theta_tilde)
  
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}


n=100
tau=0.5
lc=0.5
eps=log(rweibull(n,shape=1.5,scale=2))-log(qweibull(p = tau,shape=1.5,scale=2))

X1full=rbinom(n=n,1,prob=0.5)
X2full=rbinom(n=n,1,prob=0.5)
Xfull=cbind(X1full,X2full)
Yfulluncen=X1full+X2full+eps
Yfull=tranform.IC(Yfulluncen,l=0.75*lc,seed=1)

sum(Yfull[,2] > Yfull[,1]); mean(Yfull[,2]-Yfull[,1]);sum(Yfull[,2]-Yfull[,1]>0.1)

gamma=c(1);lengamma=length(gamma)
beta=c(0,1,1);lenbeta=length(beta)

m=2#test
m_tilde=m
theta=theta_tilde=c(1,rep(1,m))

Y=Yfull
X=Xfull

eval_f="loglik.IC.beta.general2.sigma.subfun"
STARTS=NULL


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




