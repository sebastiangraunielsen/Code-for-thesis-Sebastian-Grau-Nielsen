require(graphics)
library(nlme)
library(kerntools)
library(cubature)
library(calculus)
library(fda)
library(latex2exp)

source("Fourier_product.R")
par(mfrow=c(1,1))

# Read data
tempmin <- read.csv("Documents/data_min.csv", header=F)
tempmax <- read.csv("Documents/data_max.csv", header=F)

# Maximal temperature plot
ts.plot(tempmax, main='Daily maximal temperature in 2020 by prefecture', xlab='Day', ylab='Temperature', col='darkgrey')
lines(tempmax$V47, col='blue')
lines(tempmax$V1, col='red')
legend(200, 10, legend=c('Okinawa', 'Hokkaido'),lty=1, col=c('blue', 'red'))

# Minimal temperature plot
ts.plot(tempmin, main='Daily minimal temperature in 2020 by prefecture', xlab='Day', ylab='Temperature', col='darkgrey')
lines(tempmin$V47, col='blue')
lines(tempmin$V1, col='red')
legend(180, 5, legend=c('Okinawa', 'Hokkaido'),lty=1, col=c('blue', 'red'))

# Gaussian kernel for regular PCA
gaussian_kernel <- function(x, gamma = 0.1) {
  exp(-as.matrix(dist(x, method = "euclidean")^2) * gamma)
}

principal_gauss <- matrix(0, nrow=47, ncol=366)

for(i in 1:366){
  principal_gauss[,i] <- kPCA(gaussian_kernel(cbind(data.matrix(tempmax)[i,],data.matrix(tempmin)[i,]), gamma=0.01), center=F)$PC1

}

# Plot of principal components for Gaussian kernel (gamma=0.01)
principal_gauss[,1]
ts.plot(-1*t(principal_gauss), col='darkgrey',xlab='t',ylab='')
title(ylab=TeX(r'($\langle p_1,\Phi(x_i(t)\rangle$)'), mgp=c(2.4,1,0), cex.lab=1)

lines(-1*principal_gauss[47,], col='blue')
lines(-1*principal_gauss[1,], col='red')

# Fourier series fit using OLS
time_axis <- 1:366

# Least squares fit for minimal temperature
total_fit_1 <- lm(data.matrix(tempmin)~cos(2*pi*time_axis/366) + sin(2*pi*time_axis/366) +cos(2*2*pi*time_axis/366) + sin(2*2*pi*time_axis/366) + cos(3*2*pi*time_axis/366)
   + sin(3*2*pi*time_axis/366) + cos(4*2*pi*time_axis/366) + sin(4*2*pi*time_axis/366) + cos(5*2*pi*time_axis/366) + sin(5*2*pi*time_axis/366) + 
     cos(6*2*pi*time_axis/366) + sin(6*2*pi*time_axis/366) + cos(7*2*pi*time_axis/366) + sin(7*2*pi*time_axis/366) + cos(8*2*pi*time_axis/366) + 
     sin(8*2*pi*time_axis/366) + cos(9*2*pi*time_axis/366) + sin(9*2*pi*time_axis/366) + cos(10*2*pi*time_axis/366) + sin(10*2*pi*time_axis/366))

ts.plot(total_fit_1$fitted.values, main='Fitted Fourier series for minimal temperature', col='darkgrey')
lines(total_fit_1$fitted.values[,47], col='blue')
lines(total_fit_1$fitted.values[,1], col='red')

# Least squares fit for maximal temperature
total_fit_2 <- lm(data.matrix(tempmax)~cos(2*pi*time_axis/366) + sin(2*pi*time_axis/366) +cos(2*2*pi*time_axis/366) + sin(2*2*pi*time_axis/366) + cos(3*2*pi*time_axis/366)
                  + sin(3*2*pi*time_axis/366) + cos(4*2*pi*time_axis/366) + sin(4*2*pi*time_axis/366) + cos(5*2*pi*time_axis/366) + sin(5*2*pi*time_axis/366) + 
                    cos(6*2*pi*time_axis/366) + sin(6*2*pi*time_axis/366) + cos(7*2*pi*time_axis/366) + sin(7*2*pi*time_axis/366) + cos(8*2*pi*time_axis/366) + 
                    sin(8*2*pi*time_axis/366) + cos(9*2*pi*time_axis/366) + sin(9*2*pi*time_axis/366) + cos(10*2*pi*time_axis/366) + sin(10*2*pi*time_axis/366))

ts.plot(total_fit_2$fitted.values, main='Fitted Fourier series for maximal temperature', col='darkgrey')
lines(total_fit_2$fitted.values[,47], col='blue')
lines(total_fit_2$fitted.values[,1], col='red')

# Gaussian kernel PCA on fitted Fourier series
principal_gauss_fake <- matrix(0, nrow=47, ncol=366)

for(i in 1:366){
  principal_gauss_fake[,i] <- kPCA(gaussian_kernel(cbind(total_fit_1$fitted.values[i,],total_fit_2$fitted.values[i,]), gamma=0.01), center=F)$PC1
  
}

ts.plot(-1*t(principal_gauss_fake), col='darkgrey', xlab='t')
title(ylab=TeX(r'($\langle p_1,\Phi(x_i(t)\rangle$)'), mgp=c(2.4,1,0), cex.lab=1)
lines(-1*principal_gauss_fake[47,], col='blue')
lines(-1*principal_gauss_fake[1,], col='red')

# Fourier series function 
Fourier_evaluate <- function(x, coeffs){
  input_size <- length(x)
  n_coeffs <- length(coeffs)
  output <- numeric(input_size)
  
  for(i in 1:input_size){
    output[i] <- coeffs[1] + sum(sapply(seq(2,n_coeffs-1,2), function(n) coeffs[n]*cos((n-1)*x[i]))) + sum(sapply(seq(3, n_coeffs, 2), function(n) coeffs[n]*sin((n-1)*x[i])))
    
  }
  
  return(output)
}

# RKHM Gaussian kernel 
gaussian_kernel_RKHM <- function(x, x1, x2, y1, y2, gamma){
  exp(-gamma*((Fourier_evaluate(x, coeffs=x1) - Fourier_evaluate(x, coeffs=y1))^2
              + (Fourier_evaluate(x, coeffs=x2) - Fourier_evaluate(x, coeffs=y2))^2))
}

# Calculate L^2 projection onto Fourier series of order 10
Kernel_projection <- function(x1, x2, y1, y2, gamma=1){
  a0 <- 1/(2*pi) * pcubature(Vectorize(function(x) gaussian_kernel_RKHM(x, x1, x2, y1, y2, gamma)),lowerLimit = 0, upperLimit = 2*pi)$integral
  
  a_vec <- c(a0, rep(0, 10))
  b_vec <- numeric(10)
  coeff_vec <- c(a0, rep(0,20))
  
  for(i in 1:10){
    a_vec[i + 1] <- 1/pi*pcubature(Vectorize(function(x) gaussian_kernel_RKHM(x, x1, x2, y1, y2, gamma)*cos(i*x)),0, 2*pi)$integral
    
    b_vec[i] <- 1/pi*pcubature(Vectorize(function(x) gaussian_kernel_RKHM(x, x1, x2, y1, y2, gamma)*sin(i*x)),0, 2*pi)$integral 
  }
  
  for(i in 2:21){
    if(i%%2 == 0){
      coeff_vec[i] <- a_vec[i/2+1]
    }
    else{
      coeff_vec[i] <- b_vec[(i-1)/2]
    }
  }
  return(coeff_vec)
  
}

# Fill out Gram matrix
Gram_matrix <- array(0, dim=c(47,47,21))
counter <- 0 

# Lower triangle 
for(i in 1:47){
  for(j in 1:47){
    if(j <= i){
          
      Gram_matrix[i,j,] <- Kernel_projection(total_fit_1$coefficients[,i], total_fit_2$coefficients[,i], 
                                                                    total_fit_1$coefficients[,j], total_fit_2$coefficients[,j],gamma = 0.01)
      counter <- counter + 1
      print(counter/1128)
    }
  }
}

# Upper triangle 
for(i in 1:47){
  for(j in 1:47){
    if(j > i){
      Gram_matrix[i,j,] <- Gram_matrix[j,i,]
    }
  }
}

# Calculate coefficients of pointwise product and then remove high frequency terms
approximate_product <- function(coeffs1, coeffs2, n_terms){
  maxlen <- max(length(coeffs1), length(coeffs2))
  coeffs1 <- c(coeffs1, rep(0, maxlen - length(coeffs1)))
  coeffs2 <- c(coeffs2, rep(0, maxlen - length(coeffs2)))
  return(Fourier_product(coeffs1, coeffs2)[1:(2*n_terms+1)])
}

# Matrix multiplication using approximate product function
approximate_matrix_prod <- function(mat_1, mat_2, n_terms){
  m <- nrow(mat_1)
  n <- ncol(mat_2)
  s <- ncol(mat_1)
  k <- 2*n_terms + 1
  out_mat <- array(0, dim=c(m,n,k))
  
  for(i in 1:m){
    for(j in 1:n){
      for(l in 1:s){
        out_mat[i,j,] <- out_mat[i,j,] + approximate_product(mat_1[i,l,], mat_2[l,j,], n_terms)
      }
    }
  }
  return(out_mat)
}

# Outer product of two vectors using approximate product function 
approximate_outer_product <- function(vec_1, vec_2, n_terms){
  m <- nrow(vec_1)
  k <- 2*n_terms + 1
  out_mat <- array(0, dim=c(m,m,k))
  
  for(i in 1:m){
    for(j in 1:m){
      out_mat[i,j,] <- approximate_product(vec_1[i,,], vec_2[j,,], n_terms)
    }
  }
  
  return(out_mat)
}

approximate_inner_product <- function(vec_1, vec_2, n_terms){
  m <- nrow(vec_1)
  k <-  2*n_terms + 1
  output <- numeric(k)
  
  for(i in 1:m){
      output <- output + approximate_product(vec_1[i,,],vec_2[i,,], n_terms)
  }
  return(output)
}

dim(approximate_matrix_prod(Gram_matrix, Gram_matrix, 3))

# Computing the gradient 
gradient_compute <- function(G, c, lambda, n_terms){
  G_c <- approximate_matrix_prod(G, c, n_terms)
  G2_c <- approximate_matrix_prod(G,G_c,n_terms)
  Gc_Gc <- approximate_outer_product(G_c,G_c,n_terms)
  final <- approximate_matrix_prod(Gc_Gc, c,n_terms)
  
  return(-2*G2_c-4*lambda*G_c+4*lambda*final)
}

approximate_objective_function <- function(G, c, lambda, n_terms){
  G_c <- approximate_matrix_prod(G,c,n_terms)
  first_term <- -1*approximate_inner_product(G_c,G_c, n_terms)
  G_cc <- approximate_inner_product(G_c,c,n_terms)
  second_term <- lambda*approximate_product(G_cc,G_cc,n_terms)
  third_term <- -2*lambda * G_cc + lambda
  
  return(first_term + second_term + third_term)
}

# Do gradient descent (equation 14 in Hashimoto)
gradient_descent <- function(G, c, lambda, eta, n_terms){
  gradient_out <- gradient_compute(G, c, lambda, n_terms)
  gradient_out_padded <- array(0, dim=c(47, 1, 21))
  gradient_out_padded[,,1:(2*n_terms+1)] <- gradient_out

  return(c - eta * gradient_out_padded)
}

c_0 <- array(0, dim=c(47,1,21))
c_0[1:47,1,1] <- rep(1, 47)
c_0[1:47,1,1:7] <- t(total_fit_1$coefficients)[,1:7]

c_0[1:47,1,2] <- rep(1,47)
c_0[1:47,1,3] <- rep(1,47)
c_0[1:47,1,4] <- rep(1,47)
c_0[1:47,1,5] <- rep(1,47)

lambda <- 1e-3
eta <- 1e-4
n_iter <- 10
n_terms <- 3

# c_0[1:47,1,1] <- principal_gauss_fake[,1]
 
# 0.001 og 0.0001 virker med konstanten 1 
# 1+ cos(4*pi*x/366) fungerer med 1e-4*2 og 1e-5*2
# c_0[1:47,1,1] <- colMeans(tempmax)

objective_mat <- matrix(0, nrow=10, ncol=7)
gradient_run <- gradient_descent(Gram_matrix,c_0,lambda, eta, n_terms)
objective_mat[1,] <- approximate_objective_function(Gram_matrix,c = gradient_run,lambda,n_terms)

max(gradient_run)

for(i in 2:10){
  old <- gradient_run
  gradient_run <- gradient_descent(Gram_matrix, old, lambda, eta, n_terms)
  objective_mat[i, ] <- approximate_objective_function(Gram_matrix,c = gradient_run,lambda,n_terms)
  
}

max(gradient_run)
dim(gradient_run)

curve(Fourier_evaluate(2*pi/366*x, coeffs=objective_mat[1,]), col=rainbow(10)[1], xlim=c(0,366), ylim=c(-5e5,-0.5e5), main='Approximate loss function after each iteration',ylab=TeX(r'($f(\alpha_t)$)'))

for(i in 2:10){
  curve(Fourier_evaluate(2*pi/366*x, coeffs=objective_mat[i,]), col=rainbow(10)[i], xlim=c(0,366), add=T)
}

# plot principal components
curve(Fourier_evaluate(2*pi/366*x, coeffs=gradient_run[1,,]), xlim=c(0,366), col='darkgrey', xlab='t', ylab='', ylim=c(0,5), main='c_0=1+cos(4*pi*x/366)')
title(ylab=TeX(r'($\langle p_1,\Phi(x_i)\rangle(t)$)'), mgp=c(2.4,1,0),cex.lab=1)

for(i in 2:47){
  curve(Fourier_evaluate(2*pi/366*x, coeffs=gradient_run[i,,]), xlim=c(0,366), add=T, col='darkgrey')
  
}

curve(Fourier_evaluate(2*pi/366*x, coeffs=gradient_run[1,,]), xlim=c(0,366), col='red', add=T)
curve(Fourier_evaluate(2*pi/366*x, coeffs=gradient_run[47,,]), xlim=c(0,366), col='blue', add=T)


curve(Fourier_evaluate(2*pi*x/366, coeffs=objective_mat[1,]),ylab='',main='Approximate loss function for each iteration', xlim=c(0,366), ylim=c(-500000,-1000), col=rainbow(7)[1])
for(i in 2:7){
  curve(Fourier_evaluate(2*pi*x/366, coeffs=objective_mat[i,]), add=T, col=rainbow(7)[i])
}

