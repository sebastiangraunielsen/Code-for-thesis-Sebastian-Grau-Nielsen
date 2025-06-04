library(tidyr)
library(dplyr)
library(kerntools)
library(cubature)
library(stats)
library(diceR)
library(polynom)
library(latex2exp)

alpha <- 4

# The functions used for data generation
y1_vec <- function(x){
  exp(alpha*(x[1]-x[2]))
}

y2_vec <- function(x){
  10*x[1]*x[2]
}

y3_vec <- function(x){
  cos(10*(x[1]-x[2]))
}

# Data generation
noise <- rnorm(60, mean=0, sd=0.3)

# The 9 integrals needed 
int1 <- adaptIntegrate(y1_vec,c(0,0),c(1,1))$integral
int2 <- adaptIntegrate(y2_vec,c(0,0),c(1,1))$integral
int3 <- adaptIntegrate(y3_vec,c(0,0),c(1,1))$integral
int4 <- adaptIntegrate(function(x) y1_vec(x)*y1_vec(x), c(0,0), c(1,1))$integral
int5 <- adaptIntegrate(function(x) y1_vec(x)*y2_vec(x), c(0,0), c(1,1))$integral
int6 <- adaptIntegrate(function(x) y1_vec(x)*y3_vec(x), c(0,0), c(1,1))$integral
int7 <- adaptIntegrate(function(x) y2_vec(x)*y2_vec(x), c(0,0), c(1,1))$integral
int8 <- adaptIntegrate(function(x) y2_vec(x)*y3_vec(x), c(0,0), c(1,1))$integral
int9 <- adaptIntegrate(function(x) y3_vec(x)*y3_vec(x), c(0,0), c(1,1))$integral
intvec <- c(int1,int2,int3)
intmat <- rbind(c(int4, int5, int6), c(int5, int7,int8), c(int6, int8, int9))

which_function <- function(i,j){ 
  # Matches index to underlying function
  if(i<= 20){
    first <- 1
  }
  if(i >= 21 && i<= 40){
    first <- 2
  }
  if(i >= 41){
    first <- 3
  }
  
  if(j<= 20){
    second <- 1
  }
  if(j >= 21 && j<= 40){
    second <- 2
  }
  if(j >= 41){
    second <- 3
  }
  return(c(first, second))
}

# Polynomial multiplication using fast Fourier transform 
poly_mult <- function(a, b) {
  if (!is.numeric(a) || !is.numeric(b)) stop("Both inputs must be numeric vectors.")

  na <- length(a)
  nb <- length(b)
  n_out <- na + nb - 1  # Resulting polynomial degree + 1
  
  result <- Re(fft(fft(c(a, rep(0, nb - 1))) * fft(c(b, rep(0, na - 1))), inverse = TRUE) / n_out)
  
  return(result[1:n_out])  
}

# Multiplication of matrix of polynomials with vector of polynomials
poly_matrix_mult <- function(poly_matrix, poly_vector) {
  # Get dimensions
  N <- nrow(poly_matrix)
  M <- ncol(poly_matrix)
  D1 <- dim(poly_matrix)[3] - 1  # Max degree of polynomials in matrix
  D2 <- ncol(poly_vector) - 1    # Max degree of polynomials in vector
  
  max_degree <- D1 + D2  # Max degree after multiplication
  result <- matrix(0, nrow = N, ncol = max_degree + 1)  # Store final coefficients
  
  for (i in 1:N) {
    sum_poly <- rep(0, max_degree + 1)  # Initialize zero polynomial
    for (j in 1:M) {
      # Multiply corresponding polynomial entries and accumulate
      prod_poly <- poly_mult(poly_matrix[i, j, ], poly_vector[j, ])
      sum_poly[1:length(prod_poly)] <- sum_poly[1:length(prod_poly)] + prod_poly
    }
    result[i, ] <- sum_poly  # Store the result for row i
  }
  
  return(result)  # Returns N x (D1 + D2 + 1) matrix
}

# Computing the outer product of two polynomial vectors
poly_outer_product <- function(poly_vec1, poly_vec2) {
  # Get dimensions
  N <- nrow(poly_vec1)  # Number of polynomials in first vector
  M <- nrow(poly_vec2)  # Number of polynomials in second vector
  D1 <- ncol(poly_vec1) - 1  # Max degree in poly_vec1
  D2 <- ncol(poly_vec2) - 1  # Max degree in poly_vec2
  
  max_degree <- D1 + D2  # Maximum polynomial degree after multiplication
  result <- array(0, dim = c(N, M, max_degree + 1))  # 3D array for outer product
  
  for (i in 1:N) {
    for (j in 1:M) {
      # Multiply polynomials element-wise
      result[i, j, 1:(D1 + D2 + 1)] <- poly_mult(poly_vec1[i, ], poly_vec2[j, ])
    }
  }
  
  return(result)  # Returns N x M x (D1 + D2 + 1) array
}

# Computing the inner product of two polynomial vectors
poly_inner_product <- function(poly_vec1, poly_vec2) {
  # Ensure the vectors have the same length
  if (nrow(poly_vec1) != nrow(poly_vec2)) stop("Vectors must have the same length.")
  
  N <- nrow(poly_vec1)  # Number of polynomials
  D1 <- ncol(poly_vec1) - 1  # Max degree in poly_vec1
  D2 <- ncol(poly_vec2) - 1  # Max degree in poly_vec2
  
  max_degree <- D1 + D2  # Maximum resulting polynomial degree
  result <- rep(0, max_degree + 1)  # Initialize result polynomial (zero coefficients)
  
  for (i in 1:N) {
    # Multiply corresponding polynomials and add to result
    prod_poly <- poly_mult(poly_vec1[i, ], poly_vec2[i, ])
    result[1:length(prod_poly)] <- result[1:length(prod_poly)] + prod_poly
  }
  
  return(rbind(result))  # Returns a single polynomial (as a coefficient vector)
}


# Function for plotting polynomials given coefficient vector
poly_eval <- function(coeffs, x) {
  degree <- length(coeffs) - 1
  y <- sapply(x, function(xi) sum(coeffs * xi^(0:degree)))
  return(y)
}

# Compute the gradient given Gram matrix, current point and lambda
Gradient_coef <- function(G, alpha, lambda){
  G_alpha <- poly_matrix_mult(G, alpha)
  G2_alpha <- poly_matrix_mult(G,G_alpha)
  alpha_G_alpha <- poly_inner_product(G_alpha, alpha)
  G_alpha_polymat <- array(G_alpha, dim = c(nrow(G_alpha), 1, ncol(G_alpha)))

  Final_term <- poly_matrix_mult(G_alpha_polymat, alpha_G_alpha)
  
  numrows <- nrow(G_alpha)
  maxcol <- max(ncol(G_alpha),ncol(G2_alpha),ncol(Final_term))
  
  G_alpha <- cbind(G_alpha, matrix(0,nrow=numrows,
                                   ncol=maxcol-ncol(G_alpha)))
  G2_alpha <- cbind(G2_alpha, matrix(0,nrow=numrows,
                                     ncol=maxcol-ncol(G2_alpha)))
  Final_term <- cbind(Final_term, 
                                 matrix(0,nrow=numrows,ncol=maxcol-
                                          ncol(Final_term)))
  
  return(-2*G2_alpha-4*lambda*G_alpha+4*lambda*Final_term)
}

# Perform update using gradient descent
Grad_descent_step <- function(alpha_t, eta_t, G, lambda){
  Grad <- Gradient_coef(G, alpha_t,lambda)
  maxcol <- max(ncol(alpha_t),ncol(Grad))
  alpha_t <- cbind(alpha_t, matrix(0, nrow=nrow(alpha_t),ncol=maxcol - ncol(alpha_t)))
  Grad <- cbind(Grad, matrix(0, nrow=nrow(Grad),ncol=maxcol - ncol(Grad)))
  return(alpha_t - eta_t * Grad)
}

# The function to be minimized
objective_function <- function(alpha, G, lambda){
  G_alpha <- poly_matrix_mult(G, alpha)
  first_term <- -poly_inner_product(G_alpha, G_alpha)
  alpha_G_alpha <- poly_inner_product(G_alpha, alpha)
  second_term <- lambda*poly_mult(alpha_G_alpha,alpha_G_alpha)
  third_term <- -2*lambda*alpha_G_alpha
  third_term[1,1] <- third_term[1,1] + lambda
  
  max_col <- max(length(first_term), length(second_term), length(third_term))
  
  first_term <- cbind(rbind(first_term), matrix(0, nrow=1, ncol=max_col - length(first_term)))
  second_term <- cbind(rbind(second_term), matrix(0, nrow=1, ncol=max_col - length(second_term)))
  third_term <- cbind(rbind(third_term), matrix(0, nrow=1, ncol=max_col - length(third_term)))
  
  return(first_term + second_term + third_term)
  
}

Gramcoefficients <- array(1, dim=c(60,60,3)) # Array for second degree polynomials

# Compute the Gram matrix
for(i in 1:60){
  for(j in 1:60){
    for(k in 1:2){
      which_one <- which_function(i,j)
      if(k==1){
        Gramcoefficients[i,j,k] <- (noise[i]*noise[j] + noise[j]*intvec[which_one[1]] + 
                                                 noise[i]*intvec[which_one[2]] + intmat[which_one[1],which_one[2]])
      }
      if(k==2){
        Gramcoefficients[i,j,k] <- -1*(noise[i] + noise[j] + 
                                                            intvec[which_one[1]] + intvec[which_one[2]])
      }
    }
  }
}

Gramcoefficients[1,1,]

# Plot polynomial from Gram matrix
curve(poly_eval(Gramcoefficients[1,1,],x),0,1)

for(i in 1:60){
  for(j in 1:60){
    curve(poly_eval(Gramcoefficients[i,j,],x),0,1, add=T)
  }
}


# Run the gradient descent 
lambda <- 1e-7
eta <- 1e-8
n_iter <- 6

# First iteration 
grad_run <- Grad_descent_step(alpha_t = matrix(1,ncol=1,nrow=60),eta_t = eta,G = Gramcoefficients,lambda = lambda)

curve(poly_eval(objective_function(grad_run,Gramcoefficients,lambda = lambda), x), 0, 1,ylab=TeX(r'($f(c_t)(s)$)'),col=rainbow(9)[1], xlab='s')
# Subsequent iterations
for(i in 1:n_iter){
  old <- grad_run
  grad_run <- Grad_descent_step(old,eta,Gramcoefficients,lambda)
  new_objective <- objective_function(grad_run, Gramcoefficients, lambda)
  curve(poly_eval(new_objective, x), 0, 1, col=rainbow(9)[i], add=T)
}
legend(0.5, -2e8, legend=c('t=1', 't=2', 't=3', 't=4','t=5','t=6','t=7'),lty=1, col=rainbow(9)[1:7], cex=0.6)

# Plot the principal components
curve(poly_eval(x,coeffs=grad_run[1,]),,col='blue', xlab='t', ylab='', ylim=c(0.9,1.3))

for(i in 2:20){
  curve(poly_eval(x,coeffs=grad_run[i,]),col='blue', add=T)
}
for(i in 21:40){
  curve(poly_eval(x,coeffs=grad_run[i,]),col='red',add=T)
}
for(i in 41:60){
  curve(poly_eval(x,coeffs=grad_run[i,]),col='green',add=T)
}
legend(0.3, 1.3, legend=c('Sample set 1', 'Sample set 2', 'Sample set 3'),lty=1, col=c('blue', 'red', 'green'), cex=0.6)
title(ylab=TeX(r'($\langle p_1,\Phi(x_i)\rangle(t)$)'), mgp=c(2.4,1,0),cex.lab=1)

objective_function(grad_run, Gramcoefficients, lambda)

curve(poly_eval(objective_function(grad_run, Gramcoefficients, lambda), x),0,1)

# Regular kernel PCA using Laplacian kernel  
xvals <- seq(0, 1, length.out=11)
noisy_obs <- matrix(0, nrow=60, ncol=121)
pointpairs <- expand_grid(1:11,1:11)
indices <- cbind(pull(pointpairs[,1]),pull(pointpairs[,2]))

indices
for(i in 1:60){
  for(j in 1:121)
      if(i<= 20){
        noisy_obs[i,j] <- y1_vec(c(xvals[indices[j,1]],xvals[indices[j,2]])) + noise[i]
      }
      if(i<= 40 && i>= 21){
        noisy_obs[i,j] <- y2_vec(c(xvals[indices[j,1]],xvals[indices[j,2]])) + noise[i]
      }
      
      if(i>=41){
        noisy_obs[i,j] <- y3_vec(c(xvals[indices[j,1]],xvals[indices[j,2]])) + noise[i]
      }
}

# Define the Laplacian kernel function
laplacian_kernel <- function(x, gamma = 0.1) {
  exp(-as.matrix(dist(x, method = "euclidean")) * gamma)
}

# Find first 2 principal components and plot them. Parameter "center"
# controls centering of Gram matrix
kPCA(laplacian_kernel(noisy_obs, gamma=1/10), y=as.factor(c(rep('sample-set 1',20),rep('sample-set 2',20),rep('sample-set 3',20))),
     plot=1:2, colors=c('green', 'red', 'blue'), center = T)

