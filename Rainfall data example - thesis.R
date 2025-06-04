library(splines2)
library(ICSsmoothing)
library(purrr)
library(caret)
library(mgcv)
library(kerntools)
library(pracma)
library(latex2exp)

par(mfrow=c(1,1))
source("Documents/Poly_helpers.R") # Load functions for polynomial calculations

# Gaussian kernel 
gaussian_kernel <- function(x, gamma = 0.1) {
  exp(-as.matrix(dist(x, method = "euclidean")^2) * gamma)
}

# Read the data 
Rainfall <- read.csv("Documents/daily.txt")

# Correct spelling error
Rainfall$SEASON[Rainfall$SEASON=='SPRNG'] <- 'SPRING'
Rainfall$Data <- as.Date(Rainfall$Data, "%d/%m/%Y")
Rainfall$n_days <- 1:length(Rainfall$Data)

# Boxplots of the different seasons
boxplot(Rainfall$PREC~Rainfall$SEASON, ylab='Precipitation', xlab='Day')


plot(Rainfall$Data,Rainfall$PREC, type='l', col='grey', ylab='Precipitation', xlab='Day', main='Brazilian rainfall')

sum(Rainfall$PREC==0)/2191 # Proportion of 0s 

# Function that evaluates a polynomial at x given a coefficient vector
eval_poly <- function(coeffs, x) {
  sum(coeffs * x^(0:(length(coeffs)-1)))
}

# Functions that evaluate a spline given a matrix of polynomials and a vector of knots
spline_eval_helper <- function(x, poly_matrix, knots){
    # Find the correct interval (assumes knots are sorted)
    n_splines <- nrow(poly_matrix)
    k_index <- max(which(x >= knots)) 
  
    k_index <- min(k_index, n_splines)
    
    # Evaluate the corresponding polynomial
    return(eval_poly(poly_matrix[k_index, ], x))
}

evaluate_splines <- function(poly_matrix, knots, x_values) {
  spline_values <- sapply(x_values, function(x) spline_eval_helper(x, poly_matrix, knots))
  
  return(spline_values)
}

# Cross-validation for number of uniform knots based on MSE
cv_spline_knots <- function(x, y, max_knots, n_folds){
  
  folds <- createFolds(y[-c(1, length(y))], k = n_folds, list = TRUE, returnTrain = FALSE) 
  cv_errors <- numeric(max_knots)

  for(i in 2:(max_knots+1)){
    spline_eval <- 0
    for(j in 1:n_folds){
      currentfold <- folds[[j]] + 1
      spline_train <- cics_unif_explicit_smooth(x[-currentfold], y[-currentfold], i, plotTF = F)
      poly_matrix <- spline_train$est_spline_coeffs
      knots <- spline_train$nodes
      
      spline_eval <- spline_eval + mean((y[currentfold] - evaluate_splines(poly_matrix,knots,x[currentfold]))^2)
    }
    
    cv_errors[i - 1] <- spline_eval/n_folds
  }
  return(cics_unif_explicit_smooth(x,y,k=which.min(cv_errors)+1))
}

# Run 100-fold CV with a maximum of 70 uniformly placed knots
total_fit <- cv_spline_knots(Rainfall$n_days,y = Rainfall$PREC,70,100)

Precipitation <- Rainfall$PREC
Days <- Rainfall$n_days

# Get the polynomial coefficients and knots
total_fit_poly <- total_fit$est_spline_coeffs
total_fit_knots <- total_fit$nodes

# Determine last day of each season
change <- numeric(2191)

for(i in 2:2191){
  if(Rainfall$SEASON[i] != Rainfall$SEASON[i-1]){
    change[i-1] <- 1
  }
}

last_days <- which(change==1)
last_days

# Zoom in on plots 
plot(1:1000, evaluate_splines(total_fit_poly, total_fit_knots, 1:1000), type='l', xlab='Day', ylab='Precipitation')
for(i in 1:length(last_days)){
  abline(v=last_days[i], lty=3)
}
legend(450, 13, legend = c('Fitted spline'), lty=1)

n_seasons <- length(last_days) + 1

# Determine length of the seasons
season_length <- numeric(n_seasons)
season_length[1] <- last_days[1]

for(i in 2:length(season_length)){
  season_length[i] <- c(last_days,2191)[i] - c(last_days,2191)[i-1]
}

season_length

# Parameters for gradient descent
eta <- 1e-7
lambda <- 1e-6
n_iter <- 9

Gram_mat_standard <- array(1, dim=c(n_seasons, n_seasons, 3))

# Fill out Gram matrix
for(i in 1:n_seasons){
  for(j in 1:n_seasons){
    poly1_knot <- c(1,last_days)[i]
    poly2_knot <- c(1,last_days)[j]
    
    Gram_mat_standard[i,j,1] <-  integrate(function(x) evaluate_splines(total_fit_poly,total_fit_knots,x*season_length[i]+poly1_knot)*evaluate_splines(total_fit_poly,total_fit_knots,x*season_length[j]+poly2_knot), 0,1)$value
    Gram_mat_standard[i,j,2] <- -(integrate(function(x) evaluate_splines(total_fit_poly,total_fit_knots,season_length[i]*x+poly1_knot),0,1)$value +
                                                 integrate(function(x) evaluate_splines(total_fit_poly,total_fit_knots,season_length[j]*x+poly2_knot),0,1)$value)
  }
}

# Plot functions in Gram matrix (with standardization)
curve(poly_eval(Gram_mat_standard[1,1,],x),0,1, ylim=c(0,100),ylab='', xlab='t', main='Gram matrix functions')
title(ylab=TeX(r'($K(x_i,x_j)(t)$)'), mgp=c(2.4,1,0), cex.lab=1)
for(i in 1:n_seasons){
  for(j in 1:n_seasons){
    curve(poly_eval(Gram_mat_standard[i,j,],x),0,1, add=T)
  }
}

# First iteration 
grad_run_standard <- Grad_descent_step(matrix(1, nrow=n_seasons,ncol = 1),eta, Gram_mat_standard,lambda)
curve(poly_eval(objective_function(matrix(1, nrow=n_seasons,ncol = 1),Gram_mat_standard,lambda = lambda),x), col=rainbow(10)[1], ylab='', xlab='s', ylim=c(-5e7, -6.0e6))
title(ylab=TeX(r'($f(\alpha_t)(s)$)'), mgp=c(2.4,1,0), cex.lab=1)

# Subsequent iterations
for(i in 1:n_iter){
  old <- grad_run_standard
  grad_run_standard <- Grad_descent_step(old, eta, Gram_mat_standard, lambda)
  curve(poly_eval(objective_function(grad_run_standard,Gram_mat_standard,lambda = lambda),x), col=rainbow(10)[i], add=T)
}

col_vec <- c('red', 'green', 'red', 'blue')

# Plot principal components for the different seasons
curve(poly_eval(grad_run_standard[1,],x),0,1, ylim=c(1,4), col='red', ylab='', xlab='s')
title(ylab=TeX(r'($\langle p_1,\Phi(x_i)\rangle(s)$)'), mgp=c(2.4,1,0), cex.lab=1)
for(i in 1:(n_seasons)){
  curve(poly_eval(grad_run_standard[i,],x), col=col_vec[(i %% 4)+1], add=T)
}
legend(0.35,4.05, legend=c('Summer', 'Spring and autumn', 'Winter'),lty=1, col=c('blue', 'red', 'green'))
