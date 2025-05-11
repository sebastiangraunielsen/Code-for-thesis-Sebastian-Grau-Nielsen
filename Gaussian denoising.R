library(keras)
library(kerntools)

par(mfrow=c(1,2))

fixed_point_iter <- function(z, deltas, gamma, images){
  # Performs fixed point iteration from section 2.2
  help_vec <- numeric(length(deltas))
  help_arr <- array(0, dim=dim(images))
  
  for(i in 1:length(deltas)){
    help_vec[i] <- deltas[i] * exp(-gamma * sum((z-images[i,,])^2))
    help_arr[i,,] <- help_vec[i] * images[i,,]
  }
  return(apply(help_arr, c(2,3), sum)/sum(help_vec))
}

beta_calc <- function(x, alphas, gamma, images){
  # Calculates the betas for projection onto principal axes
  kern_vec <- numeric(dim(images)[1])
  for(i in 1:dim(images)[1]){
    kern_vec[i] <- exp(-gamma * sum((x - images[i,,])^2))
  }
  
  return(t(alphas)%*%as.matrix(kern_vec))
}

mnist <- dataset_mnist()
n_train <- 200 # Number of training images for choosing gamma
n_test <- 50 # Number of test images for choosing gamma
sd_noise <- 0.5 # Standard deviation of Gaussian noise to be applied to each image
n_pc <- 20 # Number of principal components for choosing gamma
n_iter <- 100 # Number of iterations for preimage problem

distances <- matrix(0, n_train, n_train)

x_train <- mnist$train$x[1:n_train,,]/255
x_train_noisy <- x_train + array(rnorm(n_train*28^2, 0, sd_noise), dim=c(n_train, 28, 28))
x_train_noisy <- pmax(pmin(x_train_noisy, 1), 0)

x_test <- mnist$test$x[1:n_test,,]/255
x_test_noisy <- x_test + array(rnorm(n_test*28^2, 0, sd_noise), dim=c(n_test, 28, 28))
x_test_noisy <- pmax(pmin(x_test_noisy, 1), 0)

for(i in 1:n_train){ # Compute distances for median heuristic
  for(j in 1:n_train){
    distances[i,j] <- sqrt(sum((x_train_noisy[i,,] - x_train_noisy[j,,])^2))
  }
}

gamma <- 1/(2*median(distances)^2) # median heuristic if desired

choose_gamma <- function(gamma_seq, test_pics, train_dat, n_pc, sd_noise, n_iter){
  # Function estimating the mean squared prediction error for different values of gamma  
  mean_sq <- numeric(length(gamma_seq))
  test_dims <- dim(test_pics)
  n_train <- dim(train_dat)[1]
  
  test_pics_noisy <- test_pics + array(rnorm(test_dims[1]*test_dims[2]*test_dims[3], 0, sd_noise), dim = test_dims)
  
  for(k in 1:length(gamma_seq)){
  Gram_mat <- matrix(0, n_train, n_train)
  for(i in 1:n_train){
    for(j in 1:n_train){
      Gram_mat[i,j] <- exp(-sum((train_dat[i,,] - train_dat[j,,])^2)*gamma_seq[k])
    }
  }
  
  pca_run <- kPCA(Gram_mat)
  pc_matrix <- pca_run[,1:n_pc]
  denoised_images <- test_pics_noisy # We use the noisy images as a starting value
  
  for(i in 1:n_test){
    deltas <- as.matrix(pc_matrix)%*%beta_calc(test_pics_noisy[i,,], pc_matrix, gamma_seq[k], train_dat)
    for(j in 1:n_iter){ # fixed point iteration from section 2.2 where we clip to [0,1] to always get well-defined images
      denoised_images[i,,] <-  pmax(pmin(fixed_point_iter(denoised_images[i,,], deltas, gamma_seq[k], train_dat), 1), 0)
    }
  }
  
  mean_sq[k] <- sum((denoised_images - test_pics)^2)/test_dims[1]
  print(k)
  }
  return(mean_sq)
}

gamma_choose <- choose_gamma(seq(1e-4, 1, 0.01), x_test, x_train_noisy, n_pc, sd_noise, n_iter)

# Plot mean squared prediction error as a function of gamma
par(mfrow=c(1,1))
plot(seq(1e-4, 1, 0.01)[1:100],gamma_choose[1:100],xlab='gamma',ylab='Mean squared error', type='l')
gamma <- seq(1e-4, 1, 0.01)[which.min(gamma_choose)]

# Now, we pick the best gamma and train on a larger dataset
n_train <- 1000
n_test <- 100
sd_noise <- 0.5
n_pc <- 20
n_iter <- 200

x_train <- mnist$train$x[1:n_train,,]/255
x_train_noisy <- x_train + array(rnorm(n_train*28^2, 0, sd_noise), dim=c(n_train, 28, 28))
x_train_noisy <- pmax(pmin(x_train_noisy, 1), 0)

x_test <- mnist$test$x[1:n_test,,]/255
x_test_noisy <- x_test + array(rnorm(n_test*28^2, 0, sd_noise), dim=c(n_test, 28, 28))
x_test_noisy <- pmax(pmin(x_test_noisy, 1), 0)

Gram_mat <- matrix(0, n_train, n_train)

for(i in 1:n_train){
  for(j in 1:n_train){
    Gram_mat[i,j] <- exp(-sum((x_train_noisy[i,,] - x_train_noisy[j,,])^2)*gamma)
  }
}

pca_run <- kPCA(Gram_mat)

pc_matrix <- pca_run[,1:n_pc]

denoised_images <- x_test_noisy[51:70,,]

for(i in 1:20){
  deltas <- as.matrix(pc_matrix)%*%beta_calc(x_test_noisy[i+50,,], pc_matrix, gamma, x_train_noisy)
  for(j in 1:n_iter){
    denoised_images[i,,] <-  pmax(pmin(fixed_point_iter(denoised_images[i,,], deltas, gamma, x_train_noisy), 1), 0)
  }
}

par(mfrow=c(2,3))
plot(as.raster(x_test_noisy[1+50,,]))
plot(as.raster(x_test_noisy[2+50,,]))
plot(as.raster(x_test_noisy[7+50,,]))
plot(as.raster(denoised_images[1,,]))
plot(as.raster(denoised_images[2,,]))
plot(as.raster(denoised_images[7,,]))
denoised_images


