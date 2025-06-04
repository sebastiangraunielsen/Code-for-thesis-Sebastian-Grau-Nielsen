library(rcc)

poly_mult <- function(a, b) {
  # Ensure inputs are numeric vectors
  if (!is.numeric(a) || !is.numeric(b)) stop("Both inputs must be numeric vectors.")
  
  # Compute polynomial multiplication using convolution
  na <- length(a)
  nb <- length(b)
  n_out <- na + nb - 1  # Resulting polynomial degree + 1
  
  # Use FFT for efficient polynomial multiplication
  result <- Re(fft(fft(c(a, rep(0, nb - 1))) * fft(c(b, rep(0, na - 1))), inverse = TRUE) / n_out)
  
  return(result[1:n_out])  # Extract valid coefficients
}

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


poly_eval <- function(coeffs, x) {
  degree <- length(coeffs) - 1
  y <- sapply(x, function(xi) sum(coeffs * xi^(0:degree)))
  return(y)
}

Gradient_coef <- function(G, alpha, lambda){
  G_alpha <- poly_matrix_mult(G, alpha)
  G2_alpha <- poly_matrix_mult(G,G_alpha)
  G_alpha_G_alpha <- poly_outer_product(G_alpha, G_alpha)
  
  Final_term <- poly_matrix_mult(G_alpha_G_alpha, alpha)
  
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

# Old gradient -2*G2_alpha-4*lambda*G_alpha+4*lambda*Final_term

normalize <- function(x) x / max(abs(x))

Grad_descent_step <- function(alpha_t, eta_t, G, lambda){
  Grad <- Gradient_coef(G, alpha_t,lambda)
  maxcol <- max(ncol(alpha_t),ncol(Grad))
  alpha_t <- cbind(alpha_t, matrix(0, nrow=nrow(alpha_t),ncol=maxcol - ncol(alpha_t)))
  Grad <- cbind(Grad, matrix(0, nrow=nrow(Grad),ncol=maxcol - ncol(Grad)))
  return(alpha_t - eta_t * Grad)
}

spline_mult <- function(polymat1, polymat2, common_knots){
  newcol <- ncol(polymat1) + ncol(polymat2)
  
}

objective_function <- function(alpha, G, lambda){
  G_alpha <- poly_matrix_mult(G, alpha)
  first_term <- -poly_inner_product(G_alpha, G_alpha)
  alpha_G_alpha <- poly_inner_product(G_alpha, alpha)
  second_term <- lambda*poly_mult(alpha_G_alpha,alpha_G_alpha)
  third_term <- -2*lambda*alpha_G_alpha
  third_term[1,1] <- third_term[1,1] + lambda
  
  max_col <- max(length(first_term), length(second_term), length(third_term))
  
  print(first_term)
  print(second_term)
  print(third_term)
  
  first_term <- cbind(rbind(first_term), matrix(0, nrow=1, ncol=max_col - length(first_term)))
  second_term <- cbind(rbind(second_term), matrix(0, nrow=1, ncol=max_col - length(second_term)))
  third_term <- cbind(rbind(third_term), matrix(0, nrow=1, ncol=max_col - length(third_term)))
  
  return(first_term + second_term + third_term)
  
}

