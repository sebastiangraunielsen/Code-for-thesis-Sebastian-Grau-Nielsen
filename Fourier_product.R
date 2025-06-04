library(CVXR)

pad_to_equal_length <- function(long, short) {
  short <- c(short, rep(0, length(long) - length(short)))
  return(short)
}


Fourier_product <- function(coeffs1, coeffs2){
  # Computes Fourier coefficients of a product of Fourier series
  
  a1_vec <- c(coeffs1[1], coeffs1[-1][c(T,F)])
  a2_vec <- c(coeffs2[1], coeffs2[-1][c(T,F)])
  b1_vec <- pad_to_equal_length(a1_vec, coeffs1[-1][c(F,T)])
  b2_vec <- pad_to_equal_length(a1_vec,coeffs2[-1][c(F,T)])
  N <- length(a1_vec) - 1
  c1_vec <- numeric(2*N + 1)
  c2_vec <- numeric(2*N + 1)
  c1_vec[N + 1] <- a1_vec[1]
  c2_vec[N + 1] <- a2_vec[1]
  
  for(n in (-N):N){
    if(n > 0){
      c1_vec[n + N + 1] <- 1/2 * (a1_vec[n+1] - 1i*b1_vec[n])
      c2_vec[n + N + 1] <- 1/2 * (a2_vec[n+1] - 1i*b2_vec[n])
    }
    if(n < 0){
      c1_vec[n + N + 1] <- 1/2 * (a1_vec[-n+1] + 1i*b1_vec[-n])
      c2_vec[n + N + 1] <- 1/2 * (a2_vec[-n+1] + 1i*b2_vec[-n])
    }
  }
  
  new_coeff <- convolve(c1_vec, rev(c2_vec), type='open', conj=T)
  a_out <- numeric(2*N+1)
  a_out[1] <- new_coeff[2*N+1]
  
  b_out <- numeric(2*N)
  nice_order_1 <- rev(new_coeff[1:(2*N)])
  nice_order_2 <- new_coeff[(2*N+2):(4*N+1)]
  
  for(n in 1:(2*N)){
    a_out[n+1] <- nice_order_2[n] + nice_order_1[n]
    b_out[n] <- 1i*(nice_order_1[n]-nice_order_2[n])
  }
  
  coeff_out <- numeric(4*N+1)
  coeff_out[1] <- a_out[1] 
  
  for(i in 2:(4*N+1)){
    if(i%%2==0){
      coeff_out[i] <- a_out[i/2 + 1]
    }
    else{
      coeff_out[i] <- b_out[(i-1)/2]
    }
  }
  return(round(Re(coeff_out), 5))
}
