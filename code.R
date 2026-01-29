BackSolving <- function(B, b){ # Inputs: pxp matrix B and p-length vector b
  p <- length(b) # Number of unknown parameters
  if(p==1){
    x <- b/B # If there is just one equation, the answer is trivial
    x
  }
  else{
    x <- rep(NA, p) # Create a p-length vector x
    i <- p # Start from the last element of x
    while(length(b)>=1){ # The loop will hold until b has been 
                         # reduced to a scalar
      x[i] <- b[i]/as.matrix(B)[i,i] # as.matrix is needed since, 
                                     # when B is reduced to a 1x1 matrix,
                                     # R automatically transforms it into
                                     # a scalar, in which case indicating
                                     # two dimensions does not work
      if(length(b)>1){ # b and B will keep being reduced until
                       # they both include just one value
        b <- b[1:(i-1)]-B[1:(i-1),i]*x[i]
        B <- B[1:(i-1),1:(i-1)]
        i <- i-1 # Move into the previous element of x 
      }
      else{break} # Breaking the loop is needed at the end of the process
    }
    x
  }
}

GramSchmidt <- function(A){ # Input: nxp matrix A
  U <- A # This new matrix will be transformed along the function
         # and used to store the vectors u presented in the README 
         # file and used throughout the algorithm
  n <- dim(A)[1] # Number of rows of A
  p <- dim(A)[2] # Number of columns of A
  Q <- matrix(NA, n, p) # Create NA matrix Q
  R <- matrix(0, p, p) # Create zero matrix R
  for(k in 1:p){ # Loop covers every column of Q and R sequentially
    Q[,k] <- U[,k]/sqrt(sum((U[,k])^2)) # n-length vector for column k in Q
    R[k,k] <- sum(A[,k]*Q[,k]) # kth element of R's leading diagonal
    if(k<p){ # If the number of the column is different from p, the right-hand
             # side elements of the diagonal must be filled for R
      for(i in (k+1):p){ # Loop over every element at the right-hand
                         # of the diagonal
        R[k,i] <- sum(A[,i]*Q[,k])
      }
      U[,k+1] <- A[,k+1]-(sum(A[,k+1]*Q[,1]))*Q[,1] 
      if((k+1)>2){
        for(h in 2:k){
          U[,k+1] <- U[,k+1]-(sum(A[,k+1]*Q[,h]))*Q[,h]
        }
      } # Lines 46-51 build the u vector needed in the next iteration
    }
  }
  list(A = A, Q = Q, R = R, QR = Q%*%R) # This list displays the matrix A, Q, R,
                                        # and the product of the last two
}

OLSByQRDec <- function(X,y){ # Inputs: nxp matrix X and n-length vector y
  R <- GramSchmidt(X)$R # Invoke R matrix for X using 'GramSchmidt' function
  Qtrans <- t(GramSchmidt(X)$Q) # Invoke Q matrix for X and transpose it
  Qtransy <- Qtrans%*%y # Matrix-multiply transposed-Q and y
  BackSolving(R,Qtransy) # Use back substitution to get the coefficient values
}


## Example ##

#install.packages("tidyverse") # Install the 'tidyverse' package collection
                               # in case you have not done so before

library(tidyverse) # Load 'tidyverse' 

#set.seed(45) # This randomization seed leads to the same
              # output than the one observed in the README file

data.frame(x_1 = rnorm(200,3,3), 
           x_2 = rnorm(200,3,3)) %>% 
  mutate(y = 5 + 2*x_1 + 0.5*x_2 + rnorm(length(x_1),0,1)) -> data 
# Create a simulated dataset in which y is a linear combination of the 
# randomly-selected independent variables plus a ~N(0,1) error term

X <- matrix(c(rep(1,dim(data)[1]), data$x_1, data$x_2),
            dim(data)[1],
            dim(data)[2]) # Create a matrix in which each column vector 
                          # corresponds with each of the explanatory variables
                          # plus a vector of 1s (notice y has an intercept)

y <- data$y # Vector containing the values for the dependent variable

OLSByQRDec(X,y) # After applying 'OLSByQRDec' on both matrices,
                # observe the result is almost identical to the
                # data generating process
