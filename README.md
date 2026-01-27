# Solving OLS Models Using QR Decomposition
The aim of this project is solving Ordinary Least Squares (OLS) models using QR Decomposition. For that purpose, `R` is used to code three different functions, which perform the aforementioned decomposition, solve linear systems for upper triangular matrices, and use the previous results to identify the coefficients of interest, respectively.

## Overview of the Problem
A standard result in Ordinary Least Squares (OLS) models stemming from the minimization of the Sum of Squared Residuals (SSR) is 

$$(X^TX)\beta=X^Ty,$$

where $X$ is the $n\times p$ matrix of the observed values for the indepent variables, $\beta$ is the $p$-length vector containing the unknown parameters of the OLS model, and $y$ corresponds to the $n$-length vector which includes all the observations associated with the dependent variable. We can then estimate our coefficients of interest straightforwardly by inverting $(X^TX)$, namely

$$\hat{\beta}=(X^TX)^{-1}X^Ty.$$

However, sometimes the inverting process involves difficulties, which would make solving the uninverted system of equations resulting from solving the SSR problem our preferred choice. We will carry this out by applying the QR Decomposition, a procedure that, for some given $n\times p$ matrix, finds a $n\times p$ orthogonal matrix $Q$ and a $p\times p$ upper triangular $R$ whose product of matrices yields the same output. As a result of applying this on our targeted system, we obtain

$$R\beta=Q^Ty.$$

We now can work with a simpler system, as $R$ is an upper triangular matrix.

From this overview, we can conlude the `R` code we will be using to solve our model will need two major parts: a function which performs the QR Decomposition and another one that solves linear systems with upper triangular matrices. In the following sections we will present the algorithm we will implement for each case.

## QR Decomposition

The numerical method used is the Gram-Schmidt approach. The algorithm proceeds as follows: for a $n\times p$ matrix A 
```math
\begin{pmatrix} 
a_1 & a_2 & ... & a_p
\end{pmatrix},
```
where $a_1,\ a_2,\ ...\ ,\ a_p$ correspond to the $n$-length vectors associated with each to the $p$ columns of the matrix
