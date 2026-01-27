# Solving OLS Models Using QR Decomposition
The aim of this project is solving Ordinary Least Squares (OLS) models using QR Decomposition. For that purpose, `R` is used to code three different functions, which perform the aforementioned decomposition, solve linear systems for upper triangular matrices, and use the previous results to identify the coefficients of interest, respectively.

## Overview Of The Problem
A standard result in Ordinary Least Squares (OLS) models stemming from the minimization of the Sum of Squared Residuals (SSR) is 

$$(X^TX)\beta=X^Ty,$$

where $X$ is the $n\times p$ matrix of the observed values for the indepent variables, $\beta$ is the $p$-length vector containing the unknown parameters of the OLS model, and $y$ corresponds to the $n$-length vector which includes all the observations associated with the dependent variable. We can then estimate our coefficients of interest straightforwardly by inverting $(X^TX)$, namely

$$\hat{\beta}=(X^TX)^{-1}X^Ty.$$

However, sometimes the inverting process involves difficulties, which would make solving the uninverted system of equations resulting from solving the SSR problem our preferred choice. We will carry this out by applying the QR Decomposition, a procedure that, for some given $n\times p$ matrix, finds a $n\times p$ orthogonal matrix $Q$ and a $p\times p$ upper triangular $R$ whose product of matrices yields the same output. As a result of applying this on both sides of our targeted system, we obtain

$$R\beta=Q^Ty.$$

We now can work with a simpler system, as $R$ is an upper triangular matrix.

From this overview, we can conlude the `R` code we will be using to solve our model will need two major parts: a function which performs the QR Decomposition and another one that solves linear systems with upper triangular matrices. In the following sections we will present the algorithm we will implement for each case.

### QR Decomposition

The numerical method used is the Gram-Schmidt approach. For a $n\times p$ matrix A 
```math
\begin{pmatrix} 
a_1 & a_2 & \cdots & a_p
\end{pmatrix},
```
where $a_1,\ a_2,\ \cdots\ ,\ a_p$ correspond to the $n$-length vectors associated with each of the $p$ columns of the matrix, this algorithms performs the following computations:

$$u_1=a_1,\ e1=\frac{u_1}{||u_1||},$$

$$u_2=a_2-(a_2\cdot e_1)e_1,\ e_2=\frac{u_2}{||u_2||},$$

$$u_{k+1}=a_{k+1}-(a_{k+1}\cdot e_1)e_1-\cdots -(a_{k+1}\cdot e_k)e_k,\ e_{k+1}=\frac{u_{k+1}}{||u_{k+1}||},$$

where $||\cdot ||$ is the $L_2$ norm.

The resulting QR factorization is

$$
\begin{pmatrix} 
a_1 & a_2 & \cdots & a_p
\end{pmatrix}=
\begin{pmatrix} 
e_1 & e_2 & \cdots & e_p
\end{pmatrix}
\begin{pmatrix}
a_1\cdot e_1 & a_2\cdot e_1 & \cdots & a_n\cdot e_1 \\
0       & a_2\cdot e_2 & \cdots & a_n\cdot e_2 \\
\vdots  & \vdots  & \ddots & \vdots  \\
0       & 0       & \cdots & a_n\cdot e_n
\end{pmatrix}
$$

### Solving An Upper-triangular System

We will use back substitution in order to solve this type of linear systems. For an equality $Rx=b$, where $R$ is a $p\times p$ upper-triangular matrix, $x$ is a $p$-length vector of unknowns, and $b$ is a $p$-length vector of values, the algorithm works as follows:
- If $n=1$, then $x=b/A$
- If $n>1$, then:
    - A
