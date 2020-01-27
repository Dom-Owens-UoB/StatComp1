# myOptPackage functions

#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
library(matrixcalc) #load matrixcalc for pos.def. test is.positive.definite
library(MASS) #load MASS for mvrnorm


###############################
#
#' Define function returning Rosenbrock function (a-x)^2 + b*(y-x^2)^2 parameterised by a,b along with gradient and Hessian
#'
#' @param a 
#' @param b 
#'
#' @return function  (a-x)^2 + b*(y-x^2)^2 of x and y
#' @export
#'
#' @examples
#' Rosen(1,100)
Rosen <- function(a = 1,b = 100){
  f <- deriv(expression((a-x)^2 + b*(y-x^2)^2), namevec = c('x', 'y'), function.arg=T, hessian=T) #define Rosenbrock function
  f
}



###############################
#' Gradient Descent for minimisation
#'
#' @param f1
#' @param x0
#' @param s
#' @param alpha
#' @param beta
#' @param eps
#' @param max_iter  
#'
#' @return results object with Iterations, Coordinates, Gradient, Function value at minimum. Trace coordinates and trace function values as attributes.
#' @export
#'
#' @examples
GradDesc <- function(f1,x0, s=1, alpha = 1/4, beta = 1/2, eps = 1e-2, max_iter = 50000){
  x <- x0 #set initial coordinate
  iter <- 0 #set iterations to 0
  grad <-  attr(f1(x[1], x[2]), "gradient") #compute gradient at current coordinate
  fun_val <- rep(0, max_iter) #initiate vector of function values
  coords <- matrix(NA, nrow =  max_iter, ncol = 2) #initiate vector of coordinates
  
  
  while (norm(grad, type ="2") > eps && iter < max_iter){
    iter <- iter + 1 #increment iterations
    t <- s #set t to initial stepsize
    change <- x - t* grad #set change in coords
    
    while (f1(x[1],x[2]) - f1(change[1],change[2]) < alpha * t * norm(grad, type ="2")^2) {
      t <- beta * t  #set stepsize by backtracking 
      change <- x - t* grad #set change in coords
    }
    
    x <- x-t*grad #move coords
    grad <-  attr(f1(x[1], x[2]), "gradient") #compute gradient at current coordinate
    fun_val[iter] <- f1(x[1], x[2]) #compute value of function at current coordinate
    #print(fun_val[iter])
    coords[iter,] <- t(x) #store coordinates
  }
  
  
  results <- c(Iterations = iter, #print step
               Coordinates = coords[iter,], #print optima coordinates
               Gradient = grad, #print gradient
               Value = fun_val[iter]) # print fn val
  attr(results, "coords") <- coords #set ordered coordinates as attribute
  attr(results, "fun_val") <- fun_val #set ordered function values as attribute
  return(results)
}

###############################
#' Hybrid Newton nad Gradient Method
#'
#' @param f2 function with gradient and Hessian
#' @param x0 start coordinates
#' @param s
#' @param alpha
#' @param beta
#' @param eps tolerance
#' @param max_iter maximum iterations
#'
#' @return results object with Iterations, Coordinates, Gradient, Function value at minimum. Trace coordinates and trace function values as attributes.
#' @export
#'
#' @examples
hybrid <- function(f2,x0, s=1, alpha = 1/4, beta = 1/2, eps = 1e-2, max_iter = 1000){
  x <- x0 #set initial coordinate
  iter <- 0 #set iterations to 0
  grad <-  attr(f2(x[1], x[2]), "gradient") #compute gradient at current coordinate
  hess <- matrix( attr( f2(x[1], x[2]), "hessian"), nrow = 2, ncol =2 ) #compute hessian at current coordinate
  #hval <- hessian(x) #set hessian
  fun_val <- rep(0, max_iter) #initiate vector of function values
  coords <- matrix(NA, nrow =  max_iter, ncol = 2) #initiate vector of coordinates
  
  
  if(is.positive.definite(hess)){
    d <- solve(hess, t(grad)) #if gradient is pos.def. use Newton direction
  } else {
    d <- grad #else use gradient direction
  }
  
  while (norm(grad, type ="2") > eps && iter < max_iter){
    iter <- iter + 1 #increment iterations
    t <- s #set t to initial stepsize
    change <- x - t* d #set change in coords
    
    while (f2(x[1],x[2]) - f2(change[1],change[2]) < alpha * t * grad %*% d) {
      t <- beta * t  #set stepsize by backtracking 
      change <- x - t* d #set change in coords
    }
    
    x <- x-t*d #move coords
    grad <-  attr(f1(x[1], x[2]), "gradient") #compute gradient at current coordinate
    hess <- matrix( attr( f2(x[1], x[2]), "hessian"), nrow = 2, ncol =2 ) #compute hessian at current coordinate
    fun_val[iter] <- f1(x[1], x[2]) #compute value of function at current coordinate
    #print(fun_val[iter])
    coords[iter,] <- t(x) #store coordinates
    
    if(is.positive.definite(hess)){
      d <- solve(hess, t(grad)) #if gradient is pos.def. use Newton direction
    } else {
      d <- grad #else use gradient direction
    }
  }
  
  results <- c(Iterations = iter, #print step
               Coordinates = coords[iter,], #print optima coordinates
               Gradient = grad, #print gradient
               Value = fun_val[iter]) # print fn val
  
  
  attr(results, "coords") <- coords #set ordered coordinates as attribute
  attr(results, "fun_val") <- fun_val #set ordered function values as attribute
  
  results
}

###############################
#' Simulated Annealing
#'
#' @param f function 
#' @param x0 start coordinates
#' @param alpha cooling rate
#' @param sig  perturbation variance
#' @param max_iter maximum iterations
#'
#' @return results object with Iterations, Coordinates, Function value at minimum. Trace coordinates and trace function values as attributes.
#' @export
#'
#' @examples
#define simulated annealing function
SimAnneal <- function(f,x0 =  c(0,0), alpha = 0.95, sig = 5, max_iter = 100000){
  x <- matrix(nrow = max_iter+1, ncol = 2) #initialise coordinate matrix
  x[1,] <- t(x0) #initialise x0
  Temp <- 1
  fun_val <- rep(0, max_iter)
  
  for (k in 1:max_iter) {
    Temp <- alpha * Temp #apply cooling
    delta <- mvrnorm(1, c(0,0) , Sigma = sig*diag(nrow= 2) ) 
    x_new <- x[k,] +  t(delta)  #perturb location
    dif <- f(x_new[1], x_new[2]) - f(x[k,1], x[k,2])  #calculate function change
    if(f(x_new[1], x_new[2]) < f(x[k,1], x[k,2])){
      x[k+1,] <- x_new #if decreasing, accept
    } else if (runif(1) > exp(dif/Temp) ){
      x[k+1,] <- x_new # accept according to acceptance probability
    } else { 
      x[k+1,] <- x[k,]
    }
    fun_val[k] <- f(x[k,1], x[k,2])
  }
  
  results <-  c(Iterations = k, #print step
                Coordinates = x[k,], #print optima coordinates
                #Gradient = grad, #print gradient
                Value = f(x[k,1], x[k,2]) ) # print fn val
  
  attr(results, "coords") <- coords #set ordered coordinates as attribute
  attr(results, "fun_val") <- fun_val #set ordered function values as attribute
  
  results
  
}
