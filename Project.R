# Specified number of intervals.
n <- 50

# Setting a helping constant.
h <- 2/n

# Defining a function of elastic deformation.
E <- function(x) {
  if(x <= 1) {
    return(3)
  }else{
    return(5)
  }
}

# Defining e function (Not exponential, just name "e").
e <- function(i) {
  return(function(x) {
    return(pmax(1 - abs((x-i*h)/h), 0))
  })
}

# Defining derivative of e.
eDir <- function(i) {
  return(function(x) {
    ei <- e(i)
    if(ei(x) == 0) {
      return(0)
    }
    if(x >= i*h) {
      return(-1/h)
    }
    return(1/h)
  })
}

# Defining the integrating function.
integrate <- function(f, a, b) {
  x1 <- -1/sqrt(3)
  x2 <- 1/sqrt(3)
  result <- 0
  result <- result + f((b-a)/2*x1 + (b + a)/2)
  result <- result + f((b-a)/2*x2 + (b + a)/2)
  return(result*(b - a)/2)
}

# Defining B function.
B <- function(e1, e1Dir, e2, e2Dir, lower, upper) {
  integrand <- function(x) e1Dir(x) * e2Dir(x) * E(x)
  result <- integrate(integrand, lower, upper)
  temp <- 3 * e1(0) * e2(0)
  return(result - temp)
}

# Defining L function.
L <- function(e2) {
  result <- 30 * e2(0)
  return(-result)
}

# Calculating matrix of B.
matrix <- matrix(0, ncol = n, nrow = n)
for (i in 1:n) {
  for (j in 1:n) {
    ei <- e(i - 1)
    ej <- e(j - 1)
    eiDir <- eDir(i - 1)
    ejDir <- eDir(j - 1)
    if(i != j) {
      mn <- max(0, min(i - 1, j - 1) * h)
      mx <- min(2, max(i - 1, j - 1) * h)
    } else {
      mn <- max(0, (i - 2) * h)
      mx <- min(2, (i) * h)
    }
    result <- B(ei, eiDir, ej, ejDir, mn, mx)
    matrix[i, j] <- result
  }
}

# Calculating vector of L.
vector <- rep(0, times=n)
for (i in 1:n) {
  ei <- e(i - 1)
  result <- L(ei)
  vector[i] <- result
}

# Solving system of equations.
weights <- solve(matrix, vector)

# Calculating value of u(x).
u <- function(x) {
  result <- 0
  for(i in 1:n) {
    ei <- e(i-1)
    temp <- weights[i]
    result <- result + ei(x) * temp
  }
  return(result)
}

# Creating chart
x <- seq(0, 2, by = h)
result = rep(0, times=n+1)
for (i in 1:(n+1)) {
  result[i] <- u(x[i])
}

# Displaying a chart
myMain <- paste("Chart of u(x) for n = ", n)
plot(x, result, type = "l", col = "blue", lwd = 2, xlab = "", ylab = "", main = myMain)