library("CompQuadForm")
multi_tissue_p = function(y, m) {
  # y is a vector containing the continuous outcomes;
  # m contains the levels of multiple predictors Each row represent an sample (ordered the same as y); each column represent one predictor
  # no NA is allowed in y or m
  
  # first normalize each column, and make sure no columns are duplicated
  newm = c()
  m = as.matrix(m)
  for (i in 1:ncol(m)) {
    mean_col = mean(m[,i]); sd_col = sd(m[,i])
    
    m[,i] = m[,i] -mean_col;
    if (sd_col==0) {
      next;
    } else {
      m[,i] = m[,i]/ sd_col
    }
    if (length(newm)==0) {
      newm =as.matrix(m[,i])
    } else if (not_exist(newm, m[,i])) {      
      newm =cbind(newm, m[,i])
    }
  }
  
  if (length(newm)==0) {
    return(1.0)
  }
  
  p=1.0
  try({
    p = score_test_p(y, newm)
  })
  return (p) # return the p-value of the test
}

score_test_p = function(y, G){ 
  n = length(y)
  mu = rep(mean(y), n)
  K = G %*% t(G)
  
  Q = t(y-mu) %*% K %*% (y-mu)  
  
  P0.sqrt = sd(y) * (diag(n) - matrix(rep(1,n^2),nrow=n)/n)
  
  lambdas = eigen(P0.sqrt %*% K %*% P0.sqrt)$values
  p = davies(Q, lambdas)$Qq
  return(p)
}

not_exist = function(m, vec) {
  #check if vec (vector) also exists in some columns of m
  if (length(m)==0) {
    return(TRUE)
  }
  for (i in 1:ncol(m)) {
    if (sum(m[,i]!=vec)==0) {
      return(FALSE)
    }
  }
  return (TRUE)
}

# 
# # # test case 1 (alternative model)
# n = 500
# x1 = rnorm(n); x2=x1; x3 = rnorm(n)
# y = x1+x2+x3+rnorm(n)
# G = cbind(x1,x2,x3)
# score_test_p(y,G)
# multi_tissue_p(y,G)
# 
# # 
# # # test case 2 (null model)
# x1 = rnorm(n); x2=rnorm(n); x3 = rnorm(n)
# y = as.vector(rnorm(n))
# G = as.matrix(cbind(x1,x2,x3))
# score_test_p(y,G)
# multi_tissue_p(y,G)
