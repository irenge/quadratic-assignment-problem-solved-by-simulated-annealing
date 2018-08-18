set.seed(64)

# Check if needed package installed if not install packages 

packages <- c("permute","cluster","seriation")

pckg_need <- setdiff(packages, rownames(installed.packages()))
                     
if(length(pckg_need)> 0){
  install.packages(pckg_need)
}

####################### FUNCTION NEEDED ##################################

calc_cost <- function(p,n,D,obj_value){
  for(i in 1:n){
    for(j in 1:n){
      obj_value <- obj_value + (D[p[i],p[j]]*(-abs(i-j)))
    }
  }
  obj_value
}

# Load libraries

library('permute')
library('cluster')

# load data 

data("BOD")
euclid <- dist(BOD)

#

D <- as.matrix(euclid)

n <- ncol(D)
m <- n -1
p <- sample(n)

obj_value<-0

best_value <- calc_cost(p,n,D,obj_value)

mxfail <- n*(n-1)/2
dmin <- Inf
dmax <- 0 
iteration <- 100

# Heating process 

for(k in 1:10000){
  r <- sample(1:n,1)
  s <- sample(1:m,1)
  if(s>=r){ s<-s+1}
  j <- r
  i<-s
  
  # calc delta 
  delta <- 0
  for(k in 1:n){
    if((k != i) &&(k!=j)){
      delta <- delta +((-abs(j-k)) - (-abs(i-k)))*(D[p[i],p[k]]-D[p[j],p[k]])
    }
  }
  
  d <- 2*delta
  
  ## end of delta
  
  if(d >0){
    dmin<-min(dmin,d)
    dmax <- max(dmax,d)
  }
obj_value <-obj_value +d 
p <- replace(p,c(r,s),p[c(s,r)])
}

t0 <- dmin +(dmax-dmin)/10.0
tf <- dmin 
beta <- (t0-tf)/(k*t0*tf)
fail_count <- 0
tfound <- t0
temrature <- t0 
r <- 1
s<-2

