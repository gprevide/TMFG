
detgain <- function(vid, tid, W, tri){
  k4 <- simplify2array(c(tri[tid,], vid))
  k3 <- simplify2array(c(tri[tid,]))
  det(W[k3,k3]) / det(W[k4,k4])
}

mdetgain <- function(ous, tris, W, tri){
  mapply(detgain, ous, tris, MoreArgs=list(W, tri))
}


TMFGdet <- function(W){

  N <- nrow(W)
  A <- matrix(0, nrow=N, ncol=N)
  in_v <- rep(0, N)
  tri <- matrix(0, nrow=2*N-4, ncol=3)
  separators <- matrix(0, nrow=N-4, 3)
  
# %% find 4 vertices with largest strength

  X <- W * (W > mean(W) )
  s <- rowSums(X)
  k4 <- order(s, decreasing=T)[1:4]
  in_v[1:4] <- k4
  ou_v <- setdiff(1:N, k4)
# %% build the tetrahedron with largest strength

  tri[1,] <- in_v[c(1,2,3)]
  tri[2,] <- in_v[c(2,3,4)]
  tri[3,] <- in_v[c(1,2,4)]
  tri[4,] <- in_v[c(1,3,4)]
  A[in_v, in_v] <- 1; diag(A) <- 0;
# %% build initial gain table

  gain <- matrix(-Inf, nrow=N, ncol=2*N-4)
  gain[ou_v, 1:4] <- outer(ou_v, 1:4, mdetgain,  W, tri)
  kk <- 4
  maxval <- apply(gain, 2, max)
  bestv  <- apply(gain, 2, which.max)
  for (k in 5:N) {
    tr <- which.max(maxval)
    ve <- bestv[tr]
   
    ou_v <- ou_v[ou_v != ve]
    in_v[k] <- ve
    
    # %% update adjacency matrix
    A[ve, tri[tr,]] <- 1; A[tri[tr,], ve] <-1
   
    separators[k-4, ] <- tri[tr, ]
   
    tri[kk+1,] <- simplify2array(c(tri[tr,1], tri[tr,3], ve)) # add
    tri[kk+2,] <- simplify2array(c(tri[tr,2], tri[tr,3], ve)) # add
    tri[tr  ,] <- simplify2array(c(tri[tr,1], tri[tr,2], ve)) # replace
    # %% update gain table
    
    gain[ve,   ] <- -Inf
    gain[  , tr] <- -Inf
    newt <- simplify2array(c(tr, kk+1, kk+2 ))
    gain[ou_v, newt] <- outer(ou_v, newt, mdetgain,  W, tri)
    if ( length(ou_v) > 0 ) {
      maxval[newt] <- apply(gain[, newt], 2, max)
      bestv [newt] <- apply(gain[, newt], 2, which.max)
    }
   
    kk <- kk + 2 
  }
  A <- W * as.vector (A > 0)
  # %% comutes 4-clique list
  # if nargout>3     
  #     cliques = [in_v(1:4)';separators,in_v(5:end)]; 
  # end
  cliques <- rbind(t(in_v[1:4]),
                   cbind(separators, in_v[5:N]))
  A
}

