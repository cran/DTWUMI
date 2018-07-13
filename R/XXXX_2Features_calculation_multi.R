# Local feature
.local_feature_multi <- function(q){
  out <- matrix(nrow = nrow(q), ncol = 2*ncol(q))
  out[1, ] <- q[1, ]
  for (i in 2:(nrow(q)-1))
  {
    t <- c(as.numeric(q[i, ])-as.numeric(q[(i-1), ]), as.numeric(q[i, ])-as.numeric(q[(i+1), ]))
    out[i, ] <- t
  }
  out[nrow(q), ] <- as.numeric(q[nrow(q), ])
  return(out)
}

# Global feature
.global_feature_multi <- function(q){
  out <- matrix(nrow = nrow(q), ncol = 2*ncol(q))
  out[1, ] <- q[1, ]
  out[2, ] <- c(as.numeric(q[2, ])-as.numeric(q[1, ]), as.numeric(q[2, ])-as.numeric(q[3, ]))
  for (i in 3:(nrow(q)-2))
  {
    t <- c(as.numeric(q[i, ])-as.numeric(colMeans(q[1:(i-1), ], na.rm = T)), as.numeric(q[i, ])-as.numeric(colMeans(q[(i+1):nrow(q), ], na.rm = T)))
    out[i, ] <- t
  }
  out[nrow(q)-1, ] <- c(as.numeric(q[nrow(q)-1, ])-as.numeric(q[nrow(q)-2, ]), as.numeric(q[nrow(q)-1, ])-as.numeric(q[nrow(q), ]))
  out[nrow(q), ] <- as.numeric(q[nrow(q), ])
  return(out)
}

# Distance matrix
.manhattan_dist_multi <- function(p, q){
  return(sum(abs(p-q), na.rm = T))
}

.dist_matrix_multi <- function(q, r){
  dist <- matrix(10000, nrow(q), nrow(r))
  
  for(i in 2:(nrow(q)-1)){
    for (j in 2:(nrow(r)-1)){
      t <- .manhattan_dist_multi(as.numeric(q[i, ]), as.numeric(r[j, ]))
      dist[i, j] <- t
    }
  }
  return(dist)
}