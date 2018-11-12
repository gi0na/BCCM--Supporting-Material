library(igraph)
library(BiasedUrn)
library(poweRlaw)

###########
# auxiliary functions
###########
block_to_omega <- function(B, Bsize){
  b <- nrow(B)
  B[rep(1:b, each=Bsize),rep(1:b, each=Bsize)]
}

mat2vec.ix <- function(mat, directed,
                       selfloops) {
  # Returns the indices to
  # vectorise adjacency matrices
  # removing unused entries in the
  # case of undirected or no
  # selfloops graphs
  if(nrow(mat)==ncol(mat)){
    if (!directed) {
      mat[lower.tri(mat, !selfloops)] <- NA
    } else {
      if (!selfloops)
        diag(mat) <- NA
    }
  }
  ix <- !is.na(mat)
  return(ix)
}


vec2mat <- function(vec,directed,selfloops,n){
  # Generates adjacency matrix from vector
  if(length(n)>1){
    mat <- matrix(0,n[2],n[3])
  } else{
    mat <- matrix(0,n,n)
  }
  
  idx <- mat2vec.ix(mat,directed,selfloops)
  mat[idx] <- vec
  return(mat)
}
##################

##########################
## Example 1.1
##########################

# number of blocs
b <- 5

# block matrix
B <- matrix(c(
  1 , 0.1 , 0 , 0 , 0,
  0 , 1 , 0.1 , 0 , 0,
  0 , 0 , 1 , 0.1 , 0,
  0 , 0 , 0 , 1 , 0.1,
  0.1 , 0 , 0 , 0 , 1
), nrow = 5, byrow = TRUE)

# number of vertices
n=50
# number of edges
m=500

# directionality
directed <- TRUE
selfloops <- FALSE

# size of blocks
Bsize <- n/b

# propensity matrix from block matrix
omega <- block_to_omega(B,Bsize)

# combinatorial matrix
xi <- matrix(m^2/(n*(n-1)), n, n)

# generate graph from Wallenius distribution
idx <- mat2vec.ix(xi, directed, selfloops)
adj <- vec2mat(vec = rMWNCHypergeo(nran = 1, m = xi[idx], n = m, odds = omega[idx]), directed = directed, selfloops = selfloops, n = n)
g <- graph_from_adjacency_matrix(adj)
V(g)$color <- V(g)$labels <- rep(1:b, each=Bsize)
V(g)$degree <- degree(g, mode = 'out')

write_graph(graph = g, file = 'fig3a.graphml', format = 'graphml')

# compare degree sequence obtained from planted one
set.seed(111)
gees <- replicate(n = 1000, 
                  expr = vec2mat(vec = 
                                   rMWNCHypergeo(nran = 1, m = xi[idx], n = m, odds = omega[idx]), directed = directed, selfloops = selfloops, n = n),
                  simplify = FALSE)
outdegseqs <- sapply(X = gees, FUN = function(mat) apply(mat,1,sum))
(meandegs <- apply(outdegseqs,1,mean))
cbind(planted=rep(m/n,n),empirical=meandegs)


##########################
## Example 1.2
##########################

# generate power law degree sequences
alpha <- 1.8
seed <- 125
set.seed(seed)
out_degseq <- rep(rpldis(Bsize, 1, alpha), b)
sum(out_degseq)
out_degseq[2+(0:(b-1))*Bsize] <- 2

in_degseq <- rep(sum(out_degseq[1:Bsize])/Bsize, n)

# combinatorial matrix
xi <- out_degseq %*% t(in_degseq)

# generate graph from Wallenius distribution
adj <- vec2mat(vec = rMWNCHypergeo(nran = 1, m = xi[idx], n = m, odds = omega[idx]), directed = directed, selfloops = selfloops, n = n)
g <- graph_from_adjacency_matrix(adj)
V(g)$color <- V(g)$labels <- rep(1:b, each=Bsize)
V(g)$degree <- degree(g, mode = 'out')
write_graph(graph = g, file = 'fig3b.graphml', format = 'graphml')

# compare degree sequence obtained from planted one
set.seed(111)
gees <- replicate(n = 1000, 
                  expr =  vec2mat(vec = rMWNCHypergeo(nran = 1, m = xi[idx], n = m, odds = omega[idx]), directed = directed, selfloops = selfloops, n = n),
                  simplify = FALSE)
outdegseqs <- sapply(X = gees, FUN = function(mat) apply(mat,1,sum))
(meandegs <- apply(outdegseqs,1,mean))
cbind(planted=out_degseq,empirical=meandegs)

##########################
## Example 1.3
##########################

# generate power law degree sequences
seed <- 143
set.seed(seed)
out_degseq <- rpldis(n, 1, alpha)
sum(out_degseq)
out_degseq[out_degseq==1][1:3] <- 2

in_degseq <- rep(m/n, n)

# combinatorial matrix
xi <- out_degseq %*% t(in_degseq)

# generate graph from Wallenius distribution
adj <- vec2mat(vec = rMWNCHypergeo(nran = 1, m = xi[idx], n = m, odds = omega[idx]), directed = directed, selfloops = selfloops, n = n)
g <- graph_from_adjacency_matrix(adj)
V(g)$color <- V(g)$labels <- rep(1:b, each=Bsize)
V(g)$degree <- degree(g, mode = 'out')
write_graph(graph = g, file = 'fig3c.graphml', format = 'graphml')

# compare degree sequence obtained from planted one
set.seed(111)
gees <- replicate(n = 1000, 
                  expr =  vec2mat(vec = rMWNCHypergeo(nran = 1, m = xi[idx], n = m, odds = omega[idx]), directed = directed, selfloops = selfloops, n = n),
                  simplify = FALSE)
outdegseqs <- sapply(X = gees, FUN = function(mat) apply(mat,1,sum))
(meandegs <- apply(outdegseqs,1,mean))
cbind(planted=out_degseq,empirical=meandegs)

#################
#################

##########################
## Example 2.1
##########################

# number of blocs
b <- 5

# block matrix
B <- matrix(c(
  1 , 1 , 0 , 0 , 0,
  0 , 1 , 1 , 0 , 0,
  0 , 0 , 1 , 1 , 0,
  0 , 0 , 0 , 1 , 1,
  1 , 0 , 0 , 0 , 1
), nrow = 5, byrow = TRUE)

# size of blocks
Bsize <- n/b

# propensity matrix from block matrix
omega <- block_to_omega(B,Bsize)

# generate degree sequence from power-law
alpha <- 1.8
seed <- 125
set.seed(seed)
out_degseq <- rep(rpldis(Bsize, 1, alpha), b)
sum(out_degseq)
out_degseq[2+(0:(b-1))*Bsize] <- 2

in_degseq <- rep(sum(out_degseq[1:Bsize])/Bsize, n)

# combinatorial matrix
xi <- out_degseq %*% t(in_degseq)

# generate graph from Wallenius distribution
adj <- vec2mat(vec = rMWNCHypergeo(nran = 1, m = xi[idx], n = m, odds = omega[idx]), directed = directed, selfloops = selfloops, n = n)
g <- graph_from_adjacency_matrix(adj)
V(g)$color <- V(g)$labels <- rep(1:b, each=Bsize)
V(g)$degree <- degree(g, mode = 'out')

write_graph(graph = g, file = 'fig4a.graphml', format = 'graphml')

# compare degree sequence obtained from planted one
set.seed(111)
gees <- replicate(n = 1000, 
                  expr =  vec2mat(vec = rMWNCHypergeo(nran = 1, m = xi[idx], n = m, odds = omega[idx]), directed = directed, selfloops = selfloops, n = n),
                  simplify = FALSE)
outdegseqs <- sapply(X = gees, FUN = function(mat) apply(mat,1,sum))
(meandegs <- apply(outdegseqs,1,mean))
cbind(planted=out_degseq,empirical=meandegs)

##########################
## Example 2.2
##########################

# number of blocs
b <- 5

# block matrix
B <- matrix(c(
  10  , 0.1 , 0   , 0   , 0  ,
  0   , 1   , 0.1 , 0   , 0  ,
  0   , 0   , 1   , 0.1 , 0  ,
  0   , 0   , 0   , 1   , 0.1,
  0.1 , 0   , 0   , 0   , 1
), nrow = 5, byrow = TRUE)

# propensity matrix from block matrix
omega <- block_to_omega(B,Bsize)

# generate degree sequence from power-law
alpha <- 1.8
seed <- 125
set.seed(seed)
out_degseq <- rep(rpldis(Bsize, 1, alpha), b)
sum(out_degseq)
out_degseq[2+(0:(b-1))*Bsize] <- 2

in_degseq <- rep(sum(out_degseq[1:Bsize])/Bsize, n)

# combinatorial matrix
xi <- out_degseq %*% t(in_degseq)

# generate graph from Wallenius distribution
adj <- vec2mat(vec = rMWNCHypergeo(nran = 1, m = xi[idx], n = m, odds = omega[idx]), directed = directed, selfloops = selfloops, n = n)
g <- graph_from_adjacency_matrix(adj)
V(g)$color <- V(g)$labels <- rep(1:b, each=Bsize)
V(g)$degree <- degree(g, mode = 'out')

write_graph(graph = g, file = 'fig4b.graphml', format = 'graphml')

##########################
## Example 2.3
##########################

# number of blocs
b <- 5

# block matrix
B <- matrix(c(
  1   , 0.8 , 0   , 0   , 0  ,
  0   , 1   , 0.1 , 0   , 0  ,
  0   , 0   , 1   , 0.8 , 0  ,
  0   , 0   , 0   , 1   , 0.8,
  0.1 , 0   , 0   , 0   , 1
), nrow = 5, byrow = TRUE)


# propensity matrix from block matrix
omega <- block_to_omega(B,Bsize)

# generate degree sequence from power-law
alpha <- 1.8
seed <- 125
set.seed(seed)
out_degseq <- rep(rpldis(Bsize, 1, alpha), b)
sum(out_degseq)
out_degseq[2+(0:(b-1))*Bsize] <- 2

in_degseq <- rep(sum(out_degseq[1:Bsize])/Bsize, n)

# combinatorial matrix
xi <- out_degseq %*% t(in_degseq)

# generate graph from Wallenius distribution
adj <- vec2mat(vec = rMWNCHypergeo(nran = 1, m = xi[idx], n = m, odds = omega[idx]), directed = directed, selfloops = selfloops, n = n)
g <- graph_from_adjacency_matrix(adj)
V(g)$color <- V(g)$labels <- rep(1:b, each=Bsize)
V(g)$degree <- degree(g, mode = 'out')

write_graph(graph = g, file = 'fig4c.graphml', format = 'graphml')
