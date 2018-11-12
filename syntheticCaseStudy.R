library(igraph)

# fix number of vertices and edges
n <- 500
m <- 40000

# auxiliary functions
repwrap <- function(x, each) rep(x=x, each=each)
assignws <- Vectorize(FUN = function(class, wmap) wmap$ws[wmap$class==class], vectorize.args = 'class')

# fix directionality
directed <- TRUE
selfloops <- TRUE

##########################
## Case Study 1.2
##########################

# generate list of all pairs of dyads
el <- expand.grid(from=1:n, to=1:n)

# fix seed
set.seed(123)

# generate activity sequence
degseq <- round(rexp(n = n, rate = n/m))

# block labels and assign labels to dyads
vlabels <- c(rep(1, length.out=n/2),rep(2, length.out=n/2))
classes <- as.vector(vlabels %*% t(vlabels))

# define block matrix
wblocks <- c(1,1)
wbetw <- 1
wmap <- data.frame(class=unique(classes), ws=c(wblocks[1],wbetw,wblocks[2]))

# sort activity sequence in decreasing order
degseq <- sort(degseq, decreasing = TRUE)

# assign weights to dyads for sampling
el$class <- classes
el$ws <- assignws(el$class, wmap)
sampleWeigths <- rep(degseq, n)
sampleWeigths <- sampleWeigths*el$ws/min(wmap$ws)

# sample id of dyads according to weights
idsample <- sample(x = nrow(el), size = m, replace = TRUE, prob = sampleWeigths)

# generate edgelist and graph
elg <- el[idsample,]
g <- graph_from_edgelist(as.matrix(elg[,c('from','to')]), directed = directed)
if(vcount(g)<n) g <- add_vertices(g, n-vcount(g))
V(g)$color <- vlabels
adj <- get.adjacency(g, sparse = F)


##########################
## Case Study 1.1
##########################

# generate list of all pairs of dyads
el <- expand.grid(from=1:n, to=1:n)

# fix seed
set.seed(123)

# generate activity sequence
degseq <- round(rexp(n = n, rate = n/m)) # rep(1,n) #

# block labels and assign labels to dyads
vlabels <- sample(c(rep(1, length.out=n/2),rep(2, length.out=n/2)))
classes <- as.vector(vlabels %*% t(vlabels))

# define block matrix
wblocks <- c(1,1)
wbetw <- 1
wmap <- data.frame(class=unique(classes), ws=c(wblocks[1],wbetw,wblocks[2]))

# sort activity sequence in decreasing order
degseq <- sort(degseq, decreasing = TRUE)

# assign weights to dyads for sampling
el$class <- classes
el$ws <- assignws(el$class, wmap)
sampleWeigths <- rep(degseq, n)
sampleWeigths <- sampleWeigths*el$ws/min(wmap$ws)

# sample id of dyads according to weights
idsample <- sample(x = nrow(el), size = m, replace = TRUE, prob = sampleWeigths)

# generate edgelist and graph
elg <- el[idsample,]
g <- graph_from_edgelist(as.matrix(elg[,c('from','to')]), directed = directed)
if(vcount(g)<n) g <- add_vertices(g, n-vcount(g))
V(g)$color <- vlabels
adj <- get.adjacency(g, sparse = F)


##########################
## Case Study 2.1
##########################

# generate list of all pairs of dyads
el <- expand.grid(from=1:n, to=1:n)

# fix seed
set.seed(123)

# generate activity sequence
degseq <- round(rexp(n = n, rate = n/m)) # rep(1,n) #

# block labels and assign labels to dyads
vlabels <- rep(1:2, length.out=n)
classes <- as.vector(vlabels %*% t(vlabels))

# define block matrix
wblocks <- c(3,1)
wbetw <- 0.1
wmap <- data.frame(class=unique(classes), ws=c(wblocks[1],wbetw,wblocks[2]))

# assign weights to dyads for sampling
el$class <- classes
el$ws <- assignws(el$class, wmap)
sampleWeigths <- rep(degseq, n)
sampleWeigths <- sampleWeigths*el$ws/min(wmap$ws)

# sample id of dyads according to weights
idsample <- sample(x = nrow(el), size = m, replace = TRUE, prob = sampleWeigths)

# generate edgelist and graph
elg <- el[idsample,]
g <- graph_from_edgelist(as.matrix(elg[,c('from','to')]), directed = directed)
if(vcount(g)<n) g <- add_vertices(g, n-vcount(g))
V(g)$color <- vlabels
adj <- get.adjacency(g, sparse = F)

##########################
## Case Study 2.2
##########################

# generate list of all pairs of dyads
el <- expand.grid(from=1:n, to=1:n)

# fix seed
set.seed(124)

# generate activity sequence
degseq <- round(rexp(n = n, rate = n/m)) # rep(1,n) #

# block labels and assign labels to dyads
vlabels2 <- vlabels
vlabels2[vlabels==1][1:(.25*n)] <- 3
classes <- as.vector(vlabels2 %*% t(vlabels2))

# define block matrix
wblocks <- c(1,1,3)
wbetw <- c(0.1,0.8)
wmap <- data.frame(class=unique(classes), ws=c(wblocks[3], wbetw[1], wbetw[1],wblocks[2],wbetw[2], wblocks[1]))

# assign weights to dyads for sampling
el$class <- classes
el$ws <- assignws(el$class, wmap)
sampleWeigths <- rep(degseq, n)
sampleWeigths <- sampleWeigths*el$ws/min(wmap$ws)

# sample id of dyads according to weights
idsample <- sample(x = nrow(el), size = m, replace = TRUE, prob = sampleWeigths)

# generate edgelist and graph
elg <- el[idsample,]
g <- graph_from_edgelist(as.matrix(elg[,c('from','to')]), directed = directed)
V(g)$color <- V(g)$labs <- vlabels2
write.graph(graph = g, file = 'fig5.graphml', format = 'graphml')
adj <- get.adjacency(g, sparse = F)
