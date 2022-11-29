Compactness <- function(g) {
        gra.geo <- distances(g, algorithm = 'dijkstra') ## get geodesics
        gra.rdist <- 1 / gra.geo  ## get reciprocal of geodesics
        diag(gra.rdist) <- NA   ## assign NA to diagonal
        gra.rdist[gra.rdist == Inf] <- 0 ## replace infinity with 0
          # Compactness = mean of reciprocal distances
        comp.igph <- mean(gra.rdist, na.rm=TRUE) 
        return(comp.igph)
        }
