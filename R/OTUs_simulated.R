OTUs_simulated <- function (data, nSam, nOTU, n_repeat, mu, size) 
{
  OTU <- array(0,dim<-c(nSam, nOTU, n_repeat))
  for (k in 1 : n_repeat) {
    otu_raw <- data
    otu_sum <- apply(otu_raw, 2, sum)
    otu <- otu_raw[, order(otu_sum, decreasing = T)[1:nOTU]]
    n = dim(otu)[1]
    Ni = apply(otu, 1, sum)
    N = sum(Ni)
    pi = apply(otu, 2, sum)/N
    Nc = (N - sum(Ni^2)/N)/(n-1)
    S = vector()
    
    for (j in 1:dim(otu)[2]) {
      Sj = sum(Ni*(otu[,j]/Ni-pi[j])^2)/(n-1)
      S = c(S, Sj)
    }
    
    G = vector()
    
    for (j in 1:dim(otu)[2]) {
      Gj = sum(Ni*(otu[,j]/Ni)*(1-otu[,j]/Ni))/(N-n)
      G = c(G, Gj)
    }
    
    theta = sum(S-G)/sum(S+(Nc-1)*Gj)
    otu.ids <- colnames(otu)
    p.est = pi
    
    gplus <- (1 - theta)/theta
    p.est <- p.est[otu.ids]
    g.est <- p.est * gplus
    scale2 = function(x) as.numeric(scale(x))
    comm <- matrix(0, nSam, length(g.est))
    rownames(comm) <- 1:nrow(comm)
    colnames(comm) <- names(g.est)
    comm.p <- comm
    nSeq <- rnbinom(nSam, mu = mu, size = size)
    
    for (i in 1:nSam) {
      comm.p[i, ] <- MiSPU::rdirichlet(1, g.est)[1, ]
      comm[i, ] <- rmultinom(1, nSeq[i], prob = comm.p[i, ])[, 1]
    }
    
    OTU_individual = comm[, otu.ids]
    OTU[,, k] <- OTU_individual
  }

  OTU_simulated <- apply(OTU, c(1,2), mean)
  return(list(OTU_simulated = OTU_simulated))
}