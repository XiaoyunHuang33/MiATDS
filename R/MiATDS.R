MiATDS <- function (y, otu.tab, tree = NULL, cov = NULL, model = c("gaussian", "binomial"), pow = c(1:5), comp = FALSE, CLR = FALSE, opt.ncl = 30, n.perm = 5000) {

  if (!comp) {
    otu.tab.com <- t(apply(otu.tab,1,function(x)x/sum(x)))
  }
  com.otu.tab <- otu.tab.com
  if (CLR) {
    otu.tab.com <- t(apply(com.otu.tab, 1, clr))
  }
  n <- length(y)
  p <- ncol(otu.tab.com)

  if (is.null(cov)) {
    r <- y - mean(y)
  } else {
    fit <- glm(y ~ ., family = model, data = as.data.frame(cov))
    r <- y - fitted.values(fit)
  }

  r.s <- list()

  for (j in 1:n.perm) {
    r.s[[j]] <- r[shuffle(length(r))]
  }

  n.Zs <- as.numeric(r) %*% otu.tab.com
  n.Z0s <- lapply(r.s, function(x) as.numeric(x %*% otu.tab.com))
  U0s <- lapply(apply(sapply(r.s, function(x) return(x %*% otu.tab.com)), 1, list), unlist)
  se <- sapply(U0s, sd)
  Zs <- n.Zs/se
  Z0s <- lapply(n.Z0s, function(x) x/se)
  if (class(tree) == "phylo") {
    CDis <- cophenetic(tree)
    for (j in 1:nrow(CDis)) {
      ind.t <- which(CDis[j, ] == 0)
      ind.d <- which(ind.t == j)
      ind.r <- ind.t[-ind.d]
      CDis[j, ind.r] <- min(CDis[j, -ind.t])/2
    }
    if (opt.ncl == "gmx") {
      asw <- numeric(p - 1)
      for (j in 2:(p - 1)) {
        asw[j] <- pam(CDis, j)$silinfo$avg.width
      }
      k.best <- which.max(asw)
      clust.pam <- pam(CDis, k.best, cluster.only = TRUE)
    }
    if (opt.ncl == "fmx") {
      asw <- numeric(p - 1)
      j <- 2
      asw[j] <- pam(CDis, j)$silinfo$avg.width
      while (asw[j - 1] < asw[j] & j <= p - 1) {
        j <- j + 1
        asw[j] <- pam(CDis, j)$silinfo$avg.width
      }
      k.best <- j - 1
      clust.pam <- pam(CDis, k.best, cluster.only = TRUE)
    }
    if (opt.ncl != "gmx" & opt.ncl != "fmx") {
      asw <- numeric(p - 1)
      for (j in 2:opt.ncl) {
        asw[j] <- pam(CDis, j)$silinfo$avg.width
      }
      k.best <- which.max(asw)
      clust.pam <- pam(CDis, k.best, cluster.only = TRUE)
    }
    Ws <- rep(NA, p)
    for (j in 1:p) {
      ind.1 <- match(colnames(CDis)[j], names(clust.pam))
      ind.2 <- which(clust.pam == as.numeric(clust.pam[ind.1]))
      ind.3 <- which(ind.2 == ind.1)
      ind <- ind.2[-ind.3]
      if (length(ind) == 0) {
        Ws[j] <- 1
      }
      else {
        inv.Ds <- 1/(CDis[j, ind])
        abs.Zs <- abs(Zs[ind])
        Ws[j] <- sum(inv.Ds * abs.Zs)/sum(inv.Ds) + 1
      }
    }
    Ws <- Ws/sum(Ws)
    hs <- 1
    ahc.o <- as.list(MiHC.stat(Zs = Zs, hs = hs, Ws = Ws))
    ahc.p <- lapply(apply(sapply(Z0s, function(x) MiHC.stat(Zs = x, hs = hs, Ws = Ws)), 1, list), unlist)
    pvs <- mapply(function(x, y) (length(which(x < y)) + 1)/(n.perm + 1), ahc.o, ahc.p)
    ind.pvs <- 1 - pchisq(Zs^2, df = 1)
    simes.pv <- round(min(length(ind.pvs)*ind.pvs/rank(ind.pvs)),4)
    simes.pv0s.ori <- unlist(lapply(Z0s, function(x) {xx <- 1-pchisq(x^2, df=1); min(length(xx)*xx/rank(xx))}))
    simes.pv0s <- rep(NA, n.perm)
    for (j in 1:n.perm) {
      simes.pv0s[j] <- (length(which(simes.pv0s.ori[-j] < simes.pv0s.ori[j])) + 1)/(n.perm + 1)
    }
    simes.pv <- (length(which(simes.pv0s < simes.pv)) + 1)/(n.perm + 1)
    l.hs <- length(hs)
    Tu <- min(pvs[1:l.hs], simes.pv)
    T0u <- rep(NA, n.perm)
    Tw <- min(pvs[(l.hs+1):(l.hs*2)], simes.pv)
    T0w <- rep(NA, n.perm)
    for (j in 1:n.perm) {
      T0u.o <- lapply(ahc.p[1:l.hs], function(x) x[j])
      T0u.p <- lapply(ahc.p[1:l.hs], function(x) x[-j])
      T0u[j] <- min(mapply(function(x, y) (length(which(x < y)) + 1)/(n.perm + 1), T0u.o, T0u.p))
    }
    for (j in 1:n.perm) {
      T0w.o <- lapply(ahc.p[(l.hs+1):(l.hs*2)], function(x) x[j])
      T0w.p <- lapply(ahc.p[(l.hs+1):(l.hs*2)], function(x) x[-j])
      T0w[j] <- min(mapply(function(x, y) (length(which(x < y)) + 1)/(n.perm + 1), T0w.o, T0w.p))
    }
    T0u.minp <- apply(cbind(T0u, simes.pv0s),1,min)
    T0w.minp <- apply(cbind(T0w, simes.pv0s),1,min)
    T.MiHC <- min(pvs, simes.pv)
    T0.MiHC <- apply(cbind(T0u, T0w, simes.pv0s),1,min)
    pvs <- ind.pvs
    i0 <- which(pvs < 1e-08)
    i1 <- which(pvs > 0.99999999)
    pvs[i0] <- 1e-08
    pvs[i1] <- 0.99999999
    Is <- order(order(pvs))
    otu.ids <- colnames(otu.tab)
    exp.HC <- (Is/p)/sqrt(pvs * (1 - pvs)/p)
    obs.HC <- pvs/sqrt(pvs * (1 - pvs)/p)
    pd <- abs(exp.HC-obs.HC)
    pd.weight <- apply(as.matrix(pd), 2, function(x){x/sum(pd)})
    pd.weight <- pd.weight + 1
    pd.pi.weight <- Ws*pd.weight
    total.reads <- rowSums(otu.tab)
    prop <- apply(otu.tab/total.reads, 2, scale)
    U <- as.vector(t(prop) %*% r)
    Ts.wSPU = rep(NA, length(pow))
    Ts.WSPU = rep(NA, length(pow))
    Ts.wSPU <- unlist(lapply(as.list(pow), function(x) return(sum(diag(pd.weight) %*% U^x))))
    Ts.WSPU <- unlist(lapply(as.list(pow), function(x) return(sum(diag(pd.pi.weight) %*% U^x))))
    U0 <- lapply(r.s, function(x) return(as.vector(t(prop) %*% x)))
    T0s.wSPU <- list()
    T0s.WSPU <- list()
    pvs.wSPU <- rep(NA, length(pow))
    pvs.WSPU <- rep(NA, length(pow))
    for (j in 1:length(pow)) {
      T0s.wSPU[[j]] <- T0s.e.wSPU <- sapply(U0, function(x) return(sum(diag(pd.weight) %*% x^pow[j])))
      T0s.WSPU[[j]] <- T0s.e.WSPU <- sapply(U0, function(x) return(sum(diag(pd.pi.weight) %*% x^pow[j])))
      pvs.wSPU[j] <- (sum(abs(T0s.e.wSPU) > abs(Ts.wSPU[j])) + 1)/(n.perm + 1)
      pvs.WSPU[j] <- (sum(abs(T0s.e.WSPU) > abs(Ts.WSPU[j])) + 1)/(n.perm + 1)
    }
    T.awSPU <- min(pvs.wSPU)
    T.aWSPU <- min(pvs.WSPU)
    T0.awSPU <- rep(NA, n.perm)
    T0.aWSPU <- rep(NA, n.perm)
    for (l in 1:n.perm) {
      T0s.n.wSPU <- list()
      T0s.n.WSPU <- list()
      for (m in 1:length(pow)) {
        T0s.n.wSPU[[m]] <- T0s.wSPU[[m]][-l]
        T0s.n.WSPU[[m]] <- T0s.WSPU[[m]][-l]
      }
      a.U <- as.vector(t(prop) %*% r.s[[l]])
      a.pvs.wSPU <- rep(NA, length(pow))
      a.pvs.WSPU <- rep(NA, length(pow))
      a.Ts.wSPU <- rep(NA, length(pow))
      a.Ts.WSPU <- rep(NA, length(pow))
      a.Ts.wSPU <- unlist(lapply(as.list(pow), function(x) return(sum(diag(pd.weight) %*% a.U^x))))
      a.Ts.WSPU <- unlist(lapply(as.list(pow), function(x) return(sum(diag(pd.pi.weight) %*% a.U^x))))
      a.pvs.wSPU <- unlist(mapply(function(x, y) (sum(abs(x) > abs(y)) + 1)/(n.perm + 1), T0s.n.wSPU, a.Ts.wSPU))
      a.pvs.WSPU <- unlist(mapply(function(x, y) (sum(abs(x) > abs(y)) + 1)/(n.perm + 1), T0s.n.WSPU, a.Ts.WSPU))
      T0.awSPU[l] <- min(a.pvs.wSPU)
      T0.aWSPU[l] <- min(a.pvs.WSPU)
    }
    p.awSPU <- (sum(T0.awSPU < T.awSPU) + 1)/(n.perm + 1)
    p.aWSPU <- (sum(T0.aWSPU < T.aWSPU) + 1)/(n.perm + 1)

    M.MiHC.awSPU0 <- apply(cbind(T0.MiHC, T0.awSPU), 1, min)
    M.MiHC.awSPU <- min(T.MiHC, T.awSPU)
    p.MiHC.awSPU.omnibus <- (sum(M.MiHC.awSPU0 < M.MiHC.awSPU) + 1)/(n.perm + 1)

    M.MiHC.aWSPU0 <- apply(cbind(T0.MiHC, T0.aWSPU), 1, min)
    M.MiHC.aWSPU <- min(T.MiHC, T.aWSPU)
    p.MiHC.aWSPU.omnibus <- (sum(M.MiHC.aWSPU0 < M.MiHC.aWSPU) + 1)/(n.perm + 1)

    M.MiATDS0 <- apply(cbind(T0.MiHC, T0.awSPU, T0.aWSPU), 1, min)
    M.MiATDS <- min(T.MiHC, T.awSPU, T.aWSPU)
    p.MiATDS <- (sum(M.MiATDS0 < M.MiATDS) + 1)/(n.perm + 1)

    pvs.awSPU <- c(pvs.wSPU, p.awSPU)
    pvs.aWSPU <- c(pvs.WSPU, p.aWSPU)
    omnibus.pvs <- c(p.MiHC.awSPU.omnibus, p.MiHC.aWSPU.omnibus, p.MiATDS)
    names(pvs.awSPU) <- c(paste("wSPU(", pow, ")", sep = ""), "awSPU")
    names(pvs.aWSPU) <- c(paste("WSPU(", pow, ")", sep = ""), "aWSPU")
    names(omnibus.pvs) <- c("omnibus.MiHC.awSPU", "omnibus.MiHC.aWSPU", "MiATDS")
    pd.rank <- otu.ids[order(pd, decreasing = T)]

    return(list(pd.rank = pd.rank, awSPU.pvs = pvs.awSPU, aWSPU.pvs = pvs.aWSPU, omnibus.pvs = omnibus.pvs))
  }
  if (tree != NULL) {
    hs <- 1
    ahc.o <- as.list(MiHC.stat(Zs = Zs, hs = hs, Ws = NULL))
    ahc.p <- lapply(apply(sapply(Z0s, function(x) MiHC.stat(Zs = x, hs = hs, Ws = NULL)), 1, list), unlist)
    pvs <- mapply(function(x, y) (length(which(x < y)) + 1)/(n.perm + 1), ahc.o, ahc.p)
    ind.pvs <- 1 - pchisq(Zs^2, df = 1)
    simes.pv <- min(length(ind.pvs) * ind.pvs/rank(ind.pvs))
    simes.pv0s.ori <- unlist(lapply(Z0s, function(x) {
      xx <- 1 - pchisq(x^2, df = 1)
      min(length(xx) * xx/rank(xx))
    }))
    simes.pv0s <- rep(NA, n.perm)
    for (j in 1:n.perm) {
      simes.pv0s[j] <- (length(which(simes.pv0s.ori[-j] < simes.pv0s.ori[j])) + 1)/(n.perm + 1)
    }
    simes.pv <- (length(which(simes.pv0s < simes.pv)) + 1)/(n.perm + 1)
    l.hs <- length(hs)
    Tu <- min(pvs, simes.pv)
    T0u <- rep(NA, n.perm)
    for (j in 1:n.perm) {
      T0u.o <- lapply(ahc.p, function(x) x[j])
      T0u.p <- lapply(ahc.p, function(x) x[-j])
      T0u[j] <- min(mapply(function(x, y) (length(which(x < y)) + 1)/(n.perm + 1), T0u.o, T0u.p))
    }
    T0u.minp <- apply(cbind(T0u, simes.pv0s), 1, min)
    pvs <- ind.pvs
    i0 <- which(pvs < 1e-08)
    i1 <- which(pvs > 0.99999999)
    pvs[i0] <- 1e-08
    pvs[i1] <- 0.99999999
    Is <- order(order(pvs))
    otu.ids <- colnames(otu.tab)
    exp.HC <- (Is/p)/sqrt(pvs * (1 - pvs)/p)
    obs.HC <- pvs/sqrt(pvs * (1 - pvs)/p)
    pd <- abs(exp.HC-obs.HC)
    pd.weight <- apply(as.matrix(pd),2, function(x){x/sum(pd)})
    pd.weight <- pd.weight + 1
    total.reads <- rowSums(otu.tab)
    prop <- apply(otu.tab/total.reads, 2, scale)
    U <- as.vector(t(prop) %*% r)
    Ts.wSPU = rep(NA, length(pow))
    Ts.wSPU <- unlist(lapply(as.list(pow), function(x) return(sum(diag(pd.weight) %*% U^x))))
    U0 <- lapply(r.s, function(x) return(as.vector(t(prop) %*% x)))
    T0s.wSPU <- list()
    pvs.wSPU <- rep(NA, length(pow))
    for (j in 1:length(pow)) {
      T0s.wSPU[[j]] <- T0s.e.wSPU <- sapply(U0, function(x) return(sum(diag(pd.weight) %*% x^pow[j])))
      pvs.wSPU[j] <- (sum(abs(T0s.e.wSPU) > abs(Ts.wSPU[j])) + 1)/(n.perm + 1)
    }
    T.awSPU <- min(pvs.wSPU)
    T0.awSPU <- rep(NA, n.perm)
    for (l in 1:n.perm) {
      T0s.n.wSPU <- list()
      for (m in 1:length(pow)) {
        T0s.n.wSPU[[m]] <- T0s.wSPU[[m]][-l]
      }
      a.U <- as.vector(t(prop) %*% r.s[[l]])
      a.pvs.wSPU <- rep(NA, length(pow))
      a.Ts.wSPU <- rep(NA, length(pow))
      a.Ts.wSPU <- unlist(lapply(as.list(pow), function(x) return(sum(diag(pd.weight) %*% a.U^x))))
      a.pvs.wSPU <- unlist(mapply(function(x, y) (sum(abs(x) > abs(y)) + 1)/(n.perm + 1), T0s.n.wSPU, a.Ts.wSPU))
      T0.awSPU[l] <- min(a.pvs.awSPU)
    }

    p.awSPU <- (sum(T0.awSPU < T.awSPU) + 1)/(n.perm + 1)

    M.MiATDS0 <- apply(cbind(T0u.minp, T0.aWSPU), 1, min)
    M.MiATDS <- min(Tu, T.aWSPU)
    p.MiATDS <- (sum(M.MiATDS0 < M.MiATDS) + 1)/(n.perm + 1)

    pvs.awSPU <- c(pvs.wSPU, p.awSPU)
    omnibus.pvs <- p.MiATDS
    names(pvs.awSPU) <- c(paste("wSPU(", pow, ")", sep = ""), "awSPU")
    names(omnibus.pvs) <- "MiATDS"
    pd.rank <- otu.ids[order(pd, decreasing = T)]

    return(list(pd.rank = pd.rank, awSPU.pvs = pvs.awSPU, omnibus.pvs = omnibus.pvs))
  }

}
