
sojo.phenotype <- function (sum.stat.discovery, sum.stat.validation = NULL, cor.X, 
                            v.y = 1, v.y.validation = NULL, lambda.vec = NA, standardize = T, nvar = 50) 
{
  # colnames_input <- c("SNP", "A1", "A2", "b", "se", "N")
  # colnames_lack <- setdiff(colnames_input, intersect(colnames(sum.stat.discovery), 
  #                                                    colnames_input))
  # if (length(colnames_lack) > 0) {
  #   colnames_lack <- paste(colnames_lack, collapse = ", ")
  #   stop(paste("The following columns are missing:", colnames_lack))
  # }
  # if (ncol(LD_ref) != length(snp_ref)) {
  #   stop("The SNPs in reference LD matrix and its reference allele vector does't match! Please check.")
  # }
  # rownames(LD_ref) <- colnames(LD_ref) <- names(snp_ref)
  # if (is.null(sum.stat.validation)) {
  #   snps.overlap <- intersect(sum.stat.discovery$SNP, names(snp_ref))
  #   if (length(snps.overlap) == 0) {
  #     stop("There is no overlapping SNPs between summary statistics and reference LD matrix! Please check.")
  #   }
  #   rownames(sum.stat.discovery) <- sum.stat.discovery$SNP
  #   sum.stat <- sum.stat.discovery[snps.overlap, ]
  # }
  # else {
  #   snps.overlap <- intersect(intersect(sum.stat.discovery$SNP, 
  #                                       names(snp_ref)), sum.stat.validation$SNP)
  #   if (length(snps.overlap) == 0) {
  #     stop("There is no overlapping SNPs between discovery sample, validation sample and reference sample! Please check.")
  #   }
  #   rownames(sum.stat.discovery) <- sum.stat.discovery$SNP
  #   rownames(sum.stat.validation) <- sum.stat.validation$SNP
  #   sum.stat <- sum.stat.discovery[snps.overlap, ]
  #   sum.stat.valid <- sum.stat.validation[snps.overlap, 
  #                                         ]
  # }
  # LD_mat_save <- LD_ref[snps.overlap, snps.overlap]
  # LD_mat_save[lower.tri(LD_mat_save, diag = T)] <- 0
  # LD_use <- LD_mat_save + t(LD_mat_save)
  # diag(LD_use) <- 1
  # rownames(LD_use) <- colnames(LD_use) <- snps.overlap
  # snp_ref_use <- snp_ref[snps.overlap]
  # index <- sum.stat$A2 != snp_ref_use
  # tmp <- sum.stat$A1[index]
  # sum.stat$A1[index] <- sum.stat$A2[index]
  # sum.stat$A2[index] <- tmp
  # sum.stat$b[index] <- -sum.stat$b[index]
  # betas_meta <- sum.stat$b
  # betas_se <- sum.stat$se
  # n.vec <- sum.stat$N
  p <- nrow(sum.stat.discovery)
  LD_use <- cor.X
  var.X <- sum.stat.discovery$var.x
  if (standardize == T) {
    B <- LD_use
    Xy <- sum.stat.discovery$cor.x.y * sqrt(v.y)
  }
  else {
    B <- diag(sqrt(var.X)) %*% LD_use %*% diag(sqrt(var.X))
    Xy <- sum.stat.discovery$cor.x.y * sqrt(v.y) * sqrt(var.X)
  }
  lambda.v <- c()
  beta <- numeric(p)
  lambda <- max(abs(Xy))
  beta.mat <- sA.mat <- matrix(0, p, 0)
  lambda.v <- c(lambda.v, lambda)
  j1 <- which.max(abs(Xy))
  A <- c(j1)
  nA <- (1:p)[-A]
  sA <- sign(Xy[A])
  sA.v <- numeric(p)
  sA.v[A] <- sA
  sA.mat <- cbind(sA.mat, sA.v)
  XaXa_inv <- solve(B[A, A])
  beta[A] <- XaXa_inv %*% (Xy[A] - lambda * sA)
  beta.mat <- cbind(beta.mat, beta)
  XjXa <- B[nA, A]
  while (lambda > 0 & length(A) < (nvar + 1)) {
    temp <- XjXa %*% XaXa_inv
    posi1 <- (Xy[nA] - temp %*% Xy[A])/(1 - temp %*% sA)
    nega1 <- (Xy[nA] - temp %*% Xy[A])/(-1 - temp %*% sA)
    both <- cbind(posi1, nega1)
    rownames(both) <- nA
    hit <- max(both[both < lambda - 1e-10])
    sign_j <- (-1)^(which(both == hit, arr.ind = TRUE)[2] - 
                      1)
    ind <- nA[which(both == hit, arr.ind = TRUE)[1]]
    cross_all <- (XaXa_inv %*% Xy[A])/(XaXa_inv %*% sA)
    ind_cross <- which(cross_all < lambda - 1e-10)
    if (length(ind_cross) == 0) {
      cross <- -Inf
    }
    else {
      cross <- max(cross_all[ind_cross])
    }
    ind2 <- which(cross_all == cross)
    lambda <- max(hit, cross)
    if (cross < hit) {
      A <- c(A, ind)
      sA <- c(sA, sign_j)
      sA.v <- numeric(p)
      sA.v[A] <- sA
      sA.mat <- cbind(sA.mat, sA.v)
      nA <- (1:p)[-A]
    }
    else {
      beta[A[ind2]] <- 0
      A <- A[-ind2]
      sA <- sA[-ind2]
      sA.v <- numeric(p)
      sA.v[A] <- sA
      sA.mat <- cbind(sA.mat, sA.v)
      nA <- (1:p)[-A]
    }
    if (length(A) == p) {
      lambda.v <- c(lambda.v, lambda)
      XaXa_inv <- solve(B[A, A])
      beta[A] <- XaXa_inv %*% (Xy[A] - lambda * sA)
      beta.mat <- cbind(beta.mat, beta)
      break
    }
    XaXa_inv <- solve(B[A, A])
    beta[A] <- XaXa_inv %*% (Xy[A] - lambda * sA)
    beta.mat <- cbind(beta.mat, beta)
    XjXa <- B[nA, A]
    lambda.v <- c(lambda.v, lambda)
  }
  if (standardize == T) {
    beta.mat <- diag(1/sqrt(var.X)) %*% beta.mat
  }
  if (!is.null(sum.stat.validation)) {
    r2_sum <- function(beta_est, cor.x.y, var.X, var.y) {
      cov.x.y <- cor.x.y * sqrt(var.y) * sqrt(var.X)
      cov.y_hat.y <- crossprod(cov.x.y, beta_est)
      var.y_hat <- t(sqrt(var.X) * beta_est) %*% cor.X %*%
        (sqrt(var.X) * beta_est)
      return(cov.y_hat.y^2/var.y_hat/var.y)
    }
    R2 <- numeric(ncol(beta.mat))
    for (i in 2:ncol(beta.mat)) {
      R2[i] <- r2_sum(beta_est = beta.mat[, i], cor.x.y = sum.stat.validation$cor.x.y,
                      var.X = sum.stat.validation$var.x, var.y = v.y.validation)
    }
  }
  if (is.na(lambda.vec)) {
    rownames(beta.mat) <- rownames(LD_use)
    selected.markers <- rownames(beta.mat)[A]
    if (is.null(sum.stat.validation)) {
      return(list(lambda.v = lambda.v, beta.mat = Matrix(beta.mat,
                                                         sparse = TRUE), selected.markers = selected.markers[1:nvar]))
    }
    else {
      lambda.opt <- lambda.v[which.max(R2)]
      snps.opt <- which(abs(beta.mat[, which.max(R2)]) >
                          1e-10)
      beta.opt <- beta.mat[snps.opt, which.max(R2)]
      return(list(beta.opt = beta.opt, lambda.opt = lambda.opt,
                  R2 = R2, lambda.v = lambda.v, beta.mat = Matrix(beta.mat,
                                                                  sparse = TRUE), selected.markers = selected.markers[1:nvar]))
    }
  }
  lap <- function(lambda) {
    if (lambda > max(lambda.v))
      return(numeric(p))
    if (lambda < lambda.v[length(lambda.v)] || lambda ==
        lambda.v[length(lambda.v)])
      return(XaXa_inv %*% (Xy[A] - lambda * sA))
    k <- length(which(lambda < lambda.v))
    beta <- (beta.mat[, k + 1] - beta.mat[, k])/(lambda.v[k +
                                                            1] - lambda.v[k]) * (lambda - lambda.v[k]) + beta.mat[,
                                                                                                                  k]
    return(beta)
  }
  if (min(lambda.vec) < min(lambda.v))
    stop(paste("Too many variants will be selected. Please set a larger nvar or a larger lambda."))
  bm <- matrix(0, p, length(lambda.vec))
  for (i in 1:length(lambda.vec)) {
    bm[, i] <- lap(lambda.vec[i])
  }
  rownames(bm) <- rownames(LD_use)
  return(list(lambda.v = lambda.vec, beta.mat = Matrix(bm,
                                                       sparse = TRUE)))
}