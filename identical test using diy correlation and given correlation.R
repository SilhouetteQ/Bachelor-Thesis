select_pheno <- "1807_irnt"
res <- sojo.phenotype(sum.stat.discovery = df.sojo, cor.X = xxcor, v.y = yinformat$var, nvar = 8)
predictors <- res$selected.markers
cor.x.y <- LD_phenomat[select_pheno, predictors]
cor.x <- LD_phenomat[predictors, predictors]
xinformat <-  subset(informat, phenotype != "1807_irnt", select = c(phenotype, mean, var))
row.names(xinformat) <- xinformat$phenotype
xinformat <-  xinformat[predictors, ]

cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
be <- solve(cov.x) %*% cov.x.y
lm_R2 <- c()


va <-  c("1807_irnt", res$selected.markers)
b <- be
select_data <- as.matrix(valid[, va])
na_fix <- na.omit(select_data)
p <- ncol(na_fix)
q <- nrow(na_fix)
score <- sweep(na_fix[,-1], 2, b, "*")

  lm_summary <- summary(lm(na_fix[,1]~ score[,1:8] ))
  lm_R2 <- lm_summary$r.squared
num8 <- c(b[8,1],p,q,lm_R2)
  b1 <- b
  b2 <- be

  
  
  select_pheno <- "1807_irnt"
  csv <-read.csv("G:/COURSERA/GITHUB/Prediction/44_diy.csv", sep= "", stringsAsFactors = FALSE)
  select_data <- as.matrix(indi[, c(select_pheno, res$selected.markers)])
  R2 <- cor(select_data,use = 'pairwise.complete.obs')
  cor.x.y <- phenofile[select_pheno, res$selected.markers]
  cor.x <- phenofile[res$selected.markers, res$selected.markers]
  xinformat <-  csv
  row.names(xinformat) <- xinformat$phenotype
  xinformat <-  xinformat[res$selected.markers, ]
  cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
  cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
  be <- solve(cov.x) %*% cov.x.y
  na_fix <- na.omit(select_data)
  p <- ncol(na_fix)
  q <- nrow(na_fix)
  lm_R2 <- c()
  score <- sweep(na_fix[,-1], 2, be, "*")
  lm_summary <- summary(lm(na_fix[,1]~ score[,1:8] ))
  lm_R2 <- lm_summary$r.squared
  num88 <- c(be[8,1],p,q,lm_R2)
  