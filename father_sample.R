# r2_for_sample <- function(trybd,select_var_num)
# {
# setwd("/Users/qihua/Public/BS")
# source("father_sample.R")
library(Matrix)
library(stringr)
library("data.table")
csv <- read.csv("phenotypes.both_sexes.csv", sep = ";", stringsAsFactors = FALSE)
phenofile <-readRDS("20190901_R2.rds")
source("sojo.phenotype.function.R")

csvnn <- csv[csv[,3]=="binary",]
### remove illnesses of father and sibrings
csvdata <- csvnn[!(str_detect(csvnn[,2], "Illnesses")), c(1,5,6,7,8)]
ncase <- as.numeric(csvdata[,5])
ncontrol <- as.numeric(csvdata[,4])
nnnnn1 <- as.numeric(csvdata[,2])

##### obtain mean and var of binary data
nmean <- ncase / nnnnn1
nvar <- ncontrol * ncase / (nnnnn1 * nnnnn1)
csvdata[,6] <- nmean
csvdata[,7] <- nvar
colnames(csvdata)[6] <- "mean"
colnames(csvdata)[7] <- "var"
csvb <- csvdata[,c(1,6,7)]

##### obtain mean and var of continuous data
csvc_raw <- csv[csv[,3]=="continuous_irnt",]
len <- nrow(csvc_raw)
csvc_raw[,2] <- rep(0, len)
csvc_raw[,3] <- rep(1, len)
csvc <- csvc_raw[,c(1,2,3)]
colnames(csvc) <- colnames(csvb)
csvn1 <- rbind(csvb, csvc)
csvn2 <- csvn1[csvn1[,1] != "3526_irnt", ]
csvn <- csvn2[csvn2[,1] != "4501", ]
#csvn <- csvn3[csvn3[,1] != "2139_irnt", ]

##### seek overlap between correlation matrix and phenotype table
phenotype.overlap <- intersect(csvn$phenotype, colnames(phenofile))

##### obtain correlation matrix
LD_phenomat <- phenofile[phenotype.overlap, phenotype.overlap]

##### select summary statistics concerning with and without ""
informat <-  csvn[csvn$phenotype %in% phenotype.overlap,]
yinformat <- subset(informat, phenotype == "1807_irnt", select = c(phenotype, mean, var))
xinformat <-  subset(informat, phenotype != "1807_irnt", select = c(phenotype, mean, var))

##### obtain yxcorrelation and xxcorrelation
vec <- which(colnames(LD_phenomat)== "1807_irnt")
yxcor <- LD_phenomat[colnames(LD_phenomat)== "1807_irnt", -vec]
xxcor <- LD_phenomat[-vec, -vec]


##### gain selected markers and in sample r2
df.sojo <- data.frame(var.x = xinformat$var, cor.x.y = yxcor)

res <- sojo.phenotype(sum.stat.discovery = df.sojo, cor.X = xxcor, v.y = yinformat$var, nvar = 10)
sojo.insample  <- sojo.phenotype(sum.stat.discovery =df.sojo, sum.stat.validation = df.sojo,
                               cor.X = xxcor, v.y = yinformat$var, v.y.validation = yinformat$var, nvar = 10)

##### obtain estimated beta 
vec <- which(colnames(LD_phenomat) == "1807_irnt")
cor.x.y <- LD_phenomat[colnames(LD_phenomat)=="1807_irnt", -vec]
cor.x <- LD_phenomat[-vec, -vec]
cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
colnames(cov.x) <- colnames(cor.x)
row.names(cov.x) <- row.names(cor.x)
cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)

cov.x <- cov.x[res$selected.markers,res$selected.markers]
cov.x.y <- cov.x.y[res$selected.markers]
beta <- solve(cov.x) %*% cov.x.y

##### input non White British data
tsv_data <- fread("/Users/qihua/Public/PHESANT/WAS/results/output..tsv", header=TRUE, sep='\t', data.table=FALSE)
pheno_data <- tsv_data[,-c(1,2,3)]
judge_set <- c()
for (i in 1:length(colnames(pheno_data))) {
  if(class(pheno_data[,i]) == "logical"){
    judge_set[i] = TRUE
  }else{
    judge_set[i] = FALSE
  }}
binary_set <- pheno_data[,which(judge_set)]
binary_set1 <- ifelse(binary_set == TRUE, 1, 0)
continuous_set <- pheno_data[,-which(judge_set)]
colnames(continuous_set) <- paste(colnames(continuous_set), "_irnt", sep = "")
total_set <- cbind(binary_set1, continuous_set)
n_white <- total_set[which(total_set[,"21000_1001"] == 0),]
n_whi <- n_white[ ,-which(colnames(n_white) == "22006")]
### dim(whi) = 57160 1755
vari <- c("1807_irnt", res$selected.markers)

## out_r2 
  lm_R2 <- c()
  select_data <- as.matrix(n_whi[, vari])
  na_fix <- na.omit(select_data)
    p <- ncol(na_fix)
    q <- nrow(na_fix)

  score <- sweep(na_fix[,-1], 2, beta, "*")
  for (i in 1:length(vari[-1])) {
    lm_summary <- summary(lm(na_fix[,1]~ score[,1:i] ))
    lm_R2[i] <- lm_summary$r.squared
  }
  t <- list(p-1,q,lm_R2[length(lm_R2)])
  return(t)

for (i in 1:10) {
  m <- pheno_log(i)
  var_num <- c()
  sam_size <- c()
  r2 <- c()
  var_num[i] <- m[[1]]
  sam_size <- m[[2]] 
  r2 <- m[[3]]
}
s <- list(var_num, sam_size, r2)

# pheno_log <- function(num){
#   res <- sojo.phenotype(sum.stat.discovery = df.sojo, cor.X = xxcor, v.y = yinformat$var, nvar = num)
# cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
# colnames(cov.x) <- colnames(cor.x)
# row.names(cov.x) <- row.names(cor.x)
# cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
# cov.x <- cov.x[res$selected.markers,res$selected.markers]
# cov.x.y <- cov.x.y[res$selected.markers]
# beta <- solve(cov.x) %*% cov.x.y
#   cov.x <- cov.x[res$selected.markers,res$selected.markers]
#   cov.x.y <- cov.x.y[res$selected.markers]
#   beta <- solve(cov.x) %*% cov.x.y
#   vari <- c("1807_irnt", res$selected.markers)
#   ## out_r2 
#   lm_R2 <- c()
#   select_data <- as.matrix(n_whi[, vari])
#   na_fix <- na.omit(select_data)
#   p <- ncol(na_fix)
#   q <- nrow(na_fix)
#   score <- sweep(na_fix[,-1], 2, beta, "*")
#   for (i in 1:length(vari[-1])) {
#     lm_summary <- summary(lm(na_fix[,1]~ score[,1:i] ))
#     lm_R2[i] <- lm_summary$r.squared
#   }
#   t <- list(p-1,q,lm_R2[length(lm_R2)])
#   return(t)
# }
# for (i in 1:10) {
#   m <- pheno_log(i)
#   var_num <- c()
#   sam_size <- c()
#   r2 <- c()
#   var_num[i] <- m[[1]]
#   sam_size <- m[[2]] 
#   r2 <- m[[3]]
# }
# s <- list(var_num, sam_size, r2)
# 


# ### in sample r2
# in_r2 <- function(select_var, xinformat, yinformat, cor.x.y, cor.x){
#   xinformat <- xinformat[match(select_var, xinformat$phenotype),]
#   cor.x <- cor.x[select_var,select_var]
#   cor.x.y <- cor.x.y[select_var]
#   if(length(select_var) == 1){
#     cov.x <- sqrt(xinformat$var) * cor.x * sqrt(xinformat$var)
#     cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
#     ebeta <- cov.x * cov.x.y
#     cov.y.hy2 <- cov.x.y * ebeta
#     var.hy2 <- sqrt(xinformat$var * ebeta) * cor.x * (sqrt(xinformat$var) * ebeta)
#     cor.y.hy2 <- cov.y.hy2^2/var.hy2/yinformat$var
#     return(cor.y.hy2)
#   }
#   else{
#   cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
#   cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
#   ebeta <- solve(cov.x) %*% cov.x.y
#   cov.y.hy2 <- cov.x.y %*% ebeta
#   var.hy2 <- t(sqrt(xinformat$var) * ebeta) %*% cor.x %*% (sqrt(xinformat$var) * ebeta)
#   cor.y.hy2 <- cov.y.hy2^2/var.hy2/yinformat$var
#   return(cor.y.hy2)
#   }
# }
# insample_R2 <- c()
# for (i in 1:length(select_var)) {
#   insample_R2[i] <- in_r2(select_var[1:i], xinformat, yinformat, cor.x.y, cor.x)
# }
#   
# ## out_r2 function
# out_r2 <- function(data, select_var, ebeta){
#   lm_R2 <- c()
#   pheno_names <- c()
#   select_var1 <- str_extract(select_var, "[0-9]+")
#   select_data <- as.matrix(data[, select_var1])
#   ebeta <- ebeta[select_var,]
#   na_fix <- cbind(pheno1807,select_data)
#   na_fix1 <- na.omit(na_fix)
#   score <- sweep(na_fix1[,-1], 2,ebeta, "*")
#   # anova(lm(na_fix1[,1]~., as.data.frame(score )))
#   # step(lm(na_fix1[,1]~., as.data.frame(score )))
#   for (i in 1:length(select_var)) {
#     lm_summary <- summary(lm(na_fix1[,1]~ score[,1:i] ))
#     lm_R2[i] <- lm_summary$r.squared
#   }
#   return(lm_R2)
# }
#  outsample_r2 <- out_r2(bd_select, select_var, ebeta)
#  
# info <- cbind( select_var, insample_R2, outsample_r2)
# colnames(info) <- c("var", "insample_r2", "outsample_r2")

# ### gain coeff for training and validation set data in a total way
# select_var1 <- str_extract(select_var, "[0-9]+")
# select_data <- as.matrix(bd_select[, select_var1])
# sum_coeff_list <- summary(lm(pheno1807~select_data))
# sum_coeff <- unlist(sum_coeff_list$coefficients)
# n <- length(select_var)
# cov.x1 <-cov.x[select_var,select_var]
# cov.x.y1 <- cov.x.y[select_var]
# ebeta1 <- solve(cov.x1) %*% cov.x.y1
# vt_coeff <- cbind(sum_coeff[2:(n+1)],ebeta1)
# colnames(vt_coeff) <- c("valid_coeff", "training_coeff")

# ### gain r for training and validation set data
# train_cor <- cor.x.y[select_var]
# valid_cor <- c()
# n <- length(select_var)
# for (i in 2:(n+1)) {
#   na_fix <- cbind(pheno1807,select4[,i])
#   na_fix1 <- na.omit(na_fix)
#   valid_cor[i-1] <- cor(na_fix1[,1],na_fix1[,2])
# }
# vt_r <- cbind(valid_cor,train_cor)
# colnames(vt_r) <- c("valid_r", "training_r")
# 
# total_info <- cbind(info,vt_coeff,vt_r)


# na_fix <- cbind(pheno1807,select_data)
# na_fix1 <- na.omit(na_fix)
# na_fix1 <-as.data.frame(na_fix1)
# summary(lm1 <- lm(na_fix1[,1]~.,na_fix1[,-1]) )
# slm1 <- step(lm1)
# summary(slm1)
# slm1$anova





# # ### integrate 6138 data
# without_edu <- csvn3[!(str_detect(csvn3[,1],"6138.[0-9]+")),]
# edu <- csvn3[(str_detect(csvn3[,1],"6138.[0-9]+")),]
# inte_edu <- apply(edu[,2:3], 2, mean)
# inte_edu <- c("6138",inte_edu)
# csvn <- rbind.data.frame(without_edu,inte_edu)


# install.packages("car")
# library(car)
# scatterplotMatrix(xxcor,spread=FALSE,lty.smooth=2,main="Scatter Plot Matrix")
   
# r_value <-  summary(t)
# anova(t)
# return(r_value$r.squared)
# return(t)

# r_sum1 <- function(pheno_data,beta){
#   na_fix <- cbind(pheno1807,pheno_data)
#   na_fix1 <- na.omit(na_fix)
#   score <- na_fix1[,-1] * matrix(rep(beta,nrow(na_fix1)), ncol=2 )
#  
#   
#   lm(na_fix1[,1]~ score)
#   # # r_value <-  summary(t)
#   # # anova(t)
#   # # return(r_value$r.squared)
# }
# a1 <- r_sum1(pheno_data = as.matrix(select4[,1]), beta = as.matrix(ebeta[1]))
# a2 <- r_sum1(pheno_data = as.matrix(select4[,1:2]), beta = as.matrix(ebeta[1:2]))
# b <- anova(a2, a1, test="LRT")


# b$`Pr(>Chi)`
# c <- c()
#  for (i in 1:7) {
#    a1 <- r_sum(pheno_data = as.matrix(select4[,1:i]), beta = as.matrix(ebeta[1:i]))
#  a2 <- r_sum(pheno_data = as.matrix(select4[,1:(i+1)]), beta = as.matrix(ebeta[1:(i+1)]))
#  # b <- anova(a2, a1, test="LRT")
#  # c[i] <- b$`Pr(>Chisq)`
#  
#  }




# for (i in 2 : 5) {
#   lm_R2[i] <- r_sum(pheno_data = as.matrix(select4[,1:i]), beta = as.matrix(ebeta[1:i]))
#   pheno_names[i] <- colnames(select4[i])
# }
# ebeta1 <- ebeta[select3,]
# # xymat1 <- xymat[select3]
# # insample_r2 <- c()
# # for (i in 1:7) {
# #   insample_r2[i] <- (t(ebeta1[1:i]) %*% xymat1[1:i]) / (yinformat1$var + (yinformat1$mean)^(1/2))
# # }
# r2_info <- data.frame(pheno_names, lm_R2)
# r2_info

# ### gain r2 for training and validation set data
# 
# valid_r2 <- c()
# for (i in 2:5) {
#   a <- summary(lm(pheno1807~ as.matrix(select4[,2:i])))
#   valid_r2[i] <- a$r.squared
# }


# # vt_r2 <- cbind(valid_r2, insample_r2)
# # colnames(vt_r2) <- c("valid_r2", "training_r2")
# valid_info <- data.frame(pheno_names, lm_R2, valid_r2)


# for (i in 1:7) {
#   a <- summary(lm(pheno1807~ as.matrix(select4[,1:i])))
#   valid_r2[i] <- a$r.squared
# }
# for (i in 1:7) {
#   a <- summary(valid_sum[i])
#   valid_r2[i] <- a$r.squared
# }
# cbind(insample_r2,valid_r2)

# select4 <- c()
# for (i in 5) {
#   max.loca <- which.max(lm_R2)
#   select4[i] <- max.loca
#   r2_info <- select3[,-max.loca]
# }


# sele_var <- select3[,select4]
# r2_info <-  cbind(colnames(sele_var),insample_r2,lm_R2)

# 
# for (i in 2:5) {
#   na_fix2 <-cbind(sele_var[,1:i],pheno1807)
#     na_fix3 <- na.omit(na_fix2)
#     score1 <- na_fix3[,1:i] %*% ebeta[1:i,]
#     r_value1 <-  summary(lm(na_fix3[,i+1]~ score1))
#     lm_R21[i] <- r_value1$r.squared
# 
# }


# for (i in 2:5) {
#   na_fix2 <-cbind(select3[,1:i],pheno1807)
#   na_fix3 <- na.omit(na_fix2)
#   score1 <- na_fix3[,1:i] %*% ebeta[1:i,]
#   r_value1 <-  summary(lm(na_fix3[,i+1]~ score1))
#   lm_R2[i] <- r_value1$r.squared
# }
# for (i in 6:9) {
#   na_fix2 <-cbind(select3[,5:i],pheno1807)
#   na_fix3 <- na.omit(na_fix2)
#   score1 <- na_fix3[,5:i] %*% ebeta[5:i,]
#   r_value1 <-  summary(lm(na_fix3[,i+1]~ score1))
#   lm_R2[i] <- r_value1$r.squared
# }
# for (i in 10:13) {
#   na_fix2 <-cbind(select3[,9:i],pheno1807)
#   na_fix3 <- na.omit(na_fix2)
#   score1 <- na_fix3[,9:i] %*% ebeta[9:i,]
#   r_value1 <-  summary(lm(na_fix3[,i+1]~ score1))
#   lm_R2[i] <- r_value1$r.squared
# }
# 



# ese <- sqrt(abs(matrix(rep(yinformat$var, 6), ncol = 1)- xinformat$var * ebeta^2)/(diag(xxmat) * n_non_missing))



# stad_corxy <- LD_phenomat1[colnames(LD_phenomat1)=="1807_irnt", -vec1]
# stad_corxx <- LD_phenomat1[-vec1, -vec1]
# xxcor1 <- diag(sqrt(xinformat1$var)) %*% stad_corxx %*% diag(sqrt(xinformat1$var))
# yxcor1 <- stad_corxy * sqrt(yinformat1$var) * sqrt(xinformat1$var)
# yxcov <- yxcor1 * sqrt(yinformat1$var) * sqrt(xinformat1$var)
# ycov_hat.y <- crossprod(yxcov,ebeta)
# yvar_hat <- t(sqrt(xinformat1$var *ebeta)) %*% xxcor1 %*% (sqrt(xinformat1$var) * ebeta)
# r2 <- ycov_hat.y^2/yvar_hat/yinformat1$var
# for cycle R2[i] <- r2



# ## using sojo to estimate in sample r2
# modiorisojo <- data.frame(var.x =  xinformat1$var, cor.x.y = yxcor1)
# sojo.insample  <- sojo.phenotype(sum.stat.discovery =modiorisojo, sum.stat.validation = modiorisojo,
#                                cor.X = xxcor1, v.y = yinformat1$var, v.y.validation = yinformat1$var, nvar = 11)
# 
# sojo.insample_r2 <- sojo.insample$R2[-1]
# sojo.insample_r2



# ### sojo information for training set
# sojo_form <- data.frame( xinformat1$var[select4], yxcor1[select4])
# 
# ### sojo information for validation set
# na_fix <-cbind(select3[,select4],pheno1807)
# na_fix1 <- na.omit(na_fix)
# 
# ## default 4 variables
# pheno_x_naomit <- na_fix1[,1:4]
# phone_y_naomit <- na_fix1[,5]
# x_var <- apply(pheno_x_naomit,2,var)
# 
# myfun <- function(x){
#  x <-  cor(phone_y_naomit,x)
# }
# corxy <- apply(pheno_x_naomit, 2, myfun)
# valid_form <- data.frame(x_var, corxy)
# colnames(sojo_form) <- c("var.x", "cor.x.y")
# colnames(valid_form) <- c("var.x", "cor.x.y")
# 
# corxx_selected <- xxcor1[select4,select4]
# 
# res.outsample1 <- sojo.phenotype(sum.stat.discovery = sojo_form, sum.stat.validation = valid_form,
#                                  cor.X = corxx_selected, v.y = yinformat1$var, v.y.validation =yinformat1$var, nvar = 4)
# 
# res.insample1 <- sojo.phenotype(sum.stat.discovery = sojo_form, sum.stat.validation = sojo_form,
#                                cor.X = corxx_selected, v.y = yinformat$var, v.y.validation = yinformat$var, nvar = 4)
# 
# a <- cbind(res.outsample1$selected.markers, res.outsample1$R2,res.insample1$R2)
# colnames(a) <- c("pheno", "outsample","insample")



#cbind(stand_xxcor[1:5,1:5],xxcor[1:5,1:5])
# res$selected.markers
# cbind(stand_yxcor[1:5],yxcor[1:5])

# select2 <- c(which(colnames(phenofile) == "1807_irnt"), which(colnames(phenofile) == "6145_3"),
#              which(colnames(phenofile) == "6138_1"),which(colnames(phenofile) == "6150_100"),
#             which(colnames(phenofile) == "2139_irnt"),which(colnames(phenofile) == "23099_irnt"),
#              which(colnames(phenofile) == "1787"),which(colnames(phenofile) == "2754_irnt"),
#              which(colnames(phenofile) == "6177_100"),which(colnames(phenofile) == "20002_1473"))



############################################################################################
###################################################################################
#########################################################################
# csvn1 <- rbind(csvb, csvc)                                                               
#  csvn2 <- csvn1[csvn1[,1] != "3526_irnt", ]
# # csvn2 <- csvn1[csvn1[,1] != "2754_irnt", ]
#  csvn <- csvn2[csvn2[,1] != "4501", ]

# ##### obtain mean and var of ordinal data with 3 group
# csvd_raw <- csv[csv[,3]=="ordinal", c(1,9)]
# ### obtain value matrix for each group
# csvd1_raw <- str_extract(csvd_raw[,2] ,"cat.N:.[0-9,]+.[0-9,]+.[0-9,]+")
# csvd1_na_out <- na.omit(csvd1_raw)
# csvd1_advance <- str_extract(csvd1_na_out, "([0-9]+),.([0-9]+),.([0-9]+)")
# csvd1 <- unlist(strsplit(csvd1_advance, split = ","))
# csvd1_matrx <- matrix(csvd1, ncol = 3, byrow = TRUE)
# csvd1_matrx1 <- apply(csvd1_matrx,2,as.numeric)
# ### obtain number matrix for each group
# na_select <- which(!is.na(csvd1_raw))
# csvd2_raw <- str_extract(csvd_raw[na_select,2], "order:.[0-9|]+")
# csvd2_advance <- str_extract(csvd2_raw, "[0-9]+.[0-9]+.[0-9]+")
# csvd2 <- unlist(strsplit(csvd2_advance, split = "[|]"))
# csvd2_matrx <- matrix(csvd2, ncol = 3, byrow = TRUE)
# csvd2_matrx1 <- apply(csvd2_matrx,2,as.numeric)
# ### calculate mean and variance
# case_sum <- apply(csvd1_matrx1, 1, sum)
# mean_matrix <- csvd1_matrx1 * csvd2_matrx1 / case_sum
# mean_sum <- apply(mean_matrix, 1, sum)
# var_matrix <- csvd2_matrx1 * (case_sum - csvd1_matrx1) * csvd1_matrx1 / (case_sum * case_sum)
# var_sum <- apply(var_matrix, 1, sum)
# csvd <- csvd_raw[na_select,]  
# csvd[,2] <- mean_sum
# csvd[,3] <- var_sum
# colnames(csvd) <- c("phenotype", "mean","var")
# 
# ##### obtain mean and var of ordinal data with more than 3 group
# csve_raw <- csv[csv[,3]=="ordinal", c(1,9)]
# head(csve_raw)
# csve1_raw <- str_extract(csvd_raw[,2] ,"[0-9]+\\([0-9]+\\).+[0-9]+\\([0-9]+\\)")
# 
# csve1_na_out <- na.omit(csve1_raw)
# csve1_advance1 <- str_extract_all(csve1_na_out, "[0-9]+\\(")
# csve1_advance11 <-  str_extract_all(csve1_advance1, "[0-9]+")
# 
# csve1_advance2 <- str_extract_all(csve1_na_out, "\\([0-9]+\\)")
# csve1_advance22 <-  str_extract_all(csve1_advance2, "[0-9]+")
# 
# mean_sum1 <- NA
# var_sum1 <-NA
# for (i in 1:length(csve1_advance1)) {
#   num_indi <- as.numeric(unlist(csve1_advance22[[i]]))
#   
#   value_indi <- as.numeric(unlist(csve1_advance11[[i]]))
#   mean_indi <- num_indi * value_indi / sum(num_indi)
#   mean_sum1[i] <- sum(mean_indi)
#   var_indi <- num_indi * value_indi* (sum(num_indi)-num_indi)/ (sum(num_indi)*sum(num_indi))
#   var_sum1[i] <-sum(var_indi)
# }
# csve <- csve_raw[which(!is.na(csve1_raw)),]
# csve[,2] <- mean_sum1
# csve[,3] <- var_sum1
# colnames(csve) <- c("phenotype", "mean","var")

#  csvn <- csvn2[csvn2[,1] != "4501", ]

# csvn4 <- csvn3[csvn3[,1] != "6145_3", ]
# csvn5 <- csvn4[csvn4[,1] != "6138_1", ]
# csvn6 <- csvn5[csvn5[,1] != "6138_100", ]
# csvn7 <- csvn6[csvn6[,1] != "2139_irnt", ]
# csvn8 <- csvn7[csvn7[,1] != "23099_irnt", ]
# csvn <- csvn1[csvn1[,1] != "1787", ]
# csvn <- csvn9[csvn9[,1] != "6150_100", ]
#############################################################################
#####################################################################################
#############################################################################################



#####################################################################################################
###########################################################################################
##################################################################################
### validate selected markers
# ### select the overlap part
# modified_res <-  str_extract(res$selected.markers, "[0-9]+")
# bdoverlap <- which( colnames(trybd.0.0) %in% c(modified_res,"1807"))
# bdselect1 <- trybd.0.0[ , bdoverlap]
# 
# ### remove NA value, wait for correct
# bdselect2 <- na.omit(bdselect1)
# 
# 
# ### prepare summary data
# bdmean <- apply(bdselect2, 2, mean)
# bdvar <- apply(bdselect2, 2, var)
# bdinfo <- cbind(bdmean, bdvar)
# colnames(bdinfo) <- c( "mean", "var" )
# 
# ### gain validate info xxcor, yxcor and df.sojo
# vec_valid <- which(colnames(bdselect2)== "1807")
# total_cor <- cor(bdselect2, bdselect2)
# yxcor_valid <- total_cor[vec_valid, -vec_valid]
# xxcor_valid <- total_cor[-vec_valid, -vec_valid]
# 
# 
# df.sojo_valid <- data.frame(var.x = bdvar[-vec_valid], cor.x.y = yxcor_valid)
# 
# originoverlap <- intersect( modified_res, colnames(trybd.0.0))
# orisojo <- df.sojo[res$selected.markers,]
# modiorisojo <- orisojo[-4,]

# rownames(modiorisojo) <- str_extract(row.names(modiorisojo), "[0-9]+")
# 
# orixxcor <- xxcor[res$selected.markers, res$selected.markers]
# modiorixxcor <- orixxcor[-4,-4]

# res.outsample <- sojo.phenotype(sum.stat.discovery = modiorisojo, sum.stat.validation = df.sojo_valid,
#                                 cor.X = modiorixxcor, v.y = yinformat$var, v.y.validation = yinformat$var, nvar = 20)
# res.insample <- sojo.phenotype(sum.stat.discovery = modiorisojo, sum.stat.validation = modiorisojo,
#                                cor.X = modiorixxcor, v.y = yinformat$var, v.y.validation = yinformat$var, nvar = 9)
# matrix(c(res.insample$R2, res.outsample$R2), nrow = length(res.insample$R2), byrow = FALSE, dimnames = list(originoverlap,NULL))
#############################################################################
#####################################################################################
###############################################################################################
