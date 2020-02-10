


####读variant�?
gz <- gzfile('/Volumes/ShaoYa/UKBB/BenNeale.ukb.round2/variants.tsv.bgz') ###20190730
gz_1 <- read.table(gz,header=T)
idx <- which(substr(gz_1$rsid, 1, 2) == 'rs')
gz_1 <- gz_1[idx,]
saveRDS(gz_1,'/Users/ting/pleiotropy/20190901_variant_131.rds')

####删选low MAF的SNP ，highMAF是我用来干别的事情的
highMAF_idx_5e3 <-  which(gz_1$minor_AF>5e-3)
saveRDS(highMAF_idx_5e3,'/Users/ting/pleiotropy/20190901_highMAF_idx_for_131.rds')
lowMAF_idx_5e4 <-  which(gz_1$minor_AF<5e-4)
saveRDS(lowMAF_idx_5e4,'/Users/ting/pleiotropy/20190901_lowMAF_idx_for_131.rds')
gz_1$variant->snp_name
saveRDS(snp_name,'/Users/ting/pleiotropy/20190901_right_snp_name_for_131.rds')

####因为全部太长�? 我把SNP分成 18 �? 最后再拼起�? 可以不用�?

K <- 13187546 
starts <- seq(1, K, 740000)
ends <- c(starts[-1] - 1, K)

for (core in 1:18) {
  cat('core')
  hh <- which(starts[core]<= highMAF_idx_5e3 & highMAF_idx_5e3 <=ends[core])
  highMAF_idx_in_core <- highMAF_idx_5e3[hh]
  saveRDS(highMAF_idx_in_core,paste0('/Users/ting/pleiotropy/20190901_highMAF_idx_in_core_need_mimus740000.core.',core,'.rds'))
  
  ll <- which(starts[core]<= lowMAF_idx_5e4 & lowMAF_idx_5e4 <=ends[core])
  lowMAF_idx_in_core <- lowMAF_idx_5e4[ll]
  saveRDS(lowMAF_idx_in_core,paste0('/Users/ting/pleiotropy/20190901_lowMAF_idx_in_core_need_mimus740000.core.',core,'.rds'))
}



core <- scan(n=1)
K <- 13187546 
starts <- seq(1, K, 740000)
ends <- c(starts[-1] - 1, K)


zmat <- readRDS(paste0('/Users/xia/20190806_74w_zmat.', core, '.rds'))

#pheno_1374 <- readRDS('/Users/ting/pleiotropy/20190811_1374.rds')
pheno_1511 <- readRDS('/Users/ting/pleiotropy/pheno_1511.rds')
zmat[,pheno_1511]->zmat

ll <- readRDS(paste0('/Users/ting/pleiotropy/20190901_lowMAF_idx_in_core_need_mimus740000.core.',core,'.rds'))
ll <- ll-(core-1)*740000
zmat <- zmat[ll,]
saveRDS(zmat,paste0('/Users/ting/pleiotropy/20190901_lowMAF_Rmatrix.core.',core,'.rds'))

####以上得到18�?1511列的summary stat

#
####合并18�? summary stat
core <- 1
R1 <- readRDS(paste0('/Users/ting/pleiotropy/20190901_lowMAF_Rmatrix.core.',core,'.rds'))
for (core in 2:18) {
  R <- readRDS(paste0('/Users/ting/pleiotropy/20190901_lowMAF_Rmatrix.core.',core,'.rds'))
  R1 <- rbind(R1,R)
}
saveRDS(R1,'/Users/ting/pleiotropy/20190901_lowMAF_Rmatrix.merged.rds')



####这里R1 �? 1511列summary stat matrix

####后来我直接加上sex 和age 作为新的2�? 但是那次code我没有单独存下来 只需�? age �? sex�? summary stat，把lowMAF的SNP的挑出来即可

R2 <- cor(R1,use = 'pairwise.complete.obs')


###我把correlation>0.9的都删掉 trait �?1511 变成1376�?
for (i in 1: 1699) {
  which(abs(R2[i,])>0.9 &abs(R2[i,]) <1) ->idx
  if(length(idx)>0){
    R2[-idx,-idx]->R2}
}

saveRDS(R2,'/Users/ting/pleiotropy/20190901_R2.rds')
