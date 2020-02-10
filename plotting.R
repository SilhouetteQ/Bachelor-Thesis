### Painting
setwd("G:/COURSERA/GITHUB/BS")
library(ggplot2)
library("data.table")
R1 <- read.csv("outsample_R2_1.csv", sep= "", stringsAsFactors = FALSE) 
R2 <- read.csv("outsample_R2_2.csv", sep= "", stringsAsFactors = FALSE) 
R3 <- read.csv("outsample_R2_3.csv", sep= "", stringsAsFactors = FALSE) 
R4 <- read.csv("outsample_R2_mother_2.csv", sep= "", stringsAsFactors = FALSE) 
R5 <- read.csv("outsample_R2_mother_1.csv", sep= "", stringsAsFactors = FALSE)

info <- fread("phenotypes.both_sexes.tsv", header=TRUE, sep='\t', data.table=FALSE)
row.names(info) <- info$phenotype

### 1. father top description
se <- info[R1$pheno,]
de <- se$description
de[1] <- "Death of a close relative"
de[3] <- "Do not have any education qualification"
de[4] <- "No vascular/heart problems diagnosed by doctor"
de[9] <- "No medication for cholesterol, blood pressure or diabetes"
write.csv(deta, "R1_with_description.csv",  row.names = FALSE)

### 2. mother top description
ee <- info[R4$pheno,]
de1 <- ee$description
de1[1] <- "Death of a close relative"
de1[3] <- "Do not have any education qualification"
de1[4] <- "No vascular/heart problems diagnosed by doctor"
de1[10] <- "No medication for cholesterol, blood pressure or diabetes"
de1[15] <- "is.female"
ee$description <- de1
description <- paste(letters[1:16],".",ee$description)
n <- ggplot(R4)
n+ geom_col(aes(var_num, sam_size, fill = description))+labs(x= "The number of selected variants",y = "Sample size") + scale_x_continuous(labels = letters[1:16],breaks = 1:16)
n+ geom_point(aes(var_num,r2,color = description),size = 4)+geom_smooth(aes(var_num, r2),colour = "blue",size = 0.5,method = 'loess') + labs(x= "The number of selected variants",y = "Out-of-sample R2")+ scale_y_log10()+scale_x_continuous(labels = letters[1:16],breaks = 1:16)

### 3. bar and line plot for sample size and r2 separately
# deta <- read.csv( "R1_with_description.csv", sep= "", stringsAsFactors = FALSE)
description <- paste(letters[1:10],".", de)
n <- ggplot(R1)
n+ geom_col(aes(var_num, sam_size, fill = description))+labs(x= "The number of selected variants",y = "Sample size") + scale_x_continuous(labels = letters[1:10],breaks = 1:10)
n+ geom_point(aes(var_num,r2,color = description),size = 4)+geom_smooth(aes(var_num, r2),colour = "blue",size = 0.5,method = 'loess') + labs(x= "The number of selected variants",y = "Out-of-sample R2")+ scale_x_continuous(labels = letters[1:10],breaks = 1:10)


### 4. compare sample size(or r2) between white and non-white
RR2 <- data.frame(R2, rep("Non White British",12))
colnames(RR2)[6] <- "Ancestors"
RR3 <- data.frame(R3, rep("White British",12))
colnames(RR3)[6] <- "Ancestors"
si_com <- rbind(RR2,RR3)
k <- ggplot(si_com)
k+geom_line(aes(x=var_num,y = sam_size, color = Ancestors),size = 2) + labs(x ="The number of selected variants",y = "Sample size") + scale_y_log10() + scale_x_continuous(labels =c(1:10,"sex","age"),breaks = 1:12)
k+geom_line(aes(x=var_num,y = r2, color = Ancestors),size = 2) + labs(x ="The number of selected variants",y = "Out-of-sample R2") + scale_y_log10() +scale_x_continuous(labels =c(1:10,"sex","age"),breaks = 1:12) 


### 5. compare r2 estimated by LD score and direct phenotypic correlation
q <- cbind.data.frame(R1$r2, R2$r2[-c(11,12)] )
colnames(q) <- c("LD","pheno")
w <- ggplot(q)
w + geom_smooth(aes(LD, pheno), method = lm)
q <- cbind.data.frame(R1$r2, R2$r2[-c(11,12)] )
colnames(q) <- c("LD","pheno")
w <- ggplot(q)
w + geom_smooth(aes(LD, pheno), method = lm)


### 6. 
com[1:6, 6] <- "father"
com[7:12, 6] <- "mother"
k <- ggplot(com[1:6,])
com[7:12,]

fa_de <-  paste(letters[1:6],".", de[c(1:5,8)])
k+ geom_point(aes(var_num,r2,color = fa_de),size = 4) + labs(x= "The number of selected variants",y = "Out-of-sample R2")+ scale_x_discrete(labels = letters[1:6],breaks = 1:6)

library(plotly)
p<-plot_ly(y=com$r2[1:6], x=com$var_num[1:6] , type="scatter", mode="markers+lines", name = "father")
p<-add_trace(p, y=com$r2[7:12], x=com$var_num[1:6] , type="scatter", mode="markers+lines", name = "mother")
