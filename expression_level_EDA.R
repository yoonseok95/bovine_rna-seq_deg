
################################################################################
#
# Course 1. Gene Expression visualization
#
################################################################################
# Step 1.1 Data load from RINE2020 results.
################################################################################


getwd()
setwd("E:/이지연 DEG 분석/")

library(dplyr)
data = read.csv("rpkm_results.csv", row.names=1)

colnames(data) # To check the column names of the data
View(data)

################################################################################
# Step 1.2 Extracting gene expression information from RPKM results.
################################################################################

# Extract method 1.

rpkm_table = data
dim(rpkm_table3)
head(rpkm_table)
summary(rpkm_table)
# Change column names(R에서 컬럼명 바꿀 것. csv파일 자체에서 수정하면 에러남)
colnames(rpkm_table)


###1
colnames(rpkm_table) <- c("H_2nd_12_9", "H_2nd_16_9", "H_2nd_2_9",  "H_2nd_8_9", 
                          "H_3rd_14_9", "H_3rd_16_8", "H_3rd_17_9", "H_4th_14_9",
                          "H_5th_10_9", "H_5th_9_9",  "L_2nd_18_4", "L_2nd_19_4",
                          "L_3rd_20_1", "L_3rd_4_2",  "L_3rd_6_3",  "L_3rd_7_3" ,
                          "L_3rd_8_2",  "L_4th_15_3", "L_4th_17_1", "L_5th_19_1",
                          "H_RNAlater_6_7",  "H_RNAlater_7_7",  "L_RNAlater_8_4",  "L_RNAlater_9_5")
colnames(rpkm_table)
View(rpkm_table)

###2
rpkm_table2 = subset(data, select=c("hw_H_2nd_12_htseq.count", "hw_H_2nd_16_htseq.count", "hw_H_2nd_2_htseq.count",  "hw_H_2nd_8_htseq.count", 
                                    "hw_H_3rd_14_htseq.count", "hw_H_3rd_16_htseq.count", "hw_H_3rd_17_htseq.count", "hw_H_4th_14_htseq.count",
                                    "hw_H_5th_10_htseq.count", "hw_H_5th_9_htseq.count",  "hw_L_2nd_18_htseq.count", "hw_L_2nd_19_htseq.count",
                                    "hw_L_3rd_20_htseq.count", "hw_L_3rd_4_htseq.count",  "hw_L_3rd_6_htseq.count",  "hw_L_3rd_7_htseq.count", 
                                    "hw_L_3rd_8_htseq.count",  "hw_L_4th_15_htseq.count", "hw_L_4th_17_htseq.count", "hw_L_5th_19_htseq.count"))
colnames(rpkm_table2) <- c("H_2nd_12_9", "H_2nd_16_9", "H_2nd_2_9",  "H_2nd_8_9", 
                          "H_3rd_14_9", "H_3rd_16_8", "H_3rd_17_9", "H_4th_14_9",
                          "H_5th_10_9", "H_5th_9_9",  "L_2nd_18_4", "L_2nd_19_4",
                          "L_3rd_20_1", "L_3rd_4_2",  "L_3rd_6_3",  "L_3rd_7_3" ,
                          "L_3rd_8_2",  "L_4th_15_3", "L_4th_17_1", "L_5th_19_1")
colnames(rpkm_table2)  

###3
rpkm_table3 = subset(data, select=c("hw_H_2nd_12_htseq.count", "hw_H_2nd_16_htseq.count", "hw_H_2nd_2_htseq.count",  "hw_H_2nd_8_htseq.count", 
                                                                        "hw_L_2nd_18_htseq.count", "hw_L_2nd_19_htseq.count",
                                    "hw_L_3rd_20_htseq.count", "hw_L_3rd_4_htseq.count",  "hw_L_3rd_6_htseq.count",  "hw_L_3rd_7_htseq.count", 
                                    "hw_L_3rd_8_htseq.count",  "hw_L_4th_15_htseq.count", "hw_L_4th_17_htseq.count", "hw_L_5th_19_htseq.count",
                                    "RNAlater_6_htseq.count",  "RNAlater_7_htseq.count"))

colnames(rpkm_table3) <- c("H_2nd_12_9", "H_2nd_16_9", "H_2nd_2_9",  "H_2nd_8_9", 
                                                      "L_2nd_18_4", "L_2nd_19_4",
                           "L_3rd_20_1", "L_3rd_4_2",  "L_3rd_6_3",  "L_3rd_7_3" ,
                           "L_3rd_8_2",  "L_4th_15_3", "L_4th_17_1", "L_5th_19_1","H_RNAlater_6_7",  "H_RNAlater_7_7")
colnames(rpkm_table3)  





################################################################################
# Step 1.3 Histogram
################################################################################

# When using hist(), if log(0, 2) results in "-Inf", it is excluded from the histogram.
# What number will you use to draw a histogram?
rpkm_table$H_2nd_12_9 # Checking the RPKM values
log(0.0, 2)              # Inf/-Inf (infinity) 
log(rpkm_table$H_2nd_12_9, 2)
log(rpkm_table$H_2nd_12_9+1, 2)

# Let's check the distribution of gene expression in one sample through historam.
par(mfrow=c(6,4)) # Setting the layout of the histograms

png("유전자발현량데이터EDA_histogramv1.png", width=10, height=8, units="in", res=600)

par(mfrow=c(6,4), mar=c(2, 2, 2, 2))
hist(log(rpkm_table$H_2nd_12_9, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_2nd_12_9")
hist(log(rpkm_table$H_2nd_16_9, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_2nd_16_9")
hist(log(rpkm_table$H_2nd_2_9, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_2nd_2_9")
hist(log(rpkm_table$H_2nd_8_9, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_2nd_8_9")
hist(log(rpkm_table$H_3rd_14_9, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_3rd_14_9")
hist(log(rpkm_table$H_3rd_16_8, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_3rd_16_8")
hist(log(rpkm_table$H_3rd_17_9, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_3rd_17_9")
hist(log(rpkm_table$H_4th_14_9, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_4th_14_9")
hist(log(rpkm_table$H_5th_10_9, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_5th_10_9")
hist(log(rpkm_table$H_5th_9_9, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_5th_9_9")
hist(log(rpkm_table$L_2nd_18_4, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_2nd_18_4")
hist(log(rpkm_table$L_2nd_19_4, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_2nd_19_4")
hist(log(rpkm_table$L_3rd_20_1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_3rd_20_1")
hist(log(rpkm_table$L_3rd_4_2, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_3rd_4_2")
hist(log(rpkm_table$L_3rd_6_3, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_3rd_6_3")
hist(log(rpkm_table$L_3rd_7_3, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_3rd_7_3")
hist(log(rpkm_table$L_3rd_8_2, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_3rd_8_2")
hist(log(rpkm_table$L_4th_15_3, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_4th_15_3")
hist(log(rpkm_table$L_4th_17_1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_4th_17_1")
hist(log(rpkm_table$L_5th_19_1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_5th_19_1")
hist(log(rpkm_table$H_RNAlater_6_7, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_RNAlater_6_7")
hist(log(rpkm_table$H_RNAlater_7_7, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_RNAlater_7_7")
hist(log(rpkm_table$L_RNAlater_8_4, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_RNAlater_8_4")
hist(log(rpkm_table$L_RNAlater_9_5, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_RNAlater_9_5")

dev.off()

png("유전자발현량데이터EDA_histogramv2.png", width=10, height=8, units="in", res=600)
par(mfrow=c(6,4), mar=c(2, 2, 2, 2))
hist(log(rpkm_table$H_2nd_12_9+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_2nd_12_9")
hist(log(rpkm_table$H_2nd_16_9+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_2nd_16_9")
hist(log(rpkm_table$H_2nd_2_9+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_2nd_2_9")
hist(log(rpkm_table$H_2nd_8_9+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_2nd_8_9")
hist(log(rpkm_table$H_3rd_14_9+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_3rd_14_9")
hist(log(rpkm_table$H_3rd_16_8+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_3rd_16_8")
hist(log(rpkm_table$H_3rd_17_9+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_3rd_17_9")
hist(log(rpkm_table$H_4th_14_9+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_4th_14_9")
hist(log(rpkm_table$H_5th_10_9+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_5th_10_9")
hist(log(rpkm_table$H_5th_9_9+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_5th_9_9")
hist(log(rpkm_table$L_2nd_18_4+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_2nd_18_4")
hist(log(rpkm_table$L_2nd_19_4+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_2nd_19_4")
hist(log(rpkm_table$L_3rd_20_1+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_3rd_20_1")
hist(log(rpkm_table$L_3rd_4_2+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_3rd_4_2")
hist(log(rpkm_table$L_3rd_6_3+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_3rd_6_3")
hist(log(rpkm_table$L_3rd_7_3+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_3rd_7_3")
hist(log(rpkm_table$L_3rd_8_2+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_3rd_8_2")
hist(log(rpkm_table$L_4th_15_3+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_4th_15_3")
hist(log(rpkm_table$L_4th_17_1+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_4th_17_1")
hist(log(rpkm_table$L_5th_19_1+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_5th_19_1")
hist(log(rpkm_table$H_RNAlater_6_7+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_RNAlater_6_7")
hist(log(rpkm_table$H_RNAlater_7_7+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "H_RNAlater_7_7")
hist(log(rpkm_table$L_RNAlater_8_4+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_RNAlater_8_4")
hist(log(rpkm_table$L_RNAlater_9_5+1, 2), breaks = 10, xlim = c(-25, 25), probability = T, main = "L_RNAlater_9_5")
dev.off()




par(mfrow=c(1,1)) # Resetting the layout to the default value of 1 row and 1 column.

################################################################################
# Step 1.4 Scatter plot
################################################################################

plot(rpkm_table$H_2nd_12_9, rpkm_table$L_3rd_20_1, col="black", pch =19)

# Let's examine the gene expression distribution between control and case samples 
# for each IMF score(H vs L) using scatter plots.

png("H_2nd_12_9 vs H_3rd_16_8.png", width=6, height=5, units="in", res=600)
plot(log(rpkm_table$H_2nd_12_9+1, 2), log(rpkm_table$H_3rd_16_8+1, 2), 
     col="royalblue", pch =19, main = "H_2nd_12_9 versus H_3rd_16_8")
dev.off()

#for문 사용하여 모든 조합으로 그려보기
# 'H' 그룹과 'L' 그룹의 샘플 이름 추출
H_samples <- colnames(rpkm_table2)[grep("^H_", colnames(rpkm_table))]
L_samples <- colnames(rpkm_table2)[grep("^L_", colnames(rpkm_table))]

# 모든 H 샘플과 L 샘플의 조합에 대해 산점도 그리기
for (H_sample in H_samples) {
  for (L_sample in L_samples) {
    file_name <- paste0(H_sample, "_vs_", L_sample, ".png")
    png(file_name, width=6, height=5, units="in", res=600)
    
    # 산점도 그리기
    plot(log(rpkm_table2[[H_sample]]+1, 2), log(rpkm_table2[[L_sample]]+1, 2), 
         col="royalblue", pch =19, main = paste(H_sample, "versus", L_sample))
    
    dev.off()
  }
}

#H그룹 내 샘플 비교
# 'H' 그룹의 샘플 이름 추출
H_samples <- colnames(rpkm_table)[grep("^H_", colnames(rpkm_table))]

# 모든 H 샘플 조합에 대해 산점도 그리기
for (i in 1:(length(H_samples)-1)) {
  for (j in (i+1):length(H_samples)) {
    # 파일 이름 생성
    file_name <- paste0(H_samples[i], "_vs_", H_samples[j], ".png")
    
    # 그래프 저장 시작
    png(file_name, width=6, height=5, units="in", res=600)
    
    # 산점도 그리기
    plot(log(rpkm_table[[H_samples[i]]]+1, 2), log(rpkm_table[[H_samples[j]]]+1, 2), 
         col="royalblue", pch =19, main = paste(H_samples[i], "versus", H_samples[j]))
    
    # 그래프 저장 종료
    dev.off()
  }
}




# 'H' 그룹과 'L' 그룹의 데이터 그룹으로 묶어서 진행
H_group <- rpkm_table[, grep("^H_", colnames(rpkm_table))]
L_group <- rpkm_table[, grep("^L_", colnames(rpkm_table))]
H_group_mean <- rowMeans(log(H_group+1, 2))
L_group_mean <- rowMeans(log(L_group+1, 2))

# 산점도 그리기
png("H vs L_1.png", width=6, height=5, units="in", res=600)
plot(H_group_mean, L_group_mean, col="royalblue", pch =19, main = "High IMF group versus Low IMF group(Anseong+Yeongju)")
dev.off()



H_group <- rpkm_table3[, grep("^H_", colnames(rpkm_table3))]
L_group <- rpkm_table3[, grep("^L_", colnames(rpkm_table3))]
H_group_mean <- rowMeans(log(H_group+1, 2))
L_group_mean <- rowMeans(log(L_group+1, 2))

# 산점도 그리기
png("H vs L_3.png", width=6, height=5, units="in", res=600)
plot(H_group_mean, L_group_mean, col="royalblue", pch =19, main = "High IMF group versus Low IMF group(filterinig)")
dev.off()


# Correlation plot with scatter plot !

png("correlation plot_1.png", width=12, height=12, units="in", res=600)
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

upper.panel<-function(x, y){
  points(x,y, pch = 19)
}

pairs(log(rpkm_table+1, 2), 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()

##Heatmap
library(ggplot2)
library(reshape2)

# 상관 계수 행렬을 계산
cor_matrix <- cor(rpkm_table)

# melt 함수를 사용하여 데이터를 재구조화
melted_cor_matrix <- melt(cor_matrix)


png("correlation plot_2.png", width=12, height=12, units="in", res=600)

ggplot(data = melted_cor_matrix, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "royalblue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  coord_fixed()

dev.off()


################################################################################
# Step 1.5 Box plot
################################################################################

# Simple way
boxplot(log(rpkm_table, 2), col = "aquamarine")
boxplot(log(rpkm_table+1, 2), col = "pink")

png("boxplot_rpkm0.3.png", width=12, height=8, units="in", res=600)
# Select genes with an RPKM value of 0.3 or greater in at least one sample.
rpkm_filt_table <- rpkm_table[which(apply(rpkm_table, 1, max) >= 0.3),]
dim(rpkm_table)
dim(rpkm_filt_table)
head(rpkm_filt_table)
summary(rpkm_filt_table)
png("boxplot_rpkm0.3.png", width=12, height=8, units="in", res=600)
par(mar = c(10, 10, 4, 2) + 0.1)
boxplot(log(rpkm_filt_table+1, 2), col = "pink", las=2)
dev.off()



par(mar = c(2, 2, 2, 2) + 0.1)

################################################################################
# Step 1.6 MDS
################################################################################

# 2d MDS plot
t_rpkm_table <- t(rpkm_table)

group <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,1,1,2,2)
d <- dist(log(t_rpkm_table+0.1,10), method="euclidean")
fit <- cmdscale(d, eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

png("H, L MDS plot.png", width=12, height=8, units="in", res=600)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS", type="n", 
     xlim = c(-10,10), ylim = c(-10,10))
text(x, y, labels=row.names(t_rpkm_table), cex=.9, col=group)
grid()
dev.off()

# 3d MDS plot
install.packages("rgl")
library(rgl)
fit <- cmdscale(d,eig=TRUE,k=3)
x <- fit$points[,1]
y <- fit$points[,2]
z <- fit$points[,3]
plot3d(x,y,z,xlab="Coordinate 1", ylab="Coordinate 2", main="MDS", type="n")
text3d(x,y,z,texts=row.names(t_rpkm_table), cex=.7, col=group)

