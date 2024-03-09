#https://www.genekitr.fun/plot-ora-1

install.packages("tidyverse")
library(tidyverse)

getwd()
setwd("E:/이지연 DEG 분석/")

#DAVID 데이터 불러오기
david <- read.table("GO_filtering7.txt", header = TRUE, sep = "\t") 
#GO ID와 Term 분리
david <- david %>%
  separate(Term, into = c("GO_ID", "Term"), sep = "~")
colnames(david)
head(david)

# 데이터를 'Count' 값이 큰 순서로 정렬
#sorted_bp_david <- bp_david[order(bp_david$Count, decreasing = TRUE),]

# 필요한 라이브러리 로드
library(ggplot2)
library(ggprism)

# 'Category' 열에서 'BP', 'MF', 'CC'를 선택하고 'Benjamini' 값이 0.05 미만인 것만 선택
filtered_david <- david[david$Category %in% c("GOTERM_BP_DIRECT", "GOTERM_MF_DIRECT", "GOTERM_CC_DIRECT") & david$Benjamini < 0.05,]

# 'Category'와 'Count'에 따라 정렬
# factor를 이용하여 'Category' 열의 순서를 지정
filtered_david$Category <- factor(filtered_david$Category, levels = c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT"))

# 'Category'의 이름을 변경
levels(filtered_david$Category) <- c("Biological process", "Cellular component", "Molecular function")

sorted_filtered_david <- filtered_david[order(filtered_david$Category, -filtered_david$Count),]

# ggplot2를 사용하여 barplot 생성
GO_plot <- ggplot(sorted_filtered_david, aes(x = reorder(Term, Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Term", y = "Gene count", fill = "Category") +
  scale_fill_manual(values = c("Biological process" = "#41a2ab", "Cellular component" = "#48b874", "Molecular function" = "#f59145")) +
  theme_prism(base_size = 10) + ggtitle('GO Terms')+
  theme(legend.title = element_text(), strip.text = element_blank()) + 
  facet_grid(Category ~ ., scales = "free_y", space = "free_y")

GO_plot

ggsave("GO_filtering7.1.png", dpi=600, dev='png', height=5, width=7, units="in")







#############KEGG pathway

library(patchwork)
library(tidyverse)
library(ggplot2)
library(ggprism)
#데이터 불러오기 및 다듬기
david2 <- read.table("KEGG_filtering9.txt", header = TRUE, sep = "\t") 
david2 <- david2 %>%
  separate(Term, into = c("KEGG_ID", "Term"), sep = ":")
colnames(david2)
head(david2)


# dplyr 패키지 로드
library(dplyr)

# 데이터를 Fold Enrichment 값에 따라 오름차순으로 정렬
david2 <- david2 %>% arrange(Fold.Enrichment)

# 그래프 그리기
ggplot(david2, aes(x = reorder(Term, Fold.Enrichment), y = Fold.Enrichment)) + 
  geom_point(aes(size = Count, color = PValue)) +
  theme_bw(base_size = 14) + coord_flip() + ylim(0, 6)+
  labs(x = "Term", y = "Fold Enrichment", size = "Count", color = "p-value") +
  scale_colour_gradient(limits=c(0, 0.1), low="#de3131", high = "blue") +
  theme_prism(base_size = 13)+
  ggtitle("KEGG pathway enrichment") +
  theme(plot.title = element_text(size = 18), legend.title = element_text())

# 저장
ggsave("KEGG_filtering9.png", dpi=600, dev='png', height=8, width=10, units="in")







##########top 20개 plot 그리기
# 'Fold.Enrichment'에 따라 정렬하고 상위 20개만 선택
sorted_david2 <- david2[order(-david2$Fold.Enrichment),]
top_david2 <- sorted_david2[1:20,]

#dot plot 그리기
ggplot(top_david2, aes(x = Term, y = Fold.Enrichment)) + 
  geom_point(aes(size = Count, color = PValue)) +
  theme_bw(base_size = 14) + coord_flip() + ylim(0, 6)+
  labs(x = "Term", y = "Fold Enrichment", size = "Count", color = "p-value") +
  scale_colour_gradient(limits=c(0, 0.1), low="#de3131", high = "blue") +
  theme_prism(base_size = 13)+
  ggtitle("Top of 20 pathway enrichment") +
  theme(plot.title = element_text(size = 18), legend.title = element_text())

#저장
ggsave("KEGG_filtering7 top20.png", dpi=600, dev='png', height=7, width=10, units="in")


