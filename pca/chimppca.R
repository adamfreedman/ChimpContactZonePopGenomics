library(gdsfmt)
library(SNPRelate)
library(tidyverse)
setwd("/Users/adamfreedman/Dropbox/Africa/chimps_plosgenetics/pca")


## NEUTRAL CAPTURE ARRAY SNPS
bed.fn <-"chromfix_nosingletons_Phase_3_10k_FINAL_NEUTRAL.bed" 
fam.fn <-"chromfix_nosingletons_Phase_3_10k_FINAL_NEUTRAL.fam" 
bim.fn <-"chromfix_nosingletons_Phase_3_10k_FINAL_NEUTRAL.bim" 
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "chimp_neutral_nosingletons.gds")
snpgdsSummary("chimp_neutral_nosingletons.gds")

#The file name: /Users/adamfreedman/Dropbox/Africa/chimps_plosgenetics/pca/chimp_neutral_nosingletons.gds 
#The total number of samples: 85 
#The total number of SNPs: 2674 
#SNP genotypes are stored in SNP-major mode (Sample X SNP).

neutral_snpfile<-snpgdsOpen("chimp_neutral_nosingletons.gds")
neutral_pca <- snpgdsPCA(neutral_snpfile, num.thread=1,autosome.only=FALSE)

pc.percent <- neutral_pca$varprop*100
pc.percent[1:10]

pc.percent[1:10]
#[1] 15.802281  3.873503  2.611041  2.510184  2.360212  2.304391  2.012827  1.841469
#[9]  1.806561  1.741155


neutral_tab <- data.frame(sample.id = neutral_pca$sample.id,
                  EV1 = neutral_pca$eigenvect[,1],    # the first eigenvector
                  EV2 = neutral_pca$eigenvect[,2],    # the second eigenvector
                  EV3 = neutral_pca$eigenvect[,3],
                  EV4 = neutral_pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
write.table(neutral_tab,file="chimp_neutral_nosingletons_array_PCASnprelate.tsv",sep="\t",row.names=FALSE,col.names=TRUE)

## OUTLIER PCA

bed.fn <-"chromfix_nosingletons_Phase_3_10k_FINAL_outlier.bed" 
fam.fn <-"chromfix_nosingletons_Phase_3_10k_FINAL_outlier.fam" 
bim.fn <-"chromfix_nosingletons_Phase_3_10k_FINAL_outlier.bim" 
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "chimp_outlier_nosingletons.gds")
snpgdsSummary("chimp_outlier_nosingletons.gds")

#The file name: /Users/adamfreedman/Dropbox/Africa/chimps_plosgenetics/pca/chimp_outlier_nosingletons.gds 
#The total number of samples: 85 
#The total number of SNPs: 5185 
#SNP genotypes are stored in SNP-major mode (Sample X SNP).

outlier_snpfile<-snpgdsOpen("chimp_outlier_nosingletons.gds")
outlier_pca <- snpgdsPCA(outlier_snpfile, num.thread=1,autosome.only=FALSE)

pc.percent <- outlier_pca$varprop*100
pc.percent[1:10]

#1] 18.011394  4.144386  2.758812  2.670977  2.243600  1.931111  1.860373  1.796539  1.636720  1.524807


outlier_tab <- data.frame(sample.id = outlier_pca$sample.id,
                          EV1 = outlier_pca$eigenvect[,1],    # the first eigenvector
                          EV2 = outlier_pca$eigenvect[,2],    # the second eigenvector
                          EV3 = outlier_pca$eigenvect[,3],
                          EV4 = outlier_pca$eigenvect[,4],
                          stringsAsFactors = FALSE)
write.table(outlier_tab,file="chimp_outlier_nosingletons_array_PCASnprelate.tsv",sep="\t",row.names=FALSE,col.names=TRUE)


## PCA PLOTS
pop_idfile <- read_table("allids_withpop.txt")

### NEUTRAL
neutral_tibble <- as_tibble(neutral_tab)
neutral_tibble <- left_join(neutral_tab,pop_idfile,by=c("sample.id"))
neutral_PC1v2 <- neutral_tibble %>% ggplot(aes(x=EV1,y=EV2,color=pop)) +
                 geom_point() +
                 labs(x="PC1 (15.8%)",y="PC2 (3.9%)",color="") +
  scale_color_manual(values = c(
    "PteE" = "#009900",
    "PteR"   = "mediumorchid3",
    "Ptt"        = "orange"),
    breaks=c("PteE","PteR","Ptt"),
    labels=c("PTE: ecotone","PTE: forest","PTT")) +
  theme_classic() +
  theme(
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line      = element_blank(),   # no axis lines inside
    panel.grid     = element_blank(),   # no grid
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)) +
  theme(
    legend.position = c(0.95, 0.05),        # near bottom-right inside panel
    legend.justification = c(1, 0))

ggsave(neutral_PC1v2,file="chimp_neutral_PCA.pdf",width=4,height=4)



## OUTLIER
outlier_tibble <- as_tibble(outlier_tab)
outlier_tibble <- left_join(outlier_tab,pop_idfile,by=c("sample.id"))
outlier_PC1v2 <- outlier_tibble %>% ggplot(aes(x=EV1,y=EV2,color=pop)) +
  geom_point() +
  labs(x="PC1 (18.0%)",y="PC2 (4.1%)",color="") +
  scale_color_manual(values = c(
    "PteE" = "#009900",
    "PteR"   = "mediumorchid3",
    "Ptt"        = "orange"),
    breaks=c("PteE","PteR","Ptt"),
    labels=c("PTE: ecotone","PTE: forest","PTT")) +
  theme_classic() +
  theme(
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line      = element_blank(),   # no axis lines inside
    panel.grid     = element_blank(),   # no grid
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)) +
  theme(
    legend.position = c(0.95, 0.05),        # near bottom-right inside panel
    legend.justification = c(1, 0))

ggsave(outlier_PC1v2,file="chimp_outlier_PCA.pdf",width=4,height=4)
