#!/usr/bin/Rscript

require(ggplot2)
require(ggsignif)

conservation_score_Ep <- read.table("Epromoters_conservation-score-phyloP_sum_merged.bed", header = F)
conservation_score_cp <- read.table("control_promoters_conservation-score-phyloP_sum_merged.bed", header = F)


#=================conservation-score by phyloP========================

conservation_score_df1 <- data.frame(conservation_score = conservation_score_Ep$V4,
                                     type = "Epromoter")

conservation_score_df2 <- data.frame(conservation_score = conservation_score_cp$V4,
                                     type = "control_promoter")

conservation_score_df <- rbind(conservation_score_df1, conservation_score_df2)

conservation_score_df$type <- factor(conservation_score_df$type, levels = unique(conservation_score_df$type))

pdf("241_Placental_conservation_score_phyloP_Ep_cp.pdf", width = 4.5, height = 4)

ggplot(conservation_score_df, aes(x = type, y = conservation_score, fill = type))+
  geom_violin(width = 0.9)+
  geom_boxplot(width = 0.2, outlier.shape = NA)+
  scale_y_continuous(breaks = c(-1000,0,1000,2000,3000,4000))+
  scale_x_discrete(labels = c("Epromoter", "control"))+
  scale_fill_manual(values = c("#377EB8", "gray"))+
  geom_hline(yintercept = median(conservation_score_df[conservation_score_df$type == "control_promoter", 1]), linetype = 2, color = "gray")+
  labs(title = "conservation score\n241 Zoonomia Placental Mammals", y = "conservation score per promoter")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Epromoter", "control_promoter")), y_position = 5000,
              textsize = 5, linetype = "blank")


dev.off()






