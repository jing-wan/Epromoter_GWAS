#!/usr/bin/Rscript

require(ggplot2)
require(ggsignif)

Ep_RNAPIIpeaks_num <- read.table("Epromoters_ov_RNAPII_all_peaks_count.bed", header = F)
cp_RNAPIIpeaks_num <- read.table("control_ov_RNAPII_all_peaks_count.bed", header = F)

#also consider non overlapped promoters
RNAPII_peaks_df1 <- data.frame(RNAPII_peaks_per_promoter = c(Ep_RNAPIIpeaks_num$V4, rep(0, 5743-nrow(Ep_RNAPIIpeaks_num))),
                               type = "Epromoter")
RNAPII_peaks_df2 <- data.frame(RNAPII_peaks_per_promoter = c(cp_RNAPIIpeaks_num$V4, rep(0, 5743-nrow(Ep_RNAPIIpeaks_num))),
                               type = "control")

RNAPII_peaks_df <- rbind(RNAPII_peaks_df1, RNAPII_peaks_df2)
RNAPII_peaks_df$type <- factor(RNAPII_peaks_df$type, levels = c("Epromoter", "control"))


pdf("RNAPII_peaks_all_datasets_per_promoter.pdf", width = 4.5, height = 4)

ggplot(RNAPII_peaks_df, aes(x = type, y = RNAPII_peaks_per_promoter, fill = type))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  scale_y_continuous(breaks = c(0,424,655,1000,1500))+
  scale_x_discrete(labels = c("Epromoter", "control"))+
  scale_fill_manual(values = c("#377EB8", "gray"))+
  geom_hline(yintercept = median(RNAPII_peaks_df[RNAPII_peaks_df$type == "control", 1]), linetype = 2, color = "gray")+
  # labs(title = "RNAPII peaks density (all datasets)", y = "RNAPII peaks per promoter")+
  labs(title = "RNAPII enrichment", y = "RNAPII peaks per promoter")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Epromoter", "control")), y_position = 1500,
              textsize = 4, linetype = "blank")

dev.off()



