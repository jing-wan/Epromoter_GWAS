
#=========================BRD2========================================

require(ggplot2)
require(ggsignif)


Ep_TFpeaks_num <- read.table("Epromoters_ov_BRD2_peaks_count.bed", header = F)
cp_TFpeaks_num <- read.table("control_ov_BRD2_peaks_count.bed", header = F)

#also consider non overlapped promoters
TF_peaks_df1 <- data.frame(TF_peaks_per_promoter = c(Ep_TFpeaks_num$V4, rep(0, 5743-nrow(Ep_TFpeaks_num))),
                               type = "Epromoter")
TF_peaks_df2 <- data.frame(TF_peaks_per_promoter = c(cp_TFpeaks_num$V4, rep(0, 5743-nrow(Ep_TFpeaks_num))),
                               type = "control")

TF_peaks_df <- rbind(TF_peaks_df1, TF_peaks_df2)
TF_peaks_df$type <- factor(TF_peaks_df$type, levels = c("Epromoter", "control"))


pdf("BRD2_peaks_per_promoter.pdf", width = 4.5, height = 4)

ggplot(TF_peaks_df, aes(x = type, y = TF_peaks_per_promoter, fill = type))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  scale_y_continuous(breaks = c(0,25,50,75,100))+
  scale_x_discrete(labels = c("Epromoter", "control"))+
  scale_fill_manual(values = c("#377EB8", "gray"))+
  geom_hline(yintercept = median(TF_peaks_df[TF_peaks_df$type == "control", 1]), linetype = 2, color = "gray")+
  labs(title = "BRD2 peaks density (57 datasets)", y = "BRD2 peaks per promoter")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Epromoter", "control")), y_position = 90,
              textsize = 4, linetype = "blank")

dev.off()




#=========================P300========================================


Ep_TFpeaks_num <- read.table("Epromoters_ov_P300_peaks_count.bed", header = F)
cp_TFpeaks_num <- read.table("control_ov_P300_peaks_count.bed", header = F)

#also consider non overlapped promoters
TF_peaks_df1 <- data.frame(TF_peaks_per_promoter = c(Ep_TFpeaks_num$V4, rep(0, 5743-nrow(Ep_TFpeaks_num))),
                           type = "Epromoter")
TF_peaks_df2 <- data.frame(TF_peaks_per_promoter = c(cp_TFpeaks_num$V4, rep(0, 5743-nrow(Ep_TFpeaks_num))),
                           type = "control")

TF_peaks_df <- rbind(TF_peaks_df1, TF_peaks_df2)
TF_peaks_df$type <- factor(TF_peaks_df$type, levels = c("Epromoter", "control"))


pdf("P300_peaks_per_promoter.pdf", width = 4.5, height = 4)

ggplot(TF_peaks_df, aes(x = type, y = TF_peaks_per_promoter, fill = type))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  scale_y_continuous(breaks = c(0,10,20,30,40,50))+
  scale_x_discrete(labels = c("Epromoter", "control"))+
  scale_fill_manual(values = c("#377EB8", "gray"))+
  geom_hline(yintercept = median(TF_peaks_df[TF_peaks_df$type == "control", 1]), linetype = 2, color = "gray")+
  labs(title = "P300 peaks density (43 datasets)", y = "P300 peaks per promoter")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Epromoter", "control")), y_position = 45,
              textsize = 4, linetype = "blank")

dev.off()





#=========================SMARCA4========================================


Ep_TFpeaks_num <- read.table("Epromoters_ov_SMARCA4_peaks_count.bed", header = F)
cp_TFpeaks_num <- read.table("control_ov_SMARCA4_peaks_count.bed", header = F)

#also consider non overlapped promoters
TF_peaks_df1 <- data.frame(TF_peaks_per_promoter = c(Ep_TFpeaks_num$V4, rep(0, 5743-nrow(Ep_TFpeaks_num))),
                           type = "Epromoter")
TF_peaks_df2 <- data.frame(TF_peaks_per_promoter = c(cp_TFpeaks_num$V4, rep(0, 5743-nrow(Ep_TFpeaks_num))),
                           type = "control")

TF_peaks_df <- rbind(TF_peaks_df1, TF_peaks_df2)
TF_peaks_df$type <- factor(TF_peaks_df$type, levels = c("Epromoter", "control"))


pdf("SMARCA4_peaks_per_promoter.pdf", width = 4.5, height = 4)

ggplot(TF_peaks_df, aes(x = type, y = TF_peaks_per_promoter, fill = type))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  scale_y_continuous(breaks = c(0,25,50,75,100))+
  scale_x_discrete(labels = c("Epromoter", "control"))+
  scale_fill_manual(values = c("#377EB8", "gray"))+
  geom_hline(yintercept = median(TF_peaks_df[TF_peaks_df$type == "control", 1]), linetype = 2, color = "gray")+
  labs(title = "SMARCA4 peaks density (84 datasets)", y = "SMARCA4 peaks per promoter")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Epromoter", "control")), y_position = 100,
              textsize = 4, linetype = "blank")

dev.off()





