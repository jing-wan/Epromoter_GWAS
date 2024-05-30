#!/usr/bin/Rscript

G4_Ep_CGI <- read.table("./input_data_from_Carla/covG4H1_eProm_atCGIs.gff.txt", header = F)
G4_Ep_nonCGI <- read.table("./input_data_from_Carla/covG4H1_eProm_atOtherSites.gff.txt", header = F)
G4_cp_CGI <- read.table("./input_data_from_Carla/covG4H1_prom_atCGIs.gff.txt", header = F)
G4_cp_nonCGI <- read.table("./input_data_from_Carla/covG4H1_prom_atotherSites.gff.txt", header = F)

G4_coverage_df1 <-
  rbind(data.frame(G4_coverage = G4_Ep_CGI$V12, type = "Epromoter"),
        data.frame(G4_coverage = G4_cp_CGI$V12, type = "control"))

G4_coverage_df2 <-
  rbind(data.frame(G4_coverage = G4_Ep_nonCGI$V12, type = "Epromoter"),
        data.frame(G4_coverage = G4_cp_nonCGI$V12, type = "control"))

G4_coverage_df1$type <- factor(G4_coverage_df1$type, levels = c("Epromoter", "control"))
G4_coverage_df2$type <- factor(G4_coverage_df2$type, levels = c("Epromoter", "control"))


pdf("G4 coverage in promoters with CpG island densityplot.pdf", width = 5, height = 4)

ggplot(G4_coverage_df1, aes(x = G4_coverage, fill = type))+
  geom_density(alpha = 0.8)+
  scale_x_continuous(limits = c(0,100), expand = c(0,0))+
  scale_fill_manual(values = c("#1F78B4", "gray"))+
  labs(title = "promoters with CpG island", 
       x = "G4 coverage %", y = "density")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        # axis.title.x = element_blank(),
        axis.title = element_text(size=16),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.margin = unit(c(0.5,1,0.5,0.5), "cm"))

dev.off()


pdf("G4 coverage in promoters without CpG island densityplot.pdf", width = 5, height = 4)

ggplot(G4_coverage_df2, aes(x = G4_coverage, fill = type))+
  geom_density(alpha = 0.8)+
  scale_x_continuous(limits = c(0,100), expand = c(0,0))+
  scale_fill_manual(values = c("#1F78B4", "gray"))+
  labs(title = "promoters without CpG island", 
       x = "G4 coverage %", y = "density")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        # axis.title.x = element_blank(),
        axis.title = element_text(size=16),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.margin = unit(c(0.5,1,0.5,0.5), "cm"))

dev.off()




#=======================number of G4H1 distribution===========================

setwd("C:/work/sync/工作/PhD project/Epromoter_variants.v5.1/12.Genomic_features/Epromoter_G4/")

G4_Ep_CGI <- read.table("./input_data_from_Carla/nbG4H1_ePromoters_atCGIs.gff", header = F)
G4_Ep_nonCGI <- read.table("./input_data_from_Carla/nbG4H1_ePromoters_atOtherSites.gff", header = F)
G4_cp_CGI <- read.table("./input_data_from_Carla/nbG4H1_promoters_atCGIs.gff", header = F)
G4_cp_nonCGI <- read.table("./input_data_from_Carla/nbG4H1_promoters_atOtherSites.gff", header = F)


G4_df1 <-
  rbind(data.frame(G4 = G4_Ep_CGI$V10, type = "Epromoter"),
        data.frame(G4 = G4_cp_CGI$V10, type = "control"))

G4_df2 <-
  rbind(data.frame(G4 = G4_Ep_nonCGI$V10, type = "Epromoter"),
        data.frame(G4 = G4_cp_nonCGI$V10, type = "control"))

G4_df1$type <- factor(G4_df1$type, levels = c("Epromoter", "control"))
G4_df2$type <- factor(G4_df2$type, levels = c("Epromoter", "control"))


pdf("G4 numbers in promoters with CpG island densityplot.pdf", width = 5, height = 4)

ggplot(G4_df1, aes(x = G4, fill = type))+
  geom_density(alpha = 0.8)+
  scale_x_continuous(limits = c(0,15), expand = c(0,0))+
  scale_fill_manual(values = c("#1F78B4", "gray"))+
  labs(title = "with CpG island", 
       x = "G4 numbers per promoter", y = "density")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size=18, color = "black"),
        # axis.title.x = element_blank(),
        axis.title = element_text(size=18),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        plot.margin = unit(c(0.2,0.4,0.2,0.2), "cm"))

dev.off()


pdf("G4 numbers in promoters without CpG island densityplot.pdf", width = 5, height = 4)

ggplot(G4_df2, aes(x = G4, fill = type))+
  geom_density(alpha = 0.8)+
  scale_x_continuous(limits = c(0,15), expand = c(0,0))+
  scale_fill_manual(values = c("#1F78B4", "gray"))+
  labs(title = "without CpG island", 
       x = "G4 numbers per promoter", y = "density")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size=18, color = "black"),
        # axis.title.x = element_blank(),
        axis.title = element_text(size=18),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        plot.margin = unit(c(0.2,0.4,0.2,0.2), "cm"))

dev.off()


# ks.test(G4_Ep_CGI$V10, G4_cp_CGI$V10)
# p-value = 0.01592
# ks.test(G4_Ep_nonCGI$V10, G4_cp_nonCGI$V10)
# p-value = 4.996e-14


#=======================G4 hunting score distribution===========================


setwd("C:/work/sync/工作/PhD project/Epromoter_variants.v5.1/12.Genomic_features/Epromoter_G4/")

G4_Ep_CGI <- read.table("./input_data_from_Carla/G4Hscore_eProm_CGI_column13.gff", header = F)
G4_Ep_nonCGI <- read.table("./input_data_from_Carla/G4Hscore_eProm_otherSites_column13.gff", header = F)
G4_cp_CGI <- read.table("./input_data_from_Carla/G4Hscore_prom_atCGIs_column13.gff", header = F)
G4_cp_nonCGI <- read.table("./input_data_from_Carla/G4Hscore_prom_atOtherSites_column13.gff", header = F)


G4_df1 <-
  rbind(data.frame(G4 = G4_Ep_CGI$V13, type = "Epromoter"),
        data.frame(G4 = G4_cp_CGI$V13, type = "control"))

G4_df2 <-
  rbind(data.frame(G4 = G4_Ep_nonCGI$V13, type = "Epromoter"),
        data.frame(G4 = G4_cp_nonCGI$V13, type = "control"))

G4_df1$G4 <- as.numeric(G4_df1$G4)
G4_df2$G4 <- as.numeric(G4_df2$G4)
G4_df1$type <- factor(G4_df1$type, levels = c("Epromoter", "control"))
G4_df2$type <- factor(G4_df2$type, levels = c("Epromoter", "control"))


pdf("G4 hunting score in promoters with CpG island densityplot.pdf", width = 5, height = 4)

ggplot(G4_df1, aes(x = G4, fill = type))+
  geom_density(alpha = 0.8)+
  scale_x_continuous(limits = c(1,4), expand = c(0,0))+
  scale_fill_manual(values = c("#1F78B4", "gray"))+
  labs(title = "with CpG island", 
       x = "G4 hunting score", y = "density")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size=18, color = "black"),
        # axis.title.x = element_blank(),
        axis.title = element_text(size=18),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        plot.margin = unit(c(0.2,0.4,0.2,0.2), "cm"))

dev.off()


pdf("G4 hunting score in promoters without CpG island densityplot.pdf", width = 5, height = 4)

ggplot(G4_df2, aes(x = G4, fill = type))+
  geom_density(alpha = 0.8)+
  scale_x_continuous(limits = c(1,4), expand = c(0,0))+
  scale_fill_manual(values = c("#1F78B4", "gray"))+
  labs(title = "without CpG island", 
       x = "G4 hunting score", y = "density")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size=18, color = "black"),
        # axis.title.x = element_blank(),
        axis.title = element_text(size=18),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        plot.margin = unit(c(0.2,0.4,0.2,0.2), "cm"))

dev.off()



# ks.test(as.numeric(G4_Ep_CGI$V13), as.numeric(G4_cp_CGI$V13))
# p-value = 1.782e-10
# ks.test(as.numeric(G4_Ep_nonCGI$V13), as.numeric(G4_cp_nonCGI$V13))
# p-value = 1.086e-06

