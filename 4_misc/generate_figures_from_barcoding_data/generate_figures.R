renv::init()

library(grid, quietly = TRUE)
library(Rtsne, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(gplots, quietly = TRUE)
library(egg)
library(scales)
library(dplyr)
library(data.table, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)
library(ComplexHeatmap)
library(viridis)
library(colorRamps)
library(circlize)
library(ggsci)
library(UpSetR)
library(ggpubr)
library(doBy)
library(gtable)
library(stringr)
library(reshape2)

# Retrieve all the data needed for plotting
source("retrieve_data_functions.R")
water_fall_v3 <- return_waterfall_clontracer()
cw_summary <- return_cw_summary_clontracer()
ts_sorted_v2 <- return_ts_sorted_clontracer()
cpm_ctg_dh <- return_cpm_ctg_dh()
log2fc_combo_treat <- return_log2fc_combo_treat()
log2fc_combo_treat_cc <- return_log2fc_combo_treat_cc()
clonid_traj <- return_traj_classes_javi()
counts_summary_drug_combo <- return_counts_summary_drug_combo()
sensitivities <- return_drug_sensitivities()
ltr_data <- return_ltr_data()
ltr_rep_avg <- return_ltr_rep_avg()
ltr_rep_overlap <- return_ltr_rep_overlap()
cpm_dna <- return_cpm_dna()
cpm_rna_cast_annt <- return_cpm_rna()
log2fc_rna_dna <- return_log2fc_rna_dna()
chemo_rep_overlap <- return_chemo_rep_overlap()
hcc4006_cell_counts <- return_hcc4006_cellcounts()

# Cleanup the data
water_fall_v3[timepoint == "D6", timepoint := "6 days untreated"]
water_fall_v3[timepoint == "W2", timepoint := "2 weeks EGFRi"]
water_fall_v3[timepoint == "W6", timepoint := "6 weeks EGFRi"]
water_fall_v3$condition <- factor(water_fall_v3$timepoint, c("6 days untreated","2 weeks EGFRi","6 weeks EGFRi"))



cw_summary$timepoint_replicate = paste(cw_summary$Timepoint, cw_summary$Replicate, sep="_")


# Helper functions
convert_to_cpm <- function(x){log2((x*10**6) + 1)}


# Figure S1B
# ***************************************************************************************************************************************************
water_fall_v3 %>% filter(cum_fraction <= 0.9) %>% ggplot(aes(x=cell_line, y=log2fc, fill=condition)) + 
  geom_boxplot() + 
  scale_x_discrete(limit = c("MGH707", "HCC4006", "PC9", "HCC827")) +
  theme_bw(base_size = 25) +
  theme(axis.title.y = element_text(face = "italic", size = 30), axis.title.x = element_text(face = "italic", size = 30)) +
  scale_fill_manual(breaks = c("6 days untreated","2 weeks EGFRi","6 weeks EGFRi"), values=c("#808080", "#ffb3b3", "#cc0000")) + 
  labs(x= "cell line", y = "log2fc vs 0 day untreated")
ggsave("plots/S1/S1B.pdf", width = 10, height = 7)
# ***************************************************************************************************************************************************


# Figure 1B
# ***************************************************************************************************************************************************
plot_corr<-function(cellline){
  mat_sorted_v2 = ts_sorted_v2 %>% 
    filter(cum_frac <= 0.9 & Cellline==cellline) %>% 
    as.data.table %>% 
    dcast(Barcode ~ group_replicate, value.var="Fraction", fun.aggregate = mean, fill = 0)
  
  cormat <- round(cor(mat_sorted_v2 %>% select_if(is.numeric)) ,2)
  cormat <- cormat^2
  melted_cormat <- reshape2::melt(cormat)
  
  tile_plot = ggplot(data = melted_cormat, aes(x=Var1, y=reorder(Var2, desc(Var2)), fill=value)) + 
    geom_tile(color = "white")+
    scale_fill_viridis_c(limit = c(0,1), space = "Lab", name=expression("r"^2)) +
    theme_minimal(base_size = 20) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          legend.text = element_text(size=20),
          legend.position = 'none')
  
  box_plot = cw_summary %>% 
    filter(Cellline == cellline & CumulativeWealth == 90) %>% 
    ggplot(aes(x=timepoint_replicate, y=UniqueBarcodeCount, fill=Timepoint), show.legend=F) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values=c("#C0C0C0","#808080", "#ffb3b3", "#cc0000")) +
    theme_minimal(base_size = 50) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          panel.grid.major = element_blank(),  axis.title.y = element_text(face = "italic", size = 30),
          legend.position = "none"
    ) +
    ylab("# unique barcodes\n(log10 scale)") +
    scale_y_log10(limits = c(1, 1250), labels = trans_format("log10", math_format(10^.x)))
  
  ggarrange(box_plot,tile_plot, heights = c(0.4, 2), ncol = 1, nrow = 2, align = "v", legend = "none")
  
  #cormat
}

plot_corr("HCC827") 
ggsave("plots/1/1B_replicates_clonetrace_HCC827_v2.pdf", width = 20, height = 22)

plot_corr("HCC4006")
ggsave("plots/1/1B_replicates_clonetrace_HCC4006_v2.pdf", width = 20, height = 22)

plot_corr("MGH707")
ggsave("plots/1/1B_replicates_clonetrace_MGH707_v2.pdf", width = 20, height = 22)

plot_corr("PC9")
ggsave("plots/1/1B_replicates_clonetrace_PC9_v2.pdf", width = 20, height = 22)

# ***************************************************************************************************************************************************


# Figure S1C,S1D,S1E
# ***************************************************************************************************************************************************
water_fall_v3$cl_bc = paste(water_fall_v3$cell_line, water_fall_v3$barcode, sep="_")
good_barcodes = water_fall_v3 %>% filter(cum_fraction <= 0.9) %>% select(cl_bc) %>% unique %>% as.list
a = water_fall_v3 %>% filter(cl_bc %in% good_barcodes$cl_bc) %>% as.data.table %>% dcast(cell_line+barcode ~ timepoint, value.var = "log2fc", fun.aggregate = mean)
colnames(a) = c("cl", "bc", "two_wk_tr", "six_day_unt", "six_wk_tr")
pal = brewer.pal(n = 8, name = "Accent")[c(5,6,7,8)]


ggplot(a, aes(x=six_day_unt,y=six_wk_tr,col=cl))+
  geom_point(aes(colour = factor(cl))) +
  geom_smooth(method=lm, aes(fill=cl)) +
  scale_color_manual(values=pal) +
  scale_fill_manual(values=pal) +
  facet_wrap("cl") +
  stat_cor(label.x = -7, label.y = 9, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_minimal(base_size = 25) +
  ylab("log2fc_6week_treated") +
  xlab("log2fc_6day_untreated") +
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position="none")

ggsave("plots/S1/S1C-6WT_vs_6DU.pdf", width = 10, height = 10)


ggplot(a, aes(x=six_day_unt,y=two_wk_tr,col=cl))+
  geom_point(aes(colour = factor(cl))) +
  geom_smooth(method=lm, aes(fill=cl)) +
  stat_cor(label.x = -7, label.y = 9, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  facet_wrap("cl") +
  scale_color_manual(values=pal) +
  scale_fill_manual(values=pal) +
  theme_minimal(base_size = 25) +
  ylab("log2fc_2week_treated") +
  xlab("log2fc_6day_untreated") +
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position="none")

ggsave("plots/S1/S1D-2WT_vs_6DU.pdf", width = 10, height = 10)


ggplot(a, aes(x=two_wk_tr,y=six_wk_tr,col=cl))+
  geom_point(aes(colour = factor(cl))) +
  geom_smooth(method=lm, aes(fill=cl)) +
  stat_cor(label.x = -7, label.y = 9, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  facet_wrap("cl") +
  scale_color_manual(values=pal) +
  scale_fill_manual(values=pal) +
  theme_minimal(base_size = 25) +
  ylab("log2fc_6week_treated") +
  xlab("log2fc_2week_treated") +
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position="none")

ggsave("plots/S1/S1E-6WT_vs_2WT.pdf", width = 10, height = 10)


# ***************************************************************************************************************************************************



# Figure 3C/S6E
# ***************************************************************************************************************************************************
cpm_dna_cast =  cpm_dna %>% dcast(cellline+cloneid ~ timepoint_revised, value.var="revised_cpm", fun.aggregate = mean)
cpm_dna_cast_annt =  cpm_dna_cast %>% merge(clonid_traj, by=c("cellline", "cloneid"), all.x=F)
cpm_dna_cast_annt <- return_log2fc_cpm_dna_cast_annt(cpm_dna_cast_annt)
cpm_dna_melt_annt <- return_cpm_dna_melt_annt(cpm_dna_cast_annt)


reqd_cellline = c("HCC4006", "HCC827")
#reqd_cellline = c("PC9", "MGH707")
reqd_classes = c("0", "1", "2", "3", "4", "earlier_0", "earlier_1", "later_0", "later_1")
tp_level = c("0 hr", "24 hr", "2 wk", "2 wk + 2 d")
cpm_dna_melt_annt %>% filter(cellline %in% reqd_cellline & trajectory_class %in% reqd_classes & variable %in% tp_level) %>% ggplot(aes(y=value, x=factor(variable, levels=tp_level), fill=trajectory_class)) +
  geom_bar(position=position_fill(reverse=T), stat="identity", width=1) +
  xlab("timepoint") + 
  facet_wrap("cellline") +
  ylab("trajectory class proportion") +
  theme_minimal(base_size = 30) + 
  theme(strip.text = element_text(size=35),
        axis.title.x = element_text(face = "italic", size = 25),
        axis.title.y = element_text(face = "italic", size = 25), 
        legend.text = element_text(size=25),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color="black", fill=NA)) +
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[1:10])) 

ggsave("plots/3/3C_proportion_class.pdf", width = 10, height = 7)
#ggsave("plots/S6/S6E_proportion_class.pdf", width = 10, height = 7)

# ***************************************************************************************************************************************************


# Figure 3D/S6F
# ***************************************************************************************************************************************************
reqd_classes = c("0", "1", "2", "3", "4", "earlier_0", "earlier_1", "later_0", "later_1")
tp_level = c("0 hr", "24 hr", "2 wk", "2 wk + 2 d")
totals_per_traj <- cpm_dna_melt_annt %>% filter(variable %in% tp_level) %>% group_by(cellline, trajectory_class, variable) %>% summarize(total = sum(value)) %>% as.data.table
totals_all <- cpm_dna_melt_annt %>% filter(variable %in% tp_level) %>% group_by(cellline, variable) %>% summarize(total = sum(value)) %>% as.data.table
proportion_data <- return_proportion_class(totals_per_traj, reqd_classes, tp_level)
proportion_data_melt <- return_proportion_data_melt(proportion_data)


#reqd_cellline = c("HCC4006", "HCC827")
reqd_cellline = c("PC9", "MGH707")
proportion_data_melt %>% filter(cellline %in% reqd_cellline & trajectory_class %in% reqd_classes & variable %in% tp_level) %>% 
  ggplot(aes(y=value, x=variable, group=trajectory_class, color=trajectory_class)) +
  geom_line(size=1.5) + 
  ylim(-0.75,0.75) +
  xlab("timepoint") + 
  facet_wrap("cellline") +
  ylab("trajectory class proportion change") +
  theme_minimal(base_size = 30) + 
  theme(strip.text = element_text(size=35),
        axis.title.x = element_text(face = "italic", size = 25),
        axis.title.y = element_text(face = "italic", size = 25), 
        legend.text = element_text(size=25),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color="black", fill=NA)) +
  scale_colour_manual(values = c(brewer.pal(9, "Set1")[1:10]))

#ggsave("plots/3/3D_proportion_change.pdf", width = 10, height = 7)
ggsave("plots/S6/S6F_proportion_change.pdf", width = 10, height = 7)

# ***************************************************************************************************************************************************



# Figure 3E/S6G
# ***************************************************************************************************************************************************
reqd_cellline = c("HCC4006", "HCC827")
#reqd_cellline = c("PC9", "MGH707")
cpm_rna_cast_annt %>% filter(cellline %in% reqd_cellline & trajectory_class %in% c("0", "1", "2", "3", "4")) %>% ggplot(aes(y=D14_D0_log2fc, x=factor(trajectory_class, levels=c("0", "1", "2", "3", "4")), fill=trajectory_class)) +
  geom_boxplot() +
  geom_jitter(size=2, width = 0.1, height = 0.5) +
  xlab("trajectory_class") + 
  facet_wrap("cellline") +
  ylab("difference 14 day vs 0hr") +
  #ylab("log2fc 14 day vs 0hr") +
  theme_minimal(base_size = 30) + 
  theme(strip.text = element_text(size=35),
        axis.title.x = element_text(face = "italic", size = 25),
        axis.title.y = element_text(face = "italic", size = 25), 
        legend.text = element_text(size=25),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color="black", fill=NA)) +
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[1:10])) 

ggsave("plots/3/3E_rna_log2fc.pdf", width = 15, height = 10)
#ggsave("plots/S6/S6_G_rna_log2fc.pdf", width = 15, height = 10)

# ***************************************************************************************************************************************************


# Figure S8A
# ***************************************************************************************************************************************************
ltr_rep_overlap_cid <- ltr_rep_overlap %>% merge(clonid_traj, by.x=c("cell_line", "cloneid"), by.y=c("cellline", "cloneid")) %>% unique

for(cl in ltr_rep_overlap_cid$cell_line %>% unique){
  for(gp in (ltr_rep_overlap_cid %>% filter(cell_line == cl))$group %>% unique()){
    dt <- ltr_rep_overlap_cid %>% filter(cell_line == cl, group == query_group & group == gp) %>% group_by(group, query_group, replicate, count) %>% tally %>% as.data.table
    dt$count <- as.character(dt$count)
    dt$rank <- dense_rank(dt$replicate)
    dt %>% ggplot(aes(x=rank, y=n)) + 
      geom_bar(position="fill", stat="identity", aes(fill=count)) +
      xlab("replicate number") + 
      ylab("% barcodes shared") +
      theme_minimal(base_size = 25) +
      theme(axis.title.x = element_text(face = "italic", size = 25),
            axis.title.y = element_text(face = "italic", size = 25),
            plot.title = element_text(hjust=0.5)) +
      #scale_x_discrete(labels=c("2" = "1", "3" = "2", "5" = "3")) + 
      labs(fill='shared by # replicates') +
      ggtitle(paste0(gp, "_wk")) +
      scale_fill_grey(start = 0.8, end = 0.2)
    
    ggsave(paste0("plots/S8/S8A_", cl, "_", gp, "_wk","_sharing.pdf"), width = 10, height = 10)
  }
}

# ***************************************************************************************************************************************************



# Figure 4E/S8C
# ***************************************************************************************************************************************************


ltr_rep_avg_cid <- ltr_rep_avg %>% merge(clonid_traj, by.x=c("cloneid", "cell_line"), by.y=c("cloneid", "cellline"))


cl = "HCC827"

a = ggplot(ltr_rep_avg_cid %>% filter(cell_line == cl, trajectory_class != "-" & number_replicates > 1 & group %in% c("Day_0", "EGF816_4", "EGF816_7", "EGF816-INC280_7")),
           aes(x=factor(group, level=c("Day_0", "EGF816_4", "EGF816_7", "EGF816-INC280_7")), y=median_cpm)) + 
  geom_bar(position="fill", stat="identity", aes(fill=trajectory_class)) +
  xlab("condition") + 
  ylab("% contribution") +
  theme_minimal(base_size = 45) +
  theme(axis.title.x = element_text(face = "italic", size = 45),
        axis.title.y = element_text(face = "italic", size = 45),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, size = 50),
        legend.position = "none"
  ) +
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[1:8])) +
  annotate("text", label = "69", x = "Day_0", y = 1.05, size = 15) +
  annotate("text", label = "22", x = "EGF816_4", y = 1.05, size = 15) +
  annotate("text", label = "16", x = "EGF816_7", y = 1.05, size = 15) +
  annotate("text", label = "11", x = "EGF816-INC280_7", y = 1.05, size = 15) +
  scale_x_discrete(labels=c("Day_0"="Day0","EGF816_4" = "EGF816_4wk", "EGF816_7" = "EGF816_7wk", "EGF816-INC280_7" = "EGF816-INC280_7wk")) +
  ggtitle(cl)

a
ggsave(paste0("plots/4/4E_tj_cont_", cl, "_v2.pdf"), width = 10, height = 15)


cl = "HCC4006"

a = ggplot(ltr_rep_avg_cid %>% filter(cell_line == cl, trajectory_class != "-" & number_replicates > 1), aes(x=factor(group, level=c("Day_0","EGF816_4", "EGF816_7", "EGF816_10")), y=median_cpm)) + 
  geom_bar(position="fill", stat="identity", aes(fill=trajectory_class)) +
  xlab("condition") + 
  ylab("% contribution") +
  theme_minimal(base_size = 45) +
  theme(axis.title.x = element_text(face = "italic", size = 45),
        axis.title.y = element_text(face = "italic", size = 45),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 50),
        legend.position = "none"
  ) + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[1:5])) +
  annotate("text", label = "10", x = "Day_0", y = 1.05, size = 15) +
  annotate("text", label = "8", x = "EGF816_4", y = 1.05, size = 15) +
  annotate("text", label = "8", x = "EGF816_7", y = 1.05, size = 15) +
  scale_x_discrete(labels=c("Day_0"="Day0","EGF816_4" = "EGF816_4wk", "EGF816_7" = "EGF816_7wk")) +
  ggtitle(cl)


a
ggsave("plots/S8/S8C_tj_cont_HCC4006_v2.pdf", width = 10, height = 15)

# ***************************************************************************************************************************************************


# Figure 4G/S8D
# ***************************************************************************************************************************************************

chemo_rep_overlap_cid <- chemo_rep_overlap %>% merge(clonid_traj, by.x=c("cell_line", "cloneid"), by.y=c("cellline", "cloneid"))

cl = "HCC827"

a = chemo_rep_overlap_cid %>% filter(cell_line == cl, trajectory_class != "-" & number_replicates > 1) %>% 
  ggplot(aes(x=factor(mode, level=c("Day0", "EGF816", "CarboPeme", "Doce", "EGF816-CarboPeme", "EGF816-Doce")), y=median_cpm)) + 
  geom_bar(position="fill", stat="identity", aes(fill=trajectory_class)) +
  xlab("condition") + 
  ylab("% contribution") +
  theme_minimal(base_size = 45) +
  theme(axis.title.x = element_text(face = "italic", size = 45),
        axis.title.y = element_text(face = "italic", size = 45),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 50),
        legend.position = "none"
  ) + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[1:7])) +
  annotate("text", label = "59", x = "CarboPeme", y = 1.05, size = 15) +
  annotate("text", label = "69", x = "Day0", y = 1.05, size = 15) +
  annotate("text", label = "37", x = "Doce", y = 1.05, size = 15) +
  annotate("text", label = "20", x = "EGF816", y = 1.05, size = 15) +
  annotate("text", label = "18", x = "EGF816-CarboPeme", y = 1.05, size = 15) +
  annotate("text", label = "18", x = "EGF816-Doce", y = 1.05, size = 15) +
  ggtitle("HCC827")

#leg = get_legend(a)
#as_ggplot(leg)
#ggsave(paste0("plots/4/4_legend_",cl,".pdf"))

a 
ggsave("plots/4/4G_tj_chemo_hcc827_v2.pdf", width = 10, height = 15)

hcc4006_cm = chemo_rep_overlap  %>% filter((treatment == "Cons" | is.na(treatment)) & cell_line == "HCC4006")
hcc4006_cm_cid = merge(hcc4006_cm, clonid_traj %>% filter(cellline == "HCC4006"), by="cloneid") %>% unique
hcc4006_cm_cid_cellcounts = merge(hcc4006_cm_cid %>% filter(mode != "Day0"), hcc4006_cell_counts, by=c("mode", "treatment"))
hcc4006_cm_cid_cellcounts$ctg_corrected_cpm = ((hcc4006_cm_cid_cellcounts$median_cpm * hcc4006_cm_cid_cellcounts$cellcounts) / 10 ^ 6)

a = hcc4006_cm_cid_cellcounts %>% filter(trajectory_class %in% c("0", "1", "2", "3", "4") & number_replicates > 1 
                                         & mode %in% c("EGF816", "EGF816-CarboPeme")
) %>% 
  #ggplot(aes(x=factor(mode, level=c("Day0", "EGF816", "CarboPeme", "Doce", "EGF816-CarboPeme", "EGF816-Doce")), y=median_cpm)) +
  ggplot(aes(x=factor(mode, level=c("EGF816", "CarboPeme", "EGF816-CarboPeme")), y=ctg_corrected_cpm)) +
  geom_bar(position="stack", stat="identity", aes(fill=trajectory_class)) +
  xlab("condition") + 
  ylab("cell counts") +
  theme_minimal(base_size = 45) +
  theme(axis.title.x = element_text(face = "italic", size = 45),
        axis.title.y = element_text(face = "italic", size = 45),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 50),
        legend.position = "none"
  ) + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[2:7])) 
ggtitle("HCC4006")

a
ggsave("plots/S8/S8D_cell_counts_treatments_HCC4006_v2.pdf", width = 10, height = 15)


# ***************************************************************************************************************************************************


# Figure 4H/S8E
# ***************************************************************************************************************************************************
counts_summary_drug_combo <- counts_summary_drug_combo %>% merge(clonid_traj %>% filter(cellline == "HCC827"), by="cloneid")

ggplot(counts_summary_drug_combo %>% filter(treatment %in% c("EGF816","EGF816+INC280","EGF816+Tram","EGF816+BGJ398","EGF816+Carb/Peme") & trajectory_class != "-" & !grepl("earlier", trajectory_class)), aes(y=cell_count_revised, x=factor(treatment, level = c("EGF816", "EGF816+INC280", "EGF816+BGJ398", "EGF816+Carb/Peme", "EGF816+Tram")))) + 
  geom_bar(position="stack", stat="identity", aes(fill=trajectory_class)) +
  xlab("condition") + 
  ylab("cell counts") +
  theme_minimal(base_size = 45) +
  theme(axis.title.x = element_text(face = "italic", size = 45),
        axis.title.y = element_text(face = "italic", size = 45),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") + 
  #)+
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[1:5]))#,"black",brewer.pal(9,"Set1")[7:9]))

ggsave("plots/4/4H_cell_counts_treatments_v2_HCC827.pdf", width = 10, height = 15)


a = hcc4006_cm_cid %>% filter(trajectory_class %in% c("0", "1", "2", "3", "4") & number_replicates > 1 
                              #& mode %in% c("EGF816", "EGF816-CarboPeme")
) %>% 
  ggplot(aes(x=factor(mode, level=c("Day0", "EGF816", "CarboPeme", "Doce", "EGF816-CarboPeme", "EGF816-Doce")), y=median_cpm)) +
  #ggplot(aes(x=factor(mode, level=c("EGF816", "CarboPeme", "EGF816-CarboPeme")), y=ctg_corrected_cpm)) +
  geom_bar(position="fill", stat="identity", aes(fill=trajectory_class)) +
  xlab("condition") + 
  ylab("% contribution") +
  theme_minimal(base_size = 45) +
  theme(axis.title.x = element_text(face = "italic", size = 45),
        axis.title.y = element_text(face = "italic", size = 45),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 50),
        legend.position = "none"
  ) + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[1:7])) +
  annotate("text", label = "22", x = "CarboPeme", y = 1.05, size = 15) +
  annotate("text", label = "10", x = "Day0", y = 1.05, size = 15) +
  annotate("text", label = "5", x = "Doce", y = 1.05, size = 15) +
  annotate("text", label = "7", x = "EGF816", y = 1.05, size = 15) +
  annotate("text", label = "8", x = "EGF816-CarboPeme", y = 1.05, size = 15) +
  annotate("text", label = "19", x = "EGF816-Doce", y = 1.05, size = 15) +
  ggtitle("HCC4006")


#leg = get_legend(a)
#as_ggplot(leg)
#ggsave(paste0("plots/S8/S8_legend_",cl,".pdf"))

a
ggsave("plots/S8/S8E_tj_chemo_hcc4006_v2.pdf", width = 10, height = 15)

# ***************************************************************************************************************************************************

renv::snapshot()
