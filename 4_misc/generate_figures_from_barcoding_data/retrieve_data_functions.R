library(data.table)

# Folder locations
clonTracer_folder = "/da/onc/BFx/research/krishvi7/barcoding/clontracer/2018110902_EGF816_Young_Matt/output"
expressedBarcode_folder = "/da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20190805_4cellline/output_v3_DNA_info/"
dna_experiment_folder = "/da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20191219_4cellline_DNA/output"
ltr_experiment_folder = "/da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20200415_long_term_resistance_invivo_DNA/output"
hr_experiment_folder = "/da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20200410_holiday_retreat_DNA/output"
treat_experiment_folder = "/da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20190923_HCC827_drug_DNA/output"
javi_classes_folder = "/da/onc/krishvi7/bitbucket/work-krishvi7-oneoffscripts/jupyter_notebooks/ITT/expressed_barcode"
chemo_experiment_folder = "/da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20200827_HCC827_HCC4006_Chemo/output"

# Functions to retrieve data
return_cw_summary_clontracer <- function(){
  cw_summary <- fread(paste(clonTracer_folder,"cw_summary.txt",sep="/")) 
  cw_summary$timepoint_replicate = paste(cw_summary$Timepoint, cw_summary$Replicate, sep="_") 
  cw_summary %>% return
}

return_ts_sorted_clontracer <- function(){
  fread(paste(clonTracer_folder,"tall_skinny_sorted_v2_cum_rep.txt",sep="/")) %>% return
}

return_waterfall_clontracer <- function(){
  fread(paste(clonTracer_folder,"waterfall_v3.txt",sep="/")) %>% return
}

return_cpm_ctg_dh <- function(){
  fread(paste(hr_experiment_folder, "counts_summary_annt_cloneid_ctg.txt", sep="/")) %>% return
}

return_cpm_dna <- function(){
  cpm_dna <- fread(paste(dna_experiment_folder, "counts_summary_cloneid.txt", sep="/")) 
  cpm_dna = cpm_dna %>% separate(sample, into=c("cellline", "timepoint"), sep="_", remove = F)
  cpm_dna$timepoint_revised = "DX"
  cpm_dna = cpm_dna[timepoint == "0", timepoint_revised := "0 hr"]
  cpm_dna = cpm_dna[timepoint == "1", timepoint_revised := "24 hr"]
  cpm_dna = cpm_dna[timepoint == "14", timepoint_revised := "2 wk"]
  cpm_dna = cpm_dna[timepoint == "16", timepoint_revised := "2 wk + 2 d"]
  
  cpm_dna %>% return
}


return_cpm_rna <- function(){
  cpm_rna <- fread(paste(expressedBarcode_folder, "counts_summary_annt_cloneid_revised.txt", sep="/")) 
  cpm_rna = cpm_rna %>% separate(sample, into=c("cellline", "timepoint"), sep="_", remove = F)
  
  cpm_rna$timepoint_revised = "DX"
  cpm_rna = cpm_rna[timepoint == "0", timepoint_revised := "0 hr"]
  cpm_rna = cpm_rna[timepoint == "1", timepoint_revised := "24 hr"]
  cpm_rna = cpm_rna[timepoint == "14", timepoint_revised := "2 wk"]
  cpm_rna = cpm_rna[timepoint == "16", timepoint_revised := "2 wk + 2 d"]
  cpm_rna_cast =  cpm_rna %>% dcast(cellline+cloneid ~ timepoint_revised, value.var="revised_cpm", fun.aggregate = mean, fill = 0.0)
  cpm_rna_cast_annt =  cpm_rna_cast %>% merge(clonid_traj, by=c("cellline", "cloneid"), all.x=F)
  
  cpm_rna_cast_annt$D0_minus_D0 = cpm_rna_cast_annt$'0 hr' - cpm_rna_cast_annt$'0 hr'
  cpm_rna_cast_annt$D14_minus_D0 = cpm_rna_cast_annt$'2 wk' - cpm_rna_cast_annt$'0 hr'
  cpm_rna_cast_annt$D1_minus_D0 = cpm_rna_cast_annt$'24 hr' - cpm_rna_cast_annt$'0 hr'
  cpm_rna_cast_annt$D16_minus_D0 = cpm_rna_cast_annt$'2 wk + 2 d' - cpm_rna_cast_annt$'0 hr'
  
  cpm_rna_cast_annt$D14_D0_log2fc = log2(cpm_rna_cast_annt$'2 wk' + 1) - log2(cpm_rna_cast_annt$'0 hr' + 1)
  cpm_rna_cast_annt$D1_D0_log2fc = log2(cpm_rna_cast_annt$'24 hr' + 1) - log2(cpm_rna_cast_annt$'0 hr' + 1)
  cpm_rna_cast_annt$D16_D0_log2fc = log2(cpm_rna_cast_annt$'2 wk + 2 d' + 1) - log2(cpm_rna_cast_annt$'0 hr' + 1)
  
  cpm_rna_cast_annt %>% return
}

return_log2fc_cpm_dna_cast_annt <- function(cpm_dna_cast_annt){
  cpm_dna_cast_annt$D0_minus_D0 = cpm_dna_cast_annt$'0 hr' - cpm_dna_cast_annt$'0 hr'
  cpm_dna_cast_annt$D14_minus_D0 = cpm_dna_cast_annt$'2 wk' - cpm_dna_cast_annt$'0 hr'
  cpm_dna_cast_annt$D1_minus_D0 = cpm_dna_cast_annt$'24 hr' - cpm_dna_cast_annt$'0 hr'
  cpm_dna_cast_annt$D16_minus_D0 = cpm_dna_cast_annt$'2 wk + 2 d' - cpm_dna_cast_annt$'0 hr'
  
  cpm_dna_cast_annt$D14_D0_log2fc = log2(cpm_dna_cast_annt$'2 wk' + 1) - log2(cpm_dna_cast_annt$'0 hr' + 1)
  cpm_dna_cast_annt$D1_D0_log2fc = log2(cpm_dna_cast_annt$'24 hr' + 1) - log2(cpm_dna_cast_annt$'0 hr' + 1)
  cpm_dna_cast_annt$D16_D0_log2fc = log2(cpm_dna_cast_annt$'2 wk + 2 d' + 1) - log2(cpm_dna_cast_annt$'0 hr' + 1)
  
  cpm_dna_cast_annt %>% return
}

return_cpm_dna_melt_annt <- function(cpm_dna_cast_annt){
  cpm_dna_melt_annt = cpm_dna_cast_annt %>% melt(id.vars = c("cellline", "cloneid", "trajectory_class"))
  cpm_dna_melt_annt$variable = cpm_dna_melt_annt$variable %>% as.character
  
  cpm_dna_melt_annt %>% return
}


return_proportion_class <- function(totals_per_traj, traj_class, tp_level){
  proportion_data = totals_per_traj$cellline %>% unique %>% lapply(function(cl){
    traj_class %>% lapply(function(tj){
      tp_level %>% lapply(function(tp){
        if (tj %in% (totals_per_traj %>% filter(cellline == cl))$trajectory_class %>% unique){
          sum_clones = (totals_per_traj %>% filter(trajectory_class == tj & variable == tp & cellline == cl))$total
          sum_all = (totals_all %>% filter(cellline == cl & variable == tp))$total
          v = sum_clones/sum_all
          data.table(data.frame(
            cellline = cl,
            trajectory_class = tj,
            timepoint = tp,
            sum_clones = sum_clones,
            sum_all = sum_all,
            value = v
          )) %>% return 
        }
      }) %>% rbindlist
    }) %>% rbindlist
  }) %>% rbindlist
  
  proportion_data %>% return
}

return_proportion_data_melt <- function(proportion_data){
  proportion_data_cast <- proportion_data %>% dcast(cellline + trajectory_class ~ timepoint, value.var = "value", fun.aggregate = mean)
  names(proportion_data_cast) = c("cellline", "trajectory_class", "0", "2", "2+2", "24")
  proportion_data_cast$'0 hr' = proportion_data_cast$'0' - proportion_data_cast$'0'
  proportion_data_cast$'24 hr' = proportion_data_cast$'24' - proportion_data_cast$'0'
  proportion_data_cast$'2 wk' = proportion_data_cast$'2' - proportion_data_cast$'0'
  proportion_data_cast$'2 wk + 2 d' = proportion_data_cast$'2+2' - proportion_data_cast$'0'
  proportion_data_melt = proportion_data_cast %>% melt(id.vars = c("cellline", "trajectory_class"))
  
  proportion_data_melt %>% return
}

return_log2fc_combo_treat <- function(){
  fread(paste(treat_experiment_folder,"log2fc.txt",sep="/")) %>% return
}

return_log2fc_combo_treat_cc <- function(){
  fread(paste(treat_experiment_folder,"log2fc_cellcounts.txt",sep="/")) %>% return
}

return_log2fc_rna_dna <- function(){
  fread(paste(dna_experiment_folder,"log2fc_all_rna_dna.txt",sep="/")) %>% return
}

return_counts_summary_drug_combo <- function(){
  counts_summary_drug_combo <- fread(paste(treat_experiment_folder, "counts_annt_with_cells_cloneid_rerevised.txt",sep="/"))
  
  counts_summary_drug_combo$treatment = "NA"
  counts_summary_drug_combo[grepl("Day_0", sample), treatment:="Untreated_D0"]
  counts_summary_drug_combo[grepl("Day_6", sample), treatment:="Untreated_D6"]
  counts_summary_drug_combo[grepl("EGF816_S", sample), treatment:="EGF816"]
  counts_summary_drug_combo[grepl("EGF816_INC280", sample), treatment:="EGF816+INC280"]
  counts_summary_drug_combo[grepl("EGF816_Trametinib", sample), treatment:="EGF816+Tram"]
  counts_summary_drug_combo[grepl("EGF816_BGJ398", sample), treatment:="EGF816+BGJ398"]
  counts_summary_drug_combo[grepl("Carboplatin_Pemetrexed", sample), treatment:="Carb/Peme"]
  counts_summary_drug_combo[grepl("EGF816_Carboplatin_Pemetrexed", sample), treatment:="EGF816+Carb/Peme"]
  
  counts_summary_drug_combo %>% return
}

return_traj_classes_javi <- function(){
  cloneid_traj_HCC827 =  fread(paste0(javi_classes_folder, "/CloneClassMappingHCC827.csv"))
  cloneid_traj_HCC827$cellline = "HCC827"
  cloneid_traj_PC9 =  fread(paste0(javi_classes_folder, "/CloneClassMappingPC9.csv"))
  cloneid_traj_PC9$cellline = "PC9"
  cloneid_traj_HCC4006 =  fread(paste0(javi_classes_folder, "/CloneClassMappingHCC4006.csv"))
  cloneid_traj_HCC4006$cellline = "HCC4006"
  cloneid_traj_MGH707 =  fread(paste0(javi_classes_folder, "/CloneClassMappingMGH707.csv"))
  cloneid_traj_MGH707$cellline = "MGH707"
  
  clonid_traj = cloneid_traj_HCC827 %>% rbind(cloneid_traj_PC9) %>% rbind(cloneid_traj_HCC4006) %>% rbind(cloneid_traj_MGH707)
  clonid_traj %>% filter(trajectory_class != "-") %>% select(trajectory_class, cloneid, cellline) %>% unique  %>% as.data.table %>% return
}

return_drug_sensitivities <- function(){
  fread(paste(hr_experiment_folder,"counts_summary_sensitivities_v2.txt",sep="/")) %>% return
}

return_ltr_data <- function(){
  fread(paste(ltr_experiment_folder, "counts_summary_annt_lt_cloneid.txt", sep="/")) %>% return
}

return_ltr_rep_avg <- function(){
  fread(paste(ltr_experiment_folder, "replicate_avg_with_Day0.txt", sep="/")) %>% return
}

return_ltr_rep_overlap <- function(){
  fread(paste(ltr_experiment_folder, "replicate_overlaps.txt", sep="/")) %>% return
}

return_chemo_rep_overlap <- function(){
  fread(paste(chemo_experiment_folder, "replicate_avg.txt", sep="/"),  na.strings = "NA") %>% return
}

return_hcc4006_cellcounts <- function(){
  fread(paste(chemo_experiment_folder, "hcc4006_cellcounts.txt", sep="/"),  na.strings = "NA") %>% return
}
