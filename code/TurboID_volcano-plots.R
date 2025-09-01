# DVL3 wt vs. DVL3 S/TA interactome 

# TurboID of DVL3-wt, 3 DVL3-mut (delS/T, S/TA, S/TE) and control

# Libraries ####
library(here)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggprism)

# Input ####
data <- read.csv(here("input", "ProteinGroups.csv"))

# Data processing ####
# Omit reduntant columns
data.selected <- data %>% 
  select(PG.ID, UniProt_Entry, UniProt_Entry.Name, UniProt_Gene.Names..primary., UniProt_Protein.names,
         starts_with("PG.MaxLFQ"), starts_with("LIMMA"))


# Rename columns - MaxLFQ
names(data.selected) = gsub(pattern = "TurboID.", 
                            replacement = "", 
                            x = names(data.selected)) 

patterns <- c("DVL3_delS.T_C2", "DVL3_S.T_A_B3","DVL3_S.T_E_D2","DVL3_T5.2", "STOP_clone1")
replacements <- c("DVL3.delST", "DVL3.STA", "DVL3.STE", "DVL3.wt", "ctrl")

for (i in seq_along(patterns)) {
  colnames(data.selected) <- gsub(patterns[i], replacements[i], colnames(data.selected))
}

# Rename columns - LIMMA
names(data.selected) = gsub(pattern = "LIMMA_MaxLFQ_all.samples_prec.filter100_norm_imp_", 
                            replacement = "", 
                            x = names(data.selected)) 

patterns <- c("delST.", "STA.", "STE.",  "wt.", "ctrl.")
replacements <- c("DVL3.delST_", "DVL3.STA_", "DVL3.STE_", "DVL3.wt_", "ctrl_")

for (i in seq_along(patterns)) {
  colnames(data.selected) <- gsub(patterns[i], replacements[i], colnames(data.selected))
}


# Define significantly upregulated (logFC > 1, adj.P.Value < 0.05) and downregulated (logFC < -1, adj.P.Value < 0.05) proteins
get_group <- function(logFC, adjP) {
  if (is.na(logFC) || is.na(adjP)) {
    return(NA)
  } else if (adjP > 0.05) {
    return("NS")
  } else if (logFC > 1) {
    return("Up")
  } else if (logFC < -1) {
    return("Down")
  } else {
    return("NS")
  }
}

for (colname in names(data.selected)) {
  if (grepl("_logFC$", colname)) {
    prefix <- sub("_logFC$", "", colname)
    pval_colname <- paste0(prefix, "_adj.P.Val")
    group_colname <- paste0(prefix, "_group")
    
    data.selected[[group_colname]] <- mapply(get_group, data.selected[[colname]], data.selected[[pval_colname]])
  }
}


# Omit proteins with NA in PG.MaxLFQ in all samples
cols <- grep("^PG\\.MaxLFQ", names(data.selected), value = TRUE)
cols <- cols[!grepl("_log2$", cols)]

data.selected <- data.selected[!apply(data.selected[, cols], 1, function(row) all(is.na(row))), ]

# Define DVL3-specific proteins as upregulated in at least one DVL3-wt/mut TurboID compared to TurboID ctrl
data.selected <- data.selected %>% 
  mutate(DVL3.specific = case_when(DVL3.wt_ctrl_group == "Up" ~ "yes", 
                                   DVL3.delST_ctrl_group == "Up" ~ "yes",
                                   DVL3.STA_ctrl_group == "Up" ~ "yes",
                                   DVL3.STE_ctrl_group == "Up" ~ "yes",
                                   TRUE ~ "no"))


# How many proteins were detected?
nrow(data.selected)

# How many proteins were specifically upregulated?
nrow(data.selected %>% filter(DVL3.specific == "yes"))

# How many proteins were specifically upregulated (percentage)?
nrow(data.selected %>% filter(DVL3.specific == "yes")) / nrow(data.selected) * 100

# Volcano plots #### 
# Plot function
volcano <- function(data, contrast, preys) {
  ggplot(data, aes(x = .data[[paste0(contrast, "_logFC")]], 
                   y = -log10(.data[[paste0(contrast, "_adj.P.Val")]]))) + 
    geom_point(alpha = 0.6, colour = "grey") + 
    labs(x = expression(log[2]*" fold change"),
         y = expression(-log[10]*" adj. P Value")) +
    theme_prism() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = -1, linetype = "dashed") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    # Label and enhance specific preys
    geom_label_repel(data = . %>% filter(UniProt_Gene.Names..primary. %in% preys),
                     mapping = aes(label = UniProt_Gene.Names..primary.),
                     force = 100,
                     colour = "black",
                     box.padding = 0.5) +
    geom_point(data = . %>% filter(UniProt_Gene.Names..primary. %in% preys),
               colour = "darkred")
}

# Proteins of interest
poi <- c("DVL1", "DVL2", "DVL3", "VANGL1", "VANGL2", "CSNK1E", "CSNK1D", "AXIN1", "AXIN2")
up <- c("FZD3", "FZD5", "FZD6", "ADGRA2", "RNF43", "RYK", "ZNRF3")
down <- c("CCDC88A", "CCDC88C", "TNKS", "TNKS2")
updown <- c("FZD3", "FZD5", "FZD6", "ADGRA2", "RNF43", "RYK", "ZNRF3", "CCDC88A", "CCDC88C", "TNKS", "TNKS2")

# Plots
dir.create(here("outputs"))

svg(filename = here("outputs", "DVL3.wt-vs-ctrl_poi.svg"), width = 4, height = 4)
volcano(data = data.selected,
        contrast = "DVL3.wt_ctrl",
        preys = poi) +
  xlim(-9, 10.2) +
  ylim(0, 23)
dev.off()

svg(filename = here("outputs", "DVL3.STA-vs-ctrl_poi.svg"), width = 4, height = 4)
volcano(data = data.selected,
        contrast = "DVL3.STA_ctrl",
        preys = poi) +
  xlim(-9, 10.2) +
  ylim(0, 23)
dev.off()

svg(filename = here("outputs", "DVL3.STE-vs-ctrl_poi.svg"), width = 4, height = 4)
volcano(data = data.selected,
        contrast = "DVL3.STE_ctrl",
        preys = poi) +
  xlim(-9, 10.2) +
  ylim(0, 23)
dev.off()

svg(filename = here("outputs", "DVL3.delST-vs-ctrl_poi.svg"), width = 4, height = 4)
volcano(data = data.selected,
        contrast = "DVL3.delST_ctrl",
        preys = poi) +
  xlim(-9, 10.2) +
  ylim(0, 23)
dev.off()

svg(filename = here("outputs", "DVL3.STA-vs-DVL3.wt_updown.svg"), width = 4, height = 4)
volcano(data = data.selected,
        contrast = "DVL3.STA_DVL3.wt",
        preys = updown) +
  xlim(-6.8, 8.6) +
  ylim(0, 18)
dev.off()

svg(filename = here("outputs", "DVL3.STE-vs-DVL3.wt_updown.svg"), width = 4, height = 4)
volcano(data = data.selected,
        contrast ="DVL3.STE_DVL3.wt",
        preys = updown) +
  xlim(-6.8, 8.6) +
  ylim(0, 18)
dev.off()

svg(filename = here("outputs", "DVL3.delST-vs-DVL3.wt_updown.svg"), width = 4, height = 4)
volcano(data = data.selected,
        contrast ="DVL3.delST_DVL3.wt",
        preys = updown) +
  xlim(-6.8, 8.6) +
  ylim(0, 18)
dev.off()

# Export data ####
data.export <- data.selected %>% 
  select(starts_with("UniProt"),
         DVL3.specific,
         starts_with("DVL3.wt_ctrl"),
         starts_with("DVL3.STA_ctrl"),
         starts_with("DVL3.STE_ctrl"),
         starts_with("DVL3.delST_ctrl"),
         starts_with("DVL3.STA_DVL3.wt"),
         starts_with("DVL3.STE_DVL3.wt"),
         starts_with("DVL3.delST_DVL3.wt")) %>% 
  select(-matches("_B$|_t$|_P.Value$|_delog$"))

write.csv(data.export, here("outputs", "TurboID_DVL3-STA_vs_DVL3-wt.csv"), row.names = FALSE)

