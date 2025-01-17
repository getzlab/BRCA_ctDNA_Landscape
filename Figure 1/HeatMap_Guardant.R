#load libraries
library(tidyverse)
library(ggplot2)

#Join patient classifiers with df dataset
df <- Guardant_Table
df <- unique(df)

#Select first Guardant tests and put in dataframe df
df <- df[grepl("_01", df$test_id),]
df <- unique(df)

#Prep heat map categories to all mutations (remove fusions and amplifications)
df <- df %>% 
  filter(mut_type != "amp") %>% 
  select(patient_id, gene, percentage, mut_type)

df_all <- df %>% 
  select(patient_id, gene)
df_all <- unique(df_all)
df_all <- df_all %>% 
  count(gene)
df_all <- df_all %>% 
  mutate("Total" = n) 
df_all <- df_all[,-2]

#Sort by genes and MAF
#<0.1, 0.1-0.3%, 0.4-0.6%, 0.7-1%, 1-5%, 5-10%, 10-40%, 40-60%, 60-100%
df_0.1 <- df %>% 
  filter(percentage < 0.1) %>% 
  select(patient_id, gene)
df_0.1 <- unique(df_0.1)
df_0.1 <- df_0.1 %>% 
  count(gene)
df_0.1 <- df_0.1 %>% 
  mutate("df_0.1" = n)
df_0.1 <- df_0.1[,-2]

df_0.1_0.2 <- df %>% 
  filter(percentage >= 0.1) %>% 
  filter(percentage <0.2) %>% 
  select(patient_id, gene)
df_0.1_0.2 <- unique(df_0.1_0.2)
df_0.1_0.2 <- df_0.1_0.2 %>% 
  count(gene)
df_0.1_0.2 <- df_0.1_0.2 %>% 
  mutate("df_0.1_0.2" = n)
df_0.1_0.2 <- df_0.1_0.2[,-2]

df_0.2_0.4 <- df %>% 
  filter(percentage >= 0.2) %>% 
  filter(percentage <0.4) %>% 
  select(patient_id, gene)
df_0.2_0.4 <- unique(df_0.2_0.4)
df_0.2_0.4 <- df_0.2_0.4 %>% 
  count(gene)
df_0.2_0.4 <- df_0.2_0.4 %>% 
  mutate("df_0.2_0.4" = n)
df_0.2_0.4 <- df_0.2_0.4[,-2]

df_0.4_0.6 <- df %>% 
  filter(percentage >= 0.4) %>% 
  filter(percentage <0.6) %>% 
  select(patient_id, gene)
df_0.4_0.6 <- unique(df_0.4_0.6)
df_0.4_0.6 <- df_0.4_0.6 %>% 
  count(gene)
df_0.4_0.6 <- df_0.4_0.6 %>% 
  mutate("df_0.4_0.6" = n)
df_0.4_0.6 <- df_0.4_0.6[,-2]

df_0.6_0.8 <- df %>% 
  filter(percentage >= 0.6) %>% 
  filter(percentage <0.8) %>% 
  select(patient_id, gene)
df_0.6_0.8 <- unique(df_0.6_0.8)
df_0.6_0.8 <- df_0.6_0.8 %>% 
  count(gene)
df_0.6_0.8 <- df_0.6_0.8 %>% 
  mutate("df_0.6_0.8" = n)
df_0.6_0.8 <- df_0.6_0.8[,-2]

df_1 <- df %>% 
  filter(percentage >= 0.8) %>% 
  filter(percentage <1) %>% 
  select(patient_id, gene)
df_1 <- unique(df_1)
df_1 <- df_1 %>% 
  count(gene)
df_1 <- df_1 %>% 
  mutate("df_0.8_1" = n)
df_1 <- df_1[,-2]

df_1_5 <- df %>% 
  filter(percentage >= 1) %>% 
  filter(percentage <5) %>% 
  select(patient_id, gene)
df_1_5 <- unique(df_1_5)
df_1_5 <- df_1_5 %>% 
  count(gene)
df_1_5 <- df_1_5 %>% 
  mutate("df_1_5" = n)
df_1_5 <- df_1_5[,-2]

df_5_10 <- df %>% 
  filter(percentage >= 5) %>% 
  filter(percentage <10) %>% 
  select(patient_id, gene)
df_5_10 <- unique(df_5_10)
df_5_10 <- df_5_10 %>% 
  count(gene)
df_5_10 <- df_5_10 %>% 
  mutate("df_5_10" = n) 
df_5_10 <- df_5_10[,-2]

df_10_40 <- df %>% 
  filter(percentage >= 10) %>% 
  filter(percentage <40) %>% 
  select(patient_id, gene)
df_10_40 <- unique(df_10_40)
df_10_40 <- df_10_40 %>% 
  count(gene)
df_10_40 <- df_10_40 %>% 
  mutate("df_10_40" = n)
df_10_40 <- df_10_40[,-2]

df_40_60 <- df %>% 
  filter(percentage >= 40) %>% 
  filter(percentage <60) %>% 
  select(patient_id, gene)
df_40_60 <- unique(df_40_60)
df_40_60 <- df_40_60 %>% 
  count(gene)
df_40_60 <- df_40_60 %>% 
  mutate("df_40_60" = n)
df_40_60 <- df_40_60[,-2]

df_60 <- df %>% 
  filter(percentage >= 60) %>% 
  select(patient_id, gene)
df_60 <- unique(df_60)
df_60 <- df_60 %>% 
  count(gene)
df_60 <- df_60 %>% 
  mutate("df_60_100" = n)
df_60 <- df_60[,-2]

#put all data frames into list
MAF_Combo <- list(df_all, df_0.1, df_0.1_0.2, df_0.2_0.4, df_0.4_0.6, df_0.6_0.8, df_1, df_1_5, df_5_10, df_10_40, df_40_60, df_60)

#merge all data frames together
MAF_Combo <- MAF_Combo %>% reduce(full_join, by='gene')

MAF_Combo[is.na(MAF_Combo)] <- 0


#Divide by total number of each mutation
df_0.1_Fraction <- MAF_Combo$df_0.1/MAF_Combo$Total
df_0.1_0.2_Fraction <- MAF_Combo$df_0.1_0.2/MAF_Combo$Total
df_0.2_0.4_Fraction <- MAF_Combo$df_0.2_0.4/MAF_Combo$Total
df_0.4_0.6_Fraction <- MAF_Combo$df_0.4_0.6/MAF_Combo$Total
df_0.6_0.8_Fraction <- MAF_Combo$df_0.6_0.8/MAF_Combo$Total
df_0.8_1_Fraction <- MAF_Combo$df_0.8_1/MAF_Combo$Total
df_1_5_Fraction <- MAF_Combo$df_1_5/MAF_Combo$Total
df_5_10_Fraction <- MAF_Combo$df_5_10/MAF_Combo$Total
df_10_40_Fraction <- MAF_Combo$df_10_40/MAF_Combo$Total
df_40_60_Fraction <- MAF_Combo$df_40_60/MAF_Combo$Total
df_60_100_Fraction <- MAF_Combo$df_60_100/MAF_Combo$Total

MAF_Combo$df_0.1_Fraction <- df_0.1_Fraction
MAF_Combo$df_0.1_0.2_Fraction <- df_0.1_0.2_Fraction
MAF_Combo$df_0.2_0.4_Fraction <- df_0.2_0.4_Fraction
MAF_Combo$df_0.4_0.6_Fraction <- df_0.4_0.6_Fraction
MAF_Combo$df_0.6_0.8_Fraction <- df_0.6_0.8_Fraction
MAF_Combo$df_0.8_1_Fraction <- df_0.8_1_Fraction
MAF_Combo$df_1_5_Fraction <- df_1_5_Fraction
MAF_Combo$df_5_10_Fraction <- df_5_10_Fraction
MAF_Combo$df_10_40_Fraction <- df_10_40_Fraction
MAF_Combo$df_40_60_Fraction <- df_40_60_Fraction
MAF_Combo$df_60_100_Fraction <- df_60_100_Fraction

MAF_Combo_Fraction <- MAF_Combo[ -c(2:13) ]

MAF_matrix <- as.matrix(MAF_Combo_Fraction[,-1])

rownames(MAF_matrix) <- c("AKT1","ALK", "APC", "AR", "ARAF", "ARID1A", "ATM", "BRAF", "BRCA1", "BRCA2", "CCND1", "CCND2", "CCNE1", "CDH1", "CDK12", "CDK4", "CDK6", "CDKN2A", "CTNNB1", "DDR2", "EGFR", "ERBB2", "ESR1", "EZH2", "FBXW7", "FGFR1", "FGFR2", "FGFR3", "GATA3", "GNA11", "GNAQ", "GNAS", "HNF1A", "HRAS", "IDH1", "IDH2", "JAK2", "JAK3", "KIT", "KRAS", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "MET", "MLH1", "MPL", "MTOR", "MYC", "NF1", "NFE2L2", "NOTCH1", "NPM1", "NRAS", "NTRK1", "NTRK3", "PDGFRA", "PIK3CA", "PTEN", "PTPN11", "RAF1", "RB1", "RET", "RHEB", "RHOA", "RIT1", "ROS1", "SMAD4", "SMO", "STK11", "TERT", "TP53", "TSC1", "VHL")

colnames(MAF_matrix) <- c("<0.1%", "0.1-0.2%", "0.2-0.4%", "0.4-0.6%", "0.6-0.8%", "0.8-1%", "1-5%", "5-10%", "10-40%", "40-60%", ">60%")

#Reorder matrix based on decreased percentage of patients with mutations

df_all <- df_all[order(-df_all$Total),]
MAF_matrix <- MAF_matrix[df_all$gene, ]

#Make a heat map of all mutations
#Heat map - all 1st Guardant
library("pheatmap")

pheatmap(MAF_matrix, cluster_rows = FALSE, cluster_cols = FALSE, cellwidth=5, cellheight=3, fontsize = 3, colorRampPalette(c("white", "darkblue"))(50), main = "All Patients - 1st cfDNA")

##Heat map - make a heat map of copy number changes

df_first <- Guardant_Table[grepl("_01", Guardant_Table$test_id),]
df_first <- unique(df_first)

df_CN <- df_first %>% 
  filter(mut_type == "amp") %>% 
  filter(copynumber>0) %>% 
  select(patient_id, gene, mut_type, copynumber)
df_CN <- unique(df_CN)

df_CN_total <- df_CN %>% 
  select(patient_id, gene)
df_CN_total <- unique(df_CN_total)
df_CN_total <- df_CN %>% 
  count(gene) %>% 
  mutate("Total" = n)
df_CN_total <- df_CN_total[,-2]


Amp_1_2 <- df_CN %>% 
  filter(copynumber >= 1) %>% 
  filter(copynumber < 2) %>% 
  select(patient_id, gene)
Amp_1_2 <- unique(Amp_1_2)
Amp_1_2 <- Amp_1_2 %>% 
  count(gene)
Amp_1_2 <- Amp_1_2 %>% 
  mutate("CN_1_2" = n)
Amp_1_2 <- Amp_1_2[,-2]

Amp_2_3 <- df_CN %>% 
  filter(copynumber>=2) %>% 
  filter(copynumber<3) %>% 
  select(patient_id, gene)
Amp_2_3 <- unique(Amp_2_3)
Amp_2_3 <- Amp_2_3 %>% 
  count(gene)
Amp_2_3 <- Amp_2_3 %>% 
  mutate("CN_2_3" = n)
Amp_2_3 <- Amp_2_3[,-2]

Amp_3_4 <- df_CN %>% 
  filter(copynumber>=3) %>% 
  filter(copynumber<4) %>% 
  select(patient_id, gene)
Amp_3_4 <- unique(Amp_3_4)
Amp_3_4 <- Amp_3_4 %>% 
  count(gene)
Amp_3_4 <- Amp_3_4 %>% 
  mutate("CN_3_4" = n)
Amp_3_4 <- Amp_3_4[,-2]

Amp_4_5 <- df_CN %>% 
  filter(copynumber>=4) %>% 
  filter(copynumber<5) %>% 
  select(patient_id, gene)
Amp_4_5 <- unique(Amp_4_5)
Amp_4_5 <- Amp_4_5 %>% 
  count(gene)
Amp_4_5 <- Amp_4_5 %>% 
  mutate("CN_4_5" = n)
Amp_4_5 <- Amp_4_5[,-2]

Amp_5_10 <- df_CN %>% 
  filter(copynumber>=5) %>% 
  filter(copynumber<10) %>% 
  select(patient_id, gene)
Amp_5_10 <- unique(Amp_5_10)
Amp_5_10 <- Amp_5_10 %>% 
  count(gene)
Amp_5_10 <- Amp_5_10 %>% 
  mutate("CN_5_10" = n)
Amp_5_10 <- Amp_5_10[,-2]

Amp_10_50 <- df_CN %>% 
  filter(copynumber>=10) %>% 
  filter(copynumber<50) %>% 
  select(patient_id, gene)
Amp_10_50 <- unique(Amp_10_50)
Amp_10_50 <- Amp_10_50 %>% 
  count(gene)
Amp_10_50 <- Amp_10_50 %>% 
  mutate("CN_10_50" = n)
Amp_10_50 <- Amp_10_50[,-2]

Amp_50_200 <- df_CN %>% 
  filter(copynumber>=50) %>% 
  select(patient_id, gene)
Amp_50_200 <- unique(Amp_50_200)
Amp_50_200 <- Amp_50_200 %>% 
  count(gene)
Amp_50_200 <- Amp_50_200 %>% 
  mutate("CN_50_190" = n)
Amp_50_200 <- Amp_50_200[,-2]

CN_Combo <- list(df_CN_total, Amp_1_2, Amp_2_3, Amp_3_4, Amp_4_5, Amp_5_10, Amp_10_50, Amp_50_200)

CN_Combo <- CN_Combo %>% reduce(full_join, by='gene')

CN_Combo[is.na(CN_Combo)] <- 0

CN_1_2_Fraction <- CN_Combo$CN_1_2/CN_Combo$Total
CN_2_3_Fraction <- CN_Combo$CN_2_3/CN_Combo$Total
CN_3_4_Fraction <- CN_Combo$CN_3_4/CN_Combo$Total
CN_4_5_Fraction <- CN_Combo$CN_4_5/CN_Combo$Total
CN_5_10_Fraction <- CN_Combo$CN_5_10/CN_Combo$Total
CN_10_50_Fraction <- CN_Combo$CN_10_50/CN_Combo$Total
CN_50_190_Fraction <- CN_Combo$CN_50_190/CN_Combo$Total

CN_Combo$CN_1_2_Fraction <- CN_1_2_Fraction
CN_Combo$CN_2_3_Fraction <- CN_2_3_Fraction
CN_Combo$CN_3_4_Fraction <- CN_3_4_Fraction
CN_Combo$CN_4_5_Fraction <- CN_4_5_Fraction
CN_Combo$CN_5_10_Fraction <- CN_5_10_Fraction
CN_Combo$CN_10_50_Fraction <- CN_10_50_Fraction
CN_Combo$CN_50_190_Fraction <- CN_50_190_Fraction

CN_Combo_Fraction <- CN_Combo[ -c(2:9) ]

CN_matrix <- as.matrix(CN_Combo_Fraction[,-1])

rownames(CN_matrix) <- c("AR", "BRAF", "CCND1", "CCND2", "CCNE1", "CDK4", "CDK6", "EGFR", "ERBB2", "FGFR1", "FGFR2", "KIT", "KRAS", "MET", "MYC", "PDGFRA", "PIK3CA", "RAF1")

colnames(CN_matrix) <- c("1-2", "2-3", "3-4", "4-5", "5-10%", "10-50%", "50-190%")

#Reorder matrix based on decreased percentage of patients with amplifications

df_CN_total <- df_CN_total[order(-df_CN_total$Total),]
CN_matrix <- CN_matrix[df_CN_total$gene, ]

pheatmap(CN_matrix, cluster_rows = FALSE, cluster_cols = FALSE, cellwidth=5, cellheight=3, fontsize = 3, colorRampPalette(c("white", "darkorange4"))(50))
