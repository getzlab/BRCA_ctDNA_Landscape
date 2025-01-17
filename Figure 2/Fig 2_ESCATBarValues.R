#load libraries
library(tidyverse)
library(writexl)

#Load cfDNA database
df <- guardant_table
df <- unique(df)

#Merge final subtype annotation file for final subtype annotations
df <- left_join(df, annotate, by = "patient_id")
length(unique(df$patient_id))
length(unique(df$test_id))

#Merge with oncogenic annotations (performed separately via OncoKB API - https://github.com/oncokb/oncokb-annotator)
df <- left_join(df, guardant_table_OncoKB_output, relationship = "many-to-many")
df <- unique(df)
length(unique(df$patient_id))

#Determine number of patients in each subtype based on the cfDNA-based classifier
ER <- df %>% 
  filter(ER_final == "ER+") %>% 
  select(patient_id)
ER <- unique(ER)
length(unique(ER$patient_id))

HER2 <- df %>% 
  filter(HER2_final == "HER2+") %>% 
  select(patient_id)
length(unique(HER2$patient_id))

TN <- df %>% 
  filter(TN_final == "TN") %>% 
  select(patient_id)
length(unique(TN$patient_id))

##Count patients with Level 1 (ESCAT) and non-actionable gene mutations and remove patients with ESCAT mutations in the same gene so these are displayed in separate table columns

#ESCAT 1: gBRCA1, gBRCA2, PIK3CA, ESR1, NTRK fusion, BRAF V600E, RET fusion, AKT1, PTEN (note: AKT1 and PTEN added separately given recently Level 1 after capivasertib approval)
ESCAT1_level1 <- df %>% 
  filter(HIGHEST_LEVEL == "LEVEL_1") %>% 
  filter(gene == "PIK3CA" | gene == "ESR1" | gene == "NTRK1" | gene == "NTRK3" | gene == "BRAF" | gene == "RET") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
ESCAT1_level1 <- unique(ESCAT1_level1)
length(unique(ESCAT1_level1$patient_id))

ESCAT1_cap <- df %>% 
  filter(gene == "AKT1" | gene == "PTEN") %>% 
  filter(ONCOGENIC == "Oncogenic" | ONCOGENIC == "Likely Oncogenic") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)

ESCAT1_level2 <- df %>% 
  filter(HIGHEST_LEVEL == "LEVEL_2") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
ESCAT1_level2 <- unique(ESCAT1_level2)
length(unique(ESCAT1_level2$patient_id))

ESCAT1_gBRCA <- df %>% 
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  filter(percentage >= 40) %>% 
  filter(ONCOGENIC == "Oncogenic" | ONCOGENIC == "Likely Oncogenic") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
ESCAT1_gBRCA <- unique(ESCAT1_gBRCA)
length(unique(ESCAT1_gBRCA$patient_id))

ESCAT1_ERBB2amp <- df %>% 
  filter(gene == "ERBB2" & mut_type == "amp") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
ESCAT1_ERBB2amp <- unique(ESCAT1_ERBB2amp)
length(unique(ESCAT1_ERBB2amp$patient_id))

ESCAT1 <- rbind(ESCAT1_level1, ESCAT1_level2, ESCAT1_gBRCA, ESCAT1_ERBB2amp, ESCAT1_cap)
ESCAT1 <- unique(ESCAT1)
length(unique(ESCAT1$patient_id))

ESCAT1["ESCAT"] <- "ESCAT I"

#ESCAT 2: Determine the number of ESCAT II mutations: ERBB2mut, sBRCA1, sBRCA2

ESCAT2_ERBB2mut <- df %>% 
  filter(gene == "ERBB2" & LEVEL_3A == "Neratinib") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
ESCAT2_ERBB2mut <- unique(ESCAT2_ERBB2mut)
length(unique(ESCAT2_ERBB2mut$patient_id))

ESCAT2_sBRCA <- df %>% 
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  filter(percentage < 40) %>% 
  filter(ONCOGENIC == "Oncogenic" | ONCOGENIC == "Likely Oncogenic") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
ESCAT2_sBRCA <- unique(ESCAT2_sBRCA)
length(unique(ESCAT2_sBRCA$patient_id))

ESCAT2 <- rbind(ESCAT2_ERBB2mut, ESCAT2_sBRCA)
ESCAT2 <- unique(ESCAT2)
length(unique(ESCAT2$patient_id))

ESCAT2["ESCAT"] <- "ESCAT II"

#ESCAT 3: PDGFRA, ALK fusion, IDH1, IDH2, KIT, ROS1 fusion, FGFR1/2/3, EGFR, KRAS G12C
ESCAT3_drugs <- df %>% 
  filter(LEVEL_3B == "Avapritinib,Imatinib" | LEVEL_3B == "Crizotinib,Ceritinib,Alectinib,Brigatinib,Lorlatinib" | LEVEL_3B == "Ivosidenib,Vorasidenib" | LEVEL_3B == "Enasidenib,Vorasidenib" | LEVEL_3B == "Imatinib,Sunitinib,Regorafenib,Ripretinib" | LEVEL_3B == "Erdafitinib,Pemigatinib,RLY-4008,Futibatinib" | LEVEL_3B == "Osimertinib" | LEVEL_3B == "Afatinib,Patritumab Deruxtecan,Osimertinib" | LEVEL_3B == "Erlotinib,Patritumab Deruxtecan,Erlotinib+Ramucirumab,Afatinib,Gefitinib,Osimertinib,Dacomitinib") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
ESCAT3_drugs <- unique(ESCAT3_drugs)
length(unique(ESCAT3_drugs$patient_id))

ESCAT3_KRAS <- df %>% 
  filter(gene == "KRAS" & alteration == "G12C") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
ESCAT3_KRAS <- unique(ESCAT3_KRAS)
length(unique(ESCAT3_KRAS$patient_id))

ESCAT3_ROS1 <- df %>% 
  filter(gene == "ROS" & mut_type == "fusion") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
ESCAT3_ROS1 <- unique(ESCAT3_ROS1)
length(unique(ESCAT3_ROS1$patient_id))

ESCAT3 <- rbind(ESCAT3_drugs, ESCAT3_KRAS)

ESCAT3["ESCAT"] <- "ESCAT III"

##Determine number of ESCAT 4 mutations
#All: ARID1A, ATM, ATR, CDH1, MYC, RUNX1, SF3B1, TP53
ESCAT4 <- df %>%
  filter(ONCOGENIC == "Oncogenic" | ONCOGENIC == "Likely Oncogenic") %>% 
  filter(mut_type != "synonymous") %>% 
  filter( gene == "ARID1A" | gene == "ATR" | gene == "ATM" | gene == "CDH1" | gene == "TP53" | gene == "SF3B1" | gene == "RUNX1" | gene == "MYC" | gene == "NF1" | gene == "GATA3") %>%
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)

ESCAT4 <- unique(ESCAT4)
ESCAT4["ESCAT"] <- "ESCAT IV"

#Combine ESCAT 
ESCAT_All <- rbind(ESCAT1, ESCAT2, ESCAT3, ESCAT4)
ESCAT_All <- unique(ESCAT_All)
ESCAT_consolidated <- ESCAT_All %>% 
  count(gene)

#Determine number of patients with ESCAT 1/2
ESCAT1_2 <- ESCAT_All %>% 
  filter(ESCAT == "ESCAT I" | ESCAT == "ESCAT II") %>% 
  select(patient_id, ESCAT, ER_final, HER2_final, TN_final)
ESCAT1_2 <- unique(ESCAT1_2)
length(unique(ESCAT1_2$patient_id))

#Merge ESCAT Annotations with df
df <- left_join(df, ESCAT_All, relationship =
  "many-to-many")
df <- unique(df)

###########################################

##Generate numbers to export to excel file for stacked bar graphs (only one ESCAT per pt - prioritize highest ESCAT and remove lower ESCAT values to prevent duplication) - First analysis does not account for subtype - denominator will be 11,456 given accounting only for patients with detectable cfDNA and used in cfDNA-based classifier

#All Patients
All_Bar_I <- ESCAT_All %>% 
  filter(ESCAT == "ESCAT I") %>% 
  select(patient_id, ESCAT)
All_Bar_I <- unique(All_Bar_I)  

All_No_ESCAT1 <- ESCAT_All[!(ESCAT_All$patient_id %in% All_Bar_I$patient_id),]

All_Bar_II <- All_No_ESCAT1 %>% 
  filter(ESCAT == "ESCAT II") %>% 
  select(patient_id, ESCAT)
All_Bar_II <- unique(All_Bar_II)

All_No_ESCAT_1_2 <- All_No_ESCAT1[!(All_No_ESCAT1$patient_id %in% All_Bar_II$patient_id),]

All_Bar_III <- All_No_ESCAT_1_2 %>% 
  filter(ESCAT == "ESCAT III") %>% 
  select(patient_id, ESCAT)
All_Bar_III <- unique(All_Bar_III)

All_No_ESCAT_1_3 <- All_No_ESCAT_1_2[!(All_No_ESCAT_1_2$patient_id %in% All_Bar_III$patient_id),]

All_Bar_IV <- All_No_ESCAT_1_3 %>% 
  filter(ESCAT == "ESCAT IV") %>% 
  select(patient_id, ESCAT)
All_Bar_IV <- unique(All_Bar_IV)

All_Bar_None <- All_No_ESCAT_1_3[!(All_No_ESCAT_1_3$patient_id %in% All_Bar_IV$patient_id),]
All_Bar_None <- unique(All_Bar_None)
All_Bar_None <- All_Bar_None %>% 
  filter(ESCAT == "None") %>% 
  select(patient_id, ESCAT)

Bar_all <- rbind(All_Bar_I, All_Bar_II, All_Bar_III, All_Bar_IV, All_Bar_None)
Bar_all <- unique(Bar_all)
Bar_all <- Bar_all %>% 
  count(ESCAT)
Bar_all <- Bar_all %>% 
  mutate("%" = n/11456)
write_xlsx(Bar_all,"//Path_for_file//Bar_all_2023.xlsx")

##Create excel files to create stacked bar graphs by Subtype

#ER Bar Chart 
ER_ESCAT <- ESCAT_All %>% 
  filter(ER_final == "ER+") %>% 
  select(patient_id, ESCAT)
length(unique(ER_ESCAT$patient_id))

ER_Bar_I <- ER_ESCAT %>% 
  filter(ESCAT == "ESCAT I") %>% 
  select(patient_id, ESCAT)
ER_Bar_I <- unique(ER_Bar_I)  

No_ER_ESCAT1 <- ER_ESCAT[!(ER_ESCAT$patient_id %in% ER_Bar_I$patient_id),]

ER_Bar_II <- No_ER_ESCAT1 %>% 
  filter(ESCAT == "ESCAT II") %>% 
  select(patient_id, ESCAT)
ER_Bar_II <- unique(ER_Bar_II)

No_ER_ESCAT_1_2 <- No_ER_ESCAT1[!(No_ER_ESCAT1$patient_id %in% ER_Bar_II$patient_id),]

ER_Bar_III <- No_ER_ESCAT_1_2 %>% 
  filter(ESCAT == "ESCAT III") %>% 
  select(patient_id, ESCAT)
ER_Bar_III <- unique(ER_Bar_III)

No_ER_ESCAT_1_2_3 <- No_ER_ESCAT_1_2[!(No_ER_ESCAT_1_2$patient_id %in% ER_Bar_III$patient_id),]

ER_Bar_IV <- No_ER_ESCAT_1_2_3 %>% 
  filter(ESCAT == "ESCAT IV") %>% 
  select(patient_id, ESCAT)
ER_Bar_IV <- unique(ER_Bar_IV)

ER_Bar_None <- No_ER_ESCAT_1_2_3[!(No_ER_ESCAT_1_2_3$patient_id %in% ER_Bar_IV$patient_id),]
ER_Bar_None <- ER_Bar_None %>% 
  select(patient_id, ESCAT)
ER_Bar_None <- unique(ER_Bar_None)

All_ER_Bar <- rbind(ER_Bar_I, ER_Bar_II, ER_Bar_III, ER_Bar_IV, ER_Bar_None)
All_ER_Bar <- unique(All_ER_Bar)
All_ER_Bar_short <- All_ER_Bar %>% 
  count(ESCAT)
All_ER_Bar_short <- All_ER_Bar_short %>% 
  mutate("% mut" = n/4790)
write_xlsx(All_ER_Bar_short,"//Path_for_file//ER_Bar_all_2023.xlsx")

#HER2 Bar Chart
HER2_ESCAT_all <- ESCAT_All %>% 
  filter(HER2_final == "HER2+") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final, ESCAT)
length(unique(HER2_ESCAT_all$patient_id))

#Remove ERBB2 amplification from ESCAT I to not confuse results given this is 

HER2_ESCAT <- ESCAT_All %>% 
  filter(HER2_final == "HER2+") %>% 
  select(patient_id, ESCAT, gene, alteration)
length(unique(HER2_ESCAT$patient_id))
HER2_ESCAT <- HER2_ESCAT %>% 
  filter(alteration != "AMP") %>% 
  select(patient_id, ESCAT, gene, alteration)
length(unique(HER2_ESCAT$patient_id))

HER2_Bar_I <- HER2_ESCAT %>% 
  filter(ESCAT == "ESCAT I" ) %>% 
  select(patient_id, ESCAT)
HER2_Bar_I <- unique(HER2_Bar_I)  

No_HER2_ESCAT_1 <- HER2_ESCAT[!(HER2_ESCAT$patient_id %in% HER2_Bar_I$patient_id),]

HER2_Bar_II <- No_HER2_ESCAT_1 %>% 
  filter(ESCAT == "ESCAT II") %>% 
  select(patient_id, ESCAT)
HER2_Bar_II <- unique(HER2_Bar_II)

No_HER2_ESCAT__1_2 <- No_HER2_ESCAT_1[!(No_HER2_ESCAT_1$patient_id %in% HER2_Bar_II$patient_id),]

HER2_Bar_III <- No_HER2_ESCAT__1_2 %>% 
  filter(ESCAT == "ESCAT III") %>% 
  select(patient_id, ESCAT)
HER2_Bar_III <- unique(HER2_Bar_III)

No_HER2_ESCAT__1_2_3 <- No_HER2_ESCAT__1_2[!(No_HER2_ESCAT__1_2$patient_id %in% HER2_Bar_III$patient_id),]

HER2_Bar_IV <- No_HER2_ESCAT__1_2_3 %>% 
  filter(ESCAT == "ESCAT IV") %>% 
  select(patient_id, ESCAT)
HER2_Bar_IV <- unique(HER2_Bar_IV)

HER2_Bar_None <- No_HER2_ESCAT__1_2_3[!(No_HER2_ESCAT__1_2_3$patient_id %in% HER2_Bar_IV$patient_id),]
HER2_Bar_None <- HER2_Bar_None %>% 
  select(patient_id, ESCAT)
HER2_Bar_None <- unique(HER2_Bar_None)
HER2_Bar_None["ESCAT"] <- "None"

All_HER2_Bar <- rbind(HER2_Bar_I, HER2_Bar_II, HER2_Bar_III, HER2_Bar_IV, HER2_Bar_None)
All_HER2_Bar <- unique(All_HER2_Bar)
All_HER2_Bar_short <- All_HER2_Bar %>% 
  count(ESCAT)
All_HER2_Bar_short <- All_HER2_Bar_short %>% 
  mutate("% mut" = n/491)
write_xlsx(All_HER2_Bar,"//Path_to_file//HER2_Bar_all_2023.xlsx")

#TN Bar Chart
TN_ESCAT <- ESCAT_All %>% 
  filter(TN_final == "TN") %>% 
  select(patient_id, ESCAT)
length(unique(TN_ESCAT $patient_id))

TN_Bar_I <- TN_ESCAT %>% 
  filter(ESCAT == "ESCAT I") %>% 
  select(patient_id, ESCAT)
TN_Bar_I <- unique(TN_Bar_I)  
TN_Bar_I <- TN_Bar_I

No_TN_ESCAT1 <- TN_ESCAT[!(TN_ESCAT$patient_id %in% TN_Bar_I$patient_id),]

TN_Bar_II <- No_TN_ESCAT1 %>% 
  filter(ESCAT == "ESCAT II") %>% 
  select(patient_id, ESCAT)
TN_Bar_II <- unique(TN_Bar_II)

No_TN_ESCAT_1_2 <- No_TN_ESCAT1[!(No_TN_ESCAT1$patient_id %in% TN_Bar_II$patient_id),]

TN_Bar_III <- No_TN_ESCAT_1_2 %>% 
  filter(ESCAT == "ESCAT III") %>% 
  select(patient_id, ESCAT)
TN_Bar_III <- unique(TN_Bar_III)

No_TN_ESCAT_1_2_3 <- No_TN_ESCAT_1_2[!(No_TN_ESCAT_1_2$patient_id %in% TN_Bar_III$patient_id),]

TN_Bar_IV <- No_TN_ESCAT_1_2_3 %>%
  filter(ESCAT == "ESCAT IV") %>% 
  select(patient_id, ESCAT)
TN_Bar_IV <- unique(TN_Bar_IV)

TN_Bar_None <- No_TN_ESCAT_1_2_3[!(No_TN_ESCAT_1_2_3$patient_id %in% TN_Bar_IV$patient_id),]
TN_Bar_None <- TN_Bar_None %>% 
  select(patient_id, ESCAT)
TN_Bar_None <- unique(TN_Bar_None)
TN_Bar_None["ESCAT"] <- "None"

All_TN_Bar <- rbind(TN_Bar_I, TN_Bar_II, TN_Bar_III, TN_Bar_IV, TN_Bar_None)
All_TN_Bar <- unique(All_TN_Bar)
All_TN_Bar_short <- All_TN_Bar %>% 
  count(ESCAT)
All_TN_Bar_short <- All_TN_Bar_short %>% 
  mutate("% mut" = n/1431)
write_xlsx(All_TN_Bar,"//Path_to_file//TN_Bar_all_2023.xlsx")

##Obtain values for bar graph columns in excel for top mutations and color coding for ESCAT actionability (performed in excel) - excluding synonymous variants for figure given concern over actionable and potentially actionable variants

#Identify number and percentage of patients with oncogenic or likely oncogenic alterations
Oncogenic_mut <- df %>% 
  filter(ONCOGENIC == "Oncogenic" | ONCOGENIC == "Likely Oncogenic") %>% 
  filter(mut_type != "amp") %>% 
  filter(mut_type != "synonymous") %>% 
  select(patient_id, gene)
Oncogenic_mut <- unique(Oncogenic_mut)
length(unique(Oncogenic_mut$patient_id))
Oncogenic_mut_count <- Oncogenic_mut %>% 
  count(gene)
Oncogenic_mut_count <- Oncogenic_mut_count %>% 
  mutate("Oncogenic variant (%)" = n/11456)
colnames(Oncogenic_mut_count)[which(names(Oncogenic_mut_count) == "gene")] <- "Gene"
colnames(Oncogenic_mut_count)[which(names(Oncogenic_mut_count) == "n")] <- "Oncogenic mut (n)"

ESCAT_mut <- df %>% 
  filter(ESCAT == "ESCAT I" | ESCAT == "ESCAT II" |ESCAT == "ESCAT III" |ESCAT == "ESCAT IV") %>% 
  filter(mut_type != "amp") %>% 
  filter(mut_type != "synonymous") %>% 
  select(patient_id, gene)
ESCAT_mut <- unique(ESCAT_mut)
length(unique(ESCAT_mut$patient_id))
ESCAT_mut_count <- ESCAT_mut %>% 
  count(gene)
ESCAT_mut_count <- ESCAT_mut_count %>% 
  mutate("ESCAT variant (%)" = n/11456)
colnames(ESCAT_mut_count)[which(names(ESCAT_mut_count) == "gene")] <- "Gene"
colnames(ESCAT_mut_count)[which(names(ESCAT_mut_count) == "n")] <- "ESCAT mut (n)"

Amp <- df %>% 
  filter(mut_type == "amp") %>% 
  select(patient_id, gene)
Amp <- unique(Amp)
length(unique(Amp$patient_id))
Amp_count <- Amp %>% 
  count(gene)
Amp_count <- Amp_count %>% 
  mutate("Amp (%)" = n/11456)
colnames(Amp_count)[which(names(Amp_count) == "gene")] <- "Gene"
colnames(Amp_count)[which(names(Amp_count) == "n")] <- "Amp (n)"

Total_mut <- df %>% 
  filter(mut_type != "amp") %>% 
  filter(mut_type != "synonymous") %>% 
  select(patient_id, gene)
Total_mut <- unique(Total_mut)
length(unique(Total_mut$patient_id))
Total_mut <- Total_mut %>% 
  count(gene)
colnames(Total_mut)[which(names(Total_mut) == "gene")] <- "Gene"
colnames(Total_mut)[which(names(Total_mut) == "n")] <- "Total variants"

Table1_ESCAT <- full_join(Oncogenic_mut_count, Amp_count, by = "Gene")

Table1_ESCAT <- full_join(Table1_ESCAT, ESCAT_mut_count, by = "Gene")

Table1_ESCAT <- full_join(Table1_ESCAT, Total_mut, by = "Gene")

write_xlsx(Table1_ESCAT,"//Insert_path_for_file//ESCAT_bar_values.xlsx")

##Determine number of patients with suspected germline BRCA1/2 (MAF≥40%) and somatic BRCA1/2 (<40%) - Table 1
gBRCA1 <- df %>% 
  filter(gene == "BRCA1") %>%
  filter(percentage >= 40) %>% 
  filter(ONCOGENIC == "Oncogenic" | ONCOGENIC == "Likely Oncogenic") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
gBRCA1 <- unique(gBRCA1)
length(unique(gBRCA1$patient_id))

gBRCA2 <- df %>% 
  filter(gene == "BRCA2") %>%
  filter(percentage >= 40) %>% 
  filter(ONCOGENIC == "Oncogenic" | ONCOGENIC == "Likely Oncogenic") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
gBRCA2 <- unique(gBRCA2)
length(unique(gBRCA2$patient_id))

sBRCA1 <- df %>% 
  filter(gene == "BRCA1") %>%
  filter(percentage < 40) %>% 
  filter(ONCOGENIC == "Oncogenic" | ONCOGENIC == "Likely Oncogenic") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
sBRCA1 <- unique(sBRCA1)
length(unique(sBRCA1$patient_id))

sBRCA2 <- df %>% 
  filter(gene == "BRCA2") %>%
  filter(percentage < 40) %>% 
  filter(ONCOGENIC == "Oncogenic" | ONCOGENIC == "Likely Oncogenic") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
sBRCA2 <- unique(sBRCA2)
length(unique(sBRCA2$patient_id))

#Calculate patients without oncogenic variants and ESCAT 1/2 - Table 1

gBRCA1 <- df %>% 
  filter(gene == "BRCA1") %>%
  filter(percentage >= 40) %>% 
  filter(ONCOGENIC != "Oncogenic" & ONCOGENIC != "Likely Oncogenic") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
gBRCA1 <- unique(gBRCA1)
length(unique(gBRCA1$patient_id))

gBRCA2 <- df %>% 
  filter(gene == "BRCA2") %>%
  filter(percentage >= 40) %>% 
  filter(ONCOGENIC != "Oncogenic" & ONCOGENIC != "Likely Oncogenic") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
gBRCA2 <- unique(gBRCA2)
length(unique(gBRCA2$patient_id))

sBRCA1 <- df %>% 
  filter(gene == "BRCA1") %>%
  filter(percentage < 40) %>% 
  filter(ONCOGENIC != "Oncogenic" & ONCOGENIC != "Likely Oncogenic") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
sBRCA1 <- unique(sBRCA1)
length(unique(sBRCA1$patient_id))

sBRCA2 <- df %>% 
  filter(gene == "BRCA2") %>%
  filter(percentage < 40) %>% 
  filter(ONCOGENIC != "Oncogenic" & ONCOGENIC != "Likely Oncogenic") %>% 
  select(patient_id, gene, alteration, ONCOGENIC, ER_final, HER2_final, TN_final)
sBRCA2 <- unique(sBRCA2)
length(unique(sBRCA2$patient_id))
