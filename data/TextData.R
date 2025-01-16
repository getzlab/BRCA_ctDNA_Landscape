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

##Text Values##

#Determine number of samples with detectable variants and specific variants
length(unique(df$test_id))

df_SNV <- df %>% 
  filter(mut_type == "missense" | mut_type == "nonsense" | mut_type == "synonymous" | mut_type == "splice alt" | mut_type == "VUS" )
df_SNV <- unique(df_SNV)
length(unique(df_SNV$alteration))

df_SNV_nonsyn <- df_SNV %>% 
  filter(mut_type != "synonymous")
df_SNV_nonsyn <- unique(df_SNV_nonsyn)
length(unique(df_SNV_nonsyn$alteration))

df_indel <- df %>% 
  filter(mut_type == "indel")
df_indel <- unique(df_indel)
length(unique(df_indel$cdna))

df_amp <- df %>% 
  filter(mut_type == "amp")
df_amp <- unique(df_amp)
length(unique(df_amp$gene))

df_fusion <- df %>% 
  filter(mut_type == "fusion")
df_fusion <- unique(df_fusion)
length(unique(df_fusion$alteration))

#Determine number of male and female patients
female <- df %>% 
  filter(patientgender == "Female") %>% 
  select(patient_id)
length(unique(female$patient_id))

#Calculate age statistics
df$patientage <- as.numeric(df$patientage)
summary(df$patientage)

#Determine number of serial samples
Samples_01 <- df[grepl("_01", df$test_id),]
length(unique(Samples_01$patient_id)) #n = 11,313 (11,456 total)

Samples_02 <- df[grepl("_02", df$test_id),]
length(unique(Samples_02$patient_id)) #n = 1262

Samples_03 <- df[grepl("_03", df$test_id),]
length(unique(Samples_03$patient_id)) #n = 450

Samples_04 <- df[grepl("_04", df$test_id),]
length(unique(Samples_04$patient_id)) #n = 241

Samples_05 <- df[grepl("_05", df$test_id),]
length(unique(Samples_05$patient_id)) #n = 141

Samples_06 <- df[grepl("_06", df$test_id),]
length(unique(Samples_06$patient_id)) #n = 89

Samples_07 <- df[grepl("_07", df$test_id),]
length(unique(Samples_07$patient_id)) #n = 54

Samples_08 <- df[grepl("_08", df$test_id),]
length(unique(Samples_08$patient_id)) #n = 42

Samples_09 <- df[grepl("_09", df$test_id),]
length(unique(Samples_09$patient_id)) #n = 29

Samples_10 <- df[grepl("_10", df$test_id),]
length(unique(Samples_10$patient_id)) #n = 17

Samples_11 <- df[grepl("_11", df$test_id),]
length(unique(Samples_11$patient_id)) #n = 14

Samples_12 <- df[grepl("_12", df$test_id),]
length(unique(Samples_12$patient_id)) #n = 11

Samples_13 <- df[grepl("_13", df$test_id),]
length(unique(Samples_13$patient_id)) #n = 6

Samples_14 <- df[grepl("_14", df$test_id),]
length(unique(Samples_14$patient_id)) #n = 2

Samples_15 <- df[grepl("_15", df$test_id),]
length(unique(Samples_15$patient_id)) #n = 1

#Determine top mutation, median MAF and IQR for 1st timepoint, 2nd, and 3rd (no increase in abundace)
TopMAF <- Samples_01 %>% 
  filter(mut_type != "amp") %>% 
  filter(mut_type != "fusion") %>% 
  filter(mut_type != "VUS") %>% 
  select(patient_id, gene, mut_type, percentage)
summary(TopMAF$percentage)
TopMAF <- TopMAF %>%
  arrange(desc(percentage)) %>% 
  group_by(patient_id) %>%
  slice(1)
summary(TopMAF$percentage)

TopMAF_2 <- Samples_02 %>% 
  filter(mut_type != "amp") %>% 
  filter(mut_type != "fusion") %>% 
  filter(mut_type != "VUS") %>% 
  select(patient_id, test_id, gene, mut_type, percentage)
summary(TopMAF_2$percentage)
TopMAF_2 <- TopMAF_2 %>%                          
  arrange(desc(percentage)) %>% 
  group_by(patient_id) %>%
  slice(1)
summary(TopMAF_2$percentage)

TopMAF_3 <- Samples_03 %>% 
  filter(mut_type != "amp") %>% 
  filter(mut_type != "fusion") %>% 
  filter(mut_type != "VUS") %>% 
  select(patient_id, gene, mut_type, percentage)
summary(TopMAF_3$percentage)
TopMAF_3 <- TopMAF_3 %>%                          
  arrange(desc(percentage)) %>% 
  group_by(patient_id) %>%
  slice(1)
summary(TopMAF_3$percentage)


#Determine average alterations per 1st sample and IQR
#Determine number of each type of variant in first samples [select only first samples and create dataframe called "df_1") and separately assess the combined serial samples in the larger dataframe ("df") with no repetition of same gene]
df_1 <- df[grepl("_01", df$test_id),]
length(unique(df_1$patient_id))

NumberMutations <- df_1 %>% 
  select(patient_id, gene, mut_type, alteration)
NumberMutations <- unique(NumberMutations)
sapply(NumberMutations, n_distinct)

NumberMutations <- NumberMutations %>%
  distinct(patient_id, alteration) %>%
  group_by(patient_id) %>%
  summarize("comutations" = n())
summary(NumberMutations)

#Determine number of alterations <1% AF
df_low <- df_1 %>% 
  filter(mut_type != "amp") %>% 
  filter(percentage<1) %>% 
  select(patient_id, gene, mut_type, alteration)
df_low <- unique(df_low)

#Determine percentage <1%
df_allmut <- df_1 %>% 
  filter(mut_type != "amp") %>% 
  select(patient_id, gene, mut_type, alteration)
df_allmut <- unique(df_allmut)

#length alterations <1% over all alterations (fraction of df_low/df)
30958/49675

#Determine most frequently mutated genes
df_mut <- df %>% 
  filter(mut_type != "amp") %>% 
  filter(mut_type != "fusion") %>% 
  filter(mut_type != "indel") %>% 
  select(patient_id, gene)
df_mut <- unique(df_mut)
length(unique(df_mut$patient_id))
df_mut_count <- df_mut %>% 
  count(gene)
df_mut_count <- df_mut_count %>% 
  mutate("%_mut" = (n/12827) * 100)
write_xlsx(df_mut_count,"/Insert_path_for_file//Guardant_Mutations_all.xlsx")

#Calculate number of patients with gene amplifications in condensed cohort
df_amp <- df %>%
  filter(mut_type == "amp") %>% 
  select("patient_id", "gene")
df_amp <- unique(df_amp)
df_amp <- df_amp %>% 
  count(gene)
df_amp_count <- df_amp %>% 
  mutate("%_mut" = (n/12827) * 100)
write_xlsx(df_amp_count,"/Insert_path_for_file//Guardant_amplifications_all.xlsx")

#Calculate number of patients with gene fusions in condensed cohort
df_fusion <- df %>%
  filter(mut_type == "fusion") %>% 
  select("patient_id", "gene")
df_fusion <- unique(df_fusion)
df_fusion <- df_fusion %>% 
  count(gene)
df_fusion_count <- df_fusion %>% 
  mutate("%_mut" = (n/12827) * 100)
write_xlsx(df_fusion_count,"/Insert_path_for_file//Guardant_fusions_all.xlsx")

#Calculate number of patients with gene indels in condensed cohort
df_indel <- df %>%
  filter(mut_type == "indel") %>% 
  select("patient_id", "gene")
df_indel <- unique(df_indel)
df_indel <- df_indel %>% 
  count(gene)
df_indel_count <- df_indel %>% 
  mutate("%_mut" = (n/12827) * 100)
write_xlsx(df_indel_count,"/Insert_path_for_file//Guardant_indels_all.xlsx")

#Calculate number of patients with gene amplifications >4 in "Pt_Amp" in condensed cohort
Pt_Amp_4 <- df %>%
  filter(mut_type == "amp") %>% 
  filter(copynumber>4) %>% 
  select("patient_id", "gene")
Pt_Amp_4 <- unique(Pt_Amp_4)
length(unique(Pt_Amp_4$patient_id))
