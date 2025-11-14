##### P1-Cluster-NEO-KD 
##### Author: Keyong Deng
##### Department of Clinical Epidemiology, LUMC

suppressMessages(library(Hmisc)) # describe function
suppressMessages(library(foreign)) # data load
suppressMessages(library(dplyr)) # %>% pipeline
suppressMessages(library(data.table)) # Fast Read data, fread(), fwrite()
suppressMessages(library(stringr)) # string manuiplate
suppressMessages(library(ggplot2)) # Visualization
suppressMessages(library(e1071)) # Data skewness 
suppressMessages(library(tibble))

# Read in data
rawdata <- read.dta(file = "NEO_analyse_PP104_Brainshake_mealresponse_Deng_2024-12-20.dta",
                    convert.factors = T)
head(rawdata)

# Define a function to perform the analysis for a given index, sex-based energy intake adjustment
perform_analysis <- function(index_name, data) {
  
  # Subset the data by sex
  male_data <- subset(data, sex == 1)
  female_data <- subset(data, sex == 2)
  
  # Perform linear regression for males
  lm_male <- lm(get(index_name) ~ kcal, male_data)
  male_data$residuals <- lm_male$residuals
  male_constant <- summary(lm_male)$coefficients[1]
  male_coefficient <- summary(lm_male)$coefficients[2]
  male_energy_mean <- mean(male_data$kcal, na.rm = TRUE)
  male_mean_index <- male_constant + male_coefficient * male_energy_mean
  male_data[[paste0(index_name, "_eadj")]] <- male_data$residuals + male_mean_index
  
  # Perform linear regression for females
  lm_female <- lm(get(index_name) ~ kcal, female_data)
  female_data$residuals <- lm_female$residuals
  female_constant <- summary(lm_female)$coefficients[1]
  female_coefficient <- summary(lm_female)$coefficients[2]
  female_energy_mean <- mean(female_data$kcal, na.rm = TRUE)
  female_mean_index <- female_constant + female_coefficient * female_energy_mean
  female_data[[paste0(index_name, "_eadj")]] <- female_data$residuals + female_mean_index
  
  # Merge adjusted index back into original dataset
  data <- merge(data, rbind(male_data, female_data)[c("ID", paste0(index_name, "_eadj"))], by="ID", all.x = TRUE)
  
  return(data)
}

# ID
rawdata$ID <- c(1:nrow(rawdata))
# kJ to kcal
rawdata$kcal <- rawdata$kJ/4.184 
# sex
rawdata$sex <- ifelse(rawdata$sexe == "man", 1, 2)

# Prepare the necessary dietary scores
diet_scors = c('ID', "sex",'AMED', 'hPDI', 'DII', 'EDIH', 'HEI', 'kcal')
rawdata_diet <- rawdata %>% 
  dplyr::select(all_of(diet_scors))

rawdata_AMED <- rawdata_diet[complete.cases(rawdata_diet['AMED']),]

# Energy adjusted AMED
rawdata_AMED <- perform_analysis("AMED", rawdata_AMED)
rawdata <- merge(rawdata, rawdata_AMED[c("ID", "AMED_eadj")], by="ID", all.x = TRUE)

# Energy adjusted hPDI
rawdata_hPDI <- rawdata_diet[complete.cases(rawdata_diet['hPDI']),]
rawdata_hPDI <- perform_analysis("hPDI", rawdata_hPDI)
rawdata <- merge(rawdata, rawdata_hPDI[c("ID", "hPDI_eadj")], by="ID", all.x = TRUE)

# Energy adjusted DII
rawdata_DII <- rawdata_diet[complete.cases(rawdata_diet['DII']),]
rawdata_DII <- perform_analysis("DII", rawdata_DII)
rawdata <- merge(rawdata, rawdata_DII[c("ID", "DII_eadj")], by="ID", all.x = TRUE)

# Energy adjusted EDIH
rawdata_EDIH <- rawdata_diet[complete.cases(rawdata_diet['EDIH']),]
rawdata_EDIH <- perform_analysis("EDIH", rawdata_EDIH)
rawdata <- merge(rawdata, rawdata_EDIH[c("ID", "EDIH_eadj")], by="ID", all.x = TRUE)

# Energy adjusted HEI
rawdata_HEI <- rawdata_diet[complete.cases(rawdata_diet['HEI']),]
rawdata_HEI <- perform_analysis("HEI", rawdata_HEI)
rawdata <- merge(rawdata, rawdata_HEI[c("ID", "HEI_eadj")], by="ID", all.x = TRUE)

# rename HEI_eadj to AHEI_eadj
rawdata <- rawdata %>% 
  rename(AHEI_eadj = HEI_eadj)

df_update <- rawdata
## Data Cleaning
## 1. exclude prevalent T2DM, without follow-up info, without fasting, without complete meal
table(df_update$diab_prev, useNA = "always") # 592 cases
table(df_update$eind2, useNA = "always") # 94 loss-to-follow-up
table(df_update$nuchter, useNA = "always") # fasting (30 obs)
table(df_update$devmeal, useNA = "always") # complete meal = 6211 obs

df_update <- df_update %>% 
  filter(!diab_prev %in% c("Yes"), # 6079 obs
         !eind2 %in% c("lost-to-follow-up"), # 94 obs  # 5985 obs
         !nuchter %in% c("nee"), # 23 obs without fasting meal ## 5985-5962
         devmeal %in% c("complete meal") # 5962-5552 = 410 obs without complete meal
  ) # 5552 obs
table(df_update$t2d_30d)

# Make a variable (fu_years) for time at risk variable (do this part later 2024-01)
df_update$fu_time <- ifelse(df_update$t2d_inc=="Yes", df_update$diabetes2_date-df_update$visitdd, 
                      df_update$einddatum2-df_update$visitdd)

df_update$fu_years <- df_update$fu_time /365.25 
summary(df_update$fu_years)## Median 6.7 years [5.9-7.9]

## Medications
df_update <- df_update %>%  
  mutate(medication_HeartFailure = ifelse(medC01AA5 == "yes", 1, 0),
         medication_hypertension = ifelse(medC03AA03 == "yes" | medC03BA04 == "yes" | medC07 == "yes" | medC07A == "yes" |
                                            medC07AG == "yes" | medC07B == "yes" | medC07C == "yes" | medC09C == "yes", 
                                          1, 0),
         medication_lipidlower = ifelse(medC10liplow == "yes" | medC10statin == "yes" | medC10fibrate == "yes", 1, 0),
         .before = 1) %>% 
  mutate(medication_all = ifelse(medication_HeartFailure == 1 | medication_hypertension == 1 
                                 | medication_lipidlower == 1, 1, 0)) 

table(df_update$medication_HeartFailure) # 8 Heart Failure
table(df_update$medication_hypertension) # 1138 hypertension
table(df_update$medication_lipidlower) # 633 lipid-lower medication
table(df_update$medication_all) # 1440 all medication
table(df_update$diab_prev, useNA = "always") # 5552 NO Prevalent DM
table(df_update$t2d_30d, useNA = "always") # no cases
table(df_update$t2d_inc, useNA = "always") # 281 obs with incident T2DM

## Select Metabolites 
Metabolite <- df_update %>% 
  select(c(62:519))  ## 229 metabolites per time point

## Exclude those metabolites ratios
w <- str_detect(colnames(Metabolite), "p_")
Met_remove1 <- colnames(Metabolite)[which(w == T)][1:140]

w <- str_detect(colnames(Metabolite), "FA_*")
Met_remove_2 <- colnames(Metabolite)[which(w == T)][15:30]

Met_remove_3 <- c(Met_remove1, Met_remove_2, "ApoBApoA1_1", "ApoBApoA1_3", "TGPG_1", "TGPG_3", "DAGTG_1", "DAGTG_3")
Metabolites_NoRat <- Metabolite %>% 
  select(!all_of(Met_remove_3)) # 296/2 = 148 metabolites

# Replace any value <=0 entries with an NA
Metabolites_NoRat <- Metabolites_NoRat %>% mutate_at(vars(XXLVLDLP_1:Gp_3), ~ifelse(.<=0,NA,.))

# Detect outliers
is_outlier <- function(object, nSD = nSD){
  object <- as.numeric(object)
  v_mean <- mean(object, na.rm = T)
  v_median <- median(object, na.rm = T)
  v_sd <- sd(object, na.rm = T)
  v_IQR <- IQR(object, na.rm = T)
  # r <- which(object > (v_median + nSD * v_IQR) | object < (v_median - nSD * v_IQR))
  r <- which(object > (v_mean + nSD * v_sd) | object < (v_mean - nSD * v_sd))
  return(r)
}
w <- apply(Metabolites_NoRat, 2, function(x) is_outlier(x, nSD = 5)) # list: 296 metabolites, more than 5 SD

# Replace the extreme value with NA
for (i in names(w)) {
  a <- w[[i]]
  if (!is.null(a)){
    Metabolites_NoRat[a, i] <- NA
  }
}

# Sample missing
fmis_sample = apply(Metabolites_NoRat, 1,  function(x) { sum(is.na(x)) / length(x)})
summary(fmis_sample)

# Filter samples
w_sample = which(fmis_sample > 0.4) ## exclude those individuals with more than 40% missing (232 samples)
Metabolites_NoRat <- Metabolites_NoRat[-w_sample,] # 5320 obs for metabolites_NoRatio
df_update <- df_update[-w_sample,] ## exclude those samples (5320 obs remains)

# select new variable, combined with previous dataset
select_vars <- c('weegfactor1', 'sexe', 'leeftijd', 'bmim', "glucose1", "Insuline_r1", "HBA1C", "homa1IR", "Gp_3", 
                 'AMED_eadj', 'hPDI_eadj', 'DII_eadj', 'EDIH_eadj', 'AHEI_eadj')

df_update_2 <- df_update %>% 
  select(any_of(select_vars))

# Combined with cluster information
load("./k_means_clustering.Rdata")

clusters_fasting <- kmeans_cluster_all[[1]] %>% as.data.frame() %>% 
  select(V2); table(clusters_fasting)
clusters_postprandial <- kmeans_cluster_all[[2]] %>% as.data.frame() %>% 
  select(V2); table(clusters_postprandial)
clusters_delta <- kmeans_cluster_all[[3]] %>% as.data.frame() %>% 
  select(V2); table(clusters_delta)

df3_diet <- df_update_2 %>% 
  dplyr::select(AMED_eadj, hPDI_eadj, DII_eadj, EDIH_eadj, AHEI_eadj, glucose1) %>% 
  dplyr::mutate(
         cluster_fasting = clusters_fasting$V2,
         cluster_postprandial = clusters_postprandial$V2,
         cluster_delta = clusters_delta$V2)

df3_diet$cluster_fasting <- ifelse(df3_diet$cluster_fasting == 1, 4,
                                       ifelse(df3_diet$cluster_fasting == 2, 2, 
                                              ifelse(df3_diet$cluster_fasting == 3, 1,
                                                     ifelse(df3_diet$cluster_fasting == 4, 3, NA))))

df3_diet$cluster_fasting  <- as.factor(df3_diet$cluster_fasting)

df3_diet$cluster_postprandial <- as.factor(df3_diet$cluster_postprandial)

df3_diet$cluster_delta <- ifelse(df3_diet$cluster_delta == 1, 2,
                                     ifelse(df3_diet$cluster_delta == 2, 3, 
                                            ifelse(df3_diet$cluster_delta == 3, 1,
                                                   ifelse(df3_diet$cluster_delta == 4, 4, NA))))

# Define a function to summarise data (mean, sd)
summarise_diet <- function(data, group_var, variable) {
  data %>%
    group_by({{group_var}}) %>%
    dplyr::summarise(N = n(),
              mean_mean = mean({{ variable }}, na.rm = TRUE),
              sd = sd({{ variable }}, na.rm = TRUE),
              mean_sd = sprintf("%.2f ± %.2f", mean_mean, sd)
              ) 
}

# Apply the function to each variable
score_vars <- list("AMED_eadj", "hPDI_eadj",  "DII_eadj", "EDIH_eadj", "AHEI_eadj")

# fasting_result
results_fasting <- lapply(score_vars, function(var) {
  # Use !!sym to handle string variable names
  summarise_diet(df3_diet, group_var = cluster_fasting, !!sym(var))  
})

# postprandial_result
results_postprandial <- lapply(score_vars, function(var) {
  # Use !!sym to handle string variable names
  summarise_diet(df3_diet, group_var = cluster_postprandial, !!sym(var))  
})

results_fasting_combine = do.call(rbind, results_fasting) %>% 
  mutate(diet_scores = rep(c("AMED_eadj", "hPDI_eadj",  "DII_eadj", "EDIH_eadj", "AHEI_eadj"), each = 4)) %>% 
  select(cluster_fasting, N, diet_scores, mean_sd)

# long format -> wide format
require(tidyr)
results_fasting_combine_w = results_fasting_combine %>% 
  pivot_wider(names_from = diet_scores, values_from = mean_sd) %>% 
  mutate(Status = "Fasting",
         Group = case_when(cluster_fasting == 1 ~ "Cluster1",
                           cluster_fasting == 2 ~ "Cluster2",
                           cluster_fasting == 3 ~ "Cluster3",
                           cluster_fasting == 4 ~ "Cluster4"),
         .before = N) %>% 
  select(-cluster_fasting)

results_fasting_combine_w_trans = data.table::transpose(results_fasting_combine_w, list.cols = TRUE)
# setnames for the rows using column names
rownames(results_fasting_combine_w_trans) = colnames(results_fasting_combine_w)
results_fasting_combine_w_trans

# postprandial state
results_postprandial_combine = do.call(rbind, results_postprandial) %>% 
  mutate(diet_scores = rep(c("AMED_eadj", "hPDI_eadj",  "DII_eadj", "EDIH_eadj", "AHEI_eadj"), each = 4)) %>% 
  select(cluster_postprandial, N, diet_scores, mean_sd)

results_postprandial_combine_w = results_postprandial_combine %>% 
  pivot_wider(names_from = diet_scores, values_from = mean_sd) %>% 
  mutate(Status = "Postprandial",
         Group = case_when(cluster_postprandial == 1 ~ "Cluster1",
                           cluster_postprandial == 2 ~ "Cluster2",
                           cluster_postprandial == 3 ~ "Cluster3",
                           cluster_postprandial == 4 ~ "Cluster4"),
         .before = N
         ) %>% 
  select(-cluster_postprandial)

results_postprandial_combine_w_trans = data.table::transpose(results_postprandial_combine_w, list.cols = TRUE)
rownames(results_postprandial_combine_w_trans) = colnames(results_postprandial_combine_w)
results_postprandial_combine_w_trans

# combine fasting and postprandial results
results_fasting_combine_w_trans_2 = results_fasting_combine_w_trans %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(rowname != "Status") %>% 
  mutate(Status = c(" ", " ", rep("Fasting", 5)),
         .before = 1)
results_fasting_combine_w_trans_2

results_postprandial_combine_w_trans_2 = results_postprandial_combine_w_trans %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(rowname != "Status") %>% 
  mutate(Status = c(" ", " ", rep("Postprandial", 5)),
         .before = 1)
results_postprandial_combine_w_trans_2

results_fasting_post_combine = rbind(results_fasting_combine_w_trans_2, results_postprandial_combine_w_trans_2)

# x = unlist(results_fasting_post_combine$Status)
results_fasting_post_combine = data.frame(results_fasting_post_combine, stringsAsFactors = FALSE)
results_fasting_post_combine

results_fasting_post_combine = results_fasting_post_combine[-1, ]
colnames(results_fasting_post_combine) = c("Status", "Group", "Cluster1", "Cluster2", "Cluster3", "Cluster4")
results_fasting_post_combine = results_fasting_post_combine[-7,]
results_fasting_post_combine

# export 
# fwrite(results_fasting_post_combine, file = "Results/results_fasting_post_combine_V2.csv", sep = ",", row.names = FALSE)

library(openxlsx)
write.xlsx(results_fasting_post_combine, "Results/results_fasting_post_combine_V2.xlsx")

#### Make violin plot across different cluster #####

# change from wide to long
df3_diet_l <- df3_diet %>% 
  select(-glucose1) %>% 
  tidyr::pivot_longer(cols = c('AMED_eadj', 'hPDI_eadj', 'DII_eadj', 'EDIH_eadj', 'AHEI_eadj'),
               names_to = "Diet_score",
               values_to = "Scores") %>% 
  dplyr::mutate(Diet_score = factor(Diet_score))

# color theme
pal <- c(unname(yarrr::piratepal("pony")))[c(8,5,3,2,1,9)]

dp1 <- ggplot(df3_diet_l, aes(x = cluster_fasting, y = Scores, 
                              fill = cluster_fasting)) + 
  geom_violin(trim = FALSE)+
  geom_boxplot(width = 0.1, fill="white") + 
  scale_fill_manual(values = c("1" = pal[[1]], 
                               "2" = pal[[2]], 
                               "3" = pal[[3]],
                               "4" = pal[[4]]))  +
  facet_wrap(~ Diet_score, scales = "free_y", nrow = 1) +
  xlab("Fasting") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "top"
  ) +
  guides(fill = guide_legend(title = "Fasting\nCluster"))
dp1

# save figure (SFigure12A fasting)
ggsave(dp1, filename = "Results/SFigure12A.pdf", width = 10, height = 10)

dp2 <- ggplot(df3_diet_l, aes(x = cluster_postprandial, y = Scores, 
                              fill = cluster_postprandial)) + 
  geom_violin(trim = FALSE)+
  geom_boxplot(width = 0.1, fill = "white") + theme_classic() +
  scale_fill_manual(values = c("1" = pal[[1]], "2" = pal[[2]], "3" = pal[[3]],
                               '4' = pal[[4]]))  +
  facet_wrap(~ Diet_score, scales = "free_y", nrow = 1) +
  xlab("Postprandial") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "top"
    # legend.box.margin = margin(0, 0, 0, 6)
  ) +
  guides(fill = guide_legend(title = "Postprandial\nCluster"))
dp2

# save figure (SFigure12B postprandial)
ggsave(dp2, filename = "Results/SFigure12B.pdf", width = 10, height = 10)

# combine two figures (SFigure12)
library(patchwork)
png(file = "Results/SFigure12.png", width = 12, height = 8, units = "in", res = 300)
dp1/dp2
dev.off()

##### ANOVA ANALYSIS AND TUKEY TEST #######
 # the following part is not included in the manuscript, please skip :)
 # Define a function to perform the analysis
 ANOVA_analysis <- function(status, diet_score) {

   # Dynamically create the formula for the linear model
   formula <- as.formula(paste("Scores ~", status))

   # Create the linear model
   model <- lm(formula, data = df3_diet_l[df3_diet_l$Diet_score == diet_score,])

   # Perform ANOVA
   model_aov <- aov(model)

   # Output the summary of the ANOVA
   print(summary(model_aov))

   # Perform Tukey HSD test
   tukey_test <- TukeyHSD(model_aov)

   # Return the result of Tukey HSD test
   return(tukey_test)
 }

 # ANOVA and Tukey comparision: fasting
 fasting_comparision = lapply(score_vars, function(x){
   ANOVA_analysis(status = "cluster_fasting", x)
 })

 fasting_comparision[[1]]
 fasting_comparision[[2]]
 fasting_comparision[[3]]
 fasting_comparision[[4]]
 fasting_comparision[[5]]
 
 # ANOVA and Tukey comparision: postprandial
 postprandial_comparision = lapply(score_vars, function(x){
   ANOVA_analysis(status = "cluster_postprandial", x)
 })

 postprandial_comparision[[1]]
 postprandial_comparision[[2]]
 postprandial_comparision[[3]]
 postprandial_comparision[[4]]
 postprandial_comparision[[5]]
 
save(fasting_comparision, postprandial_comparision, file = "Results/Comparision_post_hoc.Rdata")

load(file = "Results/Comparision_post_hoc.Rdata")
# score_vars <- list("AMED_eadj", "hPDI_eadj",  "DII_eadj", "EDIH_eadj", "AHEI_eadj")

for (i in 1:5) {
  fasting1 <- fasting_comparision[[i]]$cluster_fasting |> data.frame() 
  fasting1 <- round(fasting1, 1)
  fasting1$diff_ci <- paste0(fasting1$diff, " (", fasting1$lwr, ",", fasting1$upr, ")")
  print(fasting1)
}


for (i in 1:5) {
  postprandial1 <- postprandial_comparision[[i]]$cluster_postprandial |> data.frame() 
  postprandial1 <- round(postprandial1, 1)
  postprandial1$diff_ci <- paste0(postprandial1$diff, " (", postprandial1$lwr, ",", postprandial1$upr, ")")
  print(postprandial1)
}
###### Correlation between Dietary pattern score and top metabolites #####
# index <- c('AMED_eadj', 'hPDI_eadj', 'DII_eadj', 'EDIH_eadj', 'AHEI_eadj')
# Meta_vars <- c("Gln_1", "Tyr_1", "MVLDLTG_1",
#                "Glc_1", "Cit_1", "bOHBut_1", "UnSat_1",
#                "Gln_3", "bOHBut_3", "Cit_3", "SHDLC_3", 
#                "LDLD_3", "MVLDLP_3", "LLDLP_3", "Alb_3", "sexe", "leeftijd", "bmim")
# 
# cor_data <- df_update %>% 
#   dplyr::select(all_of(c(index, Meta_vars))) %>% 
#   mutate(sexe = case_when(sexe == "man" ~ 0,
#                           sexe == "vrouw" ~ 1))
# 
# cor_data_complete <- cor_data[complete.cases(cor_data),]
# 
# library(ppcor)
# library(corrplot)
# covarites <- as.matrix(cor_data_complete[, c("sexe", "leeftijd", "bmim")])
# variables_interests <- as.matrix(scale(cor_data_complete[, c(index, "Gln_1", "Tyr_1", "MVLDLTG_1",
#                                               "Glc_1", "Cit_1", "bOHBut_1", "UnSat_1",
#                                               "Gln_3", "bOHBut_3", "Cit_3", "SHDLC_3", 
#                                               "LDLD_3", "MVLDLP_3", "LLDLP_3", "Alb_3")]))
# partial_cor <- pcor(cbind(variables_interests, covarites), method = "spearman")
# print(partial_cor$estimate) 
# 
# # p_cor_df = as.data.frame(partial_cor$estimate)
# # colnames(p_cor_df) <- c(index, Meta_vars)
# 
# p_cor_df_1 = partial_cor$estimate %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Rownames") %>% 
#   filter(Rownames %in% c("Gln_1", "Tyr_1", "MVLDLTG_1",
#                          "Glc_1", "Cit_1", "bOHBut_1", "UnSat_1",
#                          "Gln_3", "bOHBut_3", "Cit_3", "SHDLC_3", 
#                          "LDLD_3", "MVLDLP_3", "LLDLP_3", "Alb_3")) %>% 
#   dplyr::select(all_of(c('Rownames', index))) %>% 
#   column_to_rownames(var = "Rownames")
# 
# p_cor_df_1 <- p_cor_df_1 %>% mutate_if(is.numeric, round, 2)
# 
# p_corp_df_1 = partial_cor$p.value %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Rownames") %>% 
#   filter(Rownames %in% c("Gln_1", "Tyr_1", "MVLDLTG_1",
#                          "Glc_1", "Cit_1", "bOHBut_1", "UnSat_1",
#                          "Gln_3", "bOHBut_3", "Cit_3", "SHDLC_3", 
#                          "LDLD_3", "MVLDLP_3", "LLDLP_3", "Alb_3")) %>% 
#   dplyr::select(all_of(c('Rownames', index))) %>% 
#   column_to_rownames(var = "Rownames")
# 
# p_cor_df_1[p_corp_df_1 >= 0.05] <- NA
# 
# # Function to format numbers for display
# format_number <- function(x) {
#   ifelse(abs(x) < 0.01, 
#          sprintf("%.2e", x), 
#          sprintf("%.2f", x))
# }
# 
# # Create text matrix for cell labels
# text_matrix <- matrix(
#   paste(
#     format_number(as.matrix(p_cor_df_1)),
#     ifelse(p_corp_df_1 < 0.001, "***",
#            ifelse(p_corp_df_1 < 0.01, "**",
#                   ifelse(p_corp_df_1 < 0.05, "*", ""))),
#     sep = "\n"
#   ),
#   nrow = nrow(p_cor_df_1)
# )
# text_matrix[p_corp_df_1 >= 0.05] <- ""
# 
# # Create the heatmap
# library(pheatmap)
# library(RColorBrewer)
# my_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
# 
# png(file = "Results/Diet_score_metabolite_cor.png", width = 8, height = 10, units = "in", res = 300)
# pheatmap(
#   p_cor_df_1,
#   color = my_colors,
#   cluster_rows = FALSE,
#   cluster_cols = FALSE,
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   # Change the color of NA cells into white
#   na_col = "white",
#   display_numbers = text_matrix,
#   fontsize_number = 8,
#   number_color = "black",
#   fontsize = 8,
#   angle_col = 45, 
#   main = "Correlation Heatmap (Only Significant Results, p < 0.05)"
# )
# dev.off()
