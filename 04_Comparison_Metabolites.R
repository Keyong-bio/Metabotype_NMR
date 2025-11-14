##### P1-Cluster-NEO-KD 
##### Author: Keyong Deng
##### Department of Clinical Epidemiology, LUMC

######## 4. Compare the metabolites levels across the groups ##########
######## - for fasting and postprandial states only ########

library(foreign) # data load
library(dplyr) # %>% pipeline
library(data.table) # Fast Read data, fread(), fwrite()
library(Hmisc) # describe function
library(stringr) # string manuiplate
library(ggplot2) # Visualization
library(e1071) # Data skewness 
library(tibble)
# install iPVs package from github
library(iPVs) # identification of principal variable
library(tidyverse)
library(ggpubr)
library(ggh4x) # change the strip color
library(scales)

rm(list = ls())

load(file = "./working_imputated_df.Rdata")
load(file = "./k_means_clustering.Rdata")

# color palette
pal <- c(unname(yarrr::piratepal("pony")))[c(8,5,3,2,1,9)]

clusters_fasting <- kmeans_cluster_all[[1]] %>% as.data.frame() %>% 
  select(V2); table(clusters_fasting)
clusters_postprandial <- kmeans_cluster_all[[2]] %>% as.data.frame() %>% 
  select(V2); table(clusters_postprandial)
clusters_delta <- kmeans_cluster_all[[3]] %>% as.data.frame() %>% 
  select(V2); table(clusters_delta)

std_error <- function(x) {
  sd(x, na.rm = T) / sqrt(length(is.na(x)))
}

# color for the facet strip (the same color as KM), change the 50% transparency
strip_x <- strip_themed(
  background_x = elem_list_rect(fill = alpha(pal, 0.5))
)

# Compare all selected metabolites across the meta-botypes (Revision) --------------------
load(file =  "./df2_baseline.Rdata")

cols_select = c(colnames(working_imputated_df[[1]]), colnames(working_imputated_df[[2]]))

# the transformed metabolites
Met_comparision1 <- bind_cols(working_imputated_df[[1]], 
                         working_imputated_df[[2]]) %>% 
  mutate(cluster_fasting = clusters_fasting$V2,
         cluster_postprandial = clusters_postprandial$V2,
         .before = 1)

# # the un-transformed metabolites
# Met_comparision2 <- df2 %>% 
#   select(all_of(cols_select)) %>% 
#   mutate(cluster_fasting = clusters_fasting$V2,
#          cluster_postprandial = clusters_postprandial$V2,
#          .before = 1)

# re-code the group
Met_comparision1$cluster_fasting <- ifelse(Met_comparision1$cluster_fasting == 1, 4,
                                      ifelse(Met_comparision1$cluster_fasting == 2, 2, 
                                             ifelse(Met_comparision1$cluster_fasting == 3, 1,
                                                    ifelse(Met_comparision1$cluster_fasting == 4, 3, NA))))

Met_comparision1$cluster_fasting <- factor(Met_comparision1$cluster_fasting, 
                                          levels = c("1", "2", "3", "4"),
                                          labels = c("Low", "Low-intermediate", 
                                                 "High-intermediate", "High"))

Met_comparision1$cluster_postprandial <- factor(Met_comparision1$cluster_postprandial, 
                                                levels = c("1", "2", "3", "4"),
                                                labels = c("Low", "Low-intermediate", 
                                                           "High-intermediate", "High"))

# wide to long format 
Met_comparision1_fasting_l <- pivot_longer(Met_comparision1 %>% select(cluster_fasting, ends_with("_1")), 
                                   cols = MVLDLTG_1:Alb_1,
                                   names_to = "Metabolites",
                                   values_to = "Value") 

Met_comparision1_postprandial_l <- pivot_longer(Met_comparision1 %>% select(cluster_postprandial, ends_with("_3")), 
                                           cols = MVLDLP_3:Alb_3,
                                           names_to = "Metabolites",
                                           values_to = "Value") 


# Met_comparision2_l_s <- Met_comparision2_l %>% group_by(cluster_fasting, Metabolites) %>% 
#   dplyr::summarise(
#     mean_val = mean(Data),
#     median_val = median(Data),
#     sd_val = sd(Data),
#     .groups = "drop"
#   ) 

p_fastingmetabolomic <- ggplot(Met_comparision1_fasting_l, 
       aes(x = cluster_fasting, 
           y = Value, 
           fill = cluster_fasting)) + 
  geom_violin(trim = FALSE)+
  geom_boxplot(width = 0.1, fill="white") + 
  scale_fill_manual(values = c("Low" = pal[[1]],
                               "Low-intermediate" = pal[[2]],
                               "High-intermediate" = pal[[3]],
                               "High" = pal[[4]]))  +
  facet_wrap(~ Metabolites, scales = "free_y", nrow = 4) +
  xlab("Fasting") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    # axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "top"
  )

pdf(file = "./Results/Metabolites_comparison_Fasting_revision.pdf", width = 10, height = 10)
p_fastingmetabolomic
dev.off()

# postprandial data
p_postprandial_metabolomic <- ggplot(Met_comparision1_postprandial_l, 
                               aes(x = cluster_postprandial, 
                                   y = Value, 
                                   fill = cluster_postprandial)) + 
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill="white") + 
  scale_fill_manual(values = c("Low" = pal[[1]],
                               "Low-intermediate" = pal[[2]],
                               "High-intermediate" = pal[[3]],
                               "High" = pal[[4]]))  +
  facet_wrap(~ Metabolites, scales = "free_y", nrow = 4) +
  xlab("Fasting") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    # axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "top"
  )

pdf(file = "./Results/Metabolites_comparison_postprandial_revision.pdf", width = 10, height = 10)
p_postprandial_metabolomic
dev.off()

# Compare the baseline characteristic for those with switching and without switching (revision) ------------
head(Met_comparision1)

sum(Met_comparision1["cluster_fasting"] == Met_comparision1["cluster_postprandial"])

# 0 as stay, 1 as switch
Met_comparision1$switch = ifelse(Met_comparision1["cluster_fasting"] == Met_comparision1["cluster_postprandial"],
                                 0, 1) 

vars_base <- c("medication_hypertension", "medication_lipidlower", "sexe", "leeftijd", "bmim",  "glucose1", "Insuline_r1",
               "HBA1C", "homa1IR", "homa1B", "choltot1", "trig1", "hdlc1", "fldl1", "choltot3", "trig3", "hdlc3", "fldl3", "TG_perc")

base_df2 <- df2 %>% 
  select(all_of(vars_base)) %>% 
  cbind(Met_comparision1[["switch"]]) %>% 
  rename(Switch = cluster_fasting) %>% 
  dplyr::mutate(Switch = as.factor(Switch))

cat_vars <- c("sexe", "Switch", "medication_hypertension", "medication_lipidlower")

pacman::p_load("gtsummary", ## summary statistics and tests
               "rstatix",  ## summary statistics and statistical tests
               "janitor", ## adding totals and percents to tables
               "scales", ## easily convert proportions to percents  
               "flextable", ## converting tables to pretty images
               "skimr", ## get overview of data
               "foreign", "dplyr", "data.table", "Hmisc", "stringr", "ggplot2", "e1071") 

table_1_1 <- base_df2 %>%
  dplyr::select(-medication_hypertension, -medication_lipidlower, dplyr::everything()) %>% 
  mutate(sexe = factor(sexe, labels = c("Female", "Male")),
         medication_hypertension = factor(medication_hypertension, levels = c(0, 1), labels = c("No", "Yes")),
         medication_lipidlower = factor(medication_lipidlower, levels = c(0, 1), labels = c("No", "Yes"))) %>% 
  tbl_summary(
    by = Switch ,
    type = all_continuous() ~ "continuous2",
    statistic = list(
      all_continuous() ~ c("{mean} ({sd})",
                           "{median} ({p25}, {p75})"),
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(sexe ~ "Gender", 
                 leeftijd ~ "Age", bmim ~ "BMI", 
                 glucose1 ~ "Fasting glucose", Insuline_r1 ~ "Fasting insuline",
                 homa1IR ~ "HOMA1-IR", homa1B ~ "HOMA1-B", choltot1 ~ "Fasting TC",
                 trig1 ~ "Fasting TG", hdlc1 ~ "Fasting HDL", fldl1 ~ "Fasting LDL",
                 trig3 ~ "Postprandial TG", hdlc3 ~ "Postprandial HDL", fldl3 ~ "Postprandial LDL",
                 choltot3 ~ "Postprandial TC",
                 TG_perc ~ "TG Percent in Liver", 
                 medication_hypertension ~ "Antihypertensive medication",
                 medication_lipidlower ~ "Lipid-lowering medication"),
    # missing_text = "(Missing)"
    missing = "no"
  ) %>% 
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% 
  # add_overall() %>%
  add_n() %>%  
  modify_header(label ~ "**Variable**")
theme_gtsummary_journal(journal = "jama")
table_1_1

as_gt(table_1_1) |> 
  gt::gtsave(filename = "Results/Stable_switch_baseline.docx")

# plot the core metabolites for the each metabotype ---------------------------
fasting_vars <- c("MVLDLTG_1", "UnSat_1", "bOHBut_1", "Tyr_1", "Glc_1", "Gln_1", "Cit_1")
Met_fasting <- bind_cols(working_imputated_df[[1]], 
                         working_imputated_df[[2]]) %>% 
  select(all_of(fasting_vars)) %>% 
  mutate(cluster_fasting = clusters_fasting$V2,
         .before = 1)

Met_fasting$cluster_fasting <- ifelse(Met_fasting$cluster_fasting == 1, 4,
                                       ifelse(Met_fasting$cluster_fasting == 2, 2, 
                                              ifelse(Met_fasting$cluster_fasting == 3, 1,
                                                     ifelse(Met_fasting$cluster_fasting == 4, 3, NA))))

Met_fasting$cluster_fasting <- factor(Met_fasting$cluster_fasting, levels = c("1", "2", "3", "4"),
                                      labels = c("Low risk", "Low-intermediate risk", 
                                                "High-intermediate risk", "High risk"))

# wide dataframe to long dataframe
Met_fasting_long <- Met_fasting %>% 
  pivot_longer(cols = MVLDLTG_1:Cit_1, names_to = "Metabolites", values_to = "Value") %>% 
  data.frame()

Met_fasting_longsum <- Met_fasting_long %>% 
  group_by(Metabolites, cluster_fasting) %>%
  summarise(mean_value = mean(Value, na.rm = T),
            median_value = median(Value, na.rm = T),
            sd_value = sd(Value, na.rm = T),
            se_value = std_error(Value),
            lci = mean_value - 1.96 * se_value,
            uci = mean_value + 1.96 * se_value,
            .groups = "drop") 

Met_fasting_longsum$highlight <- ifelse(Met_fasting_longsum$cluster_fasting == "Low risk" & Met_fasting_longsum$Metabolites %in% c("bOHBut_1", "Glc_1", "MVLDLTG_1"), "yes", 
                                        ifelse(Met_fasting_longsum$cluster_fasting == "Low-intermediate risk" & Met_fasting_longsum$Metabolites %in% c("bOHBut_1", "UnSat_1", "MVLDLTG_1"), "yes",
                                               ifelse(Met_fasting_longsum$cluster_fasting == "High-intermediate risk" & Met_fasting_longsum$Metabolites %in% c("Gln_1", "Tyr_1", "Cit_1"), "yes",
                                                      ifelse(Met_fasting_longsum$cluster_fasting == "High risk" & Met_fasting_longsum$Metabolites %in% c("UnSat_1", "MVLDLTG_1", "Gln_1"), "yes", "no"))))

Met_fasting_longsum$highlight <- factor(Met_fasting_longsum$highlight, levels = c("yes", "no"))
Met_fasting_longsum$Metabolites <- factor(Met_fasting_longsum$Metabolites, levels = c("Gln_1", "Tyr_1", "MVLDLTG_1",
                                                                                     "Glc_1", "Cit_1", "bOHBut_1", "UnSat_1"))

# SFigure11A
pdf(file = "Results/SFigure11A.pdf", width = 6, height = 12)
plot_1 <- ggplot(data = Met_fasting_longsum, aes(x = Metabolites , y = mean_value)) +
  geom_pointrange(aes(ymin = lci, ymax = uci, color = highlight),
                  position = position_dodge(0.3)) +
  geom_hline(yintercept = c(-0.5, 0, 0.5), col = "grey50", lty = "dashed") +
  scale_color_manual(values = c("red", "white")) +
  ggtitle(label = "Fasting_state") +
  xlab(NULL) +
  ylab("Standardized_Value") +
  facet_wrap2(~ cluster_fasting, strip = strip_x, ncol = 1) +
  theme_classic2() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "none",  # remove the legend 
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    #panel.background = element_blank(),  # blank background
    #panel.border = element_blank(),      # remove border
    axis.line = element_line(colour = "black", linewidth = 1),  # black axis lines
    axis.ticks = element_line(color = "black", linewidth = 1),  # black axis ticks
    panel.grid.major = element_blank(),  # remove major grid
    panel.grid.minor = element_blank(),  # remove minor grid
    plot.title = element_text(hjust = 0, size = 12)  # left the plot title if you have one
  ) 
print(plot_1)
dev.off()

# Selected postprandial metabolites
post_vars = c("Gln_3", "bOHBut_3", "Ala_3", "SHDLC_3", "LDLD_3", "MVLDLP_3", "LLDLP_3", "Alb_3")
Met_postprandial <- bind_cols(working_imputated_df[[1]], 
                              working_imputated_df[[2]]) %>% 
  select(all_of(post_vars)) %>% 
  mutate(cluster_postprandial = clusters_postprandial$V2,
         .before = 1)

Met_postprandial$cluster_postprandial <- factor(Met_postprandial$cluster_postprandial, 
                                                levels = c("1", "2", "3", "4"),
                                                labels = c("Low risk", "Low-intermediate risk", 
                                                 "High-intermediate risk", "High risk"))

Met_postprandial_long <- Met_postprandial %>% 
  pivot_longer(cols = Gln_3:Alb_3, names_to = "Metabolites", values_to = "Value") %>% 
  data.frame()

Met_postprandial_longsum <- Met_postprandial_long %>% 
  group_by(Metabolites, cluster_postprandial) %>%
  summarise(mean_value = mean(Value, na.rm = T),
            median_value = median(Value, na.rm = T),
            sd_value = sd(Value, na.rm = T),
            se_value = std_error(Value),
            lci = mean_value - 1.96 * se_value,
            uci = mean_value + 1.96 * se_value) 

Met_postprandial_longsum$highlight <- ifelse(Met_postprandial_longsum$cluster_postprandial == "Low risk" & Met_postprandial_longsum$Metabolites %in% c("Gln_3", "bOHBut_3", "MVLDLP_3"), "yes", 
                                        ifelse(Met_postprandial_longsum$cluster_postprandial == "Low-intermediate risk" & Met_postprandial_longsum$Metabolites %in% c("SHDLC_3", "LLDLP_3", "Alb_3"), "yes",
                                               ifelse(Met_postprandial_longsum$cluster_postprandial == "High-intermediate risk" & Met_postprandial_longsum$Metabolites %in% c("LDLD_3", "Gln_3", "Ala_3"), "yes",
                                                      ifelse(Met_postprandial_longsum$cluster_postprandial == "High risk" & Met_postprandial_longsum$Metabolites %in% c("SHDLC_3", "bOHBut_3", "MVLDLP_3"), "yes", "no"))))

Met_postprandial_longsum$highlight <- factor(Met_postprandial_longsum$highlight, levels = c("yes", "no"))

Met_postprandial_longsum$Metabolites <- factor(Met_postprandial_longsum$Metabolites, 
                                               levels = c("Gln_3", "MVLDLP_3", "LLDLP_3", "SHDLC_3", "LDLD_3", "Ala_3", "bOHBut_3", "Alb_3"))

# SFigure11B
pdf(file = "Results/SFigure11B.pdf", width = 6, height = 12)
plot_2 <- ggplot(data = Met_postprandial_longsum, aes(x = Metabolites , y = mean_value)) +
  geom_pointrange(aes(ymin = lci, ymax = uci, color = highlight),
                  position = position_dodge(0.3)) +
  geom_hline(yintercept = c(-0.5, 0, 0.5), col = "grey50", lty = "dashed") +
  scale_color_manual(values = c("red", "white")) +
  ggtitle(label = "Postprandial_state") +
  xlab(NULL) +
  ylab("Standardized_Value") +
  facet_wrap2(~ cluster_postprandial, strip = strip_x, ncol = 1) +
  theme_classic2() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "none",  # remove the legend 
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(colour = "black", linewidth = 1),  # black axis lines
    axis.ticks = element_line(color = "black", linewidth = 1),  # black axis ticks
    panel.grid.major = element_blank(),  # remove major grid
    panel.grid.minor = element_blank(),  # remove minor grid
    plot.title = element_text(hjust = 0, size = 12)  # left the plot title if you have one
    
  ) 
print(plot_2)
dev.off()

# combine two plots
library(cowplot)

# SFigure11
png(filename = "Results/SFigure11.png", 
    width = 8,       
    height = 12,       
    res = 500,
    units = "in")     # Using inches instead of pixels

plot_grid(plot_1, plot_2, labels = c('A', 'B'), 
          label_fontface = "bold",
          label_size = 16)

dev.off()
