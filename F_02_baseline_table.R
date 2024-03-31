########################################
#### P1-Cluster-NEO-KD##################
########################################

library(foreign) # data load
library(dplyr) # %>% pipeline
library(data.table) # Fast Read data, fread(), fwrite()
library(Hmisc) # describe function
library(stringr) # string manuiplate
library(ggplot2) # Visualization
library(e1071) # Data skewness 
library(tibble)
# ifelse("devtools" %in% rownames(installed.packages()), 
# 	 NA, install.packages("devtools"))
# devtools::install_github("hughesevoanth/iPVs")
library(iPVs) #identification of principal variable

########## Part-3 ####################

###  Baseline characteristic for the participants ####
load("Data/k_means_clustering.Rdata")

# library(amap) ## k-MEANs method
# kmeans_cluster_all <- list()
# for (i in 1:3) {
#   # i = 1
#   new_cluster_kmeans <- c()
#   for (j in 3:5) {
#     set.seed(123)
#     # j = 4
#     kmeans_cluster = amap::Kmeans(working_imputated_df[[i]], centers = j, iter.max = 200, nstart = 50)
#     new_cluster_kmeans <- cbind(new_cluster_kmeans, kmeans_cluster$cluster)
#   }
#   kmeans_cluster_all[[i]] <- new_cluster_kmeans
# }

kmeans_cluster_fasting <- kmeans_cluster_all[[1]] %>% as.data.frame() %>% 
  dplyr::rename(C3 = "V1", C4 = "V2", C5 = "V3")
table(kmeans_cluster_fasting$C3); table(kmeans_cluster_fasting$C4); table(kmeans_cluster_fasting$C5);
kmeans_cluster_fasting$ID <- 1:5320

kmeans_cluster_post <- kmeans_cluster_all[[2]] %>% as.data.frame() %>% 
  dplyr::rename(C3 = "V1", C4 = "V2", C5 = "V3")
table(kmeans_cluster_post$C3); table(kmeans_cluster_post$C4); table(kmeans_cluster_post$C5);
kmeans_cluster_post$ID <- 1:5320

kmeans_cluster_delta <- kmeans_cluster_all[[3]] %>% as.data.frame() %>% 
  dplyr::rename(C3 = "V1", C4 = "V2", C5 = "V3")
table(kmeans_cluster_delta$C3); table(kmeans_cluster_delta$C4); table(kmeans_cluster_delta$C5);
kmeans_cluster_delta$ID <- 1:5320

clusters.fasting <- kmeans_cluster_fasting %>% select(C4); table(clusters.fasting)
clusters.postprandial <- kmeans_cluster_post %>% select(C4); table(clusters.postprandial)
clusters.delta <- kmeans_cluster_delta %>% select(C4); table(clusters.delta)

### Get the ipvs variables
imputated_ipvs_vars <- list()
for (i in 1:3) {
  imputated_ipvs_vars[[i]] <- colnames(working_imputated_df[[i]])
}
imputated_ipvs_vars

###  To use PCA analysis on the transformed_states data ####
library(ggrepel)
clusters.states <- list(Fasting = clusters.fasting, 
                        Postprandial = clusters.postprandial, 
                        Delta = clusters.delta)
for (i in 1:3) {
  # i = 1
  data <- imputated_states_transformation[[i]] %>% as.data.frame()
  pca <- prcomp(data)
  pc_coords <- data.frame(cbind(pca$x, Cluster = clusters.states[[i]]))
  colnames(pc_coords)[ncol(pc_coords)] <- "Cluster"
  pc_coords$Cluster <- as.factor(pc_coords$Cluster)
  pc_variance <- pca$sdev^2/sum(pca$sdev^2)
  
  # pdf(file = "pc_variance_fasting.pdf", width = 8, height = 8)
  pdf(file = paste0("./Results_check/pc_variance_", names(clusters.states)[[i]], "_2024.pdf"), width = 8, height = 8)
  plot(pc_variance[1:30], type = "b",cex = 0.5, xlab = "Principal Component", ylab = "Proportion of Variance explained")
  abline(h = 0.01)
  dev.off()
  
  p1 <- ggplot(pc_coords, aes(x=PC1, y=PC2)) + 
    geom_point(aes(color=Cluster), size=3) + ## color = Batch
    scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c", "#41b6c4")) +
    xlab(paste("PC1 (", round(pc_variance[1]*100,2), "% variance explained)")) +
    ylab(paste("PC2 (", round(pc_variance[2]*100,2), "% variance explained)")) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ggtitle(names(imputated_states_transformation)[i])
  ggsave(filename = paste0("./Results_check/p_", i, "_pca_", names(imputated_states_transformation)[i], "PC1vsPC2.pdf"), width = 8, height = 8)
  
  ### Extract loading for the first 10 PCs
  pc_loadings <- data.frame(metabolites = rownames(pca$rotation),
                            pca$rotation[,1:10])
  pc1_loadings <- pc_loadings[order(pc_loadings$PC1, decreasing = T),]
  pc1_loadings$metabolites_rank <- 1:nrow(pc1_loadings) 
  
  pc1_loadings_plot <- ggplot(pc1_loadings, aes(x = metabolites_rank, y = PC1)) +
    geom_point() +
    geom_point(data = pc1_loadings[intersect(rownames(pc1_loadings),imputated_ipvs_vars[[i]]),], 
               aes(x = metabolites_rank, y = PC1, label = metabolites), size = 3, color = "darkred") +
    geom_label_repel(data = pc1_loadings[intersect(rownames(pc1_loadings),imputated_ipvs_vars[[i]]),],
                     aes(x = metabolites_rank, y = PC1, label = metabolites), size = 3, color = "darkred") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, size = 20),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 20))
  ggsave(filename = paste0("p_", i, "_PC1_pca_loading", names(imputated_states_transformation)[i], "_2024_2.pdf"), width = 10, height = 10)  
  
  pc2_loadings <- pc_loadings[order(pc_loadings$PC2, decreasing = T),]
  pc2_loadings$metabolites_rank <- 1:nrow(pc2_loadings) 
  pc2_loadings_plot_2 <- ggplot(pc2_loadings, aes(x = metabolites_rank, y = PC2)) +
    geom_point() +
    geom_point(data = pc2_loadings[intersect(rownames(pc2_loadings),imputated_ipvs_vars[[i]]),],
               aes(x = metabolites_rank, y = PC2, label = metabolites), size = 3, color = "darkred") +
    geom_label_repel(data = pc2_loadings[intersect(rownames(pc2_loadings),imputated_ipvs_vars[[i]]),],
                     aes(x = metabolites_rank, y = PC2, label = metabolites), size = 3, color = "darkred") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, size = 20),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 20))
  ggsave(filename = paste0("p", i, "_PC2_pca_loading", names(imputated_states_transformation)[i], "_2024_2.pdf"), width = 10, height = 10) 
  
  pc3_loadings <- pc_loadings[order(pc_loadings$PC3, decreasing = T),]
  pc3_loadings$metabolites_rank <- 1:nrow(pc3_loadings) 
  pc3_loadings_plot_2 <- ggplot(pc3_loadings, aes(x = metabolites_rank, y = PC3)) +
    geom_point() +
    geom_point(data = pc3_loadings[intersect(rownames(pc3_loadings),imputated_ipvs_vars[[i]]),], 
               aes(x = metabolites_rank, y = PC3, label = metabolites), size = 3, color = "darkred") +
    geom_label_repel(data = pc3_loadings[intersect(rownames(pc3_loadings),imputated_ipvs_vars[[i]]),],
                     aes(x = metabolites_rank, y = PC3, label = metabolites), size = 3, color = "darkred") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, size = 20),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 20))
  ggsave(filename =  paste0("p", i, "_PC3_pca_loading", names(imputated_states_transformation)[i], "_2024_2.pdf"), width = 10, height = 10) 
  
}

## PCA_variance for each states (fasting, postprandial, delta states)
for (i in 1:3) {
  # i = 1
  data <- imputated_states_transformation[[i]] %>% as.data.frame() # include all the metabolites
  pca <- prcomp(data)
  pc_variance <- pca$sdev^2/sum(pca$sdev^2) # variances 
  
  pdf(file = paste0("./Results_check/pc_variance", names(imputated_states_transformation)[i], "scree_plot.pdf"), width = 8, height = 8)
  plot(pc_variance[1:50], type = "b",cex = 0.5, xlab = "Principal Component", ylab = "Proportion of Variance explained")
  abline(h = 0.01) ## variance > 1%
  dev.off()
}


### Another part ###########
### Baseline characteristic (divided by different clusters)
load("Data/After_QC_data.Rdata") ## import df1 data
metadata_out <- cbind(df1, cluster_fasting = clusters.fasting, cluster_postprandial =  clusters.postprandial, cluster_delta = clusters.delta)
colnames(metadata_out)[521:523] <- c("cluster_fasting", "cluster_postprandial", "cluster_delta")

summary1 <- function(data, var, var1){
  data %>% 
    group_by({{var1}}) %>% 
    summarise(
      median = median({{var}}, na.rm = T),
      IQR = paste0(sprintf("%.2f", quantile({{var}}, c(0.25), na.rm = T)), "-",
                   sprintf("%.2f", quantile({{var}}, c(0.75), na.rm = T))),
      mean = mean({{var}}, na.rm = T),
      sd = sd({{var}}, na.rm = T),
      n = n(),
      n_miss = sum(is.na({{var}})),
      .groups = "drop" ## avoid the message and leave the data in an ungrouped state
    )
}

######
independent_vars <- cbind(working_imputated_df[[1]], working_imputated_df[[2]], working_imputated_df[[3]])
vars_base <- c("medication_hypertension", "medication_lipidlower", "sexe", "leeftijd", "bmim",  "glucose1", "Insuline_r1",
               "HBA1C", "homa1IR", "homa1B", "choltot1", "trig1", "hdlc1", "fldl1", "choltot3", "trig3", "hdlc3", "fldl3", "TG_perc",
               "cluster_fasting", "cluster_postprandial", "cluster_delta")

base_df1 <- metadata_out %>% dplyr::select(all_of(vars_base))

table(base_df1$sexe, useNA = "always") ## recode the sexe
cat_vars <- c("sexe", "cluster_fasting", "cluster_postprandial", "cluster_delta", "medication_hypertension", "medication_lipidlower")
base_df1 %<>% dplyr::mutate(sexe = case_when(sexe == "vrouw" ~ 0,  #### 0 denotes female
                                             sexe == "man" ~ 1)) %>% 
  mutate_at(vars(one_of(cat_vars)), funs(factor(.))) ## change those cat_vars into categorical variable
# str(base_df1)

pacman::p_load("gtsummary", ## summary statistics and tests
               "rstatix",  ## summary statistics and statistical tests
               "janitor", ## adding totals and percents to tables
               "scales", ## easily convert proportions to percents  
               "flextable", ## converting tables to pretty images
               "skimr" ## get overview of data
) 

## test the distribution and skewness of variables 
vars_exam <- c("leeftijd", "bmim",  "glucose1", "Insuline_r1","HBA1C", "homa1IR", "homa1B", "choltot1", "trig1", "hdlc1", "fldl1", "choltot3", "trig3", "hdlc3", "fldl3", "TG_perc", "cluster_fasting", "cluster_postprandial", "cluster_delta")
base_df1_exam <- base_df1  %>% select(all_of(vars_exam))
apply(base_df1_exam[,-(17:19)], 2, function(x) skewness(x, na.rm = T)) ## glucose1, Insuline_r1, HBA1C, homa1IR, homa1B, trig1, trig3, TG_perc

## str_glue function  from stringr is used to combine values from several columns into one new column.
## https://epirhandbook.com/en/descriptive-tables.html#descriptive-tables
## using select function to reorder and rename the column names
## across() function: multiple columns, when calculate the same statistics for many columns.

######. fasting states baseline #######
table1 <- base_df1 %>%
  dplyr::select(-c(cluster_postprandial, cluster_delta)) %>%  # exclude the other two columns
  dplyr::select(-medication_hypertension, -medication_lipidlower, dplyr::everything()) %>% 
  mutate(sexe = factor(sexe, labels = c("Female", "Male")),
         medication_hypertension = factor(medication_hypertension, labels = c("No", "Yes")),
         medication_lipidlower = factor(medication_lipidlower, labels = c("No", "Yes"))) %>% 
  tbl_summary(
    by = cluster_fasting ,
    type = all_continuous() ~ "continuous2",
    statistic = list(
      all_continuous() ~ c("{mean} ({sd})",
                           "{median} ({p25}, {p75})"),
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(sexe ~ "Gender", leeftijd ~ "Age", bmim ~ "BMI", 
                 glucose1 ~ "Fasting glucose", Insuline_r1 ~ "Fasting insuline",
                 homa1IR ~ "HOMA1-IR", homa1B ~ "HOMA1-B", choltot1 ~ "Fasting TC",
                 trig1 ~ "Fasting TG", hdlc1 ~ "Fasting HDL", fldl1 ~ "Fasting LDL",
                 trig3 ~ "Postprandial TG", hdlc3 ~ "Postprandial HDL", fldl3 ~ "Postprandial LDL",
                 choltot3 ~ "Postprandial TC",
                 TG_perc ~ "TG Percent in Liver", medication_hypertension ~ "Antihypertensive medication",
                 medication_lipidlower ~ "Lipid-lowering medication"),
    # missing_text = "(Missing)"
    missing = "no"
  ) %>% 
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% 
  # add_overall() %>%
  add_n() %>%  
  modify_header(label ~ "**Variable**")
theme_gtsummary_journal(journal = "jama")
table1

### postprandial state_baseline #####
table2 <- base_df1 %>%
  dplyr::select(-c(cluster_fasting, cluster_delta)) %>%  # exclude the other two columns
  dplyr::select(-medication_hypertension, -medication_lipidlower, dplyr::everything()) %>% 
  mutate(sexe = factor(sexe, labels = c("Female", "Male")),
         medication_hypertension = factor(medication_hypertension, labels = c("No", "Yes")),
         medication_lipidlower = factor(medication_lipidlower, labels = c("No", "Yes"))) %>% 
  tbl_summary(
    by = cluster_postprandial ,
    type = all_continuous() ~ "continuous2",
    statistic = list(
      all_continuous() ~ c("{mean} ({sd})",
                           "{median} ({p25}, {p75})"),
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(sexe ~ "Gender", leeftijd ~ "Age", bmim ~ "BMI", 
                 glucose1 ~ "Fasting glucose", Insuline_r1 ~ "Fasting insuline",
                 homa1IR ~ "HOMA1-IR", homa1B ~ "HOMA1-B", choltot1 ~ "Fasting TC",
                 trig1 ~ "Fasting TG", hdlc1 ~ "Fasting HDL", fldl1 ~ "Fasting LDL",
                 trig3 ~ "Postprandial TG", hdlc3 ~ "Postprandial HDL", fldl3 ~ "Postprandial LDL",
                 choltot3 ~ "Postprandial TC",
                 TG_perc ~ "TG Percent in Liver", medication_hypertension ~ "Antihypertensive medication",
                 medication_lipidlower ~ "Lipid-lowering medication"),
    # missing_text = "(Missing)"
    missing = "no"
  ) %>% 
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% 
  # add_overall() %>%
  add_n() %>%  
  modify_header(label ~ "**Variable**")
table2

####### Delta state_baseline ######
table3 <- base_df1 %>%
  dplyr::select(-c(cluster_fasting, cluster_postprandial)) %>%  # exclude the other two columns
  dplyr::select(-medication_hypertension, -medication_lipidlower, dplyr::everything()) %>% 
  mutate(sexe = factor(sexe, labels = c("Female", "Male")),
         medication_hypertension = factor(medication_hypertension, labels = c("No", "Yes")),
         medication_lipidlower = factor(medication_lipidlower, labels = c("No", "Yes"))) %>% 
  tbl_summary(
    by = cluster_delta,
    type = all_continuous() ~ "continuous2",
    statistic = list(
      all_continuous() ~ c("{mean} ({sd})",
                           "{median} ({p25}, {p75})"),
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(sexe ~ "Gender", leeftijd ~ "Age", bmim ~ "BMI", 
                 glucose1 ~ "Fasting glucose", Insuline_r1 ~ "Fasting insuline",
                 homa1IR ~ "HOMA1-IR", homa1B ~ "HOMA1-B", choltot1 ~ "Fasting TC",
                 trig1 ~ "Fasting TG", hdlc1 ~ "Fasting HDL", fldl1 ~ "Fasting LDL",
                 trig3 ~ "Postprandial TG", hdlc3 ~ "Postprandial HDL", fldl3 ~ "Postprandial LDL",
                 choltot3 ~ "Postprandial TC",
                 TG_perc ~ "TG Percent in Liver", medication_hypertension ~ "Antihypertensive medication",
                 medication_lipidlower ~ "Lipid-lowering medication"),
    # missing_text = "(Missing)"
    missing = "no"
  ) %>% 
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% 
  # add_overall() %>%
  add_n() %>%  
  modify_header(label ~ "**Variable**")
table3

## total population
table4 <- base_df1 %>%
  dplyr::select(-c(cluster_postprandial, cluster_fasting, cluster_delta)) %>%  # exclude the other two columns
  dplyr::select(-medication_hypertension, -medication_lipidlower, dplyr::everything()) %>% 
  mutate(sexe = factor(sexe, labels = c("Female", "Male")),
         medication_hypertension = factor(medication_hypertension, labels = c("No", "Yes")),
         medication_lipidlower = factor(medication_lipidlower, labels = c("No", "Yes"))) %>% 
  tbl_summary(
    # by = cluster_fasting ,
    type = all_continuous() ~ "continuous2",
    statistic = list(
      all_continuous() ~ c("{mean} ({sd})",
                           "{median} ({p25}, {p75})"),
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(sexe ~ "Gender", leeftijd ~ "Age", bmim ~ "BMI", 
                 glucose1 ~ "Fasting glucose", Insuline_r1 ~ "Fasting insuline",
                 homa1IR ~ "HOMA1-IR", homa1B ~ "HOMA1-B", choltot1 ~ "Fasting TC",
                 trig1 ~ "Fasting TG", hdlc1 ~ "Fasting HDL", fldl1 ~ "Fasting LDL",
                 trig3 ~ "Postprandial TG", hdlc3 ~ "Postprandial HDL", fldl3 ~ "Postprandial LDL",
                 choltot3 ~ "Postprandial TC",
                 TG_perc ~ "TG Percent in Liver", medication_hypertension ~ "Antihypertensive medication",
                 medication_lipidlower ~ "Lipid-lowering medication"),
    # missing_text = "(Missing)"
    missing = "no"
  ) %>% 
  # add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% 
  # add_overall() %>%
  add_n() %>%  
  modify_header(label ~ "**Variable**")
theme_gtsummary_journal(journal = "jama")
table4

############## Compare the levels across clusters ###############
### Combine the metabolites info with cluster 
clusters.states <- list(clusters.fasting, clusters.postprandial, clusters.delta)
cluster_metabolite <- list()
for (i in 1:3) {
  # i <- 1
  temp1 <- working_imputated_df[[i]]
  # temp2 <- imputated_states_min[[i]] %>%
  #   select(all_of(colnames(temp1)))
  # temp2 <- apply(temp2, 2, function(x) scale(x, center = T))
  temp2 <- clusters.states[[i]]
  temp3 <- cbind(temp1, temp2) 
  colnames(temp3)[ncol(temp3)] <- "Cluster"
  temp3 <- apply(temp3, 2, as.numeric) %>% as.data.frame()
  cluster_metabolite[[i]] <- temp3
}

## plot the line plot for selected metabolites
library(tidyr)
library(forcats)
plot_min <- list()
std_error <- function(x) {
  sd(x) / sqrt(length(x))
}
for (i in 1:3) {
  # i <- 1
  temp1 <- cluster_metabolite[[i]]
  temp1 <- pivot_longer(temp1, cols = -Cluster, names_to = "Metabolites", values_to = "Value") %>% data.frame() %>% 
    group_by(Metabolites, Cluster) %>% 
    summarise(mean_value = mean(Value),
              sd_value = sd(Value),
              se_value = std_error(Value),
              lci = mean_value - 1.96*se_value,
              uci = mean_value + 1.96*se_value) %>% 
    ungroup()
  temp1$Cluster <- as.factor(temp1$Cluster)
  plot_min[[i]] <- ggplot(data = temp1, aes(x = fct_reorder(factor(Metabolites), mean_value) , y = mean_value, color = Cluster)) +
    geom_line(aes(group = Cluster), linewidth = 1.5) + 
    geom_point() + 
    geom_errorbar(aes(ymin=lci, ymax=uci), width = 0.5)
  theme_minimal() + 
    # facet_wrap(~ Cluster, nrow = 1) +
    xlab(NULL) +
    # ylab("Standardized_Value") +
    theme(
      plot.margin = unit(c(1, 0, 1, 1), "cm"),
      axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 1, size = 12),
      panel.background = element_blank(),  # blank background
      #panel.border = element_blank(),      # remove border
      axis.line = element_line(colour = "black"),  # black axis lines
      axis.ticks = element_line(color = "black"),  # black axis ticks
      panel.grid.major = element_blank(),  # remove major grid
      panel.grid.minor = element_blank(),  # remove minor grid
      legend.position = "right",  # adjust legend position if needed
      plot.title = element_text(hjust = 0.5)  # center the plot title if you have one
    )
}

pdf(file = "Results_check/metabolites_fasting_cluster.pdf", width = 20, height = 15)
plot_min[[1]]
dev.off()

pdf(file = "Results_check/metabolites_postprandial_cluster.pdf", width = 20, height = 15)
plot_min[[2]]
dev.off()

pdf(file = "Results_check/metabolites_delta_cluster.pdf", width = 20, height = 15)
plot_min[[3]]
dev.off()

###### Radar plot #########
## Radar plot for different clusters
# library(ggpubr)

library(viridis)
#### color themes 
theme_set(theme_classic())

#### color schemes, select 8 clusters
myvc <- scale_color_viridis(discrete=T, option="D")
myvcf <- scale_fill_viridis(discrete=T, option="D")
myvp <- viridis(4, option = "D")

### Other color schemes
myvc <- scale_color_brewer(palette="Set1")
myvcf <- scale_fill_brewer(palette="Set1")

coord_radar <- function(theta = "x", start = 0, direction = 1){
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y"
  else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)}

vars <- c("sexe","leeftijd", "bmim", "glucose1", "choltot1", "HBA1C", "trig1", "fldl1", "homa1IR", "homa1B", "TG_perc", "cluster_fasting",
          "cluster_postprandial","cluster_delta")
temp <- metadata_out %>% select(all_of(vars))
temp[,2:11] <- apply(temp[,2:11], 2, function(x) scale(x)) # z-score, scale function

temp_kmeans <- temp[,c(2:14)] # select three groups
temp_kmeans <- gather(temp_kmeans, key = "name", value = "val", -(cluster_fasting:cluster_delta))
temp_kmeans$cluster_fasting <- as.factor(temp_kmeans$cluster_fasting)
temp_kmeans$cluster_postprandial <- as.factor(temp_kmeans$cluster_postprandial)
temp_kmeans$cluster_delta <- as.factor(temp_kmeans$cluster_delta)

plot_res <- function(data, var1, var_seq){
  
  library(forcats)
  library(magrittr)
  library(ggpubr)
  kmeans_res <- data %>%
    group_by({{var1}}) %>%
    mutate(meanval = sum(val, na.rm = T)) %>%
    arrange(name, desc(meanval)) %>%
    mutate(clu_ord = fct_reorder({{var1}}, meanval))
  
  results_list <- list()
  for (i in var_seq) {
    out_new <-  kmeans_res %>%
      filter({{var1}} == i) %>% 
      group_by(name) %>%
      summarise(median = median(val, na.rm = T),
                IQR    = paste0(sprintf("%.2f",quantile(val, c(0.25), na.rm = T)), "-", 
                                sprintf("%.2f", quantile(val, c(0.75), na.rm = T))),
                mean   = mean(val, na.rm = T),
                sd     = sd(val, na.rm = T)) %>%
      mutate(cluster = i)
    results_list[[i]] <- out_new
  }
  
  # append the results to the list
  results_df_kmeans <- do.call(rbind, results_list)
  results_df_kmeans$cluster <- as.factor(results_df_kmeans$cluster)
  results_df_kmeans <- results_df_kmeans %>% mutate(name = case_when(name == "glucose1" ~ "Glu",
                                                                     name == "bmim" ~ "BMI",
                                                                     name == "leeftijd" ~ "age",
                                                                     name == "choltot1" ~ "TC",
                                                                     name == "HBA1C" ~ "HBA1C",
                                                                     name == "trig1" ~ "TG",
                                                                     name == "fldl1" ~ "LDL",
                                                                     name == "homa1IR" ~ "IR",
                                                                     name == "homa1B" ~ "IS",
                                                                     name == "TG_perc" ~ "Liver(TG%)"))
  
  limits <- c(1.1*min(results_df_kmeans$mean), 1*max(results_df_kmeans$mean))
  
  #seq <- c("TC","TG", "age", "Insulin Resistance", "Insulin Sensitivity", "Liver TG%", "Glu", "LDL","BMI", "HBA1C") 

  results_df_kmeans %<>%
    arrange(cluster, name)
  results_df_kmeans$name <- as.factor(results_df_kmeans$name)
  
  p <- ggplot(results_df_kmeans, aes(x= name, y= mean, fill = cluster))+ ## ascending order in reorder function
    geom_polygon(aes(group = cluster, color = cluster ,fill = cluster),
                 alpha = 0.5, linewidth = 0.5, show.legend = TRUE) +
    geom_point(size = 1) +
    geom_path(aes(color = cluster), alpha = 0.5, linewidth = 0.5) +
    coord_radar()+
    xlab("") + ylab("") +
    # guides(color = F) +
    myvc + myvcf +
    scale_y_continuous(limits=limits) +
    # scale_x_discrete(limits = seq) + # based on specific order of variables
    theme_pubclean() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line(linetype="solid",color="black", linewidth = .1),
          panel.grid.major.y = element_line(linetype="solid",color="black", linewidth = .1),
          panel.spacing = unit(3, "lines"),
          axis.text.x = element_text(size=8,  hjust = 1, vjust = 0),
          plot.margin =  margin(3, 3, 3, 3, "cm"),
          strip.background =element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 16),
          legend.text = element_text(size = 16)) +
    facet_wrap(vars(paste("Cluster", cluster)), nrow = 1) 
  return(p)
}

fasting_Radar <- plot_res(data = temp_kmeans, var1 = cluster_fasting, var_seq = 1:4)
fasting_Radar  
ggsave(filename = "Results_check/fasting_radar_12_18_2.pdf", width = 12, height = 6)

postprandial_Radar <- plot_res(data = temp_kmeans, var1 = cluster_postprandial, var_seq = 1:4)
postprandial_Radar
ggsave(filename = "Results_check/post_radar_12_18_2.pdf", width = 12, height = 6)

delta_Radar <- plot_res(data = temp_kmeans, var1 = cluster_delta, var_seq = 1:4)
delta_Radar
ggsave(filename = "Results_check/delta_radar_12_18_2.pdf", width = 12, height = 6)

# fasting_Radar/postprandial_Radar/delta_Radar
library(cowplot)
pg <- plot_grid(fasting_Radar, postprandial_Radar, delta_Radar, 
                labels = c("Fasting", "Postprandial", "Delta"), ncol = 1)
pg
ggsave(filename = "Results_check/Radar_plot_All_min_2023_12_14_2_2.pdf", width = 25, height = 20)
