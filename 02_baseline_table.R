##### P1-Cluster-NEO-KD 
##### Author: Keyong Deng
##### Department of Clinical Epidemiology, LUMC

pacman::p_load("gtsummary", ## summary statistics and tests
               "rstatix",  ## summary statistics and statistical tests
               "janitor", ## adding totals and percents to tables
               "scales", ## easily convert proportions to percents  
               "flextable", ## converting tables to pretty images
               "skimr", ## get overview of data
               "foreign", "dplyr", "data.table", "Hmisc", "stringr", "ggplot2", "e1071") 

######  Baseline characteristic for the participants ######
load("./After_QC_Imputed_KNN_minimum.Rdata")
load("./k_means_clustering.Rdata")

kmeans_cluster_fasting <- kmeans_cluster_all[[1]] %>% as.data.frame() %>% 
  dplyr::rename(C3 = "V1", C4 = "V2", C5 = "V3")
table(kmeans_cluster_fasting$C3); table(kmeans_cluster_fasting$C4); table(kmeans_cluster_fasting$C5);

kmeans_cluster_post <- kmeans_cluster_all[[2]] %>% as.data.frame() %>% 
  dplyr::rename(C3 = "V1", C4 = "V2", C5 = "V3")
table(kmeans_cluster_post$C3); table(kmeans_cluster_post$C4); table(kmeans_cluster_post$C5);

kmeans_cluster_delta <- kmeans_cluster_all[[3]] %>% as.data.frame() %>% 
  dplyr::rename(C3 = "V1", C4 = "V2", C5 = "V3")
table(kmeans_cluster_delta$C3); table(kmeans_cluster_delta$C4); table(kmeans_cluster_delta$C5);

clusters.fasting <- kmeans_cluster_fasting %>% select(C4); table(clusters.fasting)
clusters.postprandial <- kmeans_cluster_post %>% select(C4); table(clusters.postprandial)
clusters.delta <- kmeans_cluster_delta %>% select(C4); table(clusters.delta)

# Baseline characteristic (divided by different clusters)
# import df2_baseline data
load("./df2_baseline.Rdata") 

metadata_out <- cbind(df2, cluster_fasting = clusters.fasting, cluster_postprandial =  clusters.postprandial, cluster_delta = clusters.delta)
colnames(metadata_out)[554:556] <- c("cluster_fasting", "cluster_postprandial", "cluster_delta")

###### Radar plot #########
# Radar plot for different clusters
pal <- c(unname(yarrr::piratepal("pony")))[c(8,5,3,2)]

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
# scale the variables, z-score 
temp[,2:11] <- apply(temp[,2:11], 2, function(x) scale(x)) 
temp_kmeans <- temp[,c(2:14)] # select three groups

# change to long format
temp_kmeans <- gather(temp_kmeans, key = "name", value = "val", -(cluster_fasting:cluster_delta))

factor_columns <- c("cluster_fasting", "cluster_postprandial", "cluster_delta")
temp_kmeans <- temp_kmeans %>%
  mutate(across(all_of(factor_columns), as.factor))

plot_res <- function(data, var1, var_seq, state_name){
  
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
                                                                     name == "TG_perc" ~ "Liver\nTG%"))
  
  results_df_kmeans$name <- factor(results_df_kmeans$name)
  
  limits <- c(1*min(results_df_kmeans$mean), 1.1*max(results_df_kmeans$mean))

  results_df_kmeans %<>%
    arrange(cluster, name)
  
  p <- ggplot(results_df_kmeans, aes(x= name, y= mean, fill = cluster))+ 
    geom_polygon(aes(group = cluster, color = cluster),
                 alpha = 0.5, linewidth = 0.8, show.legend = T) +
    geom_point(size = 0.5) +
    geom_path(aes(color = cluster), linewidth = 0.8) +
    coord_radar()+
    xlab("") + ylab("") +
    scale_color_manual(values=pal)+
    scale_fill_manual(values=pal)+
    scale_y_continuous(limits=limits) +
    ggtitle(state_name) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line(linetype="solid",color="grey", linewidth = .1),
          panel.grid.major.y = element_line(linetype="solid",color="grey", linewidth = .1),
          panel.spacing = unit(3, "lines"),
          axis.text.x = element_text(size=6),
          plot.margin = unit(c(0,0,0,1), "cm"),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          panel.background = element_rect(fill = "white", color = NA)
          ) +
    facet_wrap(vars(paste("Cluster", cluster)), nrow = 1, scales = "free") 
  return(p)
}

#### Recode the name of Cluster (reordering)
temp_kmeans$cluster_fasting <- ifelse(temp_kmeans$cluster_fasting == 1, 4,
                                      ifelse(temp_kmeans$cluster_fasting == 2, 2, 
                                             ifelse(temp_kmeans$cluster_fasting == 3, 1,
                                                    ifelse(temp_kmeans$cluster_fasting == 4, 3, NA))))

temp_kmeans$cluster_fasting  <- as.factor(temp_kmeans$cluster_fasting)

temp_kmeans$cluster_postprandial <- as.factor(temp_kmeans$cluster_postprandial)

temp_kmeans$cluster_delta <- ifelse(temp_kmeans$cluster_delta == 1, 2,
                                           ifelse(temp_kmeans$cluster_delta == 2, 3, 
                                                  ifelse(temp_kmeans$cluster_delta == 3, 1,
                                                         ifelse(temp_kmeans$cluster_delta == 4, 4, NA))))

temp_kmeans$cluster_delta <- as.factor(temp_kmeans$cluster_delta)

# plot the Radar plot
fasting_Radar <- plot_res(data = temp_kmeans, var1 = cluster_fasting, var_seq = 1:4, state_name = "Fasting status")

postprandial_Radar <- plot_res(data = temp_kmeans, var1 = cluster_postprandial, var_seq = 1:4, state_name = "Postprandial status")

delta_Radar <- plot_res(data = temp_kmeans, var1 = cluster_delta, var_seq = 1:4, state_name = "Delta status")

# fasting_Radar/postprandial_Radar/delta_Radar
# https://github.com/wilkelab/cowplot/blob/master/vignettes/shared_legends.Rmd#L99
library(cowplot)
pg <- plot_grid(fasting_Radar + theme(legend.position="none"), 
                postprandial_Radar + theme(legend.position="none"), 
                delta_Radar + theme(legend.position="none"), 
                labels = c("A", "B", "C"), ncol = 1, label_size = 12)

legend <- get_legend(
  # create some space to the left of the legend
  fasting_Radar + theme(legend.box.margin = margin(0, 0, 0, 12))
) 

# share the same legend for three states
pdf(file = "Results/Figure1.pdf", width = 10, height = 10)
plot_grid(pg, legend, rel_widths = c(3, .4))
dev.off()

# summary function
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
      .groups = "drop" # avoid the message and leave the data in an un-grouped state
    )
}

## recode fasting, postprandial, delta clusters
metadata_out$cluster_fasting <- as.factor(ifelse(metadata_out$cluster_fasting == 1, 4,
                                                 ifelse(metadata_out$cluster_fasting == 2, 2, 
                                                        ifelse(metadata_out$cluster_fasting == 3, 1,
                                                               ifelse(metadata_out$cluster_fasting == 4, 3, NA)))))

metadata_out$cluster_postprandial <- as.factor(metadata_out$cluster_postprandial)

metadata_out$cluster_delta <- ifelse(metadata_out$cluster_delta == 1, 2,
                                     ifelse(metadata_out$cluster_delta == 2, 3, 
                                            ifelse(metadata_out$cluster_delta == 3, 1,
                                                   ifelse(metadata_out$cluster_delta == 4, 4, NA))))

# select vars
vars_base <- c("medication_hypertension", "medication_lipidlower", "sexe", "leeftijd", "bmim",  "glucose1", "Insuline_r1",
               "HBA1C", "homa1IR", "homa1B", "choltot1", "trig1", "hdlc1", "fldl1", "choltot3", "trig3", "hdlc3", "fldl3", "TG_perc",
               "cluster_fasting", "cluster_postprandial", "cluster_delta")

base_df1 <- metadata_out %>% dplyr::select(all_of(vars_base))

cat_vars <- c("sexe", "cluster_fasting", "cluster_postprandial", "cluster_delta", "medication_hypertension", "medication_lipidlower")

base_df1 %<>% dplyr::mutate(sexe = case_when(sexe == "vrouw" ~ 0,  #### 0 denotes female
                                             sexe == "man" ~ 1)) %>% 
  mutate_at(vars(one_of(cat_vars)), funs(factor(.))) ## change those cat_vars into categorical variable

# test the distribution and skewness of variables 
vars_exam <- c("leeftijd", "bmim",  "glucose1", "Insuline_r1","HBA1C", "homa1IR", "homa1B", "choltot1", "trig1", "hdlc1", "fldl1", "choltot3", "trig3", "hdlc3", "fldl3", "TG_perc", "cluster_fasting", "cluster_postprandial", "cluster_delta")
base_df1_exam <- base_df1  %>% select(all_of(vars_exam))
apply(base_df1_exam[,-(17:19)], 2, function(x) skewness(x, na.rm = T)) 
## glucose1, Insuline_r1, HBA1C, homa1IR, homa1B, trig1, trig3, TG_perc are skewed

## str_glue function from stringr is used to combine values from several columns into one new column.
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

as_gt(table1) |> 
  gt::gtsave(filename = "Results/Stable2_fasting.docx")

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
theme_gtsummary_journal(journal = "jama")
table2

as_gt(table2) |> 
  gt::gtsave(filename = "Results/Stable3_postprandial.docx")

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
theme_gtsummary_journal(journal = "jama")
table3

as_gt(table3) |> 
  gt::gtsave(filename = "Results/Stable4_delta.docx")

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

as_gt(table4) |> 
  gt::gtsave(filename = "Results/Table1_total.docx")
