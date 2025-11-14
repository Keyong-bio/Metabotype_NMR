##### P1-Cluster-NEO-KD 
##### Author: Keyong Deng
##### Department of Clinical Epidemiology, LUMC
library(openxlsx)
library(foreign) # data load
library(dplyr) # %>% pipeline
library(data.table) # Fast Read data, fread(), fwrite()
library(Hmisc) # describe function
library(stringr) # string manuiplate
library(ggplot2) # Visualization
library(e1071) # Data skewness 
library(tibble)
library(tidyr)
library(VennDiagram) # plot venn diagram
library(iPVs) #identification of principal variable
library(pacman)
p_load(factoextra, FactorMineR)

load(file = "./After_QC_Imputed_KNN_minimum.Rdata")
load(file = "./After_QC_Imputed_KNN_minimum_Rntrans.Rdata")

############# Part2 iPVS #################

## Here, we can also use no-transformed data to do IPVs, because we apply Spearman correlation using iPVs
set.seed(123)
imputated_ipvs <- lapply(imputated_states, function(x){
  cmat_matrix <- cor(x, method = "spearman")
  result_fasting_res <- iPVs(x, 
                             cor_method = "spearman",
                             dist_method = "R",
                             hclust_meth = "complete",
                             cutheight = 0.5,
                             cmat = cmat_matrix)
})

imputated_ipvs_vars <- list()
for (i in 1:3) {
  imputated_ipvs_vars[[i]] <- imputated_ipvs[[i]]$iPV_table$PVs
}
imputated_ipvs_vars

imputated_ipvs_vars1 <- unlist(imputated_ipvs_vars) %>% as.data.frame() 
imputated_ipvs_vars1 <- rename(imputated_ipvs_vars1, selected_metabolites = .)
write.xlsx(imputated_ipvs_vars1, file = "./Results/STable1.xlsx")

# Revision ------ Variance explained by prinicipal metabolites -------------------

processed_ipv_data <- function(ipv_data, ipv_vars){
  ipv_data["iPV_table"] |> 
    data.frame() |> 
    arrange(desc(iPV_table.clustersize)) |> 
    filter(iPV_table.PVs %in% ipv_vars) |> 
    mutate(iPV_table.PVs = paste0(iPV_table.PVs, " [", iPV_table.clustersize, "]"))
}

processed_ipv_list <- lapply(1:3, function(i){
  processed_ipv_data(imputated_ipvs[[i]], imputated_ipvs_vars[[i]])
})
  
# Total variance explained
## write a quick function
est_var_exp = function(ipvs_object, wdata){
  
  cluster_variance = sapply(ipvs_object$PVresults, function(x){
    
    n = x$variable
    if(length(n) == 1){
      totalvar = var(wdata[, n ], na.rm = TRUE)
    } else {
      totalvar =  sum(apply( wdata[, n ], 2, function(y){ 
        var(y, na.rm = TRUE) 
      }))
    }
    return(totalvar)
  })
  
  ## Total cluster variance
  total_cluster_variance = sum(cluster_variance)
  
  ## explained variance: variance explained by each PV
  explained_variance = ipvs_object$iPV_table$VarExp_by_PV
  
  ## Total explained variance
  total_explained_variance = sum(cluster_variance * explained_variance)
  
  ## Proportion of total variance explained
  proportion_total_explained = total_explained_variance / total_cluster_variance
  
  return(proportion_total_explained)
  
}

x1 = est_var_exp(ipvs_object = imputated_ipvs[[1]], wdata = imputated_states[[1]])
x2 = est_var_exp(ipvs_object = imputated_ipvs[[2]], wdata = imputated_states[[2]])
x3 = est_var_exp(ipvs_object = imputated_ipvs[[3]], wdata = imputated_states[[3]])

# build a dataframe to plot
df_var <- data.frame(
  stats = c("Fasting",
            "Postprandial",
            "Delta"),
  PVs = c(16, 16, 14),
  Variance_explained = c(x1, x2, x3)
) |> 
  mutate(stats = paste0(stats, " [", PVs, "]"))

processed_ipv_list[[4]] <- df_var
names(processed_ipv_list) <- c("Fasting", "Postprandial", "Delta", "Total")

library(viridis)
p_list <- lapply(1:3, function(i){
  
  ggplot(data = processed_ipv_list[[i]], 
         aes(x = reorder(iPV_table.PVs, -iPV_table.VarExp_by_PV) , y = iPV_table.VarExp_by_PV, fill = as.factor(iPV_table.PVs))) + 
    geom_bar(stat="identity", show.legend = F) +
    geom_text(aes(label = round(iPV_table.VarExp_by_PV, 1)),
              vjust = -0.5) + # adjust vertical position of the text
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme_bw() +
    scale_fill_viridis(discrete = TRUE, option = "cividis") +
    labs(x = "", y = "Variance explained for each cluster",
         title = names(processed_ipv_list)[[i]])
})

p_list[[1]]  
p_list[[2]]    
p_list[[3]]  
p_list[[4]] <- ggplot(data = processed_ipv_list[[4]], 
       aes(x = factor(stats, levels = c("Fasting [16]", "Postprandial [16]", "Delta [14]")) , y = Variance_explained, fill = as.factor(stats))) + 
  geom_bar(stat="identity", show.legend = F) +
  geom_text(aes(label = round(Variance_explained, 1)),
            vjust = -0.5) + # adjust vertical position of the text
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_bw() +
  scale_fill_viridis(discrete = TRUE, option = "cividis") +
  labs(x = "", y = "Variance explained for each cluster",
       title = names(processed_ipv_list)[[4]])

# cowplot to combined all the figures
# library(cowplot)
# combined_plot_cow <- plot_grid(p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]],
#                                labels = "AUTO", # A, B, C, D
#                                label_size = 12,
#                                label_fontface = "bold",
#                                ncol = 2,
#                                nrow = 2,
#                                rel_widths = c(2, 2),
#                                rel_heights = c(2, 2))
# print(combined_plot_cow)
library(patchwork)
p4_spacer <-  p_list[[4]] + plot_layout(widths = c(1,1)) # the 4th will share the half of the width
combine_plot <- (p_list[[1]]|p_list[[2]])/(p_list[[3]]|p4_spacer) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold", size = 12))
print(combine_plot)

pdf(file = "./Results/Variance_explained.pdf", width = 12, height = 12)
print(combine_plot)
dev.off()

# Revision 1 ---- Venn gram -------------------
imputated_ipvs_vars1_split <- imputated_ipvs_vars1 %>% 
  separate(col = selected_metabolites,
           into = c("Metabolite", "Status"),
           sep = "_",
           remove = FALSE) %>% 
  mutate(Status = case_when(
    Status == 1 ~ "Fasting",
    Status == 3 ~ "Postprandial",
    Status == "rp" ~ "Delta"
  ))

# Intersections
A = imputated_ipvs_vars1_split %>% filter(Status == "Fasting") %>% select(Metabolite) %>% unlist()
B = imputated_ipvs_vars1_split %>% filter(Status == "Postprandial") %>% select(Metabolite) %>% unlist()
C = imputated_ipvs_vars1_split %>% filter(Status == "Delta") %>% select(Metabolite) %>% unlist()

onlyA = setdiff(A, union(B, C))
onlyB = setdiff(B, union(A, C))
onlyC = setdiff(C, union(A, B))
AB = intersect(A, B) %>% setdiff(C)
AC = intersect(A, C) %>% setdiff(B)
BC = intersect(B, C) %>% setdiff(A)
ABC = intersect(intersect(A, B), C)

# generate the venn diagram
venn.diagram(
  x = list(
    imputated_ipvs_vars1_split %>% filter(Status == "Fasting") %>% select(Metabolite) %>% unlist(),
    imputated_ipvs_vars1_split %>% filter(Status == "Postprandial") %>% select(Metabolite) %>% unlist(),
    imputated_ipvs_vars1_split %>% filter(Status == "Delta") %>% select(Metabolite) %>% unlist()
  ),
  category.names = c("Fasting", "Postprandial", "Delta"),
  filename = "./Results/Venn_revision.png",
  output = TRUE,
  imagetype = "png",
  height = 800 , 
  width = 800 , 
  resolution = 500,
  lwd = 1,
  col = c("#1794CEFF", '#6351A0FF', '#FBBB68FF'),
  fill = c(alpha("#1794CEFF",0.8), alpha('#6351A0FF',0.8), alpha('#FBBB68FF',0.8)),
  cat.col = c("#1794CEFF", '#6351A0FF', '#FBBB68FF'), # set the color of labels
  cex = 0.5,
  cat.cex = 0.3,
  cat.default.pos = "outer",
  rotation = 1
)
# add ABC into the venn diagram

# select the working imputed dataframe
working_imputated_df <- list()
states <- c("fasting_status", "postprandial_status", "delta_status")
for (i in 1:3) {
  tempdf <- imputated_states_transformation[[i]] %>% as.data.frame()
  tempvar <- imputated_ipvs[[i]]$iPV_table$PVs 
  # selected PV and name the working_imputed_df group
  working_imputated_df[[states[[i]]]] <- tempdf %>% select(all_of(tempvar))
}

# Check the variance explain for selected metabolites - Revision ---------------------------

# Run PCA 
PCA_fasting <- FactoMineR::PCA(imputated_states_transformation[[1]], 
                               scale.unit = FALSE,
                               ncp = 10, graph = TRUE)

# var_contrib <- get_pca_var(PCA_fasting)$contrib[, 1:2]
# 
# # Create a dataframe for plotting
# contrib_df <- data.frame(
#   Variable = rownames(var_contrib),
#   PC1 = var_contrib[, 1],
#   PC2 = var_contrib[, 2]
# )
# 
# # Reshape to long format for ggplot
# contrib_long <- tidyr::pivot_longer(
#   contrib_df, 
#   cols = c(PC1, PC2),
#   names_to = "Component",
#   values_to = "Contribution"
# ) %>% 
#   filter(Variable %in% colnames(working_imputated_df[[1]]))
# 
# # Create a barplot
# ggplot(contrib_long, aes(x = reorder(Variable, Contribution), y = Contribution, fill = Component)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   coord_flip() +
#   scale_fill_manual(values = c("PC1" = "#00AFBB", "PC2" = "#FC4E07")) +
#   labs(
#     title = "Contribution of Variables to PC1 and PC2",
#     subtitle = "Fasting Metabolites",
#     x = "Variables",
#     y = "Contribution (%)"
#   ) +
#   theme_minimal()

# # Visualize eigenvalues/variances (Scree plot)
# fviz_eig(PCA_fasting, addlabels = TRUE)

variance_fasting <- fviz_pca_var(PCA_fasting,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(name = colnames(working_imputated_df[[1]]))) +
  ggtitle("Fasting_Metabolites")

# Also for the other selected metabolites at postprandial
PCA_postprandial <- FactoMineR::PCA(imputated_states_transformation[[2]], 
                               scale.unit = FALSE,
                               ncp = 10, graph = TRUE)
variance_postprandial <- fviz_pca_var(PCA_postprandial,
             col.var = "contrib",
             fill.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(name = colnames(working_imputated_df[[2]])),
             gradient.vars = c(0, 1)
             ) +
  ggtitle("Postprandial_Metabolites")

# Also for the other selected metabolites at delta
PCA_delta <- FactoMineR::PCA(imputated_states_transformation[[3]], 
                                    scale.unit = FALSE,
                                    ncp = 10, graph = TRUE)

variance_delta <- fviz_pca_var(PCA_delta,
             col.var = "contrib",
             # fill.var = "contrib",
             # gradient.cols = c("blue", "white", "red"),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(name = colnames(working_imputated_df[[3]])),
             gradient.vars = c(0, 0.5)) +
  ggtitle("Delta_Metabolites")

# plot together
library(cowplot)

# SFigure11
png(filename = "Results/Revision_variance_explained.png", 
    width = 12,       
    height = 8,       
    res = 500,
    units = "in")     # Using inches instead of pixels

plot_grid(variance_fasting, variance_postprandial, variance_delta, labels = c('A', 'B', "C"), 
          label_fontface = "bold",
          label_size = 16)

dev.off()


# already exists
# save(working_imputated_df, file = "Results/working_imputated_df.Rdata")

load(file = "Results/working_imputated_df.Rdata")
# after run the 00_Preparation_revised.R, get df2, extract the sexe
# add sex to each dataframe
working_imputated_df_new <- lapply(working_imputated_df, function(x) {
  x$sexe <- df2$sexe
  return(x)
})

source("./Code/ConsensusCluster_function.R")

male_fasting = subset(working_imputated_df_new[[1]], sexe == "man")[,-ncol(working_imputated_df_new[[1]])]
female_fasting = subset(working_imputated_df_new[[1]], sexe == "vrouw")[,-ncol(working_imputated_df_new[[1]])]

res_km_male_fasting <- ConsensusClusterPlus(t(male_fasting),
                                       maxK = 8,
                                       reps = 100,
                                       pItem = 0.8, ## percent of sample to consider for clustering
                                       pFeature = 1, ## percent of features to consider
                                       plot = "png",
                                       clusterAlg = "km",
                                       title = "km_male_fasting_status",
                                       distance = "euclidean",
                                       seed = 2023,
                                       verbose = T)    ## 4 clusters?
tryCatch({
  calcICL(res_km_male_fasting, title = "km_male_fasting_status_2024_value", plot = "pdf")
}, error = function(e) {
  # Handling the error
  message("An error occurred: ", e$message)
})

res_km_female_fasting <- ConsensusClusterPlus(t(female_fasting),
                                            maxK = 8,
                                            reps = 100,
                                            pItem = 0.8, ## percent of sample to consider for clustering
                                            pFeature = 1, ## percent of features to consider
                                            plot = "png",
                                            clusterAlg = "km",
                                            title = "km_female_fasting_status",
                                            distance = "euclidean",
                                            seed = 2023,
                                            verbose = T)    ## 4 clusters?
tryCatch({
  calcICL(res_km_female_fasting, title = "km_female_fasting_status_2024_value", plot = "pdf")
}, error = function(e) {
  # Handling the error
  message("An error occurred: ", e$message)
})

## Sex difference (postprandial status) -------------------------
male_postprandial = subset(working_imputated_df_new[[2]], sexe == "man")[,-ncol(working_imputated_df_new[[2]])]
female_postprandial = subset(working_imputated_df_new[[2]], sexe == "vrouw")[,-ncol(working_imputated_df_new[[2]])]

res_km_male_postprandial <- ConsensusClusterPlus(t(male_postprandial),
                                            maxK = 8,
                                            reps = 100,
                                            pItem = 0.8, ## percent of sample to consider for clustering
                                            pFeature = 1, ## percent of features to consider
                                            plot = "png",
                                            clusterAlg = "km",
                                            title = "km_male_postprandial_status",
                                            distance = "euclidean",
                                            seed = 2023,
                                            verbose = T)    ## 4 clusters?
tryCatch({
  calcICL(res_km_male_fasting, title = "km_male_postprandial_status_2024_value", plot = "pdf")
}, error = function(e) {
  # Handling the error
  message("An error occurred: ", e$message)
})

res_km_female_postprandial <- ConsensusClusterPlus(t(female_postprandial),
                                              maxK = 8,
                                              reps = 100,
                                              pItem = 0.8, ## percent of sample to consider for clustering
                                              pFeature = 1, ## percent of features to consider
                                              plot = "png",
                                              clusterAlg = "km",
                                              title = "km_female_postprandial_status",
                                              distance = "euclidean",
                                              seed = 2023,
                                              verbose = T)    ## 4 clusters?
tryCatch({
  calcICL(res_km_female_fasting, title = "km_female_postprandial_status_2024_value", plot = "pdf")
}, error = function(e) {
  # Handling the error
  message("An error occurred: ", e$message)
})

### Scatter plots for Selected metabolite and LDL/TG ----------
# Fasting status MVLDLTG_1/IDLP_1/SHDLFC_1
# Postprandial status MVLDLP_3/LLDLP_3/SHDLC_3/LDLD_3
# Delta status MVLDLPL_rp/LLDLCE_rp/MHDLFC_rp/XLHDLC_rp/LHDLTG_rp

cols_select_fasting <- c("hdlc1", "choltot1", "trig1", "fldl1", "MVLDLTG_1", "IDLP_1", "SHDLFC_1")
cols_select_fasting <- df2[,cols_select_fasting]

cols_select_postprandial <- c("MVLDLP_3", "LLDLP_3", "SHDLC_3", "LDLD_3", "trig3", "fldl3", "trig1", "fldl1")
cols_select_postprandial <- df2[,cols_select_postprandial]

# plot the scatter plot for variables-revision -----------
# Fit linear model

plot_regression <- function(data, x_var, y_var, xlab_label = NULL, ylab_label = NULL, 
                            color = NULL, title = NULL, my_color = NULL) {
  
  # exclude the rows with NA for x and y
  complete_df <- data[complete.cases(data[[x_var]], data[[y_var]]),]
  
  # fit model
  formula <- as.formula(paste(y_var, "~", x_var))
  
  model <-  lm(formula, data = complete_df)
  summary_model <- summary(model)
  
  # get statistic
  r_squared <- round(summary_model$r.squared, 3)
  p_value <- format.pval(summary_model$coefficients[2, 4], digits = 3)
  r <- sign(model$coefficients[2]) * sqrt(r_squared)
  r <- round(r, 3)
  
  # if(is.null(title)){
  #   title <- paste("Scatter_plot:", y_var, "vs", x_var)
  # }
  
  # Create the plot
  p <- ggplot(complete_df, aes_string(x = x_var, y = y_var)) +
    geom_point(color = color) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    theme_classic() +
    # labs(title = title,
    #      x = x_var, y = y_var) +
    xlab(xlab_label) +
    ylab(ylab_label) +
    annotate("text",
             x = min(complete_df[[x_var]]) + (max(complete_df[[x_var]])-min(complete_df[[x_var]])) * 0.1,
             y = max(complete_df[[y_var]]) * 0.9,
             label = paste(paste("R2 =", r_squared),
                           paste("p-value ="), p_value),
             hjust = 0 # left aligned text
             )
  return(p)
}

p1 <- plot_regression(cols_select_fasting, x_var = "MVLDLTG_1", y_var = "trig1", 
                      color = "darkblue", xlab_label = "MVLDLTG_Fasting", ylab_label = "TG_Fasting") + ggtitle("Fasting")
p1

p2 <- plot_regression(cols_select_fasting, x_var = "IDLP_1", y_var = "fldl1", 
                      color = "olivedrab",  xlab_label = "IDLP_Fasting", ylab_label = "LDL_Fasting") + ggtitle("Fasting")
p2

p3 <- plot_regression(cols_select_postprandial, x_var = "MVLDLP_3", y_var = "trig3", 
                      color = "coral", xlab_label = "MVLDLP_Postprandial", ylab_label = "TG_Postprandial") + ggtitle("postprandial")
p3

p4 <- plot_regression(cols_select_postprandial, x_var = "LLDLP_3", y_var = "fldl3", 
                      color = "indianred", xlab_label = "LLDLP_Postprandial", ylab_label = "LDL_Postprandial") + ggtitle("postprandial")
p4

pdf("Results/Correlation_revision.pdf", width = 10, height = 10)     # Using inches instead of pixels

plot_grid(p1, p2, p3, p4, labels = c('A', 'B', "C", "D"), 
          label_fontface = "bold",
          label_size = 16)

dev.off()


#### Consensuscluster analysis to determine the Optimal cluster number #########
# detach("package:ConsensusClusterPlus", unload = T)

# Due to the heavy computation consumption, I recommend skipping this part first
# Directly check the stability of clusters

# this section for consensus_clustering, you can skip it due to longer time to run
# source("ConsensusCluster_function.R")
# res_km_fasting <- ConsensusClusterPlus(t(working_imputated_df[[1]]),
#                                        maxK = 8,
#                                        reps = 100,
#                                        pItem = 0.8, ## percent of sample to consider for clustering
#                                        pFeature = 1, ## percent of features to consider
#                                        plot = "png",
#                                        clusterAlg = "km",
#                                        title = "km_fasting_status",
#                                        distance = "euclidean",
#                                        seed = 2023,
#                                        verbose = T)    ## 4 clusters?
# tryCatch({
#   calcICL(res_km_fasting, title = "km_fasting_status_2024_value", plot = "pdf")
# }, error = function(e) {
#   # Handling the error
#   message("An error occurred: ", e$message)
# })
# 
# ## postprandial ####
# # Sys.setenv(R_QTRANS_LIMIT = 500000)  # Setting to a higher value
# res_km_post <- ConsensusClusterPlus(t(working_imputated_df[[2]]),
#                                     maxK = 8,
#                                     reps = 100,
#                                     pItem = 0.8, ## percent of sample to consider for clustering
#                                     pFeature = 1, ## percent of features to consider
#                                     plot = "png",
#                                     clusterAlg = "km",
#                                     title = "km_postprandial_status",
#                                     distance = "euclidean",
#                                     seed = 2023,
#                                     verbose = T) 
# tryCatch({
#   calcICL(res_km_post, title = "km_postprandial_status_2024_value", plot = "pdf")
# }, error = function(e) {
#   # Handling the error
#   message("An error occurred: ", e$message)
# })
# 
# ## delta #######
# # Sys.setenv(R_QTRANS_LIMIT = 500000)  # Setting to a higher value
# res_km_delta <- ConsensusClusterPlus(t(working_imputated_df[[3]]),
#                                      maxK = 8,
#                                      reps = 100,
#                                      pItem = 0.8, ## percent of sample to consider for clustering
#                                      pFeature = 1, ## percent of features to consider
#                                      plot = "png",
#                                      clusterAlg = "km",
#                                      title = "km_delta_status",
#                                      distance = "euclidean",
#                                      seed = 2023,
#                                      verbose = T) 
# tryCatch({
#   calcICL(res_km_delta, title = "km_delta_status_2024_value", plot = "pdf")
# }, error = function(e) {
#   # Handling the error
#   message("An error occurred: ", e$message)
# })


############ Check the stability of the selected metabolites (bootstrapping) ##########
rand_ipvs <- function(input_data, resample_time, resample_size_factor,  replacement = T) {
  
  # list for storing output
  out_samp <- vector(mode = "list", length = resample_time)
  
  for (i in 1:resample_time) {
    # select rows to use in new bootstrap samples
    new_rows <- sample(nrow(input_data), floor(nrow(input_data)*resample_size_factor), replace = replacement)
    boot_data <- input_data[new_rows,]
    
    # using bootstrap data to do IPVs
    cmat_delta <- cor(boot_data, method = "spearman")
    result_delta_res <- iPVs(boot_data, 
                             cor_method = "spearman",
                             dist_method = "R",
                             hclust_meth = "complete",
                             cutheight = 0.5,
                             cmat = cmat_delta)
    out_samp[[i]] <- result_delta_res$iPV_table
  }
  return(out_samp)
} 

ipvs_resample <- lapply(imputated_states, function(x){
  rand_ipvs(x, resample_time = 200, resample_size_factor = 0.75, replacement = T)
})

# save(ipvs_resample, file = "./Results/ipvs_resample.Rdata")
load("./Results/ipvs_resample.Rdata")
plot_resample <- list()
for (i in 1:3) {
  count_pvs <- NULL
  # i = 1
  ipvs_resample_1 <- ipvs_resample[[i]]
  for (j in 1:200) {
    temp <- ipvs_resample_1[[j]]$PVs
    count_pvs <- c(count_pvs, temp)
  }
  ipvs_resample_freq <- table(count_pvs) %>% as.data.frame() %>% arrange(Freq) %>% 
    dplyr::filter(count_pvs %in% c(imputated_ipvs_vars[[1]],
                                   imputated_ipvs_vars[[2]],
                                   imputated_ipvs_vars[[3]]))
  library(viridis)
  plot_resample[[i]] <- ggplot(data = ipvs_resample_freq, 
                               aes(x = reorder(count_pvs, -Freq), y = Freq, fill = as.factor(count_pvs))) + 
    geom_bar(stat="identity", show.legend = F) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme_bw() +
    scale_color_viridis(discrete = TRUE, option = "magma") +
    geom_hline(yintercept = 50, linetype=2, color = "grey") +
    geom_hline(yintercept = 100, linetype=2) +
    labs(x = "Metabolites Name", y = "Frequency for selected metabolites (Bootstrapping)",
         title = names(ipvs_resample)[i])
}

# random ipvs (SFigure3)
library(patchwork)
pdf(file = "./Results/SFigure3.pdf", width = 20, height = 6)
plot_resample[[1]] | plot_resample[[2]] | plot_resample[[3]]
dev.off()

##### Alluvial plot to see the change of the assignment ######

### K-means clustering step ####
# library(amap)
# kmeans_cluster_all <- list()
# for (i in 1:3) {
#   # set.seed(123)
#   new_cluster_kmeans <- c()
#   for (j in 3:5) {
#     # set.seed(123)
#     kmeans_cluster = amap::Kmeans(working_imputated_df[[i]], centers = j, iter.max = 200, nstart = 25)
#     new_cluster_kmeans <- cbind(new_cluster_kmeans, kmeans_cluster$cluster)
#   }
#   kmeans_cluster_all[[i]] <- new_cluster_kmeans
# }

set.seed(123)
kmeans_cluster_fasting = amap::Kmeans(working_imputated_df[[1]], 
                              centers = 4, iter.max = 200, nstart = 25)
table(kmeans_cluster_fasting$cluster)
kmeans_cluster_fasting$centers

set.seed(1234)
kmeans_cluster_postprandial = amap::Kmeans(working_imputated_df[[2]], 
                                           centers = 4, iter.max = 200, nstart = 25)
table(kmeans_cluster_postprandial$cluster)
kmeans_cluster_postprandial$centers

save(kmeans_cluster_fasting, kmeans_cluster_postprandial, file = "Mapping_file.Rdata")

## Run the cluster centers
# set.seed(1234)
# kmeans_cluster_fast <- kmeans(x = working_imputated_df[[1]], 
#                                  centers = 4, 
#                                  algorithm = "Hartigan-Wong",
#                                  iter.max = 250, 
#                                  nstart = 100)
# table(kmeans_cluster_fast$cluster)
# kmeans_cluster_fast$centers

# Manual center calculation for cluster 1
mydf <- working_imputated_df[[1]]
order <- kmeans_cluster_all[[1]][,2]
cluster_1_points <- mydf[order == 1,]
center_1 <- colMeans(cluster_1_points)

# kmeans_cluster_fast$centers
# set.seed(123456)
# kmeans_cluster_post <- kmeans(x = working_imputated_df[[2]], 
#                                  centers = 4, 
#                                  iter.max = 250, 
#                                  nstart = 50)
# table(kmeans_cluster_post$cluster)

# kmeans_cluster_post$centers
# save(kmeans_cluster_all, file = "Data/k_means_clustering.Rdata")


# Calculate distances from new points to existing centers
assign_to_clusters <- function(new_data, centers) {
  distances <- apply(new_data, 1, function(point) {
    apply(centers, 1, function(center) {
      # Euclidean distance
      sqrt(sum((point - center)^2)) 
    })
  })
  
  # Assign each point to nearest center
  apply(distances, 2, which.min)
}

# for the new data
new_clusters <- assign_to_clusters(new_dataset, original_kmeans$centers)


# for reproduce the results, you can directly loading in the k_mean_clustering.Rdata
# rm(list = ls())
load(file = "./k_means_clustering.Rdata")
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

## Re-label function
re_label <- function(old.G, new.G, alluvial.df){
  new.order <- rep(new.G, new.G)
  old.clusters <- subset(alluvial.df, G == old.G)
  new.clusters <- subset(alluvial.df, G == new.G)
  
  for (g in old.G:1) {
    ids <- subset(old.clusters, Cluster == g)$ID
    for (new.g in 1:new.G) {
      new.clusters.g <- subset(new.clusters, new.g == Cluster)
      if (nrow(subset(new.clusters.g, ID %in% ids)) > 0.5 * length(ids)) {
        new.order[new.g] <- g 
      }
    }
  }
  
  alluvial.df[alluvial.df[, "G"] == new.G, "Cluster"] <- plyr::mapvalues(alluvial.df[alluvial.df[, "G"] == new.G, "Cluster"],
                                                                         from = seq(1, new.G), new.order)
  return(alluvial.df)
}

# convert to alluvial format
alluvial.df.fast <-cbind(kmeans_cluster_fasting[, c(1,4)], G = 3)
colnames(alluvial.df.fast) <- c("Cluster", "ID", "G")
temp1 <- cbind(kmeans_cluster_fasting[, c(2,4)], G = 4)
colnames(temp1) <- c("Cluster", "ID", "G")
temp2 <- cbind(kmeans_cluster_fasting[, c(3,4)], G = 5)
colnames(temp2) <- c("Cluster", "ID", "G")

alluvial.df.fast <- rbind(alluvial.df.fast, temp1, temp2)
alluvial.df.fast$ID <- as.character(alluvial.df.fast$ID)
alluvial.df.fast$Cluster <- as.factor(alluvial.df.fast$Cluster)
alluvial.df.fast$G <- as.factor(alluvial.df.fast$G)

alluvial.df.fast <- re_label(3, 4, alluvial.df.fast)
alluvial.df.fast <- re_label(4, 5, alluvial.df.fast)

## postprandial
alluvial.df.post <-cbind(kmeans_cluster_post[, c(1,4)], G = 3)
colnames(alluvial.df.post) <- c("Cluster", "ID", "G")
temp1 <- cbind(kmeans_cluster_post[, c(2,4)], G = 4)
colnames(temp1) <- c("Cluster", "ID", "G")
temp2 <- cbind(kmeans_cluster_post[, c(3,4)], G = 5)
colnames(temp2) <- c("Cluster", "ID", "G")

alluvial.df.post <- rbind(alluvial.df.post, temp1, temp2)
alluvial.df.post$ID <- as.character(alluvial.df.post$ID)
alluvial.df.post$Cluster <- as.factor(alluvial.df.post$Cluster)
alluvial.df.post$G <- as.factor(alluvial.df.post$G)

alluvial.df.post <- re_label(3, 4, alluvial.df.post)
alluvial.df.post <- re_label(4, 5, alluvial.df.post)

## Delta
alluvial.df.delta <-cbind(kmeans_cluster_delta[, c(1,4)], G = 3)
colnames(alluvial.df.delta) <- c("Cluster", "ID", "G")
temp1 <- cbind(kmeans_cluster_delta[, c(2,4)], G = 4)
colnames(temp1) <- c("Cluster", "ID", "G")
temp2 <- cbind(kmeans_cluster_delta[, c(3,4)], G = 5)
colnames(temp2) <- c("Cluster", "ID", "G")

alluvial.df.delta <- rbind(alluvial.df.delta, temp1, temp2)
alluvial.df.delta$ID <- as.character(alluvial.df.delta$ID)
alluvial.df.delta$Cluster <- as.factor(alluvial.df.delta$Cluster)
alluvial.df.delta$G <- as.factor(alluvial.df.delta$G)

alluvial.df.delta <- re_label(3, 4, alluvial.df.delta)
alluvial.df.delta <- re_label(4, 5, alluvial.df.delta)

# install.packages("ggalluvial")
library(ggalluvial) 

### Fasting/postprandial/delta status #######
status_plot <- c("alluvial.df.fast", "alluvial.df.post", "alluvial.df.delta")
title_name <- c("Fasting_state", "Postprandial_state", "Delta_state")

alluvial_p <- function(status, title_name) {
    
    # color theme
    pal <- c(unname(yarrr::piratepal("pony")))[c(8,5,3,2,1,9)]
    
    # plot the figure, get function change characters into object
    p <- ggplot(get(status), 
                aes(x = G,
                    stratum = Cluster,
                    alluvium = ID,
                    fill = Cluster,
                    label = Cluster)) + 
      scale_x_discrete(expand = c(.1, .1)) +
      geom_flow()  +
      geom_stratum() +
      geom_text(stat = "stratum", size = 8) +
      theme_classic() +
      theme(axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            plot.title =  element_text(size = 20)) +
      theme(legend.position = "none") +
      ggtitle(title_name) +
      scale_fill_manual(values = pal) +
      xlab("Assumed number of clusters") + 
      ylab("Frequency")
    
    ggsave(filename = paste0("./Results/",status, "alluvial.png"),  # File name
           plot = p,                  # Plot object to save
           width = 8,                 # Width in inches
           height = 6,                # Height in inches
           dpi = 300,                 # Resolution in dots per inch
           bg = "white"               # Background color
    )
    
    return(p)
  
}

alluvial_p1 <- alluvial_p(status_plot[[1]], title_name[[1]])
alluvial_p2 <- alluvial_p(status_plot[[2]], title_name[[2]])
alluvial_p3 <- alluvial_p(status_plot[[3]], title_name[[3]])

library(gridExtra)
p_alluvial <- grid.arrange(alluvial_p1, alluvial_p2, alluvial_p3,
                           nrow = 2)

# SFigure7
png(file = "./Results/SFigure7.png", width = 20, height = 20, units = "in", res = 300)
grid.arrange(alluvial_p1, alluvial_p2, alluvial_p3,
             nrow = 2)
dev.off()

###### Stability evaluation ###########
# Stability of cluster analysis
library(fpc)                           
res_stable_km_all <- list()

k <- c(4, 4, 4) ## number of clusters after ConsensusCluster function
for (i in 1:3) {
  res_stable_km <- clusterboot(working_imputated_df[[i]], B = 100, bootmethod = c("boot"), 
                               clustermethod = kmeansCBI,
                               seed = 123,
                               k = k[[i]]) 
  res_stable_km_all[[i]] <- res_stable_km$bootresult %>% t() %>% data.frame() 
}

# calculate the mean
res_stable_km_all2 <- do.call(cbind,res_stable_km_all)
apply(res_stable_km_all2, 2, mean) 

# Using gather function to change wide data to long / or other ways to do
pal <- c(unname(yarrr::piratepal("pony")))[c(8,5,3,2,1,9)]

# Plot the bar plots of Jaccard index ####
states2 <- c("fasting_state", "postprandial_state", "delta_state")
Jaccard_KM <- list()
for (i in 1:3) {
  
  temp <- res_stable_km_all[[i]] %>%  dplyr::rename(Cluster1 = X1,
                                                    Cluster2 = X2,
                                                    Cluster3 = X3,
                                                    Cluster4 = X4)
  
  temp_long <- tidyr::gather(temp, Cluster, Value, Cluster1:Cluster4)
  
  # Calculate mean, sd, se, and IC
  temp_long <- temp_long %>% 
    group_by(Cluster) %>% 
    summarise(n = n(),
              mean = mean(Value),
              sd = sd(Value)) %>% 
    mutate(se = sd/sqrt(n)) %>%
    mutate(ic = se * qt((1-0.05)/2 + .5, n-1))
  
  # Barplot
  Jaccard_KM[[i]] <- ggplot(temp_long, aes(x = Cluster, y = mean, fill = Cluster)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_errorbar(aes(x = Cluster, 
                      ymin = mean-ic, 
                      ymax = mean+ic), width = 0.2, color = "black", alpha = 0.6, linewidth = 1) +
    labs(x = NULL, y = paste0("Jaccard Index_", states2[[i]])) +
    theme_bw() +
    scale_fill_manual(values = pal) +
    # scale_color_brewer(palette = "Dark2") +
    scale_y_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.75, linetype = "dashed", col = "grey") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          # Legend text size
          legend.text = element_text(size = 12)
         
    )
  ggsave(filename = paste0("./Results/", states2[[i]], "_Jaccard_index.pdf"), width = 10, height = 8)
}

# combine three plots jaccard index (SFigure8)
library(gridExtra)
library(grid)

add_label <- function(plot, label) {
  arrangeGrob(plot, top = textGrob(label, x = unit(0.1, "npc"), 
                                   y = unit(0.95, "npc"),
                                   just = c("left", "top"),
                                   # change the fontsize and fontface 
                                   gp = gpar(fontsize = 20, fontface = "bold")))
}

Jaccard1 <- Jaccard_KM[[1]] + 
  theme(legend.position = "top",
        legend.title = element_text(size = 16, face = "bold")) +
  labs(fill = expression(Cluster["fasting"])) # put the fasting at subscripts

Jaccard2 <- Jaccard_KM[[2]] + 
  theme(legend.position = "top",
        legend.title = element_text(size = 16, face = "bold")) +
  labs(fill = expression(Cluster["postprandial"]))

Jaccard3 <- Jaccard_KM[[3]] + 
  theme(legend.position = "top",
        legend.title = element_text(size = 16, face = "bold")) +
  labs(fill = expression(Cluster["delta"]))

# SFigure8
png(file = "./Results/SFigure8.png", width = 15, height = 15, units = "in", res = 300)
grid.arrange(add_label(Jaccard1, "A"), add_label(Jaccard2, "B"), add_label(Jaccard3, "C"), nrow = 2)
dev.off()
