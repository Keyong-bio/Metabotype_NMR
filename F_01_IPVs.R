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

working_imputated_df <- list()
for (i in 1:3) {
  tempdf <- imputated_states_transformation[[i]] %>% as.data.frame()
  tempvar <- imputated_ipvs[[i]]$iPV_table$PVs 
  working_imputated_df[[i]] <- tempdf %>% select(all_of(tempvar))
}

names(working_imputated_df) <- c("fasting_status", "postprandial_status", "delta_status")

for (i in 1:3) {
  write.csv(working_imputated_df[[i]], file = paste0("./Results_check/", names(working_imputated_df)[[i]],  "_iPVs_2024_0329.csv"), row.names = F)
}

### If we used the imputed_minimum data
set.seed(123)
imputated_ipvs_minimum <- lapply(imputated_states_minimum, function(x){
  cmat_matrix <- cor(x, method = "spearman")
  result_fasting_res <- iPVs(x, 
                             cor_method = "spearman",
                             dist_method = "R",
                             hclust_meth = "complete",
                             cutheight = 0.5,
                             cmat = cmat_matrix)
})

imputated_ipvs_vars_minimum <- list()
for (i in 1:3) {
  imputated_ipvs_vars_minimum[[i]] <- imputated_ipvs_minimum[[i]]$iPV_table$PVs
}
imputated_ipvs_vars_minimum

working_imputated_df_minimum <- list()
for (i in 1:3) {
  tempdf <- imputated_states_transformation_minimum[[i]] %>% as.data.frame()
  tempvar <- imputated_ipvs_minimum[[i]]$iPV_table$PVs 
  working_imputated_df_minimum[[i]] <- tempdf %>% select(all_of(tempvar))
}

names(working_imputated_df_minimum) <- c("fasting_status", "postprandial_status", "delta_status")

for (i in 1:3) {
  write.csv(working_imputated_df_minimum[[i]], file = paste0("./Results_check/", names(working_imputated_df_minimum)[[i]],  "_iPVs_minimum_2024_0329.csv"), row.names = F)
}

# working_imputated_df_minimum_2 <- do.call(cbind,working_imputated_df_minimum)


########## Check the stability of the selected metabolites (bootstrapping) ##########
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

plot_resample <- list()
for (i in 1:3) {
  count_pvs <- NULL
  # i = 1
  ipvs_resample_1 <- ipvs_resample[[i]]
  for (j in 1:200) {
    temp <- ipvs_resample_1[[j]]$PVs
    count_pvs <- c(count_pvs, temp)
  }
  ipvs_resample_freq <- table(count_pvs) %>% as.data.frame() %>% arrange(Freq) 
  library(viridis)
  plot_resample[[i]] <- ggplot(data = ipvs_resample_freq, 
                               aes(x = reorder(count_pvs, -Freq), y = Freq, fill = as.factor(count_pvs))) + 
    geom_bar(stat="identity", show.legend = F) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme_bw() +
    scale_color_viridis(discrete = TRUE, option = "magma") +
    geom_hline(yintercept = 50,linetype=2, color = "grey") +
    geom_hline(yintercept = 100,linetype=2) +
    labs(x = "Metabolites Name", y = "Frequency for selected metabolites (Bootstrapping)",
         title = names(ipvs_resample)[i])
}

library(patchwork)
pdf(file = "./Results_check/SFigure_resample_iPVs.pdf", width = 20, height = 6)
plot_resample[[1]] | plot_resample[[2]] | plot_resample[[3]]
dev.off()

### if using anther dataset (imputed_minimum) to do bootstrapping
ipvs_resample_2 <- lapply(imputated_states_minimum, function(x){
  rand_ipvs(x, resample_time = 200, resample_size_factor = 0.75, replacement = T)
})

plot_resample_2 <- list()
for (i in 1:3) {
  count_pvs <- NULL
  # i = 1
  ipvs_resample_1 <- ipvs_resample_2[[i]]
  for (j in 1:200) {
    temp <- ipvs_resample_1[[j]]$PVs
    count_pvs <- c(count_pvs, temp)
  }
  ipvs_resample_freq <- table(count_pvs) %>% as.data.frame() %>% arrange(Freq) 
  library(viridis)
  plot_resample_2[[i]] <- ggplot(data = ipvs_resample_freq, 
                               aes(x = reorder(count_pvs, -Freq), y = Freq, fill = as.factor(count_pvs))) + 
    geom_bar(stat="identity", show.legend = F) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme_bw() +
    scale_color_viridis(discrete = TRUE, option = "magma") +
    geom_hline(yintercept = 50,linetype=2, color = "grey") +
    geom_hline(yintercept = 100,linetype=2) +
    labs(x = "Metabolites Name", y = "Frequency for selected metabolites (Bootstrapping)",
         title = names(ipvs_resample)[i])
}

library(patchwork)
pdf(file = "./Results_check/SFigure_resample_iPVs_minimum.pdf", width = 20, height = 6)
plot_resample_2[[1]] | plot_resample_2[[2]] | plot_resample_2[[3]]
dev.off()

######## if include high loading metabolites based on PCA analysis (first 10 PCs) #####
######## Calculated PCA analysis ##########
library(purrr)

pc_loading_10PC <- list()
for (i in 1:3) {
  # i = 1
  data <- imputated_states_transformation[[i]] %>% as.data.frame()
  pca <- prcomp(data)
  pc_loadings <- data.frame(metabolites = rownames(pca$rotation),
                            pca$rotation[,1:10])
  pc_loadings <- pc_loadings[,-1]
  loading_10PC <- map(pc_loadings, ~ rownames(pc_loadings[order(.x, decreasing = T),])[1])
  pc_loading_10PC[[i]] <- loading_10PC
}

vars_selected_PCA <- unlist(pc_loading_10PC)
vars_selected_ipvs <- unlist(imputated_ipvs_vars)
vars_selected <- c(vars_selected_PCA, vars_selected_ipvs)

### Combined cluster information with baseline information
### 1. Fasting #####
vars_selected_fasting <- unique(vars_selected[grep("_1", vars_selected)])
fasting_raw_selected <- imputated_states_transformation[[1]] %>% as.data.frame() %>% 
  select(any_of(vars_selected_fasting))

library(corrplot)
M1 <- cor(fasting_raw_selected, method = "spearman")
pdf(file = "./Results_check/SFigure_p_cor_fasting.pdf", width = 15, height = 15)
p_cor_fasting <- corrplot(M1, col = COL2('PRGn'))
dev.off()

### 2. Postprandial ####
vars_selected_postprandial <- unique(vars_selected[grep("_3", vars_selected)])
postprandial_raw_selected <- imputated_states_transformation[[2]] %>% as.data.frame() %>% 
  select(any_of(vars_selected_postprandial))
library(corrplot)
M2 <- cor(postprandial_raw_selected, method = "spearman")
pdf(file = "./Results_check/SFig_p_cor_postprandial.pdf", width = 15, height = 15)
p_cor_postprandial <- corrplot(M2, col = COL2('PRGn'))
dev.off()

### 3. Delta ####
vars_selected_delta <- unique(vars_selected[grep("_rp", vars_selected)])
delta_raw_selected <- imputated_states_transformation[[3]] %>% as.data.frame() %>% 
  select(any_of(vars_selected_delta))
library(corrplot)
M3 <- cor(delta_raw_selected, method = "spearman")
pdf(file = "./Results_check/SFig_p_cor_delta.pdf", width = 15, height = 15)
p_cor_delta <- corrplot(M3, col = COL2('PRGn'))
dev.off()


#### If using minimum_dataset -> compare PCA and ipvs
vars_selected_ipvs_minimum <- unlist(imputated_ipvs_vars_minimum)
vars_selected_2 <- c(vars_selected_PCA, vars_selected_ipvs_minimum)

### 1. Fasting #####
vars_selected_fasting_2 <- unique(vars_selected_2[grep("_1", vars_selected_2)])
fasting_raw_selected_2 <- imputated_states_transformation_minimum[[1]] %>% as.data.frame() %>% 
  select(any_of(vars_selected_fasting_2))

library(corrplot)
M1 <- cor(fasting_raw_selected_2, method = "spearman")
pdf(file = "./Results_check/SFigure_p_cor_fasting_min.pdf", width = 15, height = 15)
p_cor_fasting_min <- corrplot(M1, col = COL2('PRGn'))
dev.off()

### 2. Postprandial ####
vars_selected_postprandial_2 <- unique(vars_selected_2[grep("_3", vars_selected_2)])
postprandial_raw_selected_2 <- imputated_states_transformation_minimum[[2]] %>% as.data.frame() %>% 
  select(any_of(vars_selected_postprandial_2))
library(corrplot)
M2 <- cor(postprandial_raw_selected_2, method = "spearman")
pdf(file = "./Results_check/SFig_p_cor_postprandial_min.pdf", width = 15, height = 15)
p_cor_postprandial_min <- corrplot(M2, col = COL2('PRGn'))
dev.off()

### 3. Delta ####
vars_selected_delta_2 <- unique(vars_selected_2[grep("_rp", vars_selected_2)])
delta_raw_selected_2 <- imputated_states_transformation_minimum[[3]] %>% as.data.frame() %>% 
  select(any_of(vars_selected_delta_2))
library(corrplot)
M3 <- cor(delta_raw_selected_2, method = "spearman")
pdf(file = "./Results_check/SFig_p_cor_delta_min.pdf", width = 15, height = 15)
p_cor_delta_min <- corrplot(M3, col = COL2('PRGn'))
dev.off()

####### if we test the fasting working_imputated_df ########
test1 <- working_imputated_df[[1]] # fasting_status_after iPVs

## umap function (don't plan to use this part: umap plot)
library(umap)
test1_umap <- umap(test1)
head(test1_umap$layout, 3)
test1_umap_df <- test1_umap$layout %>% 
  as.data.frame() %>% 
  dplyr::rename(UMAP1 = "V1",
                UMAP2 = "V2") %>% 
  mutate(ID = row_number())

plot1 <- test1_umap_df %>% ggplot(aes(x = UMAP1, y = UMAP2)) + geom_point() + labs(x = "UMAP1",
                                                                                   y = "UMAP2",
                                                                                   subtitle = "UMAP plot") + 
  theme_classic() +
  ylim(-6,6)
ggsave("./Results_check/fasting_umap_df.pdf", height = 20, width = 20)

cormatrix <- stats::cor(test1)
corrplot::corrplot(cormatrix, addCoef.col = "white") ## check the correlation after ipvs


#### Consensuscluster analysis to determine the Optimal cluster number #########
###detach("package:ConsensusClusterPlus", unload = T)
source("ConsensusCluster_function.R")

res_km_fasting <- ConsensusClusterPlus(t(working_imputated_df[[1]]),
                                       maxK = 8,
                                       reps = 100,
                                       pItem = 0.8, ## percent of sample to consider for clustering
                                       pFeature = 1, ## percent of features to consider
                                       plot = "png",
                                       clusterAlg = "km",
                                       title = "km_fasting_status_2024check",
                                       distance = "euclidean",
                                       seed = 2023,
                                       verbose = T)    ## 4 clusters?
calcICL(res_km_fasting, title = "km_fasting_status_2024_check_value", plot = "pdf")

## postprandial ####
Sys.setenv(R_QTRANS_LIMIT = 500000)  # Setting to a higher value
res_km_post <- ConsensusClusterPlus(t(working_imputated_df[[2]]),
                                    maxK = 8,
                                    reps = 100,
                                    pItem = 0.8, ## percent of sample to consider for clustering
                                    pFeature = 1, ## percent of features to consider
                                    plot = "png",
                                    clusterAlg = "km",
                                    title = "km_postprandial_status_2024check",
                                    distance = "euclidean",
                                    seed = 2023,
                                    verbose = T) 
calcICL(res_km_post, title = "km_postprandial_status_2024_check_value", plot = "pdf")

## delta #######
Sys.setenv(R_QTRANS_LIMIT = 500000)  # Setting to a higher value
res_km_delta <- ConsensusClusterPlus(t(working_imputated_df[[3]]),
                                     maxK = 8,
                                     reps = 100,
                                     pItem = 0.8, ## percent of sample to consider for clustering
                                     pFeature = 1, ## percent of features to consider
                                     plot = "png",
                                     clusterAlg = "km",
                                     title = "km_delta_status_2024check",
                                     distance = "euclidean",
                                     seed = 2023,
                                     verbose = T) 
calcICL(res_km_delta, title = "km_delta_status_2024_check_value", plot = "pdf")

#### scenario 2 #########
res_km_fasting_2 <- ConsensusClusterPlus(t(working_imputated_df_minimum[[1]]),
                                       maxK = 8,
                                       reps = 10,
                                       pItem = 0.8, ## percent of sample to consider for clustering
                                       pFeature = 1, ## percent of features to consider
                                       plot = "png",
                                       clusterAlg = "km",
                                       title = "km_fasting_status_2024check_2",
                                       distance = "euclidean",
                                       seed = 2023,
                                       verbose = T)    ## 4 clusters?
calcICL(res_km_fasting_2, title = "km_fasting_status_2024_check_value_2", plot = "pdf")

## postprandial ####
Sys.setenv(R_QTRANS_LIMIT = 500000)  # Setting to a higher value
res_km_post_2 <- ConsensusClusterPlus(t(working_imputated_df_minimum[[2]]),
                                    maxK = 8,
                                    reps = 10,
                                    pItem = 0.8, ## percent of sample to consider for clustering
                                    pFeature = 1, ## percent of features to consider
                                    plot = "png",
                                    clusterAlg = "km",
                                    title = "km_postprandial_status_2024check_2",
                                    distance = "euclidean",
                                    seed = 2023,
                                    verbose = T) 
calcICL(res_km_post_2, title = "km_postprandial_status_2024_check_value_2", plot = "pdf")

## delta #######
Sys.setenv(R_QTRANS_LIMIT = 500000)  # Setting to a higher value
res_km_delta_2 <- ConsensusClusterPlus(t(working_imputated_df_minimum[[3]]),
                                     maxK = 8,
                                     reps = 10,
                                     pItem = 0.8, ## percent of sample to consider for clustering
                                     pFeature = 1, ## percent of features to consider
                                     plot = "png",
                                     clusterAlg = "km",
                                     title = "km_delta_status_2024check2",
                                     distance = "euclidean",
                                     seed = 2023,
                                     verbose = T) 
calcICL(res_km_delta_2, title = "km_delta_status_2024_check_value_2", plot = "pdf")

###### Stability evaluation ###########
### Stability of cluster analysis
library(fpc)                           
res_stable_km_all <- list()

k <- c(4, 4, 4) ## number of clusters

for (i in 1:3) {
  res_stable_km <- clusterboot(working_imputated_df[[i]], B = 100, bootmethod = c("boot"), 
                               clustermethod = kmeansCBI,
                               seed = 1234,
                               k = k[[i]]) 
  res_stable_km_all[[i]] <- res_stable_km$bootresult %>% t() %>% data.frame() 
}
res_stable_km_all2 <- do.call(cbind,res_stable_km_all)

apply(res_stable_km_all2, 2, mean) ## calculate the mean

library(fpc)                           
res_stable_km_all <- list()

###### if using the minimum imputation #######
k <- c(4, 4, 4) ## number of clusters

for (i in 1:3) {
  res_stable_km <- clusterboot(working_imputated_df_minimum[[i]], B = 100, bootmethod = c("boot"), 
                               clustermethod = kmeansCBI,
                               seed = 1234,
                               k = k[[i]]) 
  res_stable_km_all[[i]] <- res_stable_km$bootresult %>% t() %>% data.frame() 
}
res_stable_km_all_min <- do.call(cbind,res_stable_km_all)

apply(res_stable_km_all_min, 2, mean) ## calculate the mean -> Not as stable as the KNN imputed data

## Using gather function to change wide data to long ###
## Plot the Figure of Jaccard index ####
time_name <- c("fasting_state","postprandial_state", "delta_state")    
for (i in 1:3) {
  temp <- res_stable_km_all[[i]] %>%  dplyr::rename(Cluster1 = X1,
                                                    Cluster2 = X2,
                                                    Cluster3 = X3,
                                                    Cluster4 = X4)
  temp_long <- gather(temp, Cluster, Value, Cluster1:Cluster4)
  
  ## Calculate mean, sd, se, and IC
  temp_long <- temp_long %>% 
    group_by(Cluster) %>% 
    summarise(n = n(),
              mean = mean(Value),
              sd = sd(Value)) %>% 
    mutate(se = sd/sqrt(n)) %>%
    mutate(ic = se * qt((1-0.05)/2 + .5, n-1))
  
  ## Barplot
  Jaccard_KM <- ggplot(temp_long, aes(x = Cluster, y = mean, fill = Cluster)) +
    geom_bar(stat = "identity", alpha = 0.5, width = 0.5) +
    geom_errorbar(aes(x = Cluster, 
                      ymin = mean-ic, 
                      ymax = mean+ic), width = 0.2, color = "black", alpha = 0.6, linewidth = 1) +
    labs(x = NULL, y = paste0("Jaccard Index-", time_name[[i]])) +
    theme_bw() +
    scale_color_brewer(palette = "Dark2") +
    scale_y_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.7, linetype = "dashed", col = "grey") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, size = 20),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 20)
    )
  ggsave(filename = paste0("./Results_check/", time_name[[i]], "Jaccard_index.pdf"), width = 10, height = 8)
}

##### Alluvial plot to see the change of the assignment ######
library(amap)
kmeans_cluster_all <- list()
for (i in 1:3) {
  # i = 1
  new_cluster_kmeans <- c()
  for (j in 3:5) {
    set.seed(123)
    # j = 4
    kmeans_cluster = amap::Kmeans(working_imputated_df[[i]], centers = j, iter.max = 200, nstart = 25)
    new_cluster_kmeans <- cbind(new_cluster_kmeans, kmeans_cluster$cluster)
  }
  kmeans_cluster_all[[i]] <- new_cluster_kmeans
}

## save(kmeans_cluster_all, file = "Data/k_means_clustering.Rdata")

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

library(ggalluvial) # install.packages("ggalluvial")

### Fasting status #######
p1.fast <- ggplot(alluvial.df.fast,
                  aes(x = G,
                      stratum = Cluster,
                      alluvium = ID,
                      fill = Cluster,
                      label = Cluster)) + 
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow()  +
  geom_stratum(alpha = 0.8) +
  geom_text(stat = "stratum", size = 8) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        plot.title =  element_text(size = 30)) +
  theme(legend.position = "none") +
  ggtitle("Alluvial plot of cluster membership across different numbers (NEO cohort-Fasting status)") +
  scale_fill_manual(values = c("#e3281f",
                               "#3aa534",
                               "#ff5885",
                               "#fbc926",
                               "#511e9d",
                               "#3498db")) +
  xlab("Assumed number of clusters") + 
  ylab("Frequency") 
ggsave("./Results_check/sFigure_ggalluvial_different_cluster_fasting_2_2023_12_13.pdf", width = 20, height = 20)

### Postprandial status #######
p2.post <- ggplot(alluvial.df.post,
                  aes(x = G,
                      stratum = Cluster,
                      alluvium = ID,
                      fill = Cluster,
                      label = Cluster)) + 
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow()  +
  geom_stratum(alpha = 0.8) +
  geom_text(stat = "stratum", size = 8) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        plot.title =  element_text(size = 30)) +
  theme(legend.position = "none") +
  ggtitle("Alluvial plot of cluster membership across different numbers (NEO cohort-Postprandial status)") +
  scale_fill_manual(values = c("#e3281f",
                               "#3aa534",
                               "#ff5885",
                               "#fbc926",
                               "#511e9d",
                               "#3498db")) +
  xlab("Assumed number of clusters") + 
  ylab("Frequency") 

ggsave("./Results_check/sFig_ggalluvial_different_cluster_postprandial_2_2023_12_13.pdf", width = 20, height = 20)

### Delta status #######
p3.delta <- ggplot(alluvial.df.delta,
                   aes(x = G,
                       stratum = Cluster,
                       alluvium = ID,
                       fill = Cluster,
                       label = Cluster)) + 
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow()  +
  geom_stratum(alpha = 0.8) +
  geom_text(stat = "stratum", size = 8) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        plot.title =  element_text(size = 30)) +
  theme(legend.position = "none") +
  ggtitle("Alluvial plot of cluster membership across different numbers (NEO cohort-Delta status)") +
  scale_fill_manual(values = c("#e3281f",
                               "#3aa534",
                               "#ff5885",
                               "#fbc926",
                               "#511e9d",
                               "#3498db")) +
  xlab("Assumed number of clusters") + 
  ylab("Frequency") 
ggsave("./Results_check/sFig_ggalluvial_different_cluster_delta_2_2023_12_13.pdf", width = 20, height = 20)

library(patchwork)
p_alluvial <- p1.fast | p2.post | p3.delta +
  plot_annotation(tag_levels = 'A')
pdf(file = "./Results_check/p_alluvial_2_2023_12_18_2.pdf", width = 55, height = 20)
p_alluvial ## to see the change of clusters (assigned 3, 4, 5)
dev.off()
