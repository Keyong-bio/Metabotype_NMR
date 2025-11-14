##### P1-Cluster-NEO-KD 
##### Author: Keyong Deng
##### Department of Clinical Epidemiology, LUMC

library(foreign) # data load
library(dplyr) # %>% pipeline
library(data.table) # Fast Read data, fread(), fwrite()
library(Hmisc) # describe function
library(stringr) # string manuiplate
library(ggplot2) # Visualization
library(e1071) # Data skewness 
library(tibble)
# install iPVs package from github
library(iPVs) #identification of principal variable

df <- read.dta(file = "./NEO_analyse_PP104_Brainshake_mealresponse_Deng_2024-12-20.dta")

# Data Cleaning
# 1. exclude prevalent T2DM, without follow-up info, without fasting, without complete meal
table(df$diab_prev, useNA = "always") # 592 cases
table(df$eind2, useNA = "always") # 94 loss-to-follow-up
table(df$nuchter, useNA = "always") # fasting (30 obs)
table(df$devmeal, useNA = "always") # complete meal = 6211 obs

df1 <- df %>% 
  filter(!diab_prev %in% c("Yes"), # 6079 obs
         !eind2 %in% c("lost-to-follow-up"), # 94 obs  # 5985 obs
         !nuchter %in% c("nee"), # 23 obs without fasting meal ## 5985-5962
         devmeal %in% c("complete meal") # 5962-5552 = 410 obs without complete meal
  ) # 5552 obs
# the number corresponding to the SFigure1

# Make a variable (fu_years) for time at risk variable (do this part later 2024-01)
df1$fu_time <- ifelse(df1$t2d_inc == "Yes", df1$diabetes2_date - df1$visitdd, 
                      df1$einddatum2 - df1$visitdd)

df1$fu_years <- df1$fu_time / 365.25 
summary(df1$fu_years) ## Median 6.7 years [5.9-7.9]

# Medications
df1 <- df1 %>%  
  # generate medication variable, put it as the first variable by using .before
  mutate(medication_HeartFailure = ifelse(medC01AA5 == "yes", 1, 0),
         medication_hypertension = ifelse(medC03AA03 == "yes" | medC03BA04 == "yes" | medC07 == "yes" | medC07A == "yes" |
                                            medC07AG == "yes" | medC07B == "yes" | medC07C == "yes" | medC09C == "yes", 
                                          1, 0),
         medication_lipidlower = ifelse(medC10liplow == "yes" | medC10statin == "yes" | medC10fibrate == "yes", 1, 0),
         .before = 1) %>% 
  mutate(medication_all = ifelse(medication_HeartFailure == 1 | medication_hypertension == 1 
                                 | medication_lipidlower == 1, 1, 0)) 

table(df1$medication_HeartFailure) # 8 Heart Failure
table(df1$medication_hypertension) # 1138 hypertension
table(df1$medication_lipidlower) # 633 lipid-lower medication
table(df1$medication_all) # 1440 all medication
table(df1$diab_prev, useNA = "always") # 5552 NO Prevalent DM
table(df1$t2d_30d, useNA = "always") # no cases
table(df1$t2d_inc, useNA = "always") # 281 obs with incident T2DM

# Select Metabolites info 
Metabolite <- df1 %>% 
  dplyr::select(c(61:518))  ## 229 metabolites per time point

# Exclude those metabolites ratios
w <- str_detect(colnames(Metabolite), "p_")
Met_remove1 <- colnames(Metabolite)[which(w == T)][1:140]

w <- str_detect(colnames(Metabolite), "FA_*")
Met_remove_2 <- colnames(Metabolite)[which(w == T)][15:30]

Met_remove_3 <- c(Met_remove1, Met_remove_2, "ApoBApoA1_1", "ApoBApoA1_3", "TGPG_1", "TGPG_3", "DAGTG_1", "DAGTG_3")
Metabolites_NoRat <- Metabolite %>% 
  dplyr::select(!all_of(Met_remove_3)) # 296/2 = 148 metabolites

# Replace any value <=0 entries with an NA
Metabolites_NoRat <- Metabolites_NoRat %>%
  dplyr::mutate(across(XXLVLDLP_1:Gp_3, ~ifelse(. <= 0, NA, .)))

# Detect outliers
is_outlier <- function(object, nSD = nSD){
  object <- as.numeric(object)
  v_mean <- mean(object, na.rm = T)
  v_sd <- sd(object, na.rm = T)
  r <- which(object > (v_mean + nSD * v_sd) | object < (v_mean - nSD * v_sd))
  return(r)
}
w <- apply(Metabolites_NoRat, 2, function(x) is_outlier(x, nSD = 5)) # list: 296 metabolites, more than 5 SD

# Replace the extreme value with NA
for (i in names(w)) {
  a <- w[[i]]
  if (!is.null(a)) {
    Metabolites_NoRat[a, i] <- NA
  }
}

# Sample missing
fmis_sample = apply(Metabolites_NoRat, 1,  function(x) {sum(is.na(x)) / length(x)})
summary(fmis_sample)

# Filter samples
w_sample = which(fmis_sample >= 0.4) # exclude those individuals with more than 40% missing (232 samples)

Metabolites_NoRat <- Metabolites_NoRat[-w_sample,] # 5320 obs for metabolites_NoRatio

df2 <- df1[-w_sample,] ## exclude those samples (5320 obs remains)

# already included as Rdata
save(df2, file = "./df2_baseline.Rdata")

fmis = apply(Metabolites_NoRat, 2,  function(x) {sum(is.na(x)) / length(x)})
summary(fmis)

fmis_df <- as.data.frame(fmis)
fmis_df <- fmis_df %>% rownames_to_column(var = "Meta_status")
fmis_df$Group <-  str_split(fmis_df$Meta_status, "_", simplify = T)[,2]
fmis_df$Group <- ifelse(fmis_df$Group == "1", "Fasting_state", "Postprandial_state")
fmis_df <- fmis_df[order(fmis_df$fmis),]

fmis_df  <- fmis_df %>% 
  group_by(Group) %>% 
  arrange(desc(fmis)) %>% 
  ungroup()

# plot the missingness of metabolites

p1 <- ggplot(fmis_df,aes(Meta_status,fmis)) + 
  geom_bar(stat = "identity", color = "lightblue", fill = "lightblue") +
  theme_bw() +
  labs(title = "Missing Percentage For Each Metabolites in Two States", x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 3, unit = "pt"), size = 6),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  geom_hline(yintercept  = 0.1, linetype = "dashed",linewidth = 1, col = "darkred") +
  facet_wrap(~ Group, nrow = 2, scales = "free_x")

p1

# show the missing percentage of metabolites
if (!dir.exists("Results")) {dir.create("Results")}
# SFigure2
ggsave(p1, file = "Results/SFigure2.pdf", width = 12, height = 6)

cat("- Removing feature with missingness > 10%\n")
w = which(fmis > 0.1) # 20 metabolites with missing rate > 10%
Removemetabolites <- names(w) 
Removemetabolites_fasting <- Removemetabolites[grep("_1", Removemetabolites)] # in the fasting
Removemetabolites_postprandial <- Removemetabolites[grep("_3", Removemetabolites)] # in the postprandial

# fasting
w_fasting <- grep("_1", colnames(Metabolites_NoRat))
fasting_df1 <- Metabolites_NoRat[,w_fasting] %>% 
  dplyr::select(!any_of(Removemetabolites_fasting)) # 131 Metabolites
fmis_fasting = apply(fasting_df1, 2,  function(x) { sum(is.na(x)) / length(x)}) 
summary(fmis_fasting)

# postprandial
w_postprandial <- grep("_3", colnames(Metabolites_NoRat))
postprandial_df1 <- Metabolites_NoRat[, w_postprandial] %>% 
  dplyr::select(!all_of(Removemetabolites_postprandial)) # 145 Metabolites
fmis_postprandial = apply(postprandial_df1, 2,  function(x) { sum(is.na(x)) / length(x)})
summary(fmis_postprandial)

# Delta
w_delta <- unique(str_split(Removemetabolites , "_", simplify = T)[,1]) 
RemoveMet <- c(paste0(w_delta, "_1"), paste0(w_delta, "_3"))
delta_df1 <- Metabolites_NoRat[,setdiff(colnames(Metabolites_NoRat), RemoveMet)] # 131 Metabolites

# fasting state in the delta file
w_temp <- grep("_1", colnames(delta_df1)) 
fasting <- delta_df1[,w_temp]
colnames(fasting) <- str_split(colnames(fasting) , "_", simplify = T)[,1]

postprandial <- delta_df1[,-w_temp] # postprandial state in the delta file
colnames(postprandial) <- str_split(colnames(postprandial) , "_", simplify = T)[,1]

# Make delta response variable
delta <- postprandial - fasting
colnames(delta) <- paste0(colnames(delta), "_rp")

fmis_delta = apply(delta, 2,  function(x) {sum(is.na(x)) / length(x)})
summary(fmis_delta)

# delta response variable filter 
remove_delta <- names(which(fmis_delta > 0.1))
delta_df1 <- delta %>% 
  dplyr::select(!any_of(remove_delta)) # 130 Metabolites with 5320 observations

#### Perform imputation #####
#### Using different method to do imputation ######
#### adapted from : https://github.com/krumsieklab/keto-beam-ad/blob/main/1_serum_dataprocessing.R

#### 1. KNN method
#source("KNN_Functions.R")
 
# set.seed(123)
# ## 695 missing
# fasting_df1_imp <- imputeKNN_multicore(dat = fasting_df1[,-ncol(fasting_df1)],  mc_cores = 10) 
# dim(fasting_df1_imp)
# p_missing <- unlist(lapply(fasting_df1_imp, function(x) sum(is.na(x))))/nrow(fasting_df1_imp)  
# p_missing
# 
# ## 1149 missing
# set.seed(234)
# post_df1_imp <- imputeKNN_multicore(dat = postprandial_df1[,-ncol(postprandial_df1)],  mc_cores = 10) 
# p_missing <- unlist(lapply(post_df1_imp, function(x) sum(is.na(x))))/nrow(post_df1_imp)  
# p_missing
# 
# ## 934 missing
# set.seed(345)
# delta_df1_imp <- imputeKNN_multicore(dat = delta_df1[,-ncol(delta_df1)],  mc_cores = 10) 
# p_missing <- unlist(lapply(delta_df1_imp, function(x) sum(is.na(x))))/nrow(delta_df1_imp)  
# p_missing
# 
# imputated_states <- list(fasting_status = as.data.frame(fasting_df1_imp),
#                          postprandial_status = as.data.frame(post_df1_imp),
#                          delta_status = as.data.frame(delta_df1_imp))
# save(imputated_states, file = "Data/After_QC_Imputed_KNN_minimum.Rdata")

# Load the imputated data with KNN (Please directly load the Rdata, due to the high computation consumption)
load("./After_QC_Imputed_KNN_minimum.Rdata")

rntransform = function(y, split_ties = TRUE ){
  
  # verify that x is a vector.
  if (is.vector(y) == FALSE) { stop("id_outliers() parameter y is expected to be a vector of data.") }
  # verify some variablity
  if (length(unique(y)) == 1) { stop("trait is monomorphic") }
  if (length(unique(y)) == 2) { stop("trait is binary") }
  # identify and remove the NAs
  w = (!is.na(y))
  x <- y[w]
  # empty vector of Z values
  z = rep(NA, length(w))
  # z-transform
  temp <- (x - mean(x))/sd(x)
  # define z values for all observations including the NAs
  z[w] = temp
  # rank the data
  if (split_ties == TRUE) {
    rnt <- rank(z, ties.method = "random") - 0.5
  } else {
    rnt <- rank(z) - 0.5
  }
  # insure the NAs remain so
  rnt[is.na(z)] <- NA
  
  # inverse
  rnt <- rnt/(max(rnt, na.rm = T) + 0.5)
  # quantile normalize
  rnt <- qnorm(rnt)
  return(rnt)
  
}

imputated_states_transformation <- lapply(imputated_states, function(x) apply(x, 2, rntransform))  ## using raw data for transformation

save(imputated_states_transformation, file = "./After_QC_Imputed_KNN_minimum_Rntrans.Rdata")
