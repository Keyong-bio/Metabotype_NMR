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

######## Load data ##############
setwd("/exports/clinicalepi/Keyong/Project1_Cluster")
rawdata <- read.dta(file = "NEO_analyse_PP104_Brainshake_mealresponse_Deng_2023-06-02.dta",
                    convert.factors = T)
df <- rawdata; rm(rawdata)
df$ID <- 1:6671

## Data Cleaning
## 1. exclude prevalent T2DM, without follow-up info, without fasting, without complete meal
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

# Make a variable (fu_years) for time at risk variable (do this part later 2024-01)
df1$fu_time <- ifelse(df1$t2d_inc=="Yes", df1$diabetes2_date-df1$visitdd, 
                      df1$einddatum2-df1$visitdd)

df1$fu_years <- df1$fu_time /365.25 
summary(df1$fu_years)## Median 6.7 years [5.9-7.9]

table(df1$medC10liplow, useNA = "always") # lipid-lowering medication
table(df1$medC01AA5, useNA = "always")
table(df1$medC07)

## Medications
df1 <- df1 %>%  
  mutate(medication_HeartFailure = ifelse(medC01AA5 == "yes", 1, 0),
         medication_hypertension = ifelse(medC03AA03 == "yes" | medC03BA04 == "yes" | medC07 == "yes" | medC07A == "yes" |
                                            medC07AG == "yes" | medC07B == "yes" | medC07C == "yes" | medC09C == "yes", 
                                          1, 0),
         medication_lipidlower = ifelse(medC10liplow == "yes" | medC10statin == "yes" | medC10fibrate == "yes", 1, 0),
         .before = 1) %>% 
  mutate(medication_all = ifelse(medication_HeartFailure == 1 | medication_hypertension == 1 
                                 | medication_lipidlower == 1, 1, 0)) %>% 
  mutate(Identifer = 1:5552,
         .before = 1) 
table(df1$medication_HeartFailure) # 8 Heart Failure
table(df1$medication_hypertension) # 1138 hypertension
table(df1$medication_lipidlower) # 633 lipid-lower medication
table(df1$medication_all) # 1440 all medication
table(df1$diab_prev, useNA = "always") # 5552 NO Prevalent DM
table(df1$t2d_30d, useNA = "always") # no cases
table(df1$t2d_inc, useNA = "always") # 281 obs with incident T2DM

colnames(df1)

## Select baseline variables
blvars <- c("ID","weegfactor1", "sexe", "leeftijd", "bmim",  "lengtem", "gewichtm", "glucose1", "Insuline_r1","HBA1C", "homa1IR", "homa1B", "choltot1", "trig1", "hdlc1", "fldl1","choltot3", "trig3", "hdlc3","fldl3", "medication_all", "medication_HeartFailure", "medication_hypertension", "medication_lipidlower")
bldf <- df1 %>%
  select(all_of(blvars)) # 24 variables

## Table baseline
catvars <- c("sexe","medication_all", "medication_HeartFailure", "medication_hypertension", "medication_lipidlower")
contivars <- setdiff(blvars, catvars)

## Select Metabolites 
Metabolite <- df1 %>% 
  select(c(62:519))  ## 229 metabolites per time points

## description
Nightgale_Met_describ <- Hmisc::describe(Metabolite)
Nightgale_Met_describ <- do.call(cbind, lapply(Nightgale_Met_describ ,function(x) unlist(x)[1:12]))

# the following gives a txt file in the destination that describes the data
# sink("descriptives.txt")
# describe(Metabolite)
# sink()

## Exclude those metabolites ratios
w <- str_detect(colnames(Metabolite), "p_")
Met_remove1 <- colnames(Metabolite)[which(w == T)][1:140]

w <- str_detect(colnames(Metabolite), "FA_*")
Met_remove_2 <- colnames(Metabolite)[which(w == T)][15:30]

Met_remove_3 <- c(Met_remove1, Met_remove_2, "ApoBApoA1_1", "ApoBApoA1_3", "TGPG_1", "TGPG_3", "DAGTG_1", "DAGTG_3")
Metabolites_NoRat <- Metabolite %>% 
  select(!all_of(Met_remove_3)) # 296/2 = 148 metabolites
write.csv(Met_remove_3, file = "Results/Metabolite_Ratios_list.csv", row.names = F)

## Replace any value <=0 entries with an NA
Metabolites_NoRat <- Metabolites_NoRat %>% mutate_at(vars(XXLVLDLP_1:Gp_3), ~ifelse(.<=0,NA,.))

## Detect outliers
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

## Replace the extreme value with NA
for (i in names(w)) {
  a <- w[[i]]
  if (!is.null(a)){
    Metabolites_NoRat[a, i] <- NA
  }
}

## Sample missing
fmis_sample = apply(Metabolites_NoRat, 1,  function(x) { sum(is.na(x)) / length(x)})
# head(fmis_sample)
summary(fmis_sample)

## Filter samples
w_sample = which(fmis_sample > 0.4) ## exclude those individuals with more than 40% missing (232 samples)
Metabolites_NoRat <- Metabolites_NoRat[-w_sample,] # 5320 obs for metabolites_NoRatio
df1 <- df1[-w_sample,] ## exclude those samples (5320 obs remains)

fmis = apply(Metabolites_NoRat, 2,  function(x) {sum(is.na(x)) / length(x)})
summary(fmis)

# fmis = fmis[order(fmis)]
fmis_df <- as.data.frame(fmis)
fmis_df <- fmis_df %>% rownames_to_column(var = "Meta_status")
fmis_df$Group <-  str_split(fmis_df$Meta_status, "_", simplify = T)[,2]
fmis_df$Group <- ifelse(fmis_df$Group == "1", "Fasting_state", "Postprandial_state")
fmis_df <- fmis_df[order(fmis_df$fmis),]

if(!dir.exists("Results_check")) {dir.create("Results_check")}

fmis_df  <- fmis_df %>% 
  group_by(Group) %>% 
  arrange(desc(fmis)) %>% 
  ungroup()

p1 <- ggplot(fmis_df,aes(Meta_status,fmis)) + 
  geom_bar(stat="identity", color = "lightblue", fill = "lightblue") +
  theme_bw() +
  labs(title = "Missing Percentage For Each Metabolites in Two States", x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 3, unit = "pt"), size = 6),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  geom_hline(yintercept  = 0.1, linetype = "dashed",linewidth = 1, col = "darkred") +
  facet_wrap(~ Group, nrow = 2, scales = "free_x")

## Show the missing percentage of metabolites
ggsave(p1, file = "Results_check//SFigure_Missing_rate_Metabolites_NoRat.pdf", width = 12, height = 6)

cat("- Removing feature with missingness > 10%\n")
w = which(fmis > 0.1) # 20 metabolites with missing rate > 10%
Removemetabolites <- names(w) 
Removemetabolites_fasting <- Removemetabolites[grep("_1", Removemetabolites)] # in the fasting
Removemetabolites_postprandial <- Removemetabolites[grep("_3", Removemetabolites)] # in the postprandial

## fasting
w_fasting <- grep("_1", colnames(Metabolites_NoRat))
fasting_df1 <- Metabolites_NoRat[,w_fasting] %>% 
  select(!any_of(Removemetabolites_fasting)) # 131 Metabolites
fmis_fasting = apply(fasting_df1, 2,  function(x) { sum(is.na(x)) / length(x)}) 
summary(fmis_fasting)

## postprandial
w_postprandial <- grep("_3", colnames(Metabolites_NoRat))
postprandial_df1 <- Metabolites_NoRat[, w_postprandial] %>% 
  select(!all_of(Removemetabolites_postprandial)) # 145 Metabolites
fmis_postprandial = apply(postprandial_df1, 2,  function(x) { sum(is.na(x)) / length(x)})
summary(fmis_postprandial)

## Delta
w_delta <- unique(str_split(Removemetabolites , "_", simplify = T)[,1]) 
RemoveMet <- c(paste0(w_delta, "_1"), paste0(w_delta, "_3"))
Delta_df1 <- Metabolites_NoRat[,setdiff(colnames(Metabolites_NoRat), RemoveMet)] # 262 metabolites (131 per state)

w_temp <- grep("_1", colnames(Delta_df1)) ## fasting state in the delta file
fasting <- Delta_df1[,w_temp]
colnames(fasting) <- str_split(colnames(fasting) , "_", simplify = T)[,1]

postprandial <- Delta_df1[,-w_temp] ## postprandial state in the delta file
colnames(postprandial) <- str_split(colnames(postprandial) , "_", simplify = T)[,1]

## Make delta response variable
delta <- postprandial - fasting
colnames(delta) <- paste0(colnames(delta), "_rp")

fmis_delta = apply(delta, 2,  function(x) {sum(is.na(x)) / length(x)})
summary(fmis_delta)

## delta repsonse variable filteration
remove_delta <- names(which(fmis_delta > 0.1))
delta_df1 <- delta %>% 
  select(!any_of(remove_delta)) # 130 Metabolites with 5320 observations

fasting_df1$identifier <- 1:5320
fasting_df1_completes <- fasting_df1[complete.cases(fasting_df1),] 
dim(fasting_df1_completes)# 4625 complete cases

postprandial_df1$identifier <- 1:5320
postprandial_df1_completes <- postprandial_df1[complete.cases(postprandial_df1),] 
dim(postprandial_df1_completes) # 4171 complete cases

delta_df1$identifier <- 1:5320
delta_df1_completes <- delta_df1[complete.cases(delta_df1),] 
dim(delta_df1_completes)# 4386 complete cases

#### Perform imputation #####
#### Using different method to do imputation ######
#### adapted from : https://github.com/krumsieklab/keto-beam-ad/blob/main/1_serum_dataprocessing.R
#### 1. KNN method
source("06_KNN_Functions.R")

set.seed(123)
## 695 missing
fasting_df1_imp <- imputeKNN_multicore(dat = fasting_df1[,-ncol(fasting_df1)],  mc_cores = 10) 
dim(fasting_df1_imp)
p_missing <- unlist(lapply(fasting_df1_imp, function(x) sum(is.na(x))))/nrow(fasting_df1_imp)  
p_missing

## 1149 missing
post_df1_imp <- imputeKNN_multicore(dat = postprandial_df1[,-ncol(postprandial_df1)],  mc_cores = 10) 
p_missing <- unlist(lapply(post_df1_imp, function(x) sum(is.na(x))))/nrow(post_df1_imp)  
p_missing

## 934 missing
delta_df1_imp <- imputeKNN_multicore(dat = delta_df1[,-ncol(delta_df1)],  mc_cores = 10) 
p_missing <- unlist(lapply(delta_df1_imp, function(x) sum(is.na(x))))/nrow(delta_df1_imp)  
p_missing

imputated_states <- list(fasting_status = as.data.frame(fasting_df1_imp),
                         postprandial_status = as.data.frame(post_df1_imp),
                         delta_status = as.data.frame(delta_df1_imp))

save(imputated_states, file = "Data/After_QC_Imputed_KNN_minimum.Rdata")

## Load the imputated data with KNN
load(file = "Data/After_QC_Imputed_KNN_minimum.Rdata")

### Scenarios using imputation (Previous) #####
IMPUTE_State <- function(data_1, data_2, data_3) {
  ## 1. use median value to impute missingness (fasting)
  data_median_imputation <- data_1
  for (i in 1:ncol(data_median_imputation)) {
    x <- data_median_imputation[,i]
    m <- median(x, na.rm = T)
    x[is.na(x)] <- m
    data_median_imputation[, i] <- x
  }

  ## 2. use minimum value to impute missingness
  data_minimum_imputation <- data_2
  for (i in 1:ncol(data_minimum_imputation)) {
    x <- data_minimum_imputation[,i]
    m <- min(x, na.rm = T)
    x[is.na(x)] <- m
    data_minimum_imputation[, i] <- x }

  ## 3. use half of minimum value to impute
  data_half_minimum_imputation <- data_3
  for (i in 1:ncol(data_half_minimum_imputation)) {
    data_half_minimum_imputation[, i] <- tidyr::replace_na(data_half_minimum_imputation[,i], min(data_half_minimum_imputation[, i], na.rm = T)/2)

  }

  imputed_df <- list(data_median_imputation, data_minimum_imputation, data_half_minimum_imputation)
  return(imputed_df)
}

# Do imputation (median, minimum, half_minimum)
fasting_df1_impute <- IMPUTE_State(data_1 = fasting_df1[,-ncol(fasting_df1)],
                                   data_2 = fasting_df1[,-ncol(fasting_df1)],
                                   data_3 = fasting_df1[,-ncol(fasting_df1)])
postprandial_df1_impute <- IMPUTE_State(data_1 = postprandial_df1[,-ncol(postprandial_df1)],
                                        data_2 = postprandial_df1[,-ncol(postprandial_df1)],
                                        data_3 = postprandial_df1[,-ncol(postprandial_df1)])
delta_df1_impute <- IMPUTE_State(data_1 = delta_df1[,-ncol(delta)],
                                 data_2 = delta_df1[,-ncol(delta)],
                                 data_3 = delta_df1[,-ncol(delta)])

imputated_states_minimum <- list(fasting_minimum = fasting_df1_impute[[2]],
                         postprandial_minimum = postprandial_df1_impute[[2]],
                         delta_minimum = delta_df1_impute[[2]])

## Transformation
rntransform = function( y, split_ties = TRUE ){
  ## verify that x is a vector.
  if(is.vector(y) == FALSE){ stop("id_outliers() parameter y is expected to be a vector of data.") }
  ## verify some variablity
  if (length(unique(y)) == 1){ stop("trait is monomorphic") }
  if (length(unique(y)) == 2){ stop("trait is binary") }
  
  ## identify and remove the NAs
  w = (!is.na(y))
  x <- y[w]
  
  ## empty vector of Z values
  z = rep(NA, length(w))
  ## z-transform
  temp <- (x - mean(x))/sd(x)
  ## define z values for all observations including the NAs
  z[w] = temp
  ## rank the data
  if(split_ties == TRUE){
    rnt <- rank(z, ties.method = "random") - 0.5
  } else {
    rnt <- rank(z) - 0.5
  }
  ## insure the NAs remain so
  rnt[is.na(z)] <- NA
  
  ## inverse
  rnt <- rnt/(max(rnt, na.rm = T) + 0.5)
  ## quantile normalize
  rnt <- qnorm(rnt)
  return(rnt)
}

imputated_states_transformation <- lapply(imputated_states, function(x) apply(x, 2, rntransform))  ## using raw data for transformation
imputated_states_transformation_minimum <- lapply(imputated_states_minimum, function(x) apply(x, 2, rntransform))  ## using raw data for transformation (for the imputed minimum)


save(imputated_states_transformation, imputated_states_transformation_minimum, file = "Results_check/imputated_files.Rdata")
