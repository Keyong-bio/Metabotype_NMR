##### P1-Cluster-NEO-KD 
##### Author: Keyong Deng
##### Department of Clinical Epidemiology, LUMC

rm(list = ls())

#### 03 Cox regression model & KM for different groups ######
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library(gridExtra))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(tidycmprsk))
suppressMessages(library(ggsurvfit))
suppressMessages(library(openxlsx))
suppressMessages(library(tibble))
suppressMessages(library(ggforestplot))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(survcomp))

# color palette
pal <- c(unname(yarrr::piratepal("pony")))[c(8,5,3,2,1,9)]

#############################
load("./k_means_clustering.Rdata")
load("./df2_baseline.Rdata")

clusters_fasting <- kmeans_cluster_all[[1]] %>% as.data.frame() %>% 
  dplyr::select(V2); table(clusters_fasting)
clusters_postprandial <- kmeans_cluster_all[[2]] %>% as.data.frame() %>% 
  dplyr::select(V2); table(clusters_postprandial)
clusters_delta <- kmeans_cluster_all[[3]] %>% as.data.frame() %>% 
  dplyr::select(V2); table(clusters_delta)

metadata_out <- df2 %>% 
  dplyr::select(c(1:52,"Glc_1", "Glc_3")) %>% 
  mutate(cluster_fasting = clusters_fasting$V2,
         cluster_postprandial = clusters_postprandial$V2,
         cluster_delta = clusters_delta$V2,
         .before = 1)

metadata_cox <- metadata_out

## recode fasting, postprandial, delta clusters
metadata_cox$cluster_fasting <- as.factor(ifelse(metadata_cox$cluster_fasting == 1, 4,
                                      ifelse(metadata_cox$cluster_fasting == 2, 2, 
                                             ifelse(metadata_cox$cluster_fasting == 3, 1,
                                                    ifelse(metadata_cox$cluster_fasting == 4, 3, NA)))))

metadata_cox$cluster_postprandial <- as.factor(metadata_cox$cluster_postprandial)

metadata_cox$cluster_delta <- ifelse(metadata_cox$cluster_delta == 1, 2,
                                    ifelse(metadata_cox$cluster_delta == 2, 3, 
                                           ifelse(metadata_cox$cluster_delta == 3, 1,
                                                  ifelse(metadata_cox$cluster_delta == 4, 4, NA))))

metadata_cox$cluster_delta <- as.factor(metadata_cox$cluster_delta)

table(metadata_cox$t2d_inc, metadata_cox$cluster_fasting)

## Calculate follow up time (Make a variable (fu_years) for time at risk variable)
metadata_cox$fu_time <- ifelse(metadata_cox$t2d_inc=="Yes", 
                               metadata_cox$diabetes2_date - metadata_cox$visitdd, 
                               metadata_cox$einddatum2-metadata_cox$visitdd)
metadata_cox$fu_years <- metadata_cox$fu_time /365.25

round(summary(metadata_cox$fu_years),2) # 6.77 years (median)
round(sum(metadata_cox$fu_years), 0) # 35227 PYs

metadata_cox$t2d_inc <- as.integer(metadata_cox$t2d_inc)
table(metadata_cox$t2d_inc) # 271 cases

# check the frequency of death (51 obs)
# death_df <- metadata_cox %>% dplyr::filter(eind2 == "death")

metadata_cox$status  <- metadata_cox$t2d_inc-1 # 271 cases (T2DM)

# death occurred before T2DM onset
metadata_cox$status2 <- ifelse((metadata_cox$eind2 == "death"| metadata_cox$eind2 == "death, data until GP extraction 2012/2013") & metadata_cox$t2d_inc == 1, 
                               2, metadata_cox$status)
# consider death as competing risk for T2DM 
table(metadata_cox$status2)

metadata_cox$status2 <- factor(metadata_cox$status2, levels = c(0, 1, 2), 
                               labels = c("health","t2dm","death"))

##### Incidence rate calculation #######
##### All population ###
IR_data <-  data.frame()
IR_data[1,1] <- "Total"
IR_data[1,2] <- metadata_cox %>% dplyr::summarise(N = n()) %>% dplyr::pull(N)
IR_data[1,3] <- metadata_cox %>% filter(t2d_inc== 2) %>% dplyr::summarise(New = n()) %>% dplyr::pull(New)
IR_data[1,4] <- metadata_cox %>% dplyr::summarise(PY=sum(fu_years)) %>% dplyr::pull(PY)

colnames(IR_data)[2:4] <- c("N", "New", "PY")
IR_data$IR <- IR_data$New/IR_data$PY
IR_data$LCI <- exp(log(IR_data$IR) - 1.96*(sqrt((1-IR_data$IR)/IR_data$New)))
IR_data$UCI <- exp(log(IR_data$IR) + 1.96*(sqrt((1-IR_data$IR)/IR_data$New)))
colnames(IR_data)[1] <- "Group"

# fasting status
model <- metadata_cox %>% 
  dplyr::mutate(cluster_fasting = as.factor(cluster_fasting)) %>% 
  group_by(cluster_fasting) %>% 
  dplyr::summarise(N = n())
model2<- metadata_cox %>% 
  filter(t2d_inc == 2) %>% 
  group_by(cluster_fasting) %>% 
  dplyr::summarise(New = n())
model3<- metadata_cox %>% 
  group_by(cluster_fasting) %>% 
  dplyr::summarise(PY = sum(fu_years))
N <- model[2]
New <- model2[2]
PY <- model3[2]
IR <- New/PY
LCI <- exp(log(IR)-1.96*(sqrt((1-IR)/New)))
UCI <- exp(log(IR)+1.96*(sqrt((1-IR)/New)))

IR_out <-cbind(N, New, PY, IR, LCI, UCI)
colnames(IR_out)[4:6] <- c("IR", "LCI", "UCI")
IR_out$Group <- model3$cluster_fasting
IR_out <- IR_out %>% relocate(Group, .before =N) # replace the variables
IR_out$Group <- factor(IR_out$Group, levels = c( "1", "2", "3", "4"),
                       labels = c('Cluster1', 'Cluster2', 'Cluster3', 'Cluster4'))

## postprandial
model <- metadata_cox %>% 
  group_by(cluster_postprandial) %>% 
  dplyr::summarize(N = n())
model2<- metadata_cox %>% 
  filter(t2d_inc == 2) %>% 
  group_by(cluster_postprandial) %>% 
  dplyr::summarize(New=n())
model3<- metadata_cox %>% 
  group_by(cluster_postprandial) %>% 
  dplyr::summarize(PY=sum(fu_years))
N<-model[2]
New<-model2[2]
PY<-model3[2]
IR <- New/PY
LCI <-exp(log(IR) - 1.96*(sqrt((1-IR)/New)))
UCI <-exp(log(IR) + 1.96*(sqrt((1-IR)/New)))
IR_data_post_out <-cbind(N, New, PY, IR, LCI, UCI)
IR_data_post_out$Group <- model3$cluster_postprandial
colnames(IR_data_post_out)[4:6] <- c("IR", "LCI", "UCI")
IR_data_post_out <- IR_data_post_out %>% relocate(Group, .before =N)
IR_data_post_out$Group <- factor(IR_data_post_out$Group, levels = c(1, 2, 3, 4),
                                 labels = c('Cluster1', 'Cluster2', 'Cluster3', 'Cluster4'))

## Delta state
model <- metadata_cox %>% group_by(cluster_delta) %>% dplyr::summarize(N=n())
model2<- metadata_cox %>% filter(t2d_inc==2) %>% group_by(cluster_delta) %>% dplyr::summarize(New=n())
model3<- metadata_cox %>% group_by(cluster_delta) %>% dplyr::summarize(PY=sum(fu_years))
N <- model[2]
New <- model2[2]
PY <- model3[2]
IR <- New/PY
LCI <- exp(log(IR) - 1.96*(sqrt((1-IR)/New)))
UCI <- exp(log(IR) + 1.96*(sqrt((1-IR)/New)))

IR_data_delta_out <-cbind(N, New, PY, IR, LCI, UCI)
IR_data_delta_out$Group <- model3$cluster_delta
colnames(IR_data_delta_out)[4:6] <- c("IR", "LCI", "UCI")
IR_data_delta_out <- IR_data_delta_out %>% relocate(Group, .before =N)

IR_data_delta_out$Group <- factor(IR_data_delta_out$Group, 
                                  levels = c("1", "2", "3", "4"),
                                  labels = c('Cluster1', 'Cluster2', 'Cluster3', 'Cluster4'))

IR_data_delta_out

## combine all the states 
df_IR = rbind(IR_data, IR_out, IR_data_post_out, IR_data_delta_out)
df_IR[,"State"] = c("All", rep("Fasting", 4), rep("Postprandial", 4), rep("Delta", 4))
df_IR <- df_IR %>% 
  dplyr::mutate(IR = IR * 1000,
                LCI = LCI * 1000,
                UCI = UCI * 1000)

# export to supplementary Table5
write.xlsx(df_IR, file = "./Results/STable5.xlsx")

# fix the color 
pal_new <- c(unname(yarrr::piratepal("pony")))[c(1,8,5,3,2,9)]

df_IR_plot <- df_IR %>% 
  dplyr::mutate(Group = factor(Group, 
                                  levels = c("Total", "Cluster1", "Cluster2", 
                                             "Cluster3", "Cluster4"))) %>% 
  dplyr::mutate(State = factor(State,
                                  levels = c("All", "Fasting", 
                                             "Postprandial", "Delta"))) %>% 
  ggplot(aes(x = Group, y = IR)) +
  geom_point(aes(color=Group), size = 2) +
  geom_linerange(aes(ymin = LCI, ymax = UCI, color = Group), 
                 stat = "identity",
                linewidth = 1) +
  labs(y = "T2D Cases per 1000 person year at risk") +
  scale_color_manual(values  = c("#8e8e8e", "#1794CEFF", "#6351A0FF","#F5BACFFF", "#FBBB68FF")) +
  scale_y_continuous(limits = c(0, 20)) +
  facet_grid(~State, scales = "free_x", space = "free_x") +
  # guides(color= guide_legend(override.aes = list(shape = 21))) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.grid.major = element_blank(), # no grid
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth=0.5, 
                                    colour = "black", fill=NA, 
                                    linetype=1),
        strip.background.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        # change the legend 
        legend.title = element_blank(),
        legend.box.spacing = unit(1, "pt"),
        legend.text = element_text(size = 10),
        legend.key.width = unit(10, "pt")
        ) 
df_IR_plot
ggsave("Results/Figure2.pdf", width = 12, height = 6)


# # combined all the plots
# plot = ggarrange(IR_plot1, IR_plot2, IR_plot3, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
# 
# pdf(file = "Results/IR_all.pdf", width = 8, height = 8)
# plot1 = ggdraw() +
#   draw_plot(plot, x = 0, y = 0, width = 1, height = 1) +
#   draw_plot(IR_plot0, x = .75, y = -0.015, width = 0.2, height = 0.535) +
#   draw_plot_label(label = c("A", "B", "C", "D"), size = 12,
#                   x = c(0, 0.5, 0, 0.75), y = c(1, 1, 0.5, 0.5))
# plot1
# dev.off()

################ Cox regression analysis preparation ######################
# 1. KM curves

# test_df <- metadata_cox %>% select(46:58) %>% dplyr::filter(status=="death") ## there are 59 obs were dead before get T2DM
# metadata_cox <- metadata_cox[-which(metadata_cox$status=="death"),]

metadata_cox[, 1:3] <- lapply(metadata_cox[, 1:3], 
                              function(x) factor(x, levels = c("1", "2", "3", "4"),
                                                 labels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4")
                                                 ))

# competing risks model
# 1. cause-specific hazard of a given event
# 2. Sub-distribution hazard of a given event
# cumulative incidence using 1-KM estimates is always >= cumulative incidence using competing risks methods

# see the mean follow-up time for different clusters
round(tapply(metadata_cox$fu_years, list(metadata_cox$cluster_fasting, metadata_cox$status2), mean), digits = 2) 

#### Fitting the cumulative incidence function #####
pdf(file = "Results/Figure3A_1.pdf",width = 8, height = 8)
KM_1 <- tidycmprsk::cuminc(Surv(fu_years, status2) ~ cluster_fasting, data = metadata_cox) %>% 
  ggcuminc(outcome = "t2dm", theme = list(theme_classic(), theme(legend.position = c(0.1, 0.9))), size = 2) + 
  labs(
    title = "Fasting_state",
    x = "Follow-up Years",
    y = "Cumulative incidence of T2DM"
  ) + 
  # ylim(0, 0.15) + 
  add_confidence_interval(alpha = 0.2) +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})",
                theme = list(theme_risktable_boxed()
                             ),
                stats_label = list(n.risk = "Number at Risk")) +
  scale_color_manual(values = pal, breaks = c("Cluster1", "Cluster2", "Cluster3", "Cluster4")) +
  scale_fill_manual(values = pal, breaks = c("Cluster1", "Cluster2", "Cluster3", "Cluster4")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(hjust = 1, size = 12),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position="none"
        ) + 
  annotate(geom="text", x=2, y=0.1, label="log-rank p < 0.001")
print(KM_1)
dev.off()

x1 <- tidycmprsk::cuminc(Surv(fu_years, status2) ~ cluster_fasting, data = metadata_cox)
# log_rank_test
print(x1$cmprsk$Tests)

### postprandial (Figure3B)
pdf(file = "Results/Figure3B_1.pdf",width = 8, height = 8)
KM_2 <- tidycmprsk::cuminc(Surv(fu_years, status2) ~ cluster_postprandial, data = metadata_cox) %>% 
  ggcuminc(outcome = "t2dm", theme = list(theme_classic(), theme(legend.position = c(0.1, 0.9))), size = 2) + 
  labs(
    title = "Postprandial_state",
    x = "Follow-up Years"
  ) + 
  # ylim(0, 0.15) + 
  add_confidence_interval(alpha = 0.2) +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})",
                theme = list(theme_risktable_boxed()
                             # scale_y_discrete(label = rev(c("Low risk", "Low-intermediate risk", 
                             #                    "High-intermediate risk", "High risk")))
                ),
                stats_label = list(n.risk = "Number at Risk")) +
  scale_color_manual(values = pal, breaks = c("Cluster1", "Cluster2", "Cluster3", "Cluster4")) +
  scale_fill_manual(values = pal, breaks = c("Cluster1", "Cluster2", "Cluster3", "Cluster4")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(hjust = 1, size = 12),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12),
        # plot.margin = margin(0,0,0,0),
        legend.text = element_text(size = 12),
        legend.position="none"
  ) + annotate(geom="text", x=2, y=0.1, label="log-rank p < 0.001")
print(KM_2)
dev.off()

x2 <- tidycmprsk::cuminc(Surv(fu_years, status2) ~ cluster_postprandial, data = metadata_cox)
print(x2$cmprsk$Tests)

# delta, Figure3C
pdf(file = "Results/Figure3C_1.pdf",width = 8, height = 8)
KM_3 <- tidycmprsk::cuminc(Surv(fu_years, status2) ~ cluster_delta, data = metadata_cox) %>% 
  ggcuminc(outcome = "t2dm", theme = list(theme_classic(), theme(legend.position = c(0.8, 0.9))), size = 2) + 
  labs(
    title = "Delta_state",
    x = "Follow-up Years"
    #y = "Cumulative incidence of T2DM"
  ) + 
  # ylim(0, 0.15) + 
  add_confidence_interval(alpha = 0.2) +
  # add_confidence_interval(type = "lines") +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})",
                theme = list(theme_risktable_boxed()
                             # scale_y_discrete(label = rev(c("Low risk", "Low-intermediate risk", 
                             #                    "High-intermediate risk", "High risk")))
                ),
                stats_label = list(n.risk = "Number at Risk")) +
  scale_color_manual(values = pal, breaks = c("Cluster1", "Cluster2", "Cluster3", "Cluster4")) +
  scale_fill_manual(values = pal, breaks = c("Cluster1", "Cluster2", "Cluster3", "Cluster4")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(hjust = 1, size = 12),
        axis.title.x =  element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12),
        # plot.margin = margin(0,0,0,0),
        legend.text = element_text(size = 12)
  ) + annotate(geom="text", x=2, y=0.1, label="log-rank p = 0.24")
  
print(KM_3)
dev.off()

x3 <- tidycmprsk::cuminc(Surv(fu_years, status2) ~ cluster_delta, data = metadata_cox)
print(x3$cmprsk$Tests) # log-rank test

### log-rank test for the each prandial state
survdiff(Surv(fu_years, status) ~ cluster_fasting, data = metadata_cox)
survdiff(Surv(fu_years, status) ~ cluster_postprandial, data = metadata_cox)
survdiff(Surv(fu_years, status) ~ cluster_delta, data = metadata_cox) # p value = 0.2

#### combined with the km curves together with risk tables
library(ggsurvfit)
library(patchwork)

pdf(file = "Results/Figure3_revision.pdf", width = 20, height = 8)
patchwork::wrap_plots(list(ggsurvfit::ggsurvfit_build(KM_1),
                         ggsurvfit::ggsurvfit_build(KM_2),
                         ggsurvfit::ggsurvfit_build(KM_3)), ncol = 3) +
  patchwork::plot_annotation(tag_levels = list(c('A', "", 'B', "", 'C', ""))) &
  theme(plot.tag = element_text(size = 14, face = "bold")) # change the font size to 12
dev.off()

######### Cox regression analysis fitted model #############
### without considering death as competing risk
### general cox function 
var_clusters <- c("cluster_fasting", "cluster_postprandial", "cluster_delta")

cox_cluster <- function(data){
  
  cox_res_total <- c()
  for (i in 1:3) {
    M0 <- as.formula(paste("Surv(fu_years, status)~", var_clusters[[i]]))
    M1 <- as.formula(paste("Surv(fu_years, status)~", var_clusters[[i]], "+sexe+leeftijd+bmim"))
    M2 <- as.formula(paste("Surv(fu_years, status)~", var_clusters[[i]], "+sexe+leeftijd+medication_hypertension+medication_lipidlower+bmim+Insuline_r1"))
    
    mv_fit0 <- survival::coxph(M0, data = data)
    
    ## coefficient (beta value)
    fit0_res_beta <- summary(mv_fit0)$coefficients[,1] 
    ## standard error (se value)
    fit0_res_se <- summary(mv_fit0)$coefficients[,3] 
    ## HR
    fit0_res_HR <- summary(mv_fit0)$conf.int[,1] 
    ## HR LCI
    fit0_res_LCI <- summary(mv_fit0)$conf.int[,3] 
    ## HR UCI
    fit0_res_UCI <- summary(mv_fit0)$conf.int[,4] 
    ### Check the cox assumption
    cz <- cox.zph(mv_fit0)$table["GLOBAL",][["p"]]
    
    cox_res <- cbind(fit0_res_beta, fit0_res_se, fit0_res_HR, fit0_res_LCI, fit0_res_UCI, cz) %>% as.data.frame() %>% 
      rownames_to_column() %>% 
      mutate(Model = "Model1",
             .before = 1) 
    
    cox_res <- setNames(cox_res[,1:8], c("Model","Variables", "Beta", "SE", "HR", "HR_LCI", "HR_UCI", "Global_test_assumption"))
    
    mv_fit1 <- survival::coxph(M1, data = data)
    
    fit1_res_beta <- summary(mv_fit1)$coefficients[,1] 
    fit1_res_se <- summary(mv_fit1)$coefficients[,3] 
    fit1_res_HR <- summary(mv_fit1)$conf.int[,1] 
    fit1_res_LCI <- summary(mv_fit1)$conf.int[,3] 
    fit1_res_UCI <- summary(mv_fit1)$conf.int[,4] 
    
    cz <- cox.zph(mv_fit1)$table["GLOBAL",][["p"]]
    
    cox_res1 <- cbind(fit1_res_beta, fit1_res_se, fit1_res_HR, fit1_res_LCI, fit1_res_UCI, cz) %>% 
      as.data.frame() %>% 
      rownames_to_column() %>% 
      dplyr::mutate(Model = "Model2",
             .before = 1)
    
    cox_res1 <- setNames(cox_res1[,1:8], 
                         c("Model","Variables", "Beta", "SE", "HR", "HR_LCI", "HR_UCI", "Global_test_assumption"))
    
    mv_fit2 <- survival::coxph(M2, data = data)
    
    fit2_res_beta <- summary(mv_fit2)$coefficients[,1] 
    fit2_res_se <- summary(mv_fit2)$coefficients[,3] 
    fit2_res_HR <- summary(mv_fit2)$conf.int[,1] 
    fit2_res_LCI <- summary(mv_fit2)$conf.int[,3] 
    fit2_res_UCI <- summary(mv_fit2)$conf.int[,4] 
    
    cz <- cox.zph(mv_fit2)$table["GLOBAL",][["p"]]
    
    cox_res2 <- cbind(fit2_res_beta, fit2_res_se, fit2_res_HR, fit2_res_LCI, fit2_res_UCI, cz) %>% as.data.frame() %>% 
      rownames_to_column() %>% 
      mutate(Model = "Model3",
             .before = 1)
    
    cox_res2 <- setNames(cox_res2[,1:8], c("Model","Variables", "Beta", "SE", "HR", "HR_LCI", "HR_UCI", "Global_test_assumption"))
    
    cox_res_total <- rbind(cox_res_total, rbind(cox_res, cox_res1, cox_res2))
  }
  return(cox_res_total)
}

# without sensitivity analysis
cox_res_no_sens <- cox_cluster(data = metadata_cox) 

# Check the interaction between cluster and sex (not included in the manuscript) ####
mv_fit3 <- survival::coxph(Surv(fu_years, status) ~ cluster_fasting + sexe + leeftijd + medication_hypertension + medication_lipidlower + bmim  + Insuline_r1 + cluster_delta:sexe, data = metadata_cox)
summary(mv_fit3)

mv_fit3 <- survival::coxph(Surv(fu_years, status) ~ cluster_postprandial + sexe +  leeftijd + medication_hypertension + medication_lipidlower + bmim  + Insuline_r1 + cluster_delta:sexe, data = metadata_cox)
summary(mv_fit3)

mv_fit3 <- survival::coxph(Surv(fu_years, status) ~ cluster_delta + sexe +  leeftijd + medication_hypertension + medication_lipidlower + bmim  + Insuline_r1 + cluster_delta:sexe, data = metadata_cox)
summary(mv_fit3)

######## Sensitivity analysis - Cox regression analysis #########
######### 1. Exclude those identified in the first year (Landmark analysis) ############
library(lubridate)
sum(is.na(metadata_cox$visitdd))
# year(metadata_cox$visitdd2) <- year(metadata_cox$visitdd) + 1
# sum(is.na(metadata_cox$visitdd2))
# metadata_cox[is.na(metadata_cox$visitdd2),]

# add_year function
add_year <- function(date) {
  new_date <- date + years(1)
  
  # For NA cases, check if original was Feb 29
  na_idx <- is.na(new_date)
  if(any(na_idx)) {
    orig_dates <- date[na_idx]
    
    # Check if they were Feb 29
    feb29 <- month(orig_dates) == 2 & day(orig_dates) == 29
    # Convert Feb 29 to Feb 28 of next year
    new_date[na_idx][feb29] <- ymd(paste0(
      year(orig_dates[feb29]) + 1, "-02-28"
    ))
  }
  return(new_date)
}
metadata_cox$visitdd2 <- add_year(metadata_cox$visitdd)

metadata_cox$fu_year2 <- ifelse(metadata_cox$t2d_inc=="Yes", 
                               metadata_cox$diabetes2_date - metadata_cox$visitdd2, 
                               metadata_cox$einddatum2 - metadata_cox$visitdd2)/365.25
summary(metadata_cox$fu_year2)

metadata_cox2 <- metadata_cox %>% 
  dplyr::filter(!fu_year2 <= 0) # 5268 observations

cox_cluster_2 <- function(data){
  cox_res_total <- c()
  for (i in 1:3) {
    M0 <- as.formula(paste("Surv(fu_year2, status)~", var_clusters[[i]]))
    M1 <- as.formula(paste("Surv(fu_year2, status)~", var_clusters[[i]], "+sexe+leeftijd+bmim"))
    M2 <- as.formula(paste("Surv(fu_year2, status)~", var_clusters[[i]], "+sexe+leeftijd+medication_hypertension+medication_lipidlower+bmim+Insuline_r1"))
    
    mv_fit0 <- survival::coxph(M0, data = data)
    
    ## coefficient (beta value)
    fit0_res_beta <- summary(mv_fit0)$coefficients[,1] 
    ## standard error (se value)
    fit0_res_se <- summary(mv_fit0)$coefficients[,3] 
    ## HR
    fit0_res_HR <- summary(mv_fit0)$conf.int[,1] 
    ## HR LCI
    fit0_res_LCI <- summary(mv_fit0)$conf.int[,3] 
    ## HR UCI
    fit0_res_UCI <- summary(mv_fit0)$conf.int[,4] 
    
    ### Check the cox assumption
    cz <- cox.zph(mv_fit0)$table["GLOBAL",][["p"]]
    
    cox_res <- cbind(fit0_res_beta, fit0_res_se, fit0_res_HR, fit0_res_LCI, fit0_res_UCI, cz) %>% as.data.frame() %>% 
      tibble::rownames_to_column() %>% 
      dplyr::mutate(Model = "Model1",
             .before = 1) 
    
    cox_res <- setNames(cox_res[,1:8], c("Model","Variables", "Beta", "SE", "HR", "HR_LCI", "HR_UCI", "Global_test_assumption"))
    
    mv_fit1 <- survival::coxph(M1, data = data)
    
    fit1_res_beta <- summary(mv_fit1)$coefficients[,1] 
    fit1_res_se <- summary(mv_fit1)$coefficients[,3] 
    fit1_res_HR <- summary(mv_fit1)$conf.int[,1] 
    fit1_res_LCI <- summary(mv_fit1)$conf.int[,3] 
    fit1_res_UCI <- summary(mv_fit1)$conf.int[,4] 
    
    cz <- cox.zph(mv_fit1)$table["GLOBAL",][["p"]]
    
    cox_res1 <- cbind(fit1_res_beta, fit1_res_se, fit1_res_HR, fit1_res_LCI, fit1_res_UCI, cz) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column() %>% 
      dplyr::mutate(Model = "Model2",
                    .before = 1)
    
    cox_res1 <- setNames(cox_res1[,1:8], 
                         c("Model","Variables", "Beta", "SE", "HR", "HR_LCI", "HR_UCI", "Global_test_assumption"))
    
    mv_fit2 <- survival::coxph(M2, data = data)
    
    fit2_res_beta <- summary(mv_fit2)$coefficients[,1] 
    fit2_res_se <- summary(mv_fit2)$coefficients[,3] 
    fit2_res_HR <- summary(mv_fit2)$conf.int[,1] 
    fit2_res_LCI <- summary(mv_fit2)$conf.int[,3] 
    fit2_res_UCI <- summary(mv_fit2)$conf.int[,4] 
    
    cz <- cox.zph(mv_fit2)$table["GLOBAL",][["p"]]
    
    cox_res2 <- cbind(fit2_res_beta, fit2_res_se, fit2_res_HR, fit2_res_LCI, fit2_res_UCI, cz) %>% as.data.frame() %>% 
      tibble::rownames_to_column() %>% 
      dplyr::mutate(Model = "Model3",
             .before = 1)
    
    cox_res2 <- setNames(cox_res2[,1:8], c("Model","Variables", "Beta", "SE", "HR", "HR_LCI", "HR_UCI", "Global_test_assumption"))
    
    cox_res_total <- rbind(cox_res_total, rbind(cox_res, cox_res1, cox_res2))
  }
  return(cox_res_total)
}

# with sensitivity analysis
cox_res_sens <- cox_cluster_2(data = metadata_cox2) %>% 
  filter(!Variables %in% c("sexevrouw", "leeftijd", "bmim", 
  "medication_hypertension", "medication_lipidlower", "Insuline_r1"))

library(plyr)
cox_res_sens1 <- numcolwise(round,2)(cox_res_sens)
cox_res_sens1 <- cbind(cox_res_sens[,c(1, 2)], cox_res_sens1)
detach("package:plyr", unload=TRUE) 

## export result
# write.xlsx(cox_res_no_sens, file = "Results/cox_res_all.xlsx", rowNames = F)
# write.xlsx(cox_res_sens1, file = "Results/cox_res_all_sens_revised.xlsx", rowNames = F)

######### Fine-Gray Analysis ############
#### competing risk analysis from https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html
#### Calculate the Hazard Ratio (Fine-Gray model)
fg_res_all <- c()

for (i in 1:3) {
  # i = 1
  M0 <- as.formula(paste("Surv(fu_years, status2)~", var_clusters[[i]]))
  M1 <- as.formula(paste("Surv(fu_years, status2)~", var_clusters[[i]], "+sexe+leeftijd+bmim"))
  M2 <- as.formula(paste("Surv(fu_years, status2)~", var_clusters[[i]], "+sexe+leeftijd+medication_hypertension+medication_lipidlower+bmim+Insuline_r1"))
  
  model <- crr(M0, data = metadata_cox) 
  Beta_val <- model$tidy$estimate
  se_val <- model$tidy$std.error
  HR1 <- exp(model$tidy$estimate)
  HR1_LCI <- exp(model$tidy$conf.low)
  HR1_UCI <- exp(model$tidy$conf.high)
  
  cox_fg <- cbind(model$tidy$term, Beta_val, se_val, HR1, HR1_LCI, HR1_UCI) %>% as.data.frame() %>% 
    mutate(Model = "Model1",
           .before = 1) 
  cox_fg <- setNames(cox_fg[,1:7], c("Model","Variables", "Beta_fg", "SE_fg", "HR_fg", "HR_fg_LCI", "HR_fg_UCI"))
  
  model1 <- crr(M1, data = metadata_cox) 
  Beta_val <- model1$tidy$estimate
  se_val <- model1$tidy$std.error
  HR2 <- exp(model1$tidy$estimate)
  HR2_LCI <- exp(model1$tidy$conf.low)
  HR2_UCI <- exp(model1$tidy$conf.high)
  
  cox_fg1 <- cbind(model1$tidy$term, Beta_val,se_val, HR2, HR2_LCI, HR2_UCI) %>% as.data.frame() %>% 
    mutate(Model = "Model2",
           .before = 1) 
  cox_fg1 <- setNames(cox_fg1[,1:7], c("Model","Variables", "Beta_fg", "SE_fg", "HR_fg", "HR_fg_LCI", "HR_fg_UCI"))
  
  model2 <- crr(M2, data = metadata_cox) 
  Beta_val <- model2$tidy$estimate
  se_val <- model2$tidy$std.error
  HR3 <- exp(model2$tidy$estimate)
  HR3_LCI <- exp(model2$tidy$conf.low)
  HR3_UCI <- exp(model2$tidy$conf.high)
  
  cox_fg2 <- cbind(model2$tidy$term, Beta_val, se_val, HR3, HR3_LCI, HR3_UCI) %>% as.data.frame() %>% 
    mutate(Model = "Model3",
           .before = 1) 
  cox_fg2 <- setNames(cox_fg2[,1:7], c("Model","Variables", "Beta_fg", "SE_fg", "HR_fg", "HR_fg_LCI", "HR_fg_UCI"))
  
  fg_res_all <- rbind(fg_res_all, rbind(cox_fg, cox_fg1, cox_fg2))
  
}

fg_res_all[,3:7] <- apply(fg_res_all[,3:7], 2, as.numeric)

write.xlsx(fg_res_all, file = "Results/fg_res_all.xlsx", rowNames = F)


############### Figure-for-Cox-regression ##############
library(ggforestplot)

## Cox-res-without sensitivity ###
cox_res_no_sens <- cox_cluster(data = metadata_cox) %>% 
  filter(!Variables %in% c("sexe1", "leeftijd", "bmim", 
                           "medication_hypertension", "medication_lipidlower", "Insuline_r1"))
  
cox_res_no_sens <- tidyr::separate(data = cox_res_no_sens, 
                                   col = Variables, 
                                   c("States", "Clusters"), 
                                   sep = "(?<=[a-z])(?=[A-Z])") %>% 
  dplyr::filter(!is.na(Clusters))

cox_res_no_sens$States <- factor(cox_res_no_sens$States, 
                                 levels = c("cluster_fasting", "cluster_postprandial", "cluster_delta"), 
                                 labels = c("Fasting_state", "Postprandial_state", "Delta_state"))

cox_res_sens <- tidyr::separate(data = cox_res_sens, 
                                col = Variables, c("States", "Clusters"), 
                                sep = "(?<=[a-z])(?=[A-Z])") %>% 
  dplyr::filter(!is.na(Clusters))

cox_res_sens$States <- factor(cox_res_sens$States, 
                                 levels = c("cluster_fasting", "cluster_postprandial", "cluster_delta"), 
                                 labels = c("Fasting_state", "Postprandial_state", "Delta_state"))

fg_res_all <- read.xlsx("Results/fg_res_all.xlsx") %>% 
  filter(!Variables %in% c("sexevrouw", "leeftijd", "bmim", 
                           "medication_hypertension", "medication_lipidlower", "Insuline_r1"))

fg_res_all <- tidyr::separate(data = fg_res_all, 
                              col = Variables, c("States", "Clusters"), 
                              sep = "(?<=[a-z])(?=[A-Z])") %>% 
  dplyr::filter(!is.na(Clusters))

fg_res_all$States <- factor(fg_res_all$States, 
                              levels = c("cluster_fasting", "cluster_postprandial", "cluster_delta"), 
                              labels = c("Fasting_state", "Postprandial_state", "Delta_state"))

## Forest plot for Cox regression model ####
suppressMessages(library(ggforce))

# Add the reference group
refence <- data.frame(Model = c("Model1", "Model2", "Model3"), 
                      States = rep(c("Fasting_state", "Postprandial_state", "Delta_state"), each = 3),
                      Clusters = rep(c("Cluster1", "Cluster1", "Cluster1"), each = 3),
                      Beta = 0,
                      SE = 0,
                      HR = 1,
                      HR_LCI  = 1,
                      HR_UCI = 1, Global_test_assumption = NA)

cox_res_no_sens_all <- rbind(refence, cox_res_no_sens)
cox_res_no_sens_all <- cox_res_no_sens_all[order(cox_res_no_sens_all$Model, cox_res_no_sens_all$States),]

cox_res_no_sens_all <- cox_res_no_sens_all %>% 
  dplyr::mutate(
    Model = factor(
      Model, levels = c("Model3", "Model2", "Model1")
    ),
    States = factor(
      States, 
      levels = c("Fasting_state", "Postprandial_state", "Delta_state")
    )
  ) %>% 
  dplyr::mutate(
    Clusters = factor(
      Clusters,
      levels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4")
    )
  )

# Figure4
pdf(file = "Results/Figure4.pdf", width = 10, height = 6)
library(ggforce)
p_cox <- ggforestplot::forestplot(
  df = cox_res_no_sens_all,
  name = Clusters,
  estimate = Beta,
  se = SE,
  xlab = "Hazard Ratio with 95%CI",
  colour = Model,
  logodds = T
) +
  ggforce::facet_row(
    facets = ~ States,
    scales = "free_x"
  ) +
  ggplot2::scale_color_brewer(
    palette = "Set1"
  ) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )
p_cox
dev.off()

#### Sensitivity Cox #####
cox_res_sens_all <- rbind(refence, cox_res_sens)
cox_res_sens_all <- cox_res_sens_all %>% 
  dplyr::mutate(
    Model = factor(
      Model, levels = c("Model3", "Model2", "Model1")
    ),
    States = factor(
      States, levels = c("Fasting_state", "Postprandial_state", "Delta_state")
    )
  ) %>% 
  dplyr::mutate(
    Clusters = factor(
      Clusters,
      levels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4")
    )
  )

# SFigure9
pdf(file = "Results/SFigure9.pdf", width = 10, height = 6)
library(ggforce)
p_cox <- ggforestplot::forestplot(
  df = cox_res_sens_all,
  name = Clusters,
  estimate = Beta,
  se = SE,
  xlab = "Hazard Ratio with 95%CI",
  colour = Model,
  logodds = T
) +
  ggforce::facet_row(
    facets = ~ States,
    scales = "free_x",
    # space = "free"
  ) +
  ggplot2::scale_color_brewer(
    palette = "Set1"
  ) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )
p_cox
dev.off()

### Fine-Gray model #####
fg_res_all <- read.xlsx("Results/fg_res_all.xlsx") %>% 
  filter(!Variables %in% c("sexevrouw", "leeftijd", "bmim", 
                           "medication_hypertension", "medication_lipidlower", "Insuline_r1"))

fg_res_all <- tidyr::separate(data = fg_res_all, 
                              col = Variables, c("States", "Clusters"), 
                              sep = "(?<=[a-z])(?=[A-Z])") %>% 
  dplyr::filter(!is.na(Clusters))

fg_res_all$States <- factor(fg_res_all$States, 
                            levels = c("cluster_fasting", "cluster_postprandial", "cluster_delta"), 
                            labels = c("Fasting_state", "Postprandial_state", "Delta_state"))


refence_2 <- data.frame(Model = c("Model1", "Model2", "Model3"), 
                      States = rep(c("Fasting_state", "Postprandial_state", "Delta_state"), each = 3),
                      Clusters = rep(c("Cluster1", "Cluster1", "Cluster1"), each = 3),
                      Beta_fg = 0,
                      SE_fg = 0,
                      HR_fg = 1,
                      HR_fg_LCI  = 1,
                      HR_fg_UCI = 1)

fg_res_all_2 <- rbind(refence_2, fg_res_all)

fg_res_all_2  <- fg_res_all_2  %>% 
  mutate(
    Model = factor(
      Model, levels = c("Model3", "Model2", "Model1")
    ),
    States = factor(
      States, levels = c("Fasting_state", "Postprandial_state", "Delta_state")
    )
  )

# SFigure10
pdf(file = "Results/SFigure10.pdf", width = 10, height = 6)
library(ggforce)
p_cox <- forestplot(
  df = fg_res_all_2,
  name = Clusters,
  estimate = Beta_fg,
  se = SE_fg,
  xlab = "Hazard Ratio with 95%CI",
  colour = Model,
  logodds = T
) +
  ggforce::facet_row(
    facets = ~ States,
    scales = "free_x"
  ) +
  ggplot2::scale_color_brewer(
    palette = "Set1"
  ) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(), # no grid
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )
p_cox
dev.off()

############# Added value of the metabotypes in the T2DM prediction #########
suppressMessages(library(survIDINRI))
suppressMessages(library(timeROC))
suppressMessages(library(timereg))
suppressMessages(library(nricens))
suppressMessages(library(caret))
suppressMessages(library(glmnet))
suppressMessages(library(survcomp))

# Cox regression - prediction model based on internal validation, cross-validation ------
# Define a function to create the model
create_model <- function(base_vars, addition_vars = NULL, interaction_vars = NULL){
  
  # start with base formula
  formula_str <- "~ 0"
  
  # add base variable (always included)
  if (length(base_vars) > 0) {
    formula_str <- paste(formula_str, "+", paste(base_vars, collapse = "+"))
  } 
  
  if (!is.null(interaction_vars) && length(addition_vars) > 0){
    formula_str <- paste(formula_str, "+", paste(addition_vars, collapse = "+"))
  }
  
  # add interaction items 
  if (!is.null(interaction_vars) && length(interaction_vars) > 0){
    for (i in 1:length(interaction_vars)){
      if(length(interaction_vars[[i]]) == 2) {
        formula_str <- paste(formula_str, "+", 
                             paste(interaction_vars[[i]], collapse = ":"))
      }
    }
  }
  
  # create and return the model
  formula_obj <- as.formula(formula_str)
  return(formula_obj)
  # return(model.matrix(formula_obj, data = df))
}

# for loop for analysis cross_validation ------------

# metadata_cox$sexe <- ifelse(metadata_cox$sexe  == "vrouw", 0, 1)
metadata_cox$sexe <- as.factor(metadata_cox$sexe)

set.seed(123)
k <- 10

# folds
# folds <- createFolds(metadata_cox$status, k = k, list = TRUE, returnTrain = TRUE)
# Create a grouping variable
strata_var <- interaction(metadata_cox$sexe, as.factor(metadata_cox$status))
folds <- createFolds(strata_var, k = 10, list = TRUE, returnTrain = TRUE)

# predict years
prediction_time <- c(5, 7)
time_labels <- c("5_year", "7_year")

# Create a list to store the results
result_list <- list()

for (t_index in 1:length(prediction_time)) {
  # 
  t0 <- prediction_time[t_index]
  t_label <- time_labels[t_index]
  
  # create a dataframe to store the results for each fold
  fold_results <- data.frame(
    fold = 1:k,
    c_index_base = numeric(k),
    c_index_extra1 = numeric(k),
    c_index_extra2 = numeric(k),
    c_index_extra3 = numeric(k),
    c_index_extra4 = numeric(k),
    IDI_extra1 = numeric(k),
    IDI_extra2 = numeric(k),
    IDI_extra3 = numeric(k),
    IDI_extra4 = numeric(k),
    NRI_extra1 = numeric(k),
    NRI_extra2 = numeric(k),
    NRI_extra3 = numeric(k),
    NRI_extra4 = numeric(k)
  )
  
  # perform for k fold
  for (i in 1:k) {
    # split data into two parts
    train_df <- metadata_cox[folds[[i]],] 
    test_df <- metadata_cox[-folds[[i]],]
    
    print(table(train_df$cluster_fasting))
    print(table(train_df$cluster_postprandial))
    print(table(train_df$sexe))
    print(table(test_df$cluster_fasting))
    print(table(test_df$cluster_postprandial))
    
    
    # based on different models
    base_model <- coxph(Surv(fu_years, status) ~  sexe + leeftijd + bmim, data = train_df, x = TRUE)
    extra_model1 <- coxph(Surv(fu_years, status) ~  sexe + leeftijd + bmim + cluster_fasting, data = train_df, x = TRUE)
    extra_model2 <- coxph(Surv(fu_years, status) ~  sexe + leeftijd + bmim + cluster_postprandial, data = train_df, x = TRUE)
    extra_model3 <- coxph(Surv(fu_years, status) ~  sexe + leeftijd + bmim + cluster_fasting + cluster_postprandial, data = train_df, x = TRUE)
    extra_model4 <- coxph(Surv(fu_years, status) ~  sexe + leeftijd + bmim + cluster_fasting + cluster_postprandial + cluster_fasting:cluster_postprandial, data = train_df, x = TRUE)
    
    # predict on the test data
    base_risk <- predict(base_model, test_df)
    extra_model1_risk <- predict(extra_model1, test_df)
    extra_model2_risk <- predict(extra_model2, test_df)
    extra_model3_risk <- predict(extra_model3, test_df)
    extra_model4_risk <- predict(extra_model4, test_df)
    
    # c_index
    surv_obj <- Surv(test_df$fu_years, test_df$status)
    fold_results$c_index_base[i] <- survConcordance(surv_obj ~ base_risk)$concordance
    fold_results$c_index_extra1[i] <- survConcordance(surv_obj ~ extra_model1_risk)$concordance
    fold_results$c_index_extra2[i] <- survConcordance(surv_obj ~ extra_model2_risk)$concordance
    fold_results$c_index_extra3[i] <- survConcordance(surv_obj ~ extra_model3_risk)$concordance
    fold_results$c_index_extra4[i] <- survConcordance(surv_obj ~ extra_model4_risk)$concordance
    
    # extract relevant columns for survIDINRI packages
    test_df <- test_df[complete.cases(test_df[, c("fu_years", "status", "sexe", "leeftijd", "bmim", "cluster_fasting")]),]
    surv_data <- test_df[, c("fu_years", "status")]
    
    covs0 <- model.matrix(~ 0 + sexe + leeftijd + bmim, 
                          data = test_df)[, -1]
    covs1 <- model.matrix(~ 0 + sexe + leeftijd + bmim + cluster_fasting, 
                          data = test_df)[,-1]
    covs2 <- model.matrix(~ 0 + sexe + leeftijd + bmim + cluster_postprandial, 
                          data = test_df)[,-1]
    covs3 <- model.matrix(~ 0 + sexe + leeftijd + bmim + cluster_postprandial + cluster_fasting, 
                          data = test_df)[,-1]
    covs4 <- model.matrix(~ 0 + sexe + leeftijd + bmim + cluster_postprandial + cluster_fasting + cluster_fasting:cluster_postprandial, 
                          data = test_df)[,-1]
    
    # base versus extra_model1
    result1 <- try(
      IDI.INF(indata = surv_data,
              covs0 = covs0,
              covs1 = covs1,
              t0 = t0,
              npert = 1000,
              seed1 = 2025),
      silent = TRUE
    )
    
    # base versus extra_model2
    result2 <- try(
      IDI.INF(indata = surv_data,
              covs0 = covs0,
              covs1 = covs2,
              t0 = t0,
              npert = 1000,
              seed1 = 2025),
      silent = TRUE
    )
    
    # base versus extra_model3
    result3 <- try(
      IDI.INF(indata = surv_data,
              covs0 = covs0,
              covs1 = covs3,
              t0 = t0,
              npert = 1000,
              seed1 = 2025),
      silent = TRUE
    )
    
    # base versus extra_model4
    result4 <- try(
      IDI.INF(indata = surv_data,
              covs0 = covs0,
              covs1 = covs4,
              t0 = t0,
              npert = 1000,
              seed1 = 2025),
      silent = TRUE
    )
    
    # Save IDI and NRI value 
    if(TRUE){
      summary1 <- IDI.INF.OUT(result1)
      fold_results$IDI_extra1[i] <- summary1[1,1]
      fold_results$NRI_extra1[i] <- summary1[2,1]
    }
    
    if(TRUE){
      summary2 <- IDI.INF.OUT(result2)
      fold_results$IDI_extra2[i] <- summary2[1,1]
      fold_results$NRI_extra2[i] <- summary2[2,1]
    }
    
    if(TRUE){
      summary3 <- IDI.INF.OUT(result3)
      fold_results$IDI_extra3[i] <- summary3[1,1]
      fold_results$NRI_extra3[i] <- summary3[2,1]
    }
    
    if(TRUE) {
      summary4 <- IDI.INF.OUT(result4)
      fold_results$IDI_extra4[i] <- summary4[1,1]
      fold_results$NRI_extra4[i] <- summary4[2,1]
    }
    # calculate the summary for each columns
    result_list[[t_label]] <- fold_results
  }
}

foldresult_1 <- result_list[["5_year"]]
foldresult_2 <-result_list[["7_year"]]

# calculate the confidence intervals (5 years)
foldresult_1$diff_c_extra1 <- foldresult_1$c_index_extra1 - foldresult_1$c_index_base
foldresult_1$diff_c_extra2 <- foldresult_1$c_index_extra2 - foldresult_1$c_index_base
foldresult_1$diff_c_extra3 <- foldresult_1$c_index_extra3 - foldresult_1$c_index_base
foldresult_1$diff_c_extra4 <- foldresult_1$c_index_extra4 - foldresult_1$c_index_base

# function to calculate CIs
get_mean_ci <- function(x){
  n = length(x)
  mean_diff = mean(x)
  se = sd(x)/sqrt(n)
  ci_lower <- mean_diff - 1.96*se
  ci_upper <- mean_diff + 1.96*se
  mean_CI_res <- sprintf("%.2f(%.2f,%.2f)", mean_diff, ci_lower, ci_upper)
  return(mean_CI_res)
}

# for the c-index calculation
c_extra_base_ci <- get_mean_ci(foldresult_1$c_index_base)
c_extra1_ci <- get_mean_ci(foldresult_1$c_index_extra1)
c_extra2_ci <- get_mean_ci(foldresult_1$c_index_extra2)
c_extra3_ci <- get_mean_ci(foldresult_1$c_index_extra3)
c_extra4_ci <- get_mean_ci(foldresult_1$c_index_extra4)
c_index <- cbind(c_extra_base_ci, c_extra1_ci, c_extra2_ci, c_extra3_ci, c_extra4_ci)
c_index

diff_c_extra1_ci <- get_mean_ci(foldresult_1$diff_c_extra1)
diff_c_extra2_ci <- get_mean_ci(foldresult_1$diff_c_extra2)
diff_c_extra3_ci <- get_mean_ci(foldresult_1$diff_c_extra3)
diff_c_extra4_ci <- get_mean_ci(foldresult_1$diff_c_extra4)
diff_c_ci = cbind(diff_c_extra1_ci, diff_c_extra2_ci, diff_c_extra3_ci, diff_c_extra4_ci)

# for the 5 year, IDI calculation
IDI_extra1_ci <- get_mean_ci(foldresult_1$IDI_extra1)
IDI_extra2_ci <- get_mean_ci(foldresult_1$IDI_extra2)
IDI_extra3_ci <- get_mean_ci(foldresult_1$IDI_extra3)
IDI_extra4_ci <- get_mean_ci(foldresult_1$IDI_extra4)
IDI_ci = cbind(IDI_extra1_ci, IDI_extra2_ci, IDI_extra3_ci, IDI_extra4_ci)

NRI_extra1_ci <- get_mean_ci(foldresult_1$NRI_extra1)
NRI_extra2_ci <- get_mean_ci(foldresult_1$NRI_extra2)
NRI_extra3_ci <- get_mean_ci(foldresult_1$NRI_extra3)
NRI_extra4_ci <- get_mean_ci(foldresult_1$NRI_extra4)
NRI_ci = cbind(NRI_extra1_ci, NRI_extra2_ci, NRI_extra3_ci, NRI_extra4_ci)

# for the 7 year, IDI calculation
# calculate the confidence intervals (7 years)

IDI_extra1_ci_7y <- get_mean_ci(foldresult_2$IDI_extra1)
IDI_extra2_ci_7y <- get_mean_ci(foldresult_2$IDI_extra2)
IDI_extra3_ci_7y <- get_mean_ci(foldresult_2$IDI_extra3)
IDI_extra4_ci_7y <- get_mean_ci(foldresult_2$IDI_extra4)
IDI_ci_7y = cbind(IDI_extra1_ci_7y, IDI_extra2_ci_7y, IDI_extra3_ci_7y, IDI_extra4_ci_7y)

NRI_extra1_ci_7y <- get_mean_ci(foldresult_2$NRI_extra1)
NRI_extra2_ci_7y <- get_mean_ci(foldresult_2$NRI_extra2)
NRI_extra3_ci_7y <- get_mean_ci(foldresult_2$NRI_extra3)
NRI_extra4_ci_7y <- get_mean_ci(foldresult_2$NRI_extra4)
NRI_ci_7y = cbind(NRI_extra1_ci_7y, NRI_extra2_ci_7y, NRI_extra3_ci_7y, NRI_extra4_ci_7y)

# combine all the results
metrics <- c("diff_c", "IDI_5y", "IDI_7y", "NRI_5y", "NRI_7y")
metric_objects <- list(diff_c_ci, IDI_ci, IDI_ci_7y, NRI_ci, NRI_ci_7y)
results_df_pre <- data.frame(Metric = metrics)
for (i in 1:4) {
  results_df_pre[paste0("Extra", i)] <- sapply(metric_objects, function(x) x[1,i])
}
colnames(results_df_pre)[2:5] <- c("Basic+Fasting", "Basic+Postprandial", "Basic+Fasting+Postprandial", 
                                   "Basic+Fasting*Postprandial")
results_df_pre %>% write.xlsx("./Results/Prediction_Model_revision.xlsx")

## Previous methods -------------------------
# without considering competing risk from death
# time = metadata_cox$fu_time
# event = metadata_cox$status
# 
# metadata_cox$sexe <- ifelse(metadata_cox$sexe  == "vrouw", 0, 1)
# metadata_cox$sexe <- as.factor(metadata_cox$sexe)
# 
# # basic prediction model, leeftijd, sexe, bmim, Insuline_r1, HBA1C, glucose1
# z.std = subset(metadata_cox, select = c(fu_years, status, leeftijd, sexe, bmim, HBA1C))
# z.std = z.std[complete.cases(z.std),]
# 
# # new prediction model, basic + fasting metabotypes (Insuline_r1, HBA1C, glucose1)
# z.new = subset(metadata_cox, select = c(fu_years, status, leeftijd, sexe, bmim, HBA1C, cluster_fasting))
# z.new = z.new[complete.cases(z.new),]
# 
# z.new2 = subset(metadata_cox, select = c(fu_years, status, leeftijd, sexe, bmim, HBA1C, cluster_postprandial, cluster_fasting))
# z.new2 = z.new2[complete.cases(z.new2),]
# 
# # basic model (age, sex, and bmi)
# mstd = coxph(Surv(fu_years, status) ~ sexe + leeftijd + bmim,  
#              data = z.new, 
#              x = TRUE)
# mstd_predict = predict(mstd)
# 
# tt <- concordance.index(
#   x = mstd_predict,
#   surv.time = z.new$fu_years,
#   surv.event = z.new$status,
#   method = "noether",
#   na.rm = TRUE
# )
# 
# # tt$c.index
# # tt$se
# # tt$lower
# # concordance(mstd)
# 
# # added fasting
# mnew1 = coxph(Surv(fu_years, status) ~ sexe + leeftijd + bmim + cluster_fasting,  
#              data = z.new, 
#              x = TRUE)
# mnew1_predict = predict(mnew1)
# tt_fasting <- concordance.index(
#   x = mnew1_predict,
#   surv.time = z.new$fu_years,
#   surv.event = z.new$status,
#   method = "noether",
#   na.rm = TRUE
# )
# 
# # tt_fasting$c.index
# # tt_fasting$se
# # tt_fasting$lower
# # concordance(mnew1)
# 
# library(compareC)
# compareC(z.std$fu_years, z.std$status, mstd_predict, mnew1_predict)
# 
# # added postprandial
# mnew2 = coxph(Surv(fu_years, status) ~ sexe + leeftijd + bmim + cluster_postprandial,  
#              data = z.new2, 
#              x = TRUE)
# mnew2_predict = predict(mnew2)
# 
# tt_postprandial <- concordance.index(
#   x = mnew2_predict,
#   surv.time = z.new$fu_years,
#   surv.event = z.new$status,
#   method = "noether",
#   na.rm = TRUE
# )
# # tt_postprandial$c.index
# # concordance(mnew2)
# 
# # combine fasting and postprandial metabotypes
# mnew3 = coxph(Surv(fu_years, status) ~ sexe + leeftijd + bmim + cluster_postprandial + cluster_fasting,  
#               data = z.new2, 
#               x = TRUE)
# mnew3_predict = predict(mnew3)
# 
# tt_combined <- concordance.index(
#   x = mnew3_predict,
#   surv.time = z.new$fu_years,
#   surv.event = z.new$status,
#   method = "noether",
#   na.rm = TRUE
# )
# # tt_combined$c.index
# # concordance(mnew3)
# 
# cindex.comp(tt_fasting, tt)
# cindex.comp(tt_postprandial, tt)
# cindex.comp(tt_combined, tt)
# cindex.comp(tt_combined, tt_postprandial)
# cindex.comp(tt_combined, tt_fasting)
# 
# compareC(z.std$fu_years, z.std$status, mstd_predict, mnew1_predict)
# 0.0302948 - 0.0302948/2.626426 * 1.96
# 0.0302948 + 0.0302948/2.626426 * 1.96
# 
# compareC(z.std$fu_years, z.std$status, mstd_predict, mnew2_predict)
# 0.02119123 - 0.02119123/2.128814 * 1.96
# 0.02119123 + 0.02119123/2.128814 * 1.96
# 
# compareC(z.std$fu_years, z.std$status, mstd_predict, mnew3_predict)
# 0.04155344 - 0.04155344/3.494992 * 1.96
# 0.04155344 + 0.04155344/3.494992 * 1.96
# 
# # compareC(z.std$fu_years, z.std$status, mnew1_predict, mnew3_predict)
# # compareC(z.std$fu_years, z.std$status, mnew2_predict, mnew3_predict)
# 
#### continuous IDI --------------------
# library(survIDINRI)
# metadata_cox[, c("fu_years", "status")]
# 
# metadata_cox1 = metadata_cox[complete.cases(metadata_cox[, c("fu_years", "status", "sexe", "leeftijd", "bmim", "cluster_fasting")]),]
# 
# covs0 <- model.matrix(~ 0 + sexe + leeftijd + bmim, 
#                       data = metadata_cox1)[, -1]
# covs1 <- model.matrix(~ 0 + sexe + leeftijd + bmim + cluster_fasting, 
#                       data = metadata_cox1)[,-1]
# covs2 <- model.matrix(~ 0 + sexe + leeftijd + bmim + cluster_postprandial, 
#                       data = metadata_cox1)[,-1]
# covs3 <- model.matrix(~ 0 + sexe + leeftijd + bmim + cluster_postprandial + cluster_fasting, 
#                       data = metadata_cox1)[,-1]
# 
# x1 <- IDI.INF(metadata_cox1[, c("fu_years", "status")], 
#         covs0 = covs0,
#         # compare with covs1 (fasting)
#         covs1 = covs1,
#         t0 = 5,
#         seed1 = 2025,
#         npert = 1000)
# round(IDI.INF.OUT(x1),2)
# 
# x1_1 <- IDI.INF(metadata_cox1[, c("fu_years", "status")], 
#               covs0 = covs0,
#               # compare with covs1 (fasting)
#               covs1 = covs1,
#               t0 = 7,
#               seed1 = 2025,
#               npert = 1000)
# round(IDI.INF.OUT(x1_1),2)
# 
# # m1: results of IDI
# # m2: results of continuous NRI
# # NRI 0.5342864  0.3910667307 0.6085248
# # m3: results of median improvement in risk score
# 
# x2 <- IDI.INF(metadata_cox1[, c("fu_years", "status")], 
#              covs0 = covs0,
#              # compare with covs2 (postprandial)
#              covs1 = covs2,
#              t0 = 5,
#              seed1 = 2025,
#              npert = 1000)
# round(IDI.INF.OUT(x2),2)
# 
# x2_1 <- IDI.INF(metadata_cox1[, c("fu_years", "status")], 
#               covs0 = covs0,
#               # compare with covs2 (postprandial)
#               covs1 = covs2,
#               t0 = 7,
#               seed1 = 2025,
#               npert = 1000)
# round(IDI.INF.OUT(x2_1),2)
# 
# # compare fasting and postprandial
# x3 <- IDI.INF(metadata_cox1[, c("fu_years", "status")], 
#               # covs1 (fasting)
#               covs0 = covs1,
#               # compare with covs3 (fasting + postprandial)
#               covs1 = covs3,
#               t0 = 5,
#               seed1 = 2025,
#               npert = 1000)
# round(IDI.INF.OUT(x3),2)
# 
# x3_1 <- IDI.INF(metadata_cox1[, c("fu_years", "status")], 
#               # covs1 (fasting)
#               covs0 = covs1,
#               # compare with covs3 (fasting + postprandial)
#               covs1 = covs3,
#               t0 = 7,
#               seed1 = 2025,
#               npert = 1000)
# round(IDI.INF.OUT(x3_1),2)
# 
# # IDI = sprintf("%.2f,(%.2f,%.2f)", temp[1,1], temp[1,2], temp[1,3])
# # NRI = sprintf("%.2f,(%.2f,%.2f)", temp[2,1], temp[2,2], temp[2,3])
# 
# # Make a table for the results (STable9)
# sTable9 <- data.frame(matrix(ncol = 4, nrow = 0))
# colnames(sTable9) <- c("Basic", "Basic+Fasting", "Baic+Postprandial", "Basic+ Fasting + Postprandial")
# sTable9[1,1] <- paste0(round(tt$c.index,2), "(", round(tt$lower,2), ",", round(tt$upper,2), ")")
# sTable9[1,2] <- paste0(round(tt_fasting$c.index,2), "(", round(tt_fasting$lower,2), ",", round(tt_fasting$upper,2), ")")
# sTable9[1,3] <- paste0(round(tt_postprandial$c.index,2), "(", round(tt_postprandial$lower,2), ",", round(tt_postprandial$upper,2), ")")
# sTable9[1,4] <- paste0(round(tt_combined$c.index,2), "(", round(tt_combined$lower,2), ",", round(tt_combined$upper,2), ")")
# 
# sTable9[2,1] <- "ref"
# sTable9[2,2] <- paste0(round(0.0302948, 2), "(", round(0.0302948 - 0.0302948/2.626426 * 1.96, 2), ",", 
#                        round(0.0302948 + 0.0302948/2.626426 * 1.96, 2), ")") # extract from the compareC function
# sTable9[2,3] <- paste0(round(0.02119123, 2), "(", round(0.02119123 - 0.02119123/2.128814 * 1.96, 2), ",", 
#                        round(0.02119123 - 0.02119123/2.128814 * 1.96, 2), ")")
# sTable9[2,4] <- paste0(round(0.04155344, 2), "(", round(0.04155344 - 0.04155344/3.494992 * 1.96), ",", 
#                        round(0.04155344 - 0.04155344/3.494992 * 1.96, 2), ")")
# sTable9[3,1] <- "ref"
# sTable9[3,2] <- paste0(round(IDI.INF.OUT(x1),2)[1,1], "(",round(IDI.INF.OUT(x1),2)[1,2], ",", round(IDI.INF.OUT(x1),2)[1,3], ")") # extract from the compareC function
# sTable9[3,3] <- paste0(round(IDI.INF.OUT(x2),2)[1,1], "(",round(IDI.INF.OUT(x2),2)[1,2], ",", round(IDI.INF.OUT(x2),2)[1,3], ")") # extract from the compareC function
# sTable9[3,4] <- paste0(round(IDI.INF.OUT(x3),2)[1,1], "(",round(IDI.INF.OUT(x3),2)[1,2], ",", round(IDI.INF.OUT(x3),2)[1,3], ")") # extract from the compareC function
# 
# sTable9[4,1] <- "ref"
# sTable9[4,2] <- paste0(round(IDI.INF.OUT(x1_1),2)[1,1], "(",round(IDI.INF.OUT(x1_1),2)[1,2], ",", round(IDI.INF.OUT(x1_1),2)[1,3], ")") # extract from the compareC function
# sTable9[4,3] <- paste0(round(IDI.INF.OUT(x2_1),2)[1,1], "(",round(IDI.INF.OUT(x2_1),2)[1,2], ",", round(IDI.INF.OUT(x2_1),2)[1,3], ")") # extract from the compareC function
# sTable9[4,4] <- paste0(round(IDI.INF.OUT(x3_1),2)[1,1], "(",round(IDI.INF.OUT(x3_1),2)[1,2], ",", round(IDI.INF.OUT(x3_1),2)[1,3], ")") # extract from the compareC function
# 
# sTable9[5,1] <- "ref"
# sTable9[5,2] <- paste0(round(IDI.INF.OUT(x1),2)[2,1], "(",round(IDI.INF.OUT(x1),2)[2,2], ",", round(IDI.INF.OUT(x1),2)[2,3], ")") # extract from the compareC function
# sTable9[5,3] <- paste0(round(IDI.INF.OUT(x2),2)[2,1], "(",round(IDI.INF.OUT(x2),2)[2,2], ",", round(IDI.INF.OUT(x2),2)[2,3], ")") # extract from the compareC function
# sTable9[5,4] <- paste0(round(IDI.INF.OUT(x3),2)[2,1], "(",round(IDI.INF.OUT(x3),2)[2,2], ",", round(IDI.INF.OUT(x3),2)[2,3], ")") # extract from the compareC function
# 
# sTable9[6,1] <- "ref"
# sTable9[6,2] <- paste0(round(IDI.INF.OUT(x1_1),2)[2,1], "(",round(IDI.INF.OUT(x1_1),2)[2,2], ",", round(IDI.INF.OUT(x1_1),2)[2,3], ")") # extract from the compareC function
# sTable9[6,3] <- paste0(round(IDI.INF.OUT(x2_1),2)[2,1], "(",round(IDI.INF.OUT(x2_1),2)[2,2], ",", round(IDI.INF.OUT(x2_1),2)[2,3], ")") # extract from the compareC function
# sTable9[6,4] <- paste0(round(IDI.INF.OUT(x3_1),2)[2,1], "(",round(IDI.INF.OUT(x3_1),2)[2,2], ",", round(IDI.INF.OUT(x3_1),2)[2,3], ")") # extract from the compareC function6
# 
# rownames(sTable9) <- c("c-index (95%CI)",
#                        "c-index increment (95%CI)",
#                        "5-year IDI(95%CI)",
#                        "7-year IDI(95%CI)",
#                        "5-year NRI(95%CI)",
#                        "7-year NRI(95%CI")
# 
# sTable9 <- sTable9 %>% 
#   dplyr::mutate(col = rownames(.), .before = 1)
# 
# write.xlsx(sTable9, file = "./Results/STable9.xlsx")

#### calculate time-dependent AUC value --------------
mnew_fasting_timeROC = coxph(Surv(fu_years, status) ~ sexe + leeftijd + bmim + cluster_fasting,  
                     data = metadata_cox1, 
                     x = TRUE)
metadata_cox1$mnew_fasting_timeROC_predict = predict(mnew_fasting_timeROC)

ROC.fasting <- timeROC(T = metadata_cox1$fu_years,
                  delta = metadata_cox1$status,
                  marker = as.numeric(metadata_cox1$mnew_fasting_timeROC_predict),
                  # other_markers = covs0,
                  cause=1,
                  weighting = "marginal",
                  # times = quantile(metadata_cox1$fu_years,
                  #                probs = seq(0.2,0.8,0.1)),
                  times = c(1, 3, 5, 7, 9, 11),
                  iid = TRUE)
ROC.fasting

mnew_postprandia_timeROC = coxph(Surv(fu_years, status) ~ sexe + leeftijd + bmim + cluster_postprandial,  
                             data = metadata_cox1, 
                             x = TRUE)
metadata_cox1$mnew_postprandia_timeROC_predict = predict(mnew_postprandia_timeROC)

ROC.postprandial <- timeROC(T = metadata_cox1$fu_years,
                       delta = metadata_cox1$status,
                       marker = as.numeric(metadata_cox1$mnew_postprandia_timeROC_predict),
                       # other_markers = covs0,
                       cause = 1,
                       weighting = "marginal",
                       # times = quantile(metadata_cox1$fu_years,
                       #                probs=seq(0.2,0.8,0.1)),
                       times = c(1, 3, 5, 7, 9, 11),
                       iid = TRUE)
ROC.postprandial 

mnew_timeROC = coxph(Surv(fu_years, status) ~ sexe + leeftijd + bmim + cluster_postprandial + cluster_fasting,  
              data = metadata_cox1, 
              x = TRUE)
metadata_cox1$mnew_timeROC_predict = predict(mnew_timeROC)

ROC.combined <- timeROC(T = metadata_cox1$fu_years,
                            delta = metadata_cox1$status,
                            marker = as.numeric(metadata_cox1$mnew_timeROC_predict),
                            cause = 1,
                            weighting = "marginal",
                            # times = quantile(metadata_cox1$fu_years,
                            #                  probs=seq(0.2,0.8,0.1)),
                            times = c(1, 3, 5, 7, 9, 11),
                            iid = TRUE)
ROC.combined 

# SFigure13
pdf(file = "Results/SFigure13.pdf", width = 6, height = 6)
plotAUCcurve(ROC.fasting,
             conf.int=FALSE,col="red")
plotAUCcurve(ROC.postprandial,
             conf.int=FALSE,col="blue", add = TRUE)
plotAUCcurve(ROC.combined,
             conf.int=FALSE,col="black", add = TRUE)

legend(x = 0.9, y = 1,c("Base + Fasting metabotypes","Base + Postprandial metabotypes", 
                        "Base + Fasting + Postprandial metabotypes"),
       col=c("red","blue", "black"),lwd=2,
       inset=0.02)
dev.off()
# compare(ROC.fasting, ROC.postprandial,adjusted=TRUE)


## compare the assignment for different cluster assignment -----------
head(metadata_cox)

# STable8
cross_tab <- metadata_cox %>% 
  select(cluster_fasting, cluster_postprandial) %>% 
  table(.) %>% 
  addmargins(.) %>% 
  as.data.frame() 
write.xlsx(cross_tab, file = "./Results/STable8.xlsx")

metadata_cox %>% 
  select(cluster_fasting, cluster_postprandial) %>% 
  table(.) %>% 
  prop.table(., margin = 1) * 100
