library(iPVs)


## load your dataset (mine is a flat text file tab-delimited

n = "....secret..../research/projects/fgfp/bmi_metabolon_mr/data/study_data_v1.2.Rdata"
load(n)
mydata = study_data$metabolon_data


## perform any QC that you would like
## I will filter on feature missingness
fmis = apply(mydata, 2, function(x){  sum(is.na(x))/length(x) })
w = which(fmis > 0.2)
wdata = mydata[, -w]

## NO data transformation is required
## as we will be using the non-parametric 
## Spearman's rank correlation. 

#######################################
## A quick and easy wrapper function to
## do everything for you is:
#######################################
## a series of cut heights
s = seq(0.1, 0.5, by = 0.1)
mypvs = iPVs(wdata, 
             cor_method = "spearman", 				## you can choose spearman or pearson
             dist_method = "R", 					    ## you can choose 'R' (1-abs(R)) or 'R2' (R*R)
             hclust_meth = "complete", 				## you can choose 'complete', 'average' (UPGMA), or 'mcquitty' (WPGMA)
             cutheight  = 0.5,     ## 1 or a series of cut heights
             cmat = NULL) 			                ## you can provide your own correlation matrix, which can save time.



###################################################
## Exploring the Table of PVs
###################################################
mypvs$iPV_table |> arrange( desc(clustersize) ) |> head()

mypvs$PVresults$compid_52285 |> head()


mypvs$PVresults$compid_52285 |> filter(variable != "compid_52285") |> arrange( desc(cum_vexp) )

###################################################
##
## Proportion total variance explained by PVs
##
###################################################
### Cluster Variance: total variance for each cluster
cluster_variance = sapply(mypvs$PVresults, function(x){
  
  n = x$variable
  if(length(n) == 1){
    totalvar = var(wdata[, n ], na.rm = TRUE)
  } else {
    totalvar =  sum( apply( wdata[, n ], 2, function(y){ 
      var(y, na.rm = TRUE) 
    }) )
  }
  return(totalvar)
})

## Total cluster variance
total_cluster_variance = sum(cluster_variance)

## explained variance: variance explained by each PV
explained_variance = mypvs$iPV_table$VarExp_by_PV

## Total explained variance
total_explained_variance = sum(cluster_variance * explained_variance)

## Proportion of total variance explained
proportion_total_explained = total_explained_variance / total_cluster_variance

### SUMMARY:
###     There are 411 PVs in the Metabolon Data set at a cut height of 0.5
###.    Those 411 PVs explain (proportion_total_explained * 100 = ) 45.09% of all the variance among the 1057 metabolites in the data. 




######################################################
## Looking at different tree cut heights.
## I am going to guess that the VarExp
## decreases as the cut height increases.
######################################################

s = seq(0.1, 0.8, by = 0.1)
mypvs = iPVs(wdata, 
             cor_method = "spearman", 				## you can choose spearman or pearson
             dist_method = "R", 					    ## you can choose 'R' (1-abs(R)) or 'R2' (R*R)
             hclust_meth = "complete", 				## you can choose 'complete', 'average' (UPGMA), or 'mcquitty' (WPGMA)
             cutheight  = s,     ## 1 or a series of cut heights
             cmat = NULL) 			                ## you can provide your own correlation matrix, which can save time.

## write a quick function
est_var_exp = function(ipvs_object){
  cluster_variance = sapply(ipvs_object$PVresults, function(x){
    
    n = x$variable
    if(length(n) == 1){
      totalvar = var(wdata[, n ], na.rm = TRUE)
    } else {
      totalvar =  sum( apply( wdata[, n ], 2, function(y){ 
        var(y, na.rm = TRUE) 
      }) )
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

## For each tree cut height return the proportion of variance explained by the PVs
PTVE = sapply(mypvs, est_var_exp)
PTVE = data.frame(cutheight = 1:8/10, VarExp = PTVE)

PTVE |> ggplot(aes(x = cutheight, y = VarExp)) +
  geom_point() + 
  geom_smooth(method = 'lm', formula = 'y~x') +
  geom_abline(intercept = 1, slope = -1) +
  theme_bw() +
  labs(x = "tree cut height or (1 - rho)", 
       y = "total variance explained",
       title = "Total Variance Explained by Selected Principal Variables",
       subtitle = "Metablon data")

