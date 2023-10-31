
# script to analyze Reddit bias data

# clean environment and initialize
rm(list = ls())
gc()
set.seed(123)
options(scipen = 999)

# download necessary libraries
install.packages(c("devtools","reticulate","igraph","rstudioapi","TDA"))
devtools::install_github("shaelebrown/TDApplied")
reticulate::py_install("ripser")

# set working directory to analysis directory
path <- rstudioapi::getSourceEditorContext()$path
path <- strsplit(path,split = "analysis.R")[[1]][[1]]
setwd(path)

# verify that python can be found
reticulate::py_config()
reticulate::py_available()

# load TDApplied
library(TDApplied)

# read in four datasets
race_df <- read.csv("../data/race_df.csv")
race_df$X <- NULL
gender_df <- read.csv("../data/gender_df.csv")
gender_df$X <- NULL
orientation_df <- read.csv("../data/orientation_df.csv")
orientation_df$X <- NULL
religion_df <- read.csv("../data/religion_df.csv")
religion_df$X <- NULL
all_df <- do.call(rbind,list(race_df[sample(1:nrow(race_df),size = round(0.25*nrow(race_df))),],
                             gender_df[sample(1:nrow(gender_df),size = round(0.25*nrow(gender_df))),],
                             orientation_df[sample(1:nrow(orientation_df),size = round(0.25*nrow(orientation_df))),],
                             religion_df[sample(1:nrow(religion_df),size = round(0.25*nrow(religion_df))),]))
all_df <- read.csv("../data/all_df.csv")

# compute persistence diagrams with boostrapping
ripser <- import_ripser()
race_PH <- bootstrap_persistence_thresholds(race_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.2,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T,alpha = 0.1)
race_PH <- readRDS("race_PH.rds")
gender_PH <- bootstrap_persistence_thresholds(gender_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.3,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T,alpha = 0.1)
gender_PH <- readRDS("gender_PH.rds")
orientation_PH <- bootstrap_persistence_thresholds(orientation_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.2,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T,alpha = 0.1)
orientation_PH <- readRDS("orientation_PH.rds")
religion_PH <- bootstrap_persistence_thresholds(religion_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.2,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T,alpha = 0.1)
religion_PH <- readRDS("religion_PH.rds")
all_PH <- bootstrap_persistence_thresholds(all_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.45,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T,alpha = 0.1)
all_PH <- readRDS("all_PH.rds")

# plot diagrams with and without thresholds
plot_diagram(D = race_PH$diag,title = "Race",max_radius = 1.3)
plot_diagram(D = race_PH$diag,title = "Race with thresholds",thresholds = race_PH,max_radius = 1.3)
plot_diagram(D = gender_PH$diag,title = "Gender",max_radius = 1.3)
plot_diagram(D = gender_PH$diag,title = "Gender with thresholds",thresholds = gender_PH,max_radius = 1.3)
plot_diagram(D = orientation_PH$diag,title = "Orientation",max_radius = 1.3)
plot_diagram(D = orientation_PH$diag,title = "Orientation with thresholds",thresholds = orientation_PH,max_radius = 1.3)
plot_diagram(D = religion_PH$diag,title = "Religion",max_radius = 1.3)
plot_diagram(D = religion_PH$diag,title = "Religion with thresholds",thresholds = religion_PH,max_radius = 1.3)
plot_diagram(D = all_PH$diag,title = "All data")
plot_diagram(D = all_PH$diag,title = "All data with thresholds",thresholds = all_PH)

# how many significant features existed in each dimension?
table(race_PH$subsetted_diag$dimension)
table(gender_PH$subsetted_diag$dimension)
table(orientation_PH$subsetted_diag$dimension)
table(religion_PH$subsetted_diag$dimension)
table(all_PH$subsetted_diag$dimension)

# one 2-sphere in race df and one loop in religion df
# most clusters in orientation df then gender df
# all data had two clusters and one loop, what were these features?
D_all <- as.matrix(dist(all_df))
all_reps <- TDA::ripsDiag(X = D_all,maxdimension = 1,maxscale = 1.45,dist = "arbitrary",library = "dionysus",location = T)
saveRDS(all_reps,"all_reps.rds")
all_reps <- readRDS("all_reps.rds")
H0_death_1 <- all_PH$subsetted_diag$death[[1]]
H0_death_2 <- all_PH$subsetted_diag$death[[2]]
H1_birth <- all_PH$subsetted_diag$birth[[3]]
H1_death <- all_PH$subsetted_diag$death[[3]]
H0_ind_1 <- 144
H0_ind_2 <- 142
H1_ind <- 563
H0_1 <- unique(as.vector(all_reps$cycleLocation[[H0_ind_1]]))
H0_2 <- unique(as.vector(all_reps$cycleLocation[[H0_ind_2]]))
H1 <- unique(as.vector(all_reps$cycleLocation[[H1_ind]]))

# plot using VR graphs
g <- vr_graphs(X = all_df,eps = c(0.99*H0_death_1,0.99*H0_death_2,H1_birth))
biases <- c(rep("blue",round(0.25*nrow(race_df))),
            rep("red",round(0.25*nrow(gender_df))),
            rep("green",round(0.25*nrow(orientation_df))),
            rep("yellow",round(0.25*nrow(religion_df))))
# race is blue, gender is red, orientation is green, religion is yellow
plot_vr_graph(g,eps = H1_birth,vertex_labels = F,component_of = H1[[1]],cols = biases,title = "Orientation loop") 
# orientation loop, which was not considered significant just in the orientation dataset
plot_vr_graph(g,eps = H1_birth,vertex_labels = F,cols = biases,title = "All data at loop scale") 
# all data at the same scale
plot_vr_graph(g,eps = 0.99*H0_death_1,vertex_labels = F,cols = biases,title = "All data at cluster death scale")
# all groups are separable at this scale! but mostly gender
# how separable are the classes by neighbors?
edges = g$graphs[[1]]$graph
neighborhood_selectivity <- unlist(lapply(X = g$vertices,FUN = function(X){
  
  X <- as.numeric(X)
  neighbors <- as.numeric(unique(c(edges[which(edges[,1] == X),2L],edges[which(edges[,2] == X),1L])))
  if(length(neighbors) == 0)
  {
    return(0)
  }
  bias <- biases[[X]]
  neighbor_bias <- biases[neighbors]
  return(length(which(neighbor_bias == bias))/length(neighbor_bias))
  
}))
length(which(neighborhood_selectivity < 0.25))/length(neighborhood_selectivity) # about 0.4 percent
length(which(neighborhood_selectivity < 0.5))/length(neighborhood_selectivity) # about 2 percent
length(which(neighborhood_selectivity == 1))/length(neighborhood_selectivity) # about 28 percent
# so accuracy of this geometric model would be about 99.6% on this data for majority vote
# or 98% on 50% vote

# what about in the remainder of the data?
res <- lapply(X = list(race_df,gender_df,orientation_df,religion_df),FUN = function(X){
  
  df <- as.data.frame(t(setdiff(as.data.frame(t(X)),as.data.frame(t(all_df)))))
  return(list(df,nrow(df)))
  
})
test_df <- do.call(rbind,lapply(res,"[[",1)) # missing a couple rows due to rounding
test_biases <- c(rep("blue",res[[1]][[2]]),rep("red",res[[2]][[2]]),rep("green",res[[3]][[2]]),rep("yellow",res[[4]][[2]]))
g_test <- vr_graphs(X = test_df,eps = c(0.99*H0_death_1))
edges_test = g_test$graphs[[1]]$graph
neighborhood_selectivity_test <- unlist(lapply(X = 1:length(g_test$vertices),FUN = function(X){
  
  Y <- g_test$vertices[[X]]
  neighbors <- unique(c(edges_test[which(edges_test[,1] == Y),2L],edges_test[which(edges_test[,2] == Y),1L]))
  if(length(neighbors) == 0)
  {
    return(0)
  }
  neighbors <- match(neighbors,g_test$vertices)
  bias <- test_biases[[X]]
  neighbor_bias <- test_biases[neighbors]
  return(length(which(neighbor_bias == bias))/length(neighbor_bias))
  
}))
length(which(neighborhood_selectivity_test < 0.25))/length(neighborhood_selectivity_test) # about 0.3 percent
length(which(neighborhood_selectivity_test < 0.5))/length(neighborhood_selectivity_test) # about 1.7 percent
length(which(neighborhood_selectivity_test == 1))/length(neighborhood_selectivity_test) # about 20 percent
# so about 99.7% accurate on test set by majority vote, or 98.3% on 50% vote.
