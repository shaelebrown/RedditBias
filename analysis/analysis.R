
# script to analyze Reddit bias data

# clean environment and initialize
rm(list = ls())
gc()
set.seed(123)
options(scipen = 999)

# download necessary libraries
install.packages(c("devtools","reticulate","igraph","rstudioapi"))
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

# compute persistence diagrams with boostrapping
ripser <- import_ripser()
race_PH <- bootstrap_persistence_thresholds(race_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.2,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T)
gender_PH <- bootstrap_persistence_thresholds(gender_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.3,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T)
orientation_PH <- bootstrap_persistence_thresholds(orientation_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.2,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T)
religion_PH <- bootstrap_persistence_thresholds(religion_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.2,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T)



